module operators_mod

  use const_mod
  use mesh_mod
  use static_mod
  use state_mod
  use tend_mod
  use parallel_mod
  use namelist_mod
  use log_mod
  use pv_mod
  use ke_mod
  use reduce_mod

  implicit none

  private

  public operators_prepare
  public nonlinear_coriolis_operator
  public energy_gradient_operator
  public mass_flux_divergence_operator
  public calc_mass_on_edge   ! Shold be called when mass is updated.
  public calc_mass_on_vertex ! Same as above.

contains

  subroutine operators_prepare(state)

    type(state_type), intent(inout) :: state

    call calc_normal_mass_flux(state)
    call calc_tangent_mass_flux(state)
    call calc_pv_on_vertex(state)
    call calc_ke_on_cell(state)

  end subroutine operators_prepare

  subroutine calc_mass_on_edge(state)

    type(state_type), intent(inout) :: state

    integer i, j

    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        state%m_lon(i,j) = (state%mesh%lon_edge_left_area (j) * state%gd(i,  j) + &
                            state%mesh%lon_edge_right_area(j) * state%gd(i+1,j)   &
                           ) / state%mesh%lon_edge_area(j) / g
      end do
    end do

    do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        state%m_lat(i,j) = (state%mesh%lat_edge_up_area  (j) * state%gd(i,j  ) + &
                            state%mesh%lat_edge_down_area(j) * state%gd(i,j-1)   &
                           ) / state%mesh%lat_edge_area(j) / g
#else
        state%m_lat(i,j) = (state%mesh%lat_edge_up_area  (j) * state%gd(i,j+1) + &
                            state%mesh%lat_edge_down_area(j) * state%gd(i,j  )   &
                           ) / state%mesh%lat_edge_area(j) / g
#endif
      end do
    end do

  end subroutine calc_mass_on_edge

  subroutine calc_mass_on_vertex(state)

    type(state_type), intent(inout) :: state

    integer i, j
    real(r8) pole

    do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        state%m_vtx(i,j) = ((state%gd(i,j-1) + state%gd(i+1,j-1)) * state%mesh%subcell_area(2,j-1) + &
                            (state%gd(i,j  ) + state%gd(i+1,j  )) * state%mesh%subcell_area(1,j  )   &
                           ) / state%mesh%vertex_area(j) / g
#else
        state%m_vtx(i,j) = ((state%gd(i,j  ) + state%gd(i+1,j  )) * state%mesh%subcell_area(2,j  ) + &
                            (state%gd(i,j+1) + state%gd(i+1,j+1)) * state%mesh%subcell_area(1,j+1)   &
                           ) / state%mesh%vertex_area(j) / g
#endif
      end do
    end do
#ifdef STAGGER_V_ON_POLE
    if (state%mesh%has_south_pole()) then
      j = state%mesh%half_lat_start_idx
      pole = 0.0_r8
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        pole = pole + state%gd(i,j)
      end do
      call parallel_zonal_sum(pole)
      pole = pole / state%mesh%num_half_lon / g
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        state%m_vtx(i,j) = pole
      end do
    end if
    if (state%mesh%has_north_pole()) then
      j = state%mesh%half_lat_end_idx
      pole = 0.0_r8
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        pole = pole + state%gd(i,j-1)
      end do
      call parallel_zonal_sum(pole)
      pole = pole / state%mesh%num_half_lon / g
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        state%m_vtx(i,j) = pole
      end do
    end if
#endif

  end subroutine calc_mass_on_vertex

  subroutine calc_normal_mass_flux(state)

    type(state_type), intent(inout) :: state

    integer i, j

    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        state%mf_lon_n(i,j) = state%m_lon(i,j) * state%u(i,j)
      end do
    end do

    do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        state%mf_lat_n(i,j) = state%m_lat(i,j) * state%v(i,j)
      end do
    end do

    call parallel_fill_halo(state%mesh, state%mf_lon_n)
    call parallel_fill_halo(state%mesh, state%mf_lat_n)

  end subroutine calc_normal_mass_flux

  subroutine calc_tangent_mass_flux(state)

    type(state_type), intent(inout), target :: state

    integer i, j

    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        state%mf_lon_t(i,j) = state%mesh%full_tangent_wgt(1,j) * (state%mf_lat_n(i,j  ) + state%mf_lat_n(i+1,j  )) + &
                              state%mesh%full_tangent_wgt(2,j) * (state%mf_lat_n(i,j+1) + state%mf_lat_n(i+1,j+1))
#else
        state%mf_lon_t(i,j) = state%mesh%full_tangent_wgt(1,j) * (state%mf_lat_n(i,j-1) + state%mf_lat_n(i+1,j-1)) + &
                              state%mesh%full_tangent_wgt(2,j) * (state%mf_lat_n(i,j  ) + state%mf_lat_n(i+1,j  ))
#endif
      end do
    end do

    do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        state%mf_lat_t(i,j) = state%mesh%half_tangent_wgt(1,j) * (state%mf_lon_n(i-1,j-1) + state%mf_lon_n(i,j-1)) + &
                              state%mesh%half_tangent_wgt(2,j) * (state%mf_lon_n(i-1,j  ) + state%mf_lon_n(i,j  ))
#else
        state%mf_lat_t(i,j) = state%mesh%half_tangent_wgt(1,j) * (state%mf_lon_n(i-1,j  ) + state%mf_lon_n(i,j  )) + &
                              state%mesh%half_tangent_wgt(2,j) * (state%mf_lon_n(i-1,j+1) + state%mf_lon_n(i,j+1))
#endif
      end do
    end do

    call parallel_fill_halo(state%mesh, state%mf_lat_t)
    call parallel_fill_halo(state%mesh, state%mf_lon_t)

  end subroutine calc_tangent_mass_flux

  subroutine nonlinear_coriolis_operator(state, tend, dt)

    type(state_type), intent(inout) :: state
    type(tend_type ), intent(inout) :: tend
    real(r8)        , intent(in   ) :: dt

    integer i, j

    select case (pv_scheme)
    case (1)
      call calc_pv_on_edge_midpoint(state)
    case (2)
      call calc_pv_on_edge_upwind(state)
    case (3)
      call calc_pv_on_edge_apvm(state, dt)
    case (4)
      call calc_pv_on_edge_scale_aware_apvm(state)
    case default
      call log_error('Unknown PV scheme!')
    end select

    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        tend%qhv(i,j) = state%mesh%full_tangent_wgt(1,j) * (                      &
          state%mf_lat_n(i  ,j  ) * (state%pv_lon(i,j) + state%pv_lat(i  ,j)  ) + &
          state%mf_lat_n(i+1,j  ) * (state%pv_lon(i,j) + state%pv_lat(i+1,j)  )   &
        ) * 0.5_r8
#else
        tend%qhv(i,j) = state%mesh%full_tangent_wgt(1,j) * (                      &
          state%mf_lat_n(i  ,j-1) * (state%pv_lon(i,j) + state%pv_lat(i  ,j-1)) + &
          state%mf_lat_n(i+1,j-1) * (state%pv_lon(i,j) + state%pv_lat(i+1,j-1))   &
        ) * 0.5_r8
#endif
      end do
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        tend%qhv(i,j) = tend%qhv(i,j) + state%mesh%full_tangent_wgt(2,j) * (      &
          state%mf_lat_n(i  ,j+1) * (state%pv_lon(i,j) + state%pv_lat(i  ,j+1)) + &
          state%mf_lat_n(i+1,j+1) * (state%pv_lon(i,j) + state%pv_lat(i+1,j+1))   &
        ) * 0.5_r8
#else
        tend%qhv(i,j) = tend%qhv(i,j) + state%mesh%full_tangent_wgt(2,j) * (      &
          state%mf_lat_n(i  ,j  ) * (state%pv_lon(i,j) + state%pv_lat(i  ,j  )) + &
          state%mf_lat_n(i+1,j  ) * (state%pv_lon(i,j) + state%pv_lat(i+1,j  ))   &
        ) * 0.5_r8
#endif
      end do
    end do

    do j = state%mesh%half_lat_start_idx, state%mesh%half_lat_end_idx
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        tend%qhu(i,j) = state%mesh%half_tangent_wgt(1,j) * (                      &
          state%mf_lon_n(i-1,j-1) * (state%pv_lat(i,j) + state%pv_lon(i-1,j-1)) + &
          state%mf_lon_n(i  ,j-1) * (state%pv_lat(i,j) + state%pv_lon(i  ,j-1))   &
        ) * 0.5_r8
#else
        tend%qhu(i,j) = state%mesh%half_tangent_wgt(1,j) * (                      &
          state%mf_lon_n(i-1,j  ) * (state%pv_lat(i,j) + state%pv_lon(i-1,j  )) + &
          state%mf_lon_n(i  ,j  ) * (state%pv_lat(i,j) + state%pv_lon(i  ,j  ))   &
        ) * 0.5_r8
#endif
      end do
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        tend%qhu(i,j) = tend%qhu(i,j) + state%mesh%half_tangent_wgt(2,j) * (      &
          state%mf_lon_n(i-1,j  ) * (state%pv_lat(i,j) + state%pv_lon(i-1,j  )) + &
          state%mf_lon_n(i  ,j  ) * (state%pv_lat(i,j) + state%pv_lon(i  ,j  ))   &
        ) * 0.5_r8
#else
        tend%qhu(i,j) = tend%qhu(i,j) + state%mesh%half_tangent_wgt(2,j) * (      &
          state%mf_lon_n(i-1,j+1) * (state%pv_lat(i,j) + state%pv_lon(i-1,j+1)) + &
          state%mf_lon_n(i  ,j+1) * (state%pv_lat(i,j) + state%pv_lon(i  ,j+1))   &
        ) * 0.5_r8
#endif
      end do
    end do

  end subroutine nonlinear_coriolis_operator

  subroutine energy_gradient_operator(static, state, tend)

    type(static_type), intent(in   ) :: static
    type(state_type ), intent(inout) :: state
    type(tend_type  ), intent(inout) :: tend

    integer i, j, move

    ! --------------------------------------------------------------------------
    !                           Zonal geopotential energy gradient
    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      if (reduced_full_mesh(j)%reduce_factor > 0) then
        tend%dpedlon(:,j) = 0.0_r8
        do move = 1, reduced_full_mesh(j)%reduce_factor
          do i = reduced_full_mesh(j)%full_lon_start_idx, reduced_full_mesh(j)%full_lon_end_idx
            reduced_full_tend(j)%dpedlon(i) = (                                               &
              reduced_full_state (j)%gd (i+1,0,move) - reduced_full_state (j)%gd (i,0,move) + &
              reduced_full_static(j)%ghs(i+1,0,move) - reduced_full_static(j)%ghs(i,0,move)   &
            ) / reduced_full_mesh(j)%de_lon(0)
          end do
          call reduce_append_array(move, reduced_full_mesh(j)  , &
                                   reduced_full_tend(j)%dpedlon, &
                                   state%mesh, tend%dpedlon(:,j))
        end do
        call parallel_overlay_inner_halo(state%mesh, tend%dpedlon(:,j), left_halo=.true.)
      else
        do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
          tend%dpedlon(i,j) = (                   &
            state %gd (i+1,j) - state %gd (i,j) + &
            static%ghs(i+1,j) - static%ghs(i,j)   &
          ) / state%mesh%de_lon(j)
        end do
      end if
    end do

    ! --------------------------------------------------------------------------
    !                       Zonal kinetic energy gradient
    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        tend%dkedlon(i,j) = (state%ke_cell(i+1,j) - state%ke_cell(i,j)) / state%mesh%de_lon(j)
      end do
    end do

    ! --------------------------------------------------------------------------
    !                         Meridional energy gradient
    do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        tend%dpedlat(i,j) = (state %gd (i,j) - state %gd (i,j-1) + &
                             static%ghs(i,j) - static%ghs(i,j-1)   &
                            ) / state%mesh%de_lat(j)
        tend%dkedlat(i,j) = (state%ke_cell(i,j) - state%ke_cell(i,j-1)) / state%mesh%de_lat(j)
#else
        tend%dpedlat(i,j) = (state %gd (i,j+1) - state %gd (i,j) + &
                             static%ghs(i,j+1) - static%ghs(i,j)   &
                            ) / state%mesh%de_lat(j)
        tend%dkedlat(i,j) = (state%ke_cell(i,j+1) - state%ke_cell(i,j)) / state%mesh%de_lat(j)
#endif
      end do
    end do

  end subroutine energy_gradient_operator

  subroutine mass_flux_divergence_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j, move
    real(r8) pole

    ! --------------------------------------------------------------------------
    !                       Zonal mass flux divergence
    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      if (reduced_full_mesh(j)%reduce_factor > 0) then
        tend%mf_div(:,j) = 0.0_r8
        do move = 1, reduced_full_mesh(j)%reduce_factor
          do i = reduced_full_mesh(j)%full_lon_start_idx, reduced_full_mesh(j)%full_lon_end_idx
            reduced_full_tend(j)%mf_div_lon(i) = (                                                  &
              reduced_full_state(j)%mf_lon_n(i,0,move) - reduced_full_state(j)%mf_lon_n(i-1,0,move) &
            ) * reduced_full_mesh (j)%le_lon(0) / reduced_full_mesh (j)%cell_area(0)
          end do
          call reduce_append_array(move, reduced_full_mesh(j),      &
                                   reduced_full_tend(j)%mf_div_lon, &
                                   state%mesh, tend%mf_div(:,j))
        end do
        call parallel_overlay_inner_halo(state%mesh, tend%mf_div(:,j), left_halo=.true.)
      else
        do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
          tend%mf_div(i,j) = (                          &
            state%mf_lon_n(i,j) - state%mf_lon_n(i-1,j) &
          ) * mesh %le_lon(j) / mesh %cell_area(j)
        end do
      end if
    end do

    ! --------------------------------------------------------------------------
    !                    Meridional mass flux divergence
    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        tend%mf_div(i,j) = tend%mf_div(i,j) + (                                  &
#ifdef STAGGER_V_ON_POLE
                            (state%mf_lat_n(i,j+1) * mesh%le_lat(j+1) - &
                             state%mf_lat_n(i,j  ) * mesh%le_lat(j  ))  &
#else
                            (state%mf_lat_n(i,j  ) * mesh%le_lat(j  ) - &
                             state%mf_lat_n(i,j-1) * mesh%le_lat(j-1))  &
#endif
                           ) / mesh%cell_area(j)
      end do
    end do

#ifndef STAGGER_V_ON_POLE
    if (state%mesh%has_south_pole()) then
      j = state%mesh%full_lat_start_idx
      pole = 0.0d0
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        pole = pole + state%mf_lat_n(i,j)
      end do
      call parallel_zonal_sum(pole)
      pole = pole * mesh%le_lat(j) / state%mesh%num_full_lon / mesh%cell_area(j)
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        tend%mf_div(i,j) = pole
      end do
    end if
    if (state%mesh%has_north_pole()) then
      j = state%mesh%full_lat_end_idx
      pole = 0.0d0
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        pole = pole - state%mf_lat_n(i,j-1)
      end do
      call parallel_zonal_sum(pole)
      pole = pole * mesh%le_lat(j-1) / state%mesh%num_full_lon / mesh%cell_area(j)
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        tend%mf_div(i,j) = pole
      end do
    end if
#endif

  end subroutine mass_flux_divergence_operator

end module operators_mod
