module pv_mod

  use const_mod
  use mesh_mod
  use allocator_mod
  use namelist_mod
  use parallel_mod
  use state_mod
  use tend_mod

  implicit none

  private

  public calc_pv_on_vertex
  public calc_pv_on_edge_midpoint
  public calc_pv_on_edge_upwind
  public calc_pv_on_edge_apvm
  public calc_pv_on_edge_scale_aware_apvm

contains

  subroutine calc_pv_on_vertex(state)

    type(state_type), intent(inout) :: state

    integer i, j
    real(r8) pole

    do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        state%pv(i,j) = ((state%u(i  ,j-1) * state%mesh%cell_lon_dist(j-1) - &
                          state%u(i  ,j  ) * state%mesh%cell_lon_dist(j  ) + &
                          state%v(i+1,j  ) * state%mesh%cell_lat_dist(j  ) - &
                          state%v(i  ,j  ) * state%mesh%cell_lat_dist(j  )   &
                         ) / state%mesh%vertex_area(j) + state%mesh%half_f(j)    &
                        ) / state%m_vtx(i,j)
#else
        state%pv(i,j) = ((state%u(i  ,j  ) * state%mesh%cell_lon_dist(j  ) - &
                          state%u(i  ,j+1) * state%mesh%cell_lon_dist(j+1) + &
                          state%v(i+1,j  ) * state%mesh%cell_lat_dist(j  ) - &
                          state%v(i  ,j  ) * state%mesh%cell_lat_dist(j  )   &
                         ) / state%mesh%vertex_area(j) + state%mesh%half_f(j)    &
                        ) / state%m_vtx(i,j)
#endif
      end do
    end do
#ifdef STAGGER_V_ON_POLE
    if (state%mesh%has_south_pole()) then
      j = state%mesh%half_lat_start_idx
      pole = 0.0_r8
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        pole = pole - state%u(i,j) * state%mesh%cell_lon_dist(j)
      end do
      call parallel_zonal_sum(pole)
      pole = pole / state%mesh%num_half_lon / state%mesh%vertex_area(j)
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        state%pv(i,j) = (pole + state%mesh%half_f(j)) / state%m_vtx(i,j)
      end do
    end if
    if (state%mesh%has_north_pole()) then
      j = state%mesh%half_lat_end_idx
      pole = 0.0_r8
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        pole = pole + state%u(i,j-1) * state%mesh%cell_lon_dist(j-1)
      end do
      call parallel_zonal_sum(pole)
      pole = pole / state%mesh%num_half_lon / state%mesh%vertex_area(j)
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        state%pv(i,j) = (pole + state%mesh%half_f(j)) / state%m_vtx(i,j)
      end do
    end if
#else
    if (pv_pole_stokes) then
      ! Special treatment of vorticity around Poles
      if (state%mesh%has_south_pole()) then
        j = state%mesh%half_lat_start_idx
        pole = 0.0_r8
        do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
          pole = pole - state%u(i,j+1) * state%mesh%cell_lon_dist(j+1)
        end do
        call parallel_zonal_sum(pole)
        pole = pole / state%mesh%num_half_lon / state%mesh%vertex_area(j)
        do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
          state%pv(i,j) = (pole + state%mesh%half_f(j)) / state%m_vtx(i,j)
        end do
      end if
      if (state%mesh%has_north_pole()) then
        j = state%mesh%half_lat_end_idx
        pole = 0.0_r8
        do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
          pole = pole + state%u(i,j) * state%mesh%cell_lon_dist(j)
        end do
        call parallel_zonal_sum(pole)
        pole = pole / state%mesh%num_half_lon / state%mesh%vertex_area(j)
        do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
          state%pv(i,j) = (pole + state%mesh%half_f(j)) / state%m_vtx(i,j)
        end do
      end if
    end if
#endif
    call parallel_fill_halo(state%mesh, state%pv)

  end subroutine calc_pv_on_vertex

  subroutine calc_dpv_on_edge(state)

    type(state_type), intent(inout) :: state

    integer i, j

    ! Tangent pv difference
    do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        state%dpv_lat_t(i,j) = state%pv(i,j) - state%pv(i-1,j)
      end do
    end do
    call parallel_fill_halo(state%mesh, state%dpv_lat_t)

    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        state%dpv_lon_t(i,j) = state%pv(i,j+1) - state%pv(i,j)
#else
        state%dpv_lon_t(i,j) = state%pv(i,j) - state%pv(i,j-1)
#endif
      end do
    end do
    call parallel_fill_halo(state%mesh, state%dpv_lon_t)

    ! Normal pv difference
    do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        state%dpv_lon_n(i,j) = 0.25_r8 * (state%dpv_lat_t(i-1,j-1) + state%dpv_lat_t(i,j-1) + &
                                          state%dpv_lat_t(i-1,j  ) + state%dpv_lat_t(i,j  ))
#else
        state%dpv_lon_n(i,j) = 0.25_r8 * (state%dpv_lat_t(i-1,j  ) + state%dpv_lat_t(i,j  ) + &
                                          state%dpv_lat_t(i-1,j+1) + state%dpv_lat_t(i,j+1))
#endif
      end do
    end do

    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        state%dpv_lat_n(i,j) = 0.25_r8 * (state%dpv_lon_t(i,j  ) + state%dpv_lon_t(i+1,j  ) + &
                                          state%dpv_lon_t(i,j+1) + state%dpv_lon_t(i+1,j+1))
#else
        state%dpv_lat_n(i,j) = 0.25_r8 * (state%dpv_lon_t(i,j-1) + state%dpv_lon_t(i+1,j-1) + &
                                          state%dpv_lon_t(i,j  ) + state%dpv_lon_t(i+1,j  ))
#endif
      end do
    end do

  end subroutine calc_dpv_on_edge

  subroutine calc_pv_on_edge_midpoint(state)

    type(state_type), intent(inout) :: state

    integer i, j

    do j = state%mesh%half_lat_start_idx, state%mesh%half_lat_end_idx
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        state%pv_lat(i,j) = 0.5_r8 * (state%pv(i-1,j) + state%pv(i,j))
      end do 
    end do 
      
    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        state%pv_lon(i,j) = 0.5_r8 * (state%pv(i,j) + state%pv(i,j+1))
#else
        state%pv_lon(i,j) = 0.5_r8 * (state%pv(i,j) + state%pv(i,j-1))
#endif
      end do 
    end do 

    call parallel_fill_halo(state%mesh, state%pv_lon)
    call parallel_fill_halo(state%mesh, state%pv_lat)

  end subroutine calc_pv_on_edge_midpoint

  subroutine calc_pv_on_edge_upwind(state)

    type(state_type), intent(inout) :: state

    real(r8), parameter :: beta0 = 1.0_r8
    real(r8), parameter :: dpv0 = 1.0e-9_r8
    real(r8) dpv, beta
    integer i, j

    call calc_dpv_on_edge(state)

    do j = state%mesh%half_lat_start_idx, state%mesh%half_lat_end_idx
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        ! beta = state%mesh%half_upwind_beta(j)
        dpv = min(abs(state%dpv_lat_t(i,j)), abs(state%dpv_lat_n(i,j)))
        beta = beta0 * exp(-(dpv0 / dpv)**2)
        state%pvc_lat(i,j) = beta * 0.5_r8 * &
          (state%dpv_lat_t(i,j) * sign(1.0_r8, state%mf_lat_t(i,j) + &
           state%dpv_lat_n(i,j) * sign(1.0_r8, state%v(i,j))))
        state%pv_lat(i,j) = 0.5_r8 * (state%pv(i,j) + state%pv(i-1,j)) - state%pvc_lat(i,j)
      end do
    end do

    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole 
     do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        ! beta = state%mesh%full_upwind_beta(j)
        dpv = min(abs(state%dpv_lon_t(i,j)), abs(state%dpv_lon_n(i,j)))
        beta = beta0 * exp(-(dpv0 / dpv)**2)
        state%pvc_lon(i,j) = beta * 0.5_r8 * &
          (state%dpv_lon_t(i,j) * sign(1.0_r8, state%mf_lon_t(i,j) + &
           state%dpv_lon_n(i,j) * sign(1.0_r8, state%u(i,j))))
        state%pv_lon(i,j) = 0.5_r8 * (state%pv(i,j+1) + state%pv(i,j)) - state%pvc_lon(i,j)
      end do
    end do

    call parallel_fill_halo(state%mesh, state%pv_lon)
    call parallel_fill_halo(state%mesh, state%pv_lat)

  end subroutine calc_pv_on_edge_upwind

  subroutine calc_pv_on_edge_apvm(state, dt)

    type(state_type), intent(inout) :: state
    real(r8)        , intent(in   ) :: dt

    real(r8) un, vn, ut, vt, le, de
    integer i, j

    call calc_dpv_on_edge(state)

    do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
      le = state%mesh%vertex_lon_dist(j)
      de = state%mesh%cell_lat_dist(j)
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        ut = state%mf_lat_t(i,j) / state%m_lat(i,j)
        vn = state%v(i,j)
        state%pvc_lat(i,j) = 0.5_r8 * (ut * state%dpv_lat_t(i,j) / le + vn * state%dpv_lat_n(i,j) / de) * dt
        state%pv_lat(i,j) = 0.5_r8 * (state%pv(i,j) + state%pv(i-1,j)) - state%pvc_lat(i,j)
      end do
    end do
#ifdef STAGGER_V_ON_POLE
    state%pv_lat(:,state%mesh%half_lat_start_idx) = state%pv(:,state%mesh%half_lat_start_idx)
    state%pv_lat(:,state%mesh%half_lat_end_idx  ) = state%pv(:,state%mesh%half_lat_end_idx  )
#endif

    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      le = state%mesh%vertex_lat_dist(j)
      de = state%mesh%cell_lon_dist(j)
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        un = state%u(i,j)
        vt = state%mf_lon_t(i,j) / state%m_lon(i,j)
        state%pvc_lon(i,j) = 0.5_r8 * (un * state%dpv_lon_n(i,j) / de + vt * state%dpv_lon_t(i,j) / le) * dt
#ifdef STAGGER_V_ON_POLE
        state%pv_lon(i,j) = 0.5_r8 * (state%pv(i,j+1) + state%pv(i,j)) - state%pvc_lon(i,j)
#else
        state%pv_lon(i,j) = 0.5_r8 * (state%pv(i,j-1) + state%pv(i,j)) - state%pvc_lon(i,j)
#endif
      end do
    end do

    call parallel_fill_halo(state%mesh, state%pv_lon)
    call parallel_fill_halo(state%mesh, state%pv_lat)

  end subroutine calc_pv_on_edge_apvm

  subroutine calc_pv_on_edge_scale_aware_apvm(state)

    type(state_type), intent(inout) :: state

    real(r8), parameter :: alpha = 0.0013_r8
    real(r8) un, vn, ut, vt, le, de, ke, h, pv_adv
    integer i, j

    call calc_dpv_on_edge(state)

    ke = state%total_ke**(-3.0_r8 / 4.0_r8)
    h  = state%total_m / state%mesh%total_area / g

    do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
      le = state%mesh%vertex_lon_dist(j)
      de = state%mesh%cell_lat_dist(j)
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        ut = state%mf_lat_t(i,j) / state%m_lat(i,j)
        vn = state%v(i,j)
        pv_adv = ut * state%dpv_lat_t(i,j) / le + vn * state%dpv_lat_n(i,j) / de
        state%pvc_lat(i,j) = alpha * ke * h * de**3 * abs(pv_adv) * pv_adv
        state%pv_lat (i,j) = 0.5_r8 * (state%pv(i,j) + state%pv(i-1,j)) - state%pvc_lat(i,j)
      end do
    end do
#ifdef STAGGER_V_ON_POLE
    state%pv_lat(:,state%mesh%half_lon_start_idx) = state%pv(:,state%mesh%half_lon_start_idx)
    state%pv_lat(:,state%mesh%half_lon_end_idx  ) = state%pv(:,state%mesh%half_lon_end_idx  )
#endif

    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      le = state%mesh%vertex_lat_dist(j)
      de = state%mesh%cell_lon_dist(j)
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        un = state%u(i,j)
        vt = state%mf_lon_t(i,j) / state%m_lon(i,j)
        pv_adv = un * state%dpv_lon_n(i,j) / de + vt * state%dpv_lon_t(i,j) / le
        state%pvc_lon(i,j) = alpha * ke * h * de**3 * abs(pv_adv) * pv_adv
#ifdef STAGGER_V_ON_POLE
        state%pv_lon (i,j) = 0.5_r8 * (state%pv(i,j+1) + state%pv(i,j)) - state%pvc_lon(i,j)
#else
        state%pv_lon (i,j) = 0.5_r8 * (state%pv(i,j-1) + state%pv(i,j)) - state%pvc_lon(i,j)
#endif
      end do
    end do

    call parallel_fill_halo(state%mesh, state%pv_lon)
    call parallel_fill_halo(state%mesh, state%pv_lat)

  end subroutine calc_pv_on_edge_scale_aware_apvm

end module pv_mod