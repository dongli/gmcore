module ke_mod

  use const_mod
  use mesh_mod
  use allocator_mod
  use parallel_mod
  use state_mod
  use tend_mod

  implicit none

  private

  public calc_ke_on_cell

contains

  subroutine calc_ke_on_cell(state)

    type(state_type), intent(inout) :: state

    integer i, j
    real(r8) pole

    do j = state%mesh%full_lat_start_idx, state%mesh%full_lat_end_idx
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        state%ke_cell(i,j) = (state%mesh%lon_edge_right_area(j  ) * state%u(i-1,j  )**2 + &
                              state%mesh%lon_edge_left_area (j  ) * state%u(i  ,j  )**2 + &
#ifdef STAGGER_V_ON_POLE
                              state%mesh%lat_edge_up_area   (j  ) * state%v(i  ,j  )**2 + &
                              state%mesh%lat_edge_down_area (j+1) * state%v(i  ,j+1)**2   &
#else
                              state%mesh%lat_edge_up_area   (j-1) * state%v(i  ,j-1)**2 + &
                              state%mesh%lat_edge_down_area (j  ) * state%v(i  ,j  )**2   &
#endif
                             ) / state%mesh%cell_area(j)
      end do
    end do
#ifndef STAGGER_V_ON_POLE
    if (state%mesh%has_south_pole()) then
      j = state%mesh%full_lat_start_idx
      pole = 0.0d0
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        pole = pole + state%mesh%lat_edge_down_area(j) * state%v(i,j)**2
      end do
      call parallel_zonal_sum(pole)
      pole = pole / state%mesh%num_full_lon / state%mesh%cell_area(j)
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        state%ke_cell(i,j) = pole
      end do
    end if
    if (state%mesh%has_north_pole()) then
      j = state%mesh%full_lat_end_idx
      pole = 0.0d0
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        pole = pole + state%mesh%lat_edge_up_area(j-1) * state%v(i,j-1)**2
      end do
      call parallel_zonal_sum(pole)
      pole = pole / state%mesh%num_full_lon / state%mesh%cell_area(j)
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        state%ke_cell(i,j) = pole
      end do
    end if
#endif

    call parallel_fill_halo(state%mesh, state%ke_cell)

  end subroutine calc_ke_on_cell

end module ke_mod