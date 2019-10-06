module debug_mod

  use const_mod
  use log_mod
  use mesh_mod
  use static_mod
  use state_mod
  use tend_mod

  private

  public debug_check_space_operators
  public debug_print_min_max

contains

  subroutine debug_check_space_operators(static, state, tend)

    type(static_type), intent(in) :: static
    type(state_type ), intent(in) :: state
    type(tend_type  ), intent(in) :: tend

    integer i, j
    type(mesh_type), pointer :: mesh
    real(r8) ip_coriolis
    real(r8) ip_energe_gradient
    real(r8) ip_mf_div

    mesh => state%mesh
    ip_coriolis = 0.0d0
    ip_energe_gradient = 0.0d0
    ip_mf_div = 0.0d0

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        ip_coriolis = ip_coriolis + tend%qhv(i,j) * state%mf_lon_n(i,j) * mesh%lon_edge_area(j) / radius**2
        ip_energe_gradient = ip_energe_gradient + tend%dEdlon(i,j) * state%mf_lon_n(i,j) * mesh%lon_edge_area(j)
      end do
    end do

    do j = mesh%half_lat_start_idx, mesh%half_lat_end_idx
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        ip_coriolis = ip_coriolis + tend%qhu(i,j) * state%mf_lat_n(i,j) * mesh%lat_edge_area(j) / radius**2
        ip_energe_gradient = ip_energe_gradient + tend%dEdlat(i,j) * state%mf_lat_n(i,j) * mesh%lat_edge_area(j)
      end do
    end do

    do j = mesh%full_lat_start_idx, mesh%full_lat_end_idx
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        ip_energe_gradient = ip_energe_gradient + tend%mf_div(i,j) * (state%gd(i,j) + static%ghs(i,j)) * mesh%cell_area(j)
        ip_mf_div = ip_mf_div + tend%mf_div(i,j) * mesh%cell_area(j)
      end do
    end do

    call log_add_diag('coriolis', ip_coriolis / radius**2)
    call log_add_diag('energy_gradient', ip_energe_gradient / radius**2)
    call log_add_diag('mf_div', ip_mf_div / radius**2)

  end subroutine debug_check_space_operators

  subroutine debug_print_min_max(array)

    real(r8), intent(in) :: array(:,:)

    print *, minval(array), maxval(array)

  end subroutine debug_print_min_max

end module debug_mod