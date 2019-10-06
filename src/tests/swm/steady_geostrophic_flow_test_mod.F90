module steady_geostrophic_flow_test_mod

  use flogger
  use const_mod
  use parallel_mod
  use mesh_mod
  use state_mod
  use static_mod

  implicit none

  private

  public steady_geostrophic_flow_test_set_initial_condition

  real(r8), parameter :: alpha = 0.0_r8
  real(r8), parameter :: u0 = 2.0_r8 * pi * radius / (12.0_r8 * 86400.0_r8)
  real(r8), parameter :: gd0 = 2.94e4_r8 ! m2 s-2

contains

  subroutine steady_geostrophic_flow_test_set_initial_condition()

    real(r8) cos_lat, sin_lat, cos_lon, sin_lon, cos_alpha, sin_alpha
    integer i, j

    call log_notice('Use steady geostrophic flow initial condition.')

    cos_alpha = cos(alpha)
    sin_alpha = sin(alpha)
    static%ghs(:,:) = 0.0

    do j = states(1)%mesh%full_lat_start_idx, states(1)%mesh%full_lat_end_idx
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = states(1)%mesh%half_lon_start_idx, states(1)%mesh%half_lon_end_idx
        cos_lon = mesh%half_cos_lon(i)
        states(1)%u(i,j) = u0 * (cos_lat * cos_alpha + cos_lon * sin_lat * sin_alpha)
      end do
    end do

    call parallel_fill_halo(states(1)%mesh, states(1)%u)

    do j = states(1)%mesh%half_lat_start_idx, states(1)%mesh%half_lat_end_idx
      do i = states(1)%mesh%full_lon_start_idx, states(1)%mesh%full_lon_end_idx
        sin_lon = mesh%full_cos_lon(i)
        states(1)%v(i,j) = - u0 * sin_lon * sin_alpha
      end do
    end do

    call parallel_fill_halo(states(1)%mesh, states(1)%v)

    do j = states(1)%mesh%full_lat_start_idx, states(1)%mesh%full_lat_end_idx
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = states(1)%mesh%full_lon_start_idx, states(1)%mesh%full_lon_end_idx
        cos_lon = mesh%full_cos_lon(i)
        states(1)%gd(i,j) = gd0 - (radius * omega * u0 + u0**2 * 0.5) * (sin_lat * cos_alpha - cos_lon * cos_lat * sin_alpha)**2
      end do
    end do

    call parallel_fill_halo(states(1)%mesh, states(1)%gd)

  end subroutine steady_geostrophic_flow_test_set_initial_condition

end module steady_geostrophic_flow_test_mod
