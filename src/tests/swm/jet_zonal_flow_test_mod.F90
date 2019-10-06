module jet_zonal_flow_test_mod

  use flogger
  use const_mod
  use parallel_mod
  use mesh_mod
  use state_mod
  use static_mod

  implicit none

  private

  public jet_zonal_flow_test_set_initial_condition

  real(r8), parameter :: u_max = 80.0_r8
  real(r8), parameter :: lat0 = pi / 7.0_r8
  real(r8), parameter :: lat1 = pi / 2.0_r8 - lat0
  real(r8), parameter :: en = exp(-4.0_r8 / (lat1 - lat0)**2_r8)
  real(r8), parameter :: gh0 = g * 1.0e4_r8
  real(r8), parameter :: ghd = g * 120_r8
  real(r8), parameter :: lat2 = pi / 4.0_r8
  real(r8), parameter :: alpha = 1.0_r8 / 3.0_r8
  real(r8), parameter :: beta = 1.0_r8 / 15.0_r8

contains

  subroutine jet_zonal_flow_test_set_initial_condition()

    integer i, j, neval, ierr
    real(r8) abserr, lon
    type(mesh_type), pointer :: mesh

    call log_notice('Use jet zonal flow initial condition.')

    mesh => states(1)%mesh

    static%ghs = 0.0_r8

    do j = mesh%full_lat_start_idx, mesh%full_lat_end_idx
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        states(1)%u(i,j) = u_function(mesh%full_lat(j))
      end do
    end do

    call parallel_fill_halo(mesh, states(1)%u)

    states(1)%v = 0.0_r8

    do j = mesh%full_lat_start_idx, mesh%full_lat_end_idx
      if (j == mesh%full_lat_start_idx) then
        states(1)%gd(0,j) = gh0
      else
        call qags(gh_integrand, -0.5*pi, mesh%full_lat(j), 1.0e-10, 1.0e-3, states(1)%gd(0,j), abserr, neval, ierr)
        if (ierr /= 0) then
          call log_error('Failed to calculate integration at (' // to_string(i) // ',' // to_string(j) // ')!')
        end if
        states(1)%gd(0,j) = gh0 - states(1)%gd(0,j)
      end if
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        states(1)%gd(i,j) = states(1)%gd(0,j)
        ! Add perturbation.
        states(1)%gd(i,j) = states(1)%gd(i,j) + ghd * &
          cos(mesh%full_lat(j)) * &
          exp(-(merge(mesh%full_lon(i) - 2.0_r8 * pi, mesh%full_lon(i), mesh%full_lon(i) > pi)  / alpha)**2) * &
          exp(-((lat2 - mesh%full_lat(j)) / beta)**2)
      end do
    end do

    call parallel_fill_halo(mesh, states(1)%gd)

  end subroutine jet_zonal_flow_test_set_initial_condition

  real(r8) function gh_integrand(lat) result(res)

    real(r8), intent(in) :: lat

    real(r8) u, f

    u = u_function(lat)
    f = 2 * omega * sin(lat)
    res = radius * u * (f + tan(lat) / radius * u)

  end function gh_integrand

  real(r8) function u_function(lat) result(res)

    real(r8), intent(in) :: lat

    if (lat <= lat0 .or. lat >= lat1) then
      res = 0.0_r8
    else
      res = u_max / en * exp(1 / (lat - lat0) / (lat - lat1))
    end if

  end function u_function

end module jet_zonal_flow_test_mod
