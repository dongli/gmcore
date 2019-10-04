module namelist_mod

  use time_mod, start_time => start_time_array, end_time => end_time_array

  implicit none

  character(256) :: case_desc = 'N/A'
  character(256) :: case_name = 'N/A'
  character(30) :: test_case = 'N/A'
  character(30) :: history_interval(1) = 'N/A'

  integer num_lon
  integer num_lat

  character(30) :: tangent_wgt_scheme = 'classic'

  integer :: pv_scheme = 2
  logical :: pv_pole_stokes = .true.

  integer :: fast_cycles = 5
  character(30) :: split_scheme = ''
  character(30) :: time_scheme = 'predict_correct'

  integer :: reduce_factors(20) = 0

  namelist /gmcore_swm_control/ &
    case_name,                  &
    test_case,                  &
    case_desc,                  &
    num_lon,                    &
    num_lat,                    &
    start_time,                 &
    end_time,                   &
    dt_in_seconds,              &
    run_hours,                  &
    run_days,                   &
    history_interval,           &
    tangent_wgt_scheme,         &
    pv_scheme,                  &
    pv_pole_stokes,             &
    split_scheme,               &
    reduce_factors

contains

  subroutine parse_namelist(file_path)

    character(*), intent(in) :: file_path

    open(10, file=file_path, status='old')
    read(10, nml=gmcore_swm_control)
    close(10)

  end subroutine parse_namelist

end module namelist_mod
