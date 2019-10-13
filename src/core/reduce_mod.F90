module reduce_mod

  use flogger
  use const_mod
  use namelist_mod
  use sphere_geometry_mod
  use mesh_mod
  use static_mod
  use state_mod
  use parallel_mod

  implicit none

  private

  public reduced_mesh_type
  public reduced_state_type
  public reduce_init
  public reduce_run
  public reduce_append_array
  public reduce_final

  public reduced_full_mesh
  public reduced_full_static
  public reduced_full_state
  public reduced_full_tend

  type reduced_mesh_type
    integer reduce_factor
    integer halo_width
    integer num_full_lon
    integer num_half_lon
    integer full_lon_start_idx
    integer full_lon_end_idx
    integer half_lon_start_idx
    integer half_lon_end_idx
    integer full_lon_lb
    integer full_lon_ub
    integer half_lon_lb
    integer half_lon_ub
    real(r8) dlon
    real(r8), allocatable, dimension(  :) :: full_lon
    real(r8), allocatable, dimension(  :) :: half_lon
    real(r8), allocatable, dimension(  :) :: full_lat
    real(r8), allocatable, dimension(  :) :: half_lat
    real(r8), allocatable, dimension(  :) :: full_cos_lat
    real(r8), allocatable, dimension(  :) :: half_cos_lat
    real(r8), allocatable, dimension(  :) :: full_sin_lat
    real(r8), allocatable, dimension(  :) :: half_sin_lat
    real(r8), allocatable, dimension(  :) :: cell_area
    real(r8), allocatable, dimension(  :) :: lon_edge_area
    real(r8), allocatable, dimension(  :) :: lon_edge_left_area
    real(r8), allocatable, dimension(  :) :: lon_edge_right_area
    real(r8), allocatable, dimension(  :) :: lat_edge_area
    real(r8), allocatable, dimension(  :) :: lat_edge_up_area
    real(r8), allocatable, dimension(  :) :: lat_edge_down_area
    real(r8), allocatable, dimension(  :) :: vertex_area
    real(r8), allocatable, dimension(:,:) :: subcell_area
    real(r8), allocatable, dimension(  :) :: de_lon
    real(r8), allocatable, dimension(  :) :: de_lat
    real(r8), allocatable, dimension(  :) :: le_lat
    real(r8), allocatable, dimension(  :) :: le_lon
    real(r8), allocatable, dimension(:,:) :: full_tangent_wgt
    real(r8), allocatable, dimension(:,:) :: half_tangent_wgt
    real(r8), allocatable, dimension(  :) :: full_f
    real(r8), allocatable, dimension(  :) :: half_f
  contains
    final :: reduced_mesh_final
  end type reduced_mesh_type

  type(reduced_mesh_type), allocatable :: reduced_full_mesh(:)

  type reduced_static_type
    real(r8), allocatable, dimension(:,:,:) :: ghs
  contains
    final :: reduced_static_final
  end type reduced_static_type

  type(reduced_static_type), allocatable :: reduced_full_static(:)

  type reduced_state_type
    real(r8), allocatable, dimension(:,:,:) :: u
    real(r8), allocatable, dimension(:,:,:) :: v
    real(r8), allocatable, dimension(:,:,:) :: gd
    real(r8), allocatable, dimension(:,:,:) :: pv
    real(r8), allocatable, dimension(:,:,:) :: pv_lon
    real(r8), allocatable, dimension(:,:,:) :: pv_lat
    real(r8), allocatable, dimension(:,:,:) :: dpv_lon_t
    real(r8), allocatable, dimension(:,:,:) :: dpv_lat_t
    real(r8), allocatable, dimension(:,:,:) :: dpv_lon_n
    real(r8), allocatable, dimension(:,:,:) :: dpv_lat_n
    real(r8), allocatable, dimension(:,:,:) :: m_lon
    real(r8), allocatable, dimension(:,:,:) :: m_lat
    real(r8), allocatable, dimension(:,:,:) :: mf_lon_n
    real(r8), allocatable, dimension(:,:,:) :: mf_lat_n
    real(r8), allocatable, dimension(:,:,:) :: mf_lon_t
    real(r8), allocatable, dimension(:,:,:) :: mf_lat_t
  contains
    final :: reduced_state_final
  end type reduced_state_type

  type(reduced_state_type), allocatable :: reduced_full_state(:)

  type reduced_tend_type
    real(r8), allocatable, dimension(:) :: qhv
    real(r8), allocatable, dimension(:) :: qhu
    real(r8), allocatable, dimension(:) :: mf_div_lon
    real(r8), allocatable, dimension(:) :: dpedlon
  contains
    final :: reduced_tend_final
  end type reduced_tend_type

  type(reduced_tend_type), allocatable :: reduced_full_tend(:)

  interface
    subroutine reduce_sub_interface(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)
      import mesh_type, state_type, reduced_mesh_type, reduced_state_type, r8
      integer, intent(in) :: j
      integer, intent(in) :: buf_j
      integer, intent(in) :: move
      type(mesh_type), intent(in) :: raw_mesh
      type(state_type), intent(in) :: raw_state
      type(reduced_mesh_type), intent(in) :: reduced_mesh
      type(reduced_state_type), intent(inout) :: reduced_state
      real(r8), intent(in) :: dt
    end subroutine reduce_sub_interface
  end interface

contains

  subroutine reduce_init()

    integer j, full_j

    allocate(reduced_full_mesh  (mesh%full_lat_lb:mesh%full_lat_ub))
    allocate(reduced_full_static(mesh%full_lat_start_idx:mesh%full_lat_end_idx))
    allocate(reduced_full_state (mesh%full_lat_start_idx:mesh%full_lat_end_idx))
    allocate(reduced_full_tend  (mesh%full_lat_start_idx:mesh%full_lat_end_idx))

    do j = 1, size(reduce_factors)
      if (reduce_factors(j) == 0) exit
      if (mod(mesh%num_full_lon, reduce_factors(j)) /= 0) then
        call log_error('Zonal reduce factor ' // to_string(reduce_factors(j)) // ' cannot divide zonal grid number ' // to_string(mesh%num_full_lon) // '!')
      end if
      if (mesh%has_south_pole()) then
#ifdef STAGGER_V_ON_POLE
        full_j = mesh%full_lat_start_idx+j-1
#else
        full_j = mesh%full_lat_start_idx+j
#endif
        call reduce_full_mesh(reduce_factors(j), full_j, mesh, reduced_full_mesh(full_j))
      end if
      if (mesh%has_north_pole()) then
#ifdef STAGGER_V_ON_POLE
        full_j = mesh%full_lat_end_idx-j+1
#else
        full_j = mesh%full_lat_end_idx-j
#endif
        call reduce_full_mesh(reduce_factors(j), full_j, mesh, reduced_full_mesh(full_j))
      end if
    end do

    do j = mesh%full_lat_start_idx, mesh%full_lat_end_idx
      if (reduced_full_mesh(j)%reduce_factor > 0) then
        call allocate_reduced_full_static(reduced_full_mesh(j), reduced_full_static(j))
        call reduce_full_static(j, mesh, static, reduced_full_mesh(j), reduced_full_static(j))
        call allocate_reduced_full_state(reduced_full_mesh(j), reduced_full_state(j))
        call allocate_reduced_full_tend(reduced_full_mesh(j), reduced_full_tend(j))
      end if
    end do

  end subroutine reduce_init

  subroutine reduce_run(state, dt)

    type(state_type), intent(in) :: state
    real(r8), intent(in) :: dt

    integer j

    do j = state%mesh%full_lat_start_idx, state%mesh%full_lat_end_idx
      if (reduced_full_mesh(j)%reduce_factor > 0) then
        call reduce_full_state(j, state%mesh, state, reduced_full_mesh(j), reduced_full_state(j), dt)
      end if
    end do

  end subroutine reduce_run

  subroutine reduce_full_mesh(reduce_factor, j, raw_mesh, reduced_mesh)

    integer, intent(in) :: reduce_factor
    integer, intent(in) :: j
    type(mesh_type), intent(in) :: raw_mesh
    type(reduced_mesh_type), intent(inout) :: reduced_mesh

    real(r8) x(3), y(3), z(3)
    integer i, buf_j

    reduced_mesh%reduce_factor      = reduce_factor
    reduced_mesh%halo_width         = raw_mesh%halo_width
    reduced_mesh%num_full_lon       = raw_mesh%num_full_lon / reduce_factor
    reduced_mesh%num_half_lon       = raw_mesh%num_half_lon / reduce_factor
    reduced_mesh%full_lon_start_idx = raw_mesh%full_lon_start_idx
    reduced_mesh%full_lon_end_idx   = raw_mesh%full_lon_start_idx + reduced_mesh%num_full_lon - 1
    reduced_mesh%half_lon_start_idx = raw_mesh%half_lon_start_idx
    reduced_mesh%half_lon_end_idx   = raw_mesh%half_lon_start_idx + reduced_mesh%num_half_lon - 1
    reduced_mesh%full_lon_lb        = reduced_mesh%full_lon_start_idx - raw_mesh%halo_width
    reduced_mesh%full_lon_ub        = reduced_mesh%full_lon_end_idx + raw_mesh%halo_width
    reduced_mesh%half_lon_lb        = reduced_mesh%half_lon_start_idx - raw_mesh%halo_width
    reduced_mesh%half_lon_ub        = reduced_mesh%half_lon_end_idx + raw_mesh%halo_width
    reduced_mesh%dlon               = raw_mesh%dlon * reduce_factor

    allocate(reduced_mesh%full_lon(reduced_mesh%num_full_lon))
    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
      reduced_mesh%full_lon(i) = raw_mesh%full_lon(raw_mesh%full_lon_start_idx) + (i - 1) * reduced_mesh%dlon
    end do
    allocate(reduced_mesh%half_lon(reduced_mesh%num_half_lon))
    do i = reduced_mesh%half_lon_start_idx, reduced_mesh%half_lon_end_idx
      reduced_mesh%half_lon(i) = reduced_mesh%full_lon(reduced_mesh%full_lon_start_idx) + (i - 0.5) * reduced_mesh%dlon
    end do

    allocate(reduced_mesh%full_lat           (  -2:2)); reduced_mesh%full_lat            = inf
    allocate(reduced_mesh%half_lat           (  -2:2)); reduced_mesh%half_lat            = inf
    allocate(reduced_mesh%full_cos_lat       (  -2:2)); reduced_mesh%full_cos_lat        = inf
    allocate(reduced_mesh%half_cos_lat       (  -2:2)); reduced_mesh%half_cos_lat        = inf
    allocate(reduced_mesh%full_sin_lat       (  -2:2)); reduced_mesh%full_sin_lat        = inf
    allocate(reduced_mesh%half_sin_lat       (  -2:2)); reduced_mesh%half_sin_lat        = inf
    allocate(reduced_mesh%cell_area          (  -1:1)); reduced_mesh%cell_area           = inf
    allocate(reduced_mesh%lon_edge_area      (  -1:1)); reduced_mesh%lon_edge_area       = inf
    allocate(reduced_mesh%lon_edge_left_area (  -1:1)); reduced_mesh%lon_edge_left_area  = inf
    allocate(reduced_mesh%lon_edge_right_area(  -1:1)); reduced_mesh%lon_edge_right_area = inf
    allocate(reduced_mesh%lat_edge_area      (  -1:1)); reduced_mesh%lat_edge_area       = inf
    allocate(reduced_mesh%lat_edge_up_area   (  -1:1)); reduced_mesh%lat_edge_up_area    = inf
    allocate(reduced_mesh%lat_edge_down_area (  -1:1)); reduced_mesh%lat_edge_down_area  = inf
    allocate(reduced_mesh%vertex_area        (  -1:1)); reduced_mesh%vertex_area         = inf
    allocate(reduced_mesh%subcell_area       (2,-1:1)); reduced_mesh%subcell_area        = inf
    allocate(reduced_mesh%de_lon             (  -1:1)); reduced_mesh%de_lon              = inf
    allocate(reduced_mesh%de_lat             (  -1:1)); reduced_mesh%de_lat              = inf
    allocate(reduced_mesh%le_lat             (  -2:2)); reduced_mesh%le_lat              = inf
    allocate(reduced_mesh%le_lon             (  -2:2)); reduced_mesh%le_lon              = inf
    allocate(reduced_mesh%full_tangent_wgt   (2,-1:1)); reduced_mesh%full_tangent_wgt    = inf
    allocate(reduced_mesh%half_tangent_wgt   (2,-1:1)); reduced_mesh%full_tangent_wgt    = inf
    allocate(reduced_mesh%full_f             (  -1:1)); reduced_mesh%full_f              = inf
    allocate(reduced_mesh%half_f             (  -1:1)); reduced_mesh%half_f              = inf
    reduced_mesh%full_lat     = raw_mesh%full_lat    (j+lbound(reduced_mesh%full_lat    , 1):j+ubound(reduced_mesh%full_lat    , 1))
    reduced_mesh%half_lat     = raw_mesh%half_lat    (j+lbound(reduced_mesh%half_lat    , 1):j+ubound(reduced_mesh%half_lat    , 1))
    reduced_mesh%full_cos_lat = raw_mesh%full_cos_lat(j+lbound(reduced_mesh%full_cos_lat, 1):j+ubound(reduced_mesh%full_cos_lat, 1))
    reduced_mesh%half_cos_lat = raw_mesh%half_cos_lat(j+lbound(reduced_mesh%half_cos_lat, 1):j+ubound(reduced_mesh%half_cos_lat, 1))
    reduced_mesh%full_sin_lat = raw_mesh%full_sin_lat(j+lbound(reduced_mesh%full_sin_lat, 1):j+ubound(reduced_mesh%full_sin_lat, 1))
    reduced_mesh%half_sin_lat = raw_mesh%half_sin_lat(j+lbound(reduced_mesh%half_sin_lat, 1):j+ubound(reduced_mesh%half_sin_lat, 1))
    reduced_mesh%full_f       = raw_mesh%full_f      (j+lbound(reduced_mesh%full_f      , 1):j+ubound(reduced_mesh%full_f      , 1))
    reduced_mesh%half_f       = raw_mesh%half_f      (j+lbound(reduced_mesh%half_f      , 1):j+ubound(reduced_mesh%half_f      , 1))

#ifdef STAGGER_V_ON_POLE
    ! Cell area
    do buf_j = lbound(reduced_mesh%cell_area, 1), ubound(reduced_mesh%cell_area, 1)
      if (j + buf_j >= max(raw_mesh%num_full_lat, raw_mesh%num_half_lat) .or. j + buf_j < 1) cycle
      if (buf_j > 0) then
        reduced_mesh%subcell_area(1,buf_j) = radius**2 * 0.5_r8 * reduced_mesh%dlon * (reduced_mesh%full_sin_lat(buf_j) - reduced_mesh%half_sin_lat(buf_j))
      end if
      if (buf_j < 1) then
        reduced_mesh%subcell_area(2,buf_j) = radius**2 * 0.5_r8 * reduced_mesh%dlon * (reduced_mesh%half_sin_lat(buf_j+1) - reduced_mesh%full_sin_lat(buf_j))
      end if
      reduced_mesh%cell_area(buf_j) = radius**2 * reduced_mesh%dlon * (reduced_mesh%half_sin_lat(buf_j+1) - reduced_mesh%half_sin_lat(buf_j))
      call cartesian_transform(reduced_mesh%full_lon(1), reduced_mesh%full_lat(buf_j  ), x(1), y(1), z(1))
      call cartesian_transform(reduced_mesh%half_lon(1), reduced_mesh%half_lat(buf_j  ), x(2), y(2), z(2))
      call cartesian_transform(reduced_mesh%half_lon(1), reduced_mesh%half_lat(buf_j+1), x(3), y(3), z(3))
      reduced_mesh%lon_edge_left_area(buf_j) = calc_area(x, y, z)
      reduced_mesh%lon_edge_right_area(buf_j) = reduced_mesh%lon_edge_left_area(buf_j)
      reduced_mesh%lon_edge_area(buf_j) = reduced_mesh%lon_edge_left_area(buf_j) + reduced_mesh%lon_edge_right_area(buf_j)
    end do
    ! Vertex area
    do buf_j = lbound(reduced_mesh%vertex_area, 1), ubound(reduced_mesh%vertex_area, 1)
      if (j + buf_j < 1) cycle
      if (reduced_mesh%full_sin_lat(buf_j-1) == inf) then
        reduced_mesh%vertex_area(buf_j) = radius**2 * reduced_mesh%dlon * (reduced_mesh%full_sin_lat(buf_j) + 1)
        cycle
      else if (reduced_mesh%full_sin_lat(buf_j) == inf) then
        reduced_mesh%vertex_area(buf_j) = radius**2 * reduced_mesh%dlon * (1 - reduced_mesh%full_sin_lat(buf_j-1))
        cycle
      else
        reduced_mesh%vertex_area(buf_j) = radius**2 * reduced_mesh%dlon * (reduced_mesh%full_sin_lat(buf_j) - reduced_mesh%full_sin_lat(buf_j-1))
      end if
      call cartesian_transform(reduced_mesh%full_lon(2), reduced_mesh%full_lat(buf_j  ), x(1), y(1), z(1))
      call cartesian_transform(reduced_mesh%half_lon(1), reduced_mesh%half_lat(buf_j  ), x(2), y(2), z(2))
      call cartesian_transform(reduced_mesh%half_lon(2), reduced_mesh%half_lat(buf_j  ), x(3), y(3), z(3))
      reduced_mesh%lat_edge_up_area(buf_j) = calc_area_with_last_small_arc(x, y, z)
      call cartesian_transform(reduced_mesh%full_lon(2), reduced_mesh%full_lat(buf_j-1), x(1), y(1), z(1))
      call cartesian_transform(reduced_mesh%half_lon(2), reduced_mesh%half_lat(buf_j  ), x(2), y(2), z(2))
      call cartesian_transform(reduced_mesh%half_lon(1), reduced_mesh%half_lat(buf_j  ), x(3), y(3), z(3))
      reduced_mesh%lat_edge_down_area(buf_j) = calc_area_with_last_small_arc(x, y, z)
      reduced_mesh%lat_edge_area(buf_j) = reduced_mesh%lat_edge_up_area(buf_j) + reduced_mesh%lat_edge_down_area(buf_j)
    end do
#else
#endif
    ! Edge lengths and cell distances
    do buf_j = lbound(reduced_mesh%le_lat, 1), ubound(reduced_mesh%le_lat, 1)
      if (raw_mesh%le_lat(j+buf_j) /= inf) then
        reduced_mesh%le_lat(buf_j) = raw_mesh%le_lat(j+buf_j) * reduce_factor
      end if
    end do
    do buf_j = lbound(reduced_mesh%de_lat, 1), ubound(reduced_mesh%de_lat, 1)
      if (reduced_mesh%le_lat(buf_j) /= inf) then
        reduced_mesh%de_lat(buf_j) = 2.0_r8 * reduced_mesh%lat_edge_area(buf_j) / reduced_mesh%le_lat(buf_j)
      end if
    end do
    do buf_j = lbound(reduced_mesh%le_lon, 1), ubound(reduced_mesh%le_lon, 1)
      if (raw_mesh%le_lon(j+buf_j) /= inf) then
        reduced_mesh%le_lon(buf_j) = raw_mesh%le_lon(j+buf_j)
      end if
    end do
    do buf_j = lbound(reduced_mesh%de_lon, 1), ubound(reduced_mesh%de_lon, 1)
      if (reduced_mesh%le_lon(buf_j) /= inf) then
        reduced_mesh%de_lon(buf_j) = 2.0_r8 * reduced_mesh%lon_edge_area(buf_j) / reduced_mesh%le_lon(buf_j)
      end if
    end do

#ifdef STAGGER_V_ON_POLE
    do buf_j = lbound(reduced_mesh%full_tangent_wgt, 2), ubound(reduced_mesh%full_tangent_wgt, 2)
      if (reduced_mesh%le_lat(buf_j) /= inf .and. reduced_mesh%de_lon(buf_j) /= inf) then
        reduced_mesh%full_tangent_wgt(1,buf_j) = reduced_mesh%le_lat(buf_j  ) / reduced_mesh%de_lon(buf_j) * 0.25_r8
      end if
      if (reduced_mesh%le_lat(buf_j+1) /= inf .and. reduced_mesh%de_lon(buf_j) /= inf) then
        reduced_mesh%full_tangent_wgt(2,buf_j) = reduced_mesh%le_lat(buf_j+1) / reduced_mesh%de_lon(buf_j) * 0.25_r8
      end if
    end do
#else
    do buf_j = lbound(reduced_mesh%full_tangent_wgt, 2), ubound(reduced_mesh%full_tangent_wgt, 2)
      if (reduced_mesh%le_lat(buf_j-1) /= inf .and. reduced_mesh%de_lon(buf_j) /= inf) then
        reduced_mesh%full_tangent_wgt(1,buf_j) = reduced_mesh%le_lat(buf_j-1) / reduced_mesh%de_lon(buf_j) * 0.25_r8
      end if
      if (reduced_mesh%le_lat(buf_j) /= inf .and. reduced_mesh%de_lon(buf_j) /= inf) then
        reduced_mesh%full_tangent_wgt(2,buf_j) = reduced_mesh%le_lat(buf_j  ) / reduced_mesh%de_lon(buf_j) * 0.25_r8
      end if
    end do
#endif
#ifdef STAGGER_V_ON_POLE
    do buf_j = lbound(reduced_mesh%half_tangent_wgt, 2), ubound(reduced_mesh%half_tangent_wgt, 2)
      if (reduced_mesh%le_lon(buf_j-1) /= inf .and. reduced_mesh%de_lat(buf_j) /= inf) then
        reduced_mesh%half_tangent_wgt(1,buf_j) = reduced_mesh%le_lon(buf_j-1) / reduced_mesh%de_lat(buf_j) * 0.25_r8
      end if
      if (reduced_mesh%le_lon(buf_j  ) /= inf .and. reduced_mesh%de_lat(buf_j) /= inf) then
        reduced_mesh%half_tangent_wgt(2,buf_j) = reduced_mesh%le_lon(buf_j  ) / reduced_mesh%de_lat(buf_j) * 0.25_r8
      end if
    end do
#else
    do buf_j = lbound(reduced_mesh%half_tangent_wgt, 2), ubound(reduced_mesh%half_tangent_wgt, 2)
      if (reduced_mesh%le_lon(buf_j  ) /= inf .and. reduced_mesh%de_lat(buf_j) /= inf) then
        reduced_mesh%half_tangent_wgt(1,buf_j) = reduced_mesh%le_lon(buf_j  ) / reduced_mesh%de_lat(buf_j) * 0.25_r8
      end if
      if (reduced_mesh%le_lon(buf_j+1) /= inf .and. reduced_mesh%de_lat(buf_j) /= inf) then
        reduced_mesh%half_tangent_wgt(2,buf_j) = reduced_mesh%le_lon(buf_j+1) / reduced_mesh%de_lat(buf_j) * 0.25_r8
      end if
    end do
#endif

  end subroutine reduce_full_mesh

  subroutine allocate_reduced_full_state(reduced_mesh, reduced_state)

    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state

#ifdef STAGGER_V_ON_POLE
    allocate(reduced_state%u        (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%v        (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub, 0:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%gd       (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%pv       (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-1:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%pv_lon   (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub, 0:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%pv_lat   (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub, 0:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lon_t(reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lon_n(reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub, 0:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lat_t(reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub, 0:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lat_n(reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub, 0:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%m_lon    (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%m_lat    (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub, 0:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%mf_lon_n (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%mf_lon_t (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub, 0:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%mf_lat_n (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub, 0:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%mf_lat_t (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub, 0:1,reduced_mesh%reduce_factor))
#else
#endif

  end subroutine allocate_reduced_full_state

  subroutine reduce_full_state(j, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    reduced_state%pv       (:,:,:) = inf
    reduced_state%pv_lon   (:,:,:) = inf
    reduced_state%pv_lat   (:,:,:) = inf
    reduced_state%dpv_lon_t(:,:,:) = inf
    reduced_state%dpv_lat_t(:,:,:) = inf
    reduced_state%dpv_lon_n(:,:,:) = inf
    reduced_state%dpv_lat_n(:,:,:) = inf
    reduced_state%m_lon    (:,:,:) = inf
    reduced_state%m_lat    (:,:,:) = inf
    reduced_state%mf_lon_n (:,:,:) = inf
    reduced_state%mf_lat_n (:,:,:) = inf
    reduced_state%mf_lon_t (:,:,:) = inf
    reduced_state%mf_lat_t (:,:,:) = inf

    call apply_reduce(lbound(reduced_state%u        , 2), ubound(reduced_state%u        , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_u          , dt)
    call apply_reduce(lbound(reduced_state%v        , 2), ubound(reduced_state%v        , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_v          , dt)
    call apply_reduce(lbound(reduced_state%gd       , 2), ubound(reduced_state%gd       , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_gd         , dt)
    call apply_reduce(lbound(reduced_state%pv       , 2), ubound(reduced_state%pv       , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_pv         , dt)
    call apply_reduce(lbound(reduced_state%m_lon    , 2), ubound(reduced_state%m_lon    , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_m_lon      , dt)
    call apply_reduce(lbound(reduced_state%mf_lon_n , 2), ubound(reduced_state%mf_lon_n , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_mf_lon_n   , dt)
    call apply_reduce(lbound(reduced_state%mf_lat_t , 2), ubound(reduced_state%mf_lat_t , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_mf_lat_t   , dt)
    call apply_reduce(lbound(reduced_state%m_lat    , 2), ubound(reduced_state%m_lat    , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_m_lat      , dt)
    call apply_reduce(lbound(reduced_state%mf_lat_n , 2), ubound(reduced_state%mf_lat_n , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_mf_lat_n   , dt)
    call apply_reduce(lbound(reduced_state%mf_lon_t , 2), ubound(reduced_state%mf_lon_t , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_mf_lon_t   , dt)
    call apply_reduce(lbound(reduced_state%dpv_lon_t, 2), ubound(reduced_state%dpv_lon_t, 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_dpv_lon_t  , dt)
    call apply_reduce(lbound(reduced_state%dpv_lat_n, 2), ubound(reduced_state%dpv_lat_n, 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_dpv_lat_n  , dt)
    call apply_reduce(lbound(reduced_state%dpv_lat_t, 2), ubound(reduced_state%dpv_lat_t, 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_dpv_lat_t  , dt)
    call apply_reduce(lbound(reduced_state%dpv_lon_n, 2), ubound(reduced_state%dpv_lon_n, 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_dpv_lon_n  , dt)
    call apply_reduce(lbound(reduced_state%pv_lon   , 2), ubound(reduced_state%pv_lon   , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_pv_lon_apvm, dt)
    call apply_reduce(lbound(reduced_state%pv_lat   , 2), ubound(reduced_state%pv_lat   , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_pv_lat_apvm, dt)

  end subroutine reduce_full_state

  subroutine reduce_ghs(j, buf_j, move, raw_mesh, raw_static, reduced_mesh, reduced_static)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(static_type), intent(in) :: raw_static
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(inout) :: reduced_static

    integer raw_i, i

    raw_i = move
    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
      reduced_static%ghs(i,buf_j,move) = sum(raw_static%ghs(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j))
      raw_i = raw_i + reduced_mesh%reduce_factor
    end do
    reduced_static%ghs(:,buf_j, move) = reduced_static%ghs(:,buf_j, move) / reduced_mesh%reduce_factor
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_static%ghs(:,buf_j, move))

  end subroutine reduce_ghs

  subroutine apply_reduce(buf_lb, buf_ub, j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_sub, dt)

    integer, intent(in) :: buf_lb
    integer, intent(in) :: buf_ub
    integer, intent(in) :: j
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    procedure(reduce_sub_interface) reduce_sub
    real(r8), intent(in) :: dt

    integer buf_j, move

    do move = 1, reduced_mesh%reduce_factor
      do buf_j = buf_lb, buf_ub
        call reduce_sub(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)
      end do
    end do

  end subroutine apply_reduce

  subroutine reduce_u(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer raw_i, i

    raw_i = move
    do i = reduced_mesh%half_lon_start_idx, reduced_mesh%half_lon_end_idx
      reduced_state%u(i,buf_j,move) = sum(raw_state%u(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j))
      raw_i = raw_i + reduced_mesh%reduce_factor
    end do
    reduced_state%u(:,buf_j,move) = reduced_state%u(:,buf_j,move) / reduced_mesh%reduce_factor
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%u(:,buf_j,move))

  end subroutine reduce_u

  subroutine reduce_v(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer raw_i, i

    raw_i = move
    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
      reduced_state%v(i,buf_j,move) = sum(raw_state%v(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j))
      raw_i = raw_i + reduced_mesh%reduce_factor
    end do
    reduced_state%v(:,buf_j,move) = reduced_state%v(:,buf_j,move) / reduced_mesh%reduce_factor
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%v(:,buf_j,move))

  end subroutine reduce_v

  subroutine reduce_gd(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer raw_i, i

    raw_i = move
    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
      reduced_state%gd(i,buf_j,move) = sum(raw_state%gd(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j))
      raw_i = raw_i + reduced_mesh%reduce_factor
    end do
    reduced_state%gd(:,buf_j,move) = reduced_state%gd(:,buf_j,move) / reduced_mesh%reduce_factor
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%gd(:,buf_j,move))

  end subroutine reduce_gd

  subroutine reduce_pv(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer raw_i, i

    raw_i = move
    do i = reduced_mesh%half_lon_start_idx, reduced_mesh%half_lon_end_idx
      reduced_state%pv(i,buf_j,move) = sum(raw_state%pv(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j))
      raw_i = raw_i + reduced_mesh%reduce_factor
    end do
    reduced_state%pv(:,buf_j,move) = reduced_state%pv(:,buf_j,move) / reduced_mesh%reduce_factor
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%pv(:,buf_j,move))

  end subroutine reduce_pv

  subroutine reduce_m_lon(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i

    do i = reduced_mesh%half_lon_start_idx, reduced_mesh%half_lon_end_idx
      reduced_state%m_lon(i,buf_j,move) = (                                          &
        reduced_mesh%lon_edge_left_area (buf_j) * reduced_state%gd(i  ,buf_j,move) + &
        reduced_mesh%lon_edge_right_area(buf_j) * reduced_state%gd(i+1,buf_j,move)   &
      ) / reduced_mesh%lon_edge_area(buf_j) / g
    end do

  end subroutine reduce_m_lon

  subroutine reduce_m_lat(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i

    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
#ifdef STAGGER_V_ON_POLE
      reduced_state%m_lat(i,buf_j,move) = (                                         &
        reduced_mesh%lat_edge_up_area  (buf_j) * reduced_state%gd(i,buf_j  ,move) + &
        reduced_mesh%lat_edge_down_area(buf_j) * reduced_state%gd(i,buf_j-1,move)   &
      ) / reduced_mesh%lat_edge_area(buf_j) / g
#else
#endif
    end do

  end subroutine reduce_m_lat

  subroutine reduce_mf_lon_n(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i

    do i = reduced_mesh%half_lon_start_idx, reduced_mesh%half_lon_end_idx
      reduced_state%mf_lon_n(i,buf_j,move) = reduced_state%m_lon(i,buf_j,move) * reduced_state%u(i,buf_j,move)
    end do
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%mf_lon_n(:,buf_j,move))

  end subroutine reduce_mf_lon_n

  subroutine reduce_mf_lat_n(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i

    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
      reduced_state%mf_lat_n(i,buf_j,move) = reduced_state%m_lat(i,buf_j,move) * reduced_state%v(i,buf_j,move)
    end do
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%mf_lat_n(:,buf_j,move))

  end subroutine reduce_mf_lat_n

  subroutine reduce_mf_lon_t(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i

    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
#ifdef STAGGER_V_ON_POLE
      reduced_state%mf_lon_t(i,buf_j,move) =                                                                                           &
        reduced_mesh%full_tangent_wgt(1,buf_j) * (reduced_state%mf_lat_n(i,buf_j  ,move) + reduced_state%mf_lat_n(i+1,buf_j  ,move)) + &
        reduced_mesh%full_tangent_wgt(2,buf_j) * (reduced_state%mf_lat_n(i,buf_j+1,move) + reduced_state%mf_lat_n(i+1,buf_j+1,move))
#else
#endif
    end do
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%mf_lon_t(:,buf_j,move))

  end subroutine reduce_mf_lon_t

  subroutine reduce_mf_lat_t(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i

    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
#ifdef STAGGER_V_ON_POLE
      reduced_state%mf_lat_t(i,buf_j,move) =                                                                                           &
        reduced_mesh%half_tangent_wgt(1,buf_j) * (reduced_state%mf_lon_n(i-1,buf_j-1,move) + reduced_state%mf_lon_n(i,buf_j-1,move)) + &
        reduced_mesh%half_tangent_wgt(2,buf_j) * (reduced_state%mf_lon_n(i-1,buf_j  ,move) + reduced_state%mf_lon_n(i,buf_j  ,move))
#else
#endif
    end do
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%mf_lat_t(:,buf_j,move))

  end subroutine reduce_mf_lat_t

  subroutine reduce_dpv_lon_t(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i

    do i = reduced_mesh%half_lon_start_idx, reduced_mesh%half_lon_end_idx
#ifdef STAGGER_V_ON_POLE
      reduced_state%dpv_lon_t(i,buf_j,move) = reduced_state%pv(i,buf_j+1,move) - reduced_state%pv(i,buf_j  ,move)
#else
      reduced_state%dpv_lon_t(i,buf_j,move) = reduced_state%pv(i,buf_j  ,move) - reduced_state%pv(i,buf_j-1,move)
#endif
    end do
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%dpv_lon_t(:,buf_j,move))

  end subroutine reduce_dpv_lon_t

  subroutine reduce_dpv_lat_t(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i

    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
      reduced_state%dpv_lat_t(i,buf_j,move) = reduced_state%pv(i+1,buf_j,move) - reduced_state%pv(i,buf_j,move)
    end do
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%dpv_lat_t(:,buf_j,move))

  end subroutine reduce_dpv_lat_t

  subroutine reduce_dpv_lon_n(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i

    do i = reduced_mesh%half_lon_start_idx, reduced_mesh%half_lon_end_idx
#ifdef STAGGER_V_ON_POLE
      reduced_state%dpv_lon_n(i,buf_j,move) = 0.25_r8 * ( &
        reduced_state%dpv_lat_t(i  ,buf_j  ,move) +       &
        reduced_state%dpv_lat_t(i+1,buf_j  ,move) +       &
        reduced_state%dpv_lat_t(i  ,buf_j+1,move) +       &
        reduced_state%dpv_lat_t(i+1,buf_j+1,move)         &
      )
#else
      reduced_state%dpv_lon_n(i,buf_j,move) = 0.25_r8 * ( &
        reduced_state%dpv_lat_t(i  ,buf_j-1,move) +       &
        reduced_state%dpv_lat_t(i+1,buf_j-1,move) +       &
        reduced_state%dpv_lat_t(i  ,buf_j  ,move) +       &
        reduced_state%dpv_lat_t(i+1,buf_j  ,move)         &
      )
#endif
    end do

  end subroutine reduce_dpv_lon_n

  subroutine reduce_dpv_lat_n(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i

    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
#ifdef STAGGER_V_ON_POLE
      reduced_state%dpv_lat_n(i,buf_j,move) = 0.25_r8 * ( &
        reduced_state%dpv_lon_t(i-1,buf_j-1,move) +       &
        reduced_state%dpv_lon_t(i  ,buf_j-1,move) +       &
        reduced_state%dpv_lon_t(i-1,buf_j  ,move) +       &
        reduced_state%dpv_lon_t(i  ,buf_j  ,move)         &
      )
#else
      reduced_state%dpv_lat_n(i,buf_j,move) = 0.25_r8 * ( &
        reduced_state%dpv_lon_t(i-1,buf_j  ,move) +       &
        reduced_state%dpv_lon_t(i  ,buf_j  ,move) +       &
        reduced_state%dpv_lon_t(i-1,buf_j+1,move) +       &
        reduced_state%dpv_lon_t(i  ,buf_j+1,move)         &
      )
#endif
    end do

  end subroutine reduce_dpv_lat_n

  subroutine reduce_pv_lon_apvm(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    real(r8) un, vn, ut, vt, le, de
    integer i

    le = reduced_mesh%le_lon(buf_j)
    de = reduced_mesh%de_lon(buf_j)
    if (le == inf .or. de == inf) return
    do i = reduced_mesh%half_lon_start_idx, reduced_mesh%half_lon_end_idx
      un = reduced_state%u(i,buf_j,move)
      vt = reduced_state%mf_lon_t(i,buf_j,move) / reduced_state%m_lon(i,buf_j,move)
#ifdef STAGGER_V_ON_POLE
      reduced_state%pv_lon(i,buf_j,move) = 0.5_r8 * (     &
        reduced_state%pv(i,buf_j+1,move) +                &
        reduced_state%pv(i,buf_j  ,move)                  &
      ) - 0.5_r8 * (                                      &
        un * reduced_state%dpv_lon_n(i,buf_j,move) / de + &
        vt * reduced_state%dpv_lon_t(i,buf_j,move) / le   &
      ) * dt
#else
      reduced_state%pv_lon(i,buf_j,move) = 0.5_r8 * (     &
        reduced_state%pv(i,buf_j-1,move) +                &
        reduced_state%pv(i,buf_j  ,move)                  &
      ) - 0.5_r8 * (                                      &
        un * reduced_state%dpv_lon_n(i,buf_j,move) / de + &
        vt * reduced_state%dpv_lon_t(i,buf_j,move) / le   &
      ) * dt
#endif
    end do

  end subroutine reduce_pv_lon_apvm

  subroutine reduce_pv_lat_apvm(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    real(r8) un, vn, ut, vt, le, de
    integer i

    le = reduced_mesh%le_lat(buf_j)
    de = reduced_mesh%de_lat(buf_j)
    if (le == inf .or. de == inf) return
    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
      ut = reduced_state%mf_lat_t(i,buf_j,move) / reduced_state%m_lat(i,buf_j,move)
      vn = reduced_state%v(i,buf_j,move)
      reduced_state%pv_lat(i,buf_j,move) = 0.5_r8 * (     &
        reduced_state%pv(i  ,buf_j,move) +                &
        reduced_state%pv(i-1,buf_j,move)                  &
      ) - 0.5_r8 * (                                      &
        ut * reduced_state%dpv_lat_t(i,buf_j,move) / le + &
        vn * reduced_state%dpv_lat_n(i,buf_j,move) / de   &
      ) * dt
    end do

  end subroutine reduce_pv_lat_apvm

  subroutine reduce_append_array(move, reduced_mesh, reduced_array, raw_mesh, raw_array)

    integer, intent(in) :: move
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    real(r8), intent(in) :: reduced_array(reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub)
    type(mesh_type), intent(in) :: raw_mesh
    real(r8), intent(inout) :: raw_array(raw_mesh%full_lon_lb:raw_mesh%full_lon_ub)

    integer i, raw_i

    raw_i = raw_mesh%full_lon_start_idx + move - 1
    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
      raw_array(raw_i:raw_i+reduced_mesh%reduce_factor-1) = raw_array(raw_i:raw_i+reduced_mesh%reduce_factor-1) + reduced_array(i) / reduced_mesh%reduce_factor
      raw_i = raw_i + reduced_mesh%reduce_factor
    end do

  end subroutine reduce_append_array

  subroutine allocate_reduced_full_static(reduced_mesh, reduced_static)

    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(inout) :: reduced_static

    allocate(reduced_static%ghs(reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub, 0:0,reduced_mesh%reduce_factor))

  end subroutine allocate_reduced_full_static

  subroutine allocate_reduced_full_tend(reduced_mesh, reduced_tend)

    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_tend_type), intent(inout) :: reduced_tend

    allocate(reduced_tend%qhu       (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub))
    allocate(reduced_tend%qhv       (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub))
    allocate(reduced_tend%mf_div_lon(reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub))
    allocate(reduced_tend%dpedlon   (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub))

  end subroutine allocate_reduced_full_tend

  subroutine reduce_full_static(j, raw_mesh, raw_static, reduced_mesh, reduced_static)

    integer, intent(in) :: j
    type(mesh_type), intent(in) :: raw_mesh
    type(static_type), intent(in) :: raw_static
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(inout) :: reduced_static

    integer buf_j, move

    do move = 1, reduced_mesh%reduce_factor
      do buf_j = lbound(reduced_static%ghs, 2), ubound(reduced_static%ghs, 2)
        call reduce_ghs(j, buf_j, move, raw_mesh, raw_static, reduced_mesh, reduced_static)
      end do
    end do

  end subroutine reduce_full_static

  subroutine reduced_mesh_final(this)

    type(reduced_mesh_type), intent(inout) :: this

    if (allocated(this%full_lon)) deallocate(this%full_lon)
    if (allocated(this%half_lon)) deallocate(this%half_lon)

  end subroutine reduced_mesh_final

  subroutine reduced_static_final(this)

    type(reduced_static_type), intent(inout) :: this

    if (allocated(this%ghs)) deallocate(this%ghs)

  end subroutine reduced_static_final

  subroutine reduced_state_final(this)

    type(reduced_state_type), intent(inout) :: this

    if (allocated(this%u        )) deallocate(this%u        )
    if (allocated(this%v        )) deallocate(this%v        )
    if (allocated(this%gd       )) deallocate(this%gd       )
    if (allocated(this%pv       )) deallocate(this%pv       )
    if (allocated(this%pv_lon   )) deallocate(this%pv_lon   )
    if (allocated(this%pv_lat   )) deallocate(this%pv_lat   )
    if (allocated(this%dpv_lon_n)) deallocate(this%dpv_lon_n)
    if (allocated(this%dpv_lat_n)) deallocate(this%dpv_lat_n)
    if (allocated(this%dpv_lon_t)) deallocate(this%dpv_lon_t)
    if (allocated(this%dpv_lat_t)) deallocate(this%dpv_lat_t)
    if (allocated(this%m_lon    )) deallocate(this%m_lon    )
    if (allocated(this%m_lat    )) deallocate(this%m_lat    )
    if (allocated(this%mf_lon_n )) deallocate(this%mf_lon_n )
    if (allocated(this%mf_lat_n )) deallocate(this%mf_lat_n )
    if (allocated(this%mf_lon_t )) deallocate(this%mf_lon_t )
    if (allocated(this%mf_lat_t )) deallocate(this%mf_lat_t )

  end subroutine reduced_state_final

  subroutine reduced_tend_final(this)

    type(reduced_tend_type), intent(inout) :: this

    if (allocated(this%qhv       )) deallocate(this%qhv       )
    if (allocated(this%qhu       )) deallocate(this%qhu       )
    if (allocated(this%mf_div_lon)) deallocate(this%mf_div_lon)
    if (allocated(this%dpedlon   )) deallocate(this%dpedlon   )

  end subroutine reduced_tend_final

  subroutine reduce_final()

    if (allocated(reduced_full_mesh  )) deallocate(reduced_full_mesh  )
    if (allocated(reduced_full_static)) deallocate(reduced_full_static)
    if (allocated(reduced_full_state )) deallocate(reduced_full_state )
    if (allocated(reduced_full_tend  )) deallocate(reduced_full_tend  )

  end subroutine reduce_final

end module reduce_mod
