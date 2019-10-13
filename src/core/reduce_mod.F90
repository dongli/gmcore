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
    real(r8), allocatable :: full_lon(:)
    real(r8), allocatable :: half_lon(:)
    real(r8), dimension(  -2:2) :: full_lat
    real(r8), dimension(  -2:2) :: half_lat
    real(r8), dimension(  -1:1) :: full_cos_lat
    real(r8), dimension(  -1:1) :: half_cos_lat
    real(r8), dimension(  -2:2) :: full_sin_lat
    real(r8), dimension(  -2:2) :: half_sin_lat
    real(r8), dimension(  -1:1) :: cell_area
    real(r8), dimension(  -1:1) :: lon_edge_area
    real(r8), dimension(  -1:1) :: lon_edge_left_area
    real(r8), dimension(  -1:1) :: lon_edge_right_area
    real(r8), dimension(  -1:1) :: lat_edge_area
    real(r8), dimension(  -1:1) :: lat_edge_up_area
    real(r8), dimension(  -1:1) :: lat_edge_down_area
    real(r8), dimension(  -1:1) :: vertex_area
    real(r8), dimension(2,-1:1) :: subcell_area
    real(r8), dimension(  -1:1) :: de_lon
    real(r8), dimension(  -1:1) :: de_lat
    real(r8), dimension(  -1:1) :: le_lat
    real(r8), dimension(  -1:1) :: le_lon
    real(r8), dimension(2,-1:1) :: full_tangent_wgt
    real(r8), dimension(2,-1:1) :: half_tangent_wgt
    real(r8), dimension(  -1:1) :: full_f
    real(r8), dimension(  -1:1) :: half_f
  contains
    final :: reduced_mesh_final
  end type reduced_mesh_type

  type(reduced_mesh_type), allocatable :: reduced_full_mesh(:)
  type(reduced_mesh_type), allocatable :: reduced_half_mesh(:)

  type reduced_static_type
    real(r8), allocatable, dimension(:,:,:) :: ghs
  contains
    final :: reduced_static_final
  end type reduced_static_type

  type(reduced_static_type), allocatable :: reduced_full_static(:)
  type(reduced_static_type), allocatable :: reduced_half_static(:)

  type reduced_state_type
    ! Prognostic variables
    real(r8), allocatable, dimension(:,:,:) :: u
    real(r8), allocatable, dimension(:,:,:) :: v
    real(r8), allocatable, dimension(:,:,:) :: gd
    ! Diagnostic variables (FIXME: Should we only keep diagnostic variables on center lat?)
    real(r8), allocatable, dimension(:,:,:) :: pv
    real(r8), allocatable, dimension(:,:,:) :: pv_lon
    real(r8), allocatable, dimension(:,:,:) :: pv_lat
    real(r8), allocatable, dimension(:,:,:) :: dpv_lon_t
    real(r8), allocatable, dimension(:,:,:) :: dpv_lat_t
    real(r8), allocatable, dimension(:,:,:) :: dpv_lon_n
    real(r8), allocatable, dimension(:,:,:) :: dpv_lat_n
    real(r8), allocatable, dimension(:,:,:) :: m_vtx
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
  type(reduced_state_type), allocatable :: reduced_half_state(:)

  type reduced_tend_type
    real(r8), allocatable, dimension(:) :: mf_div_lon
    real(r8), allocatable, dimension(:) :: dpedlon
    real(r8), allocatable, dimension(:) :: dkedlon
  contains
    final :: reduced_tend_final
  end type reduced_tend_type

  type(reduced_tend_type), allocatable :: reduced_full_tend(:)

contains

  subroutine reduce_init()

    integer j, full_j, half_j

    ! --------------------------------------------------------------------------
    !                                  MESH
    if (allocated(reduced_full_mesh)) deallocate(reduced_full_mesh)
    if (allocated(reduced_half_mesh)) deallocate(reduced_half_mesh)
    allocate(reduced_full_mesh(mesh%full_lat_start_idx:mesh%full_lat_end_idx))
    allocate(reduced_half_mesh(mesh%half_lat_start_idx:mesh%half_lat_end_idx))

    do j = 1, size(reduce_factors)
      if (reduce_factors(j) == 0) exit
      if (mod(mesh%num_full_lon, reduce_factors(j)) /= 0) then
        call log_error('Zonal reduce factor ' // to_string(reduce_factors(j)) // ' cannot divide zonal grid number ' // to_string(mesh%num_full_lon) // '!')
      end if
      if (mesh%has_south_pole()) then
#ifdef STAGGER_V_ON_POLE
        full_j = mesh%full_lat_start_idx+j-1
        half_j = mesh%half_lat_start_idx+j
#else
        full_j = mesh%full_lat_start_idx+j
        half_j = mesh%half_lat_start_idx+j-1
#endif
        call reduce_mesh(reduce_factors(j), full_j, mesh, reduced_full_mesh(full_j))
        call reduce_mesh(reduce_factors(j), half_j, mesh, reduced_half_mesh(half_j))
      end if
      if (mesh%has_north_pole()) then
#ifdef STAGGER_V_ON_POLE
        full_j = mesh%full_lat_end_idx-j+1
        half_j = mesh%half_lat_end_idx-j
#else
        full_j = mesh%full_lat_end_idx-j
        half_j = mesh%half_lat_end_idx-j+1
#endif
        call reduce_mesh(reduce_factors(j), full_j, mesh, reduced_full_mesh(full_j))
        call reduce_mesh(reduce_factors(j), half_j, mesh, reduced_half_mesh(half_j))
      end if
    end do

    ! --------------------------------------------------------------------------
    !                               STATIC
    if (allocated(reduced_full_static)) deallocate(reduced_full_static)
    if (allocated(reduced_half_static)) deallocate(reduced_half_static)
    allocate(reduced_full_static(mesh%full_lat_start_idx:mesh%full_lat_end_idx))
    allocate(reduced_half_static(mesh%half_lat_start_idx:mesh%half_lat_end_idx))

    do j = 1, size(reduce_factors)
      if (reduce_factors(j) == 0) exit
      if (mesh%has_south_pole()) then
#ifdef STAGGER_V_ON_POLE
        full_j = mesh%full_lat_start_idx+j-1
        half_j = mesh%half_lat_start_idx+j
#else
        full_j = mesh%full_lat_start_idx+j
        half_j = mesh%half_lat_start_idx+j-1
#endif
        call allocate_reduced_static(reduced_full_mesh(full_j), reduced_full_static(full_j))
        call allocate_reduced_static(reduced_half_mesh(half_j), reduced_half_static(half_j))
      end if
      if (mesh%has_north_pole()) then
#ifdef STAGGER_V_ON_POLE
        full_j = mesh%full_lat_end_idx-j+1
        half_j = mesh%half_lat_end_idx-j
#else
        full_j = mesh%full_lat_end_idx-j
        half_j = mesh%half_lat_end_idx-j+1
#endif
        call allocate_reduced_static(reduced_full_mesh(full_j), reduced_full_static(full_j))
        call allocate_reduced_static(reduced_half_mesh(half_j), reduced_half_static(half_j))
      end if
    end do
    do j = mesh%full_lat_start_idx, mesh%full_lat_end_idx
      call reduce_static(j, mesh, static, reduced_full_mesh(j), reduced_full_static(j))
    end do
    do j = mesh%half_lat_start_idx, mesh%half_lat_end_idx
      call reduce_static(j, mesh, static, reduced_half_mesh(j), reduced_half_static(j))
    end do

    ! --------------------------------------------------------------------------
    !                               STATE
    if (allocated(reduced_full_state)) deallocate(reduced_full_state)
    if (allocated(reduced_half_state)) deallocate(reduced_half_state)
    allocate(reduced_full_state(mesh%full_lat_start_idx:mesh%full_lat_end_idx))
    allocate(reduced_half_state(mesh%half_lat_start_idx:mesh%half_lat_end_idx))

    do j = 1, size(reduce_factors)
      if (reduce_factors(j) == 0) exit
      if (mesh%has_south_pole()) then
#ifdef STAGGER_V_ON_POLE
        full_j = mesh%full_lat_start_idx+j-1
        half_j = mesh%half_lat_start_idx+j
#else
        full_j = mesh%full_lat_start_idx+j
        half_j = mesh%half_lat_start_idx+j-1
#endif
        call allocate_reduced_state(reduced_full_mesh(full_j), reduced_full_state(full_j))
        call allocate_reduced_state(reduced_half_mesh(half_j), reduced_half_state(half_j))
      end if
      if (mesh%has_north_pole()) then
#ifdef STAGGER_V_ON_POLE
        full_j = mesh%full_lat_end_idx-j+1
        half_j = mesh%half_lat_end_idx-j
#else
        full_j = mesh%full_lat_end_idx-j
        half_j = mesh%half_lat_end_idx-j+1
#endif
        call allocate_reduced_state(reduced_full_mesh(full_j), reduced_full_state(full_j))
        call allocate_reduced_state(reduced_half_mesh(half_j), reduced_half_state(half_j))
      end if
    end do

    ! --------------------------------------------------------------------------
    !                               TENDENCY
    if (allocated(reduced_full_tend)) deallocate(reduced_full_tend)
    allocate(reduced_full_tend(mesh%full_lat_start_idx:mesh%full_lat_end_idx))

    do j = 1, size(reduce_factors)
      if (reduce_factors(j) == 0) exit
      if (mesh%has_south_pole()) then
#ifdef STAGGER_V_ON_POLE
        full_j = mesh%full_lat_start_idx+j-1
        half_j = mesh%half_lat_start_idx+j
#else
        full_j = mesh%full_lat_start_idx+j
        half_j = mesh%half_lat_start_idx+j-1
#endif
        call allocate_reduced_tend(reduced_full_mesh(full_j), reduced_full_tend(full_j))
      end if
      if (mesh%has_north_pole()) then
#ifdef STAGGER_V_ON_POLE
        full_j = mesh%full_lat_end_idx-j+1
        half_j = mesh%half_lat_end_idx-j
#else
        full_j = mesh%full_lat_end_idx-j
        half_j = mesh%half_lat_end_idx-j+1
#endif
        call allocate_reduced_tend(reduced_full_mesh(full_j), reduced_full_tend(full_j))
      end if
    end do

  end subroutine reduce_init

  subroutine reduce_run(state, dt)

    type(state_type), intent(in) :: state
    real(r8), intent(in) :: dt

    integer j


    do j = state%mesh%full_lat_start_idx, state%mesh%full_lat_end_idx
      call reduce_state(j, state%mesh, state, reduced_full_mesh(j), reduced_full_state(j), dt)
    end do
    do j = state%mesh%half_lat_start_idx, state%mesh%half_lat_end_idx
      call reduce_state(j, state%mesh, state, reduced_half_mesh(j), reduced_half_state(j), dt)
    end do

  end subroutine reduce_run

  subroutine reduce_mesh(reduce_factor, j, raw_mesh, reduced_mesh)

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
    allocate(reduced_mesh%half_lon(reduced_mesh%num_half_lon))
    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
      reduced_mesh%full_lon(i) = raw_mesh%full_lon(raw_mesh%full_lon_start_idx) + (i - 1) * reduced_mesh%dlon
    end do
    do i = reduced_mesh%half_lon_start_idx, reduced_mesh%half_lon_end_idx
      reduced_mesh%half_lon(i) = reduced_mesh%full_lon(reduced_mesh%full_lon_start_idx) + (i - 0.5) * reduced_mesh%dlon
    end do
    reduced_mesh%full_lat         = raw_mesh%full_lat        (  j-2:j+2)
    reduced_mesh%half_lat         = raw_mesh%half_lat        (  j-2:j+2)
    reduced_mesh%full_sin_lat     = raw_mesh%full_sin_lat    (  j-2:j+2)
    reduced_mesh%half_sin_lat     = raw_mesh%half_sin_lat    (  j-2:j+2)
    reduced_mesh%full_cos_lat     = raw_mesh%full_cos_lat    (  j-1:j+1)
    reduced_mesh%half_cos_lat     = raw_mesh%half_cos_lat    (  j-1:j+1)
    reduced_mesh%de_lon    = raw_mesh%de_lon   (  j-1:j+1) * reduce_factor
    reduced_mesh%de_lat    = raw_mesh%de_lat   (  j-1:j+1)
    reduced_mesh%le_lat  = raw_mesh%le_lat (  j-1:j+1) * reduce_factor
    reduced_mesh%le_lon  = raw_mesh%le_lon (  j-1:j+1)
    reduced_mesh%full_tangent_wgt = raw_mesh%full_tangent_wgt(:,j-1:j+1)
    reduced_mesh%half_tangent_wgt = raw_mesh%half_tangent_wgt(:,j-1:j+1)
    reduced_mesh%full_f           = raw_mesh%full_f          (  j-1:j+1)
    reduced_mesh%half_f           = raw_mesh%half_f          (  j-1:j+1)

#ifdef STAGGER_V_ON_POLE
    ! Cell area
    do buf_j = -1, 1
      if (j + buf_j >= max(raw_mesh%num_full_lat, raw_mesh%num_half_lat) .or. j + buf_j < 1) cycle
      reduced_mesh%cell_area(buf_j) = radius**2 * reduced_mesh%dlon * (reduced_mesh%half_sin_lat(buf_j+1) - reduced_mesh%half_sin_lat(buf_j))
      reduced_mesh%subcell_area(1,buf_j) = radius**2 * 0.5d0 * reduced_mesh%dlon * (reduced_mesh%full_sin_lat(buf_j) - reduced_mesh%half_sin_lat(buf_j))
      reduced_mesh%subcell_area(2,buf_j) = radius**2 * 0.5d0 * reduced_mesh%dlon * (reduced_mesh%half_sin_lat(buf_j+1) - reduced_mesh%full_sin_lat(buf_j))
      call cartesian_transform(reduced_mesh%full_lon(1), reduced_mesh%full_lat(buf_j  ), x(1), y(1), z(1))
      call cartesian_transform(reduced_mesh%half_lon(1), reduced_mesh%half_lat(buf_j  ), x(2), y(2), z(2))
      call cartesian_transform(reduced_mesh%half_lon(1), reduced_mesh%half_lat(buf_j+1), x(3), y(3), z(3))
      reduced_mesh%lon_edge_left_area(buf_j) = calc_area(x, y, z)
      reduced_mesh%lon_edge_right_area(buf_j) = reduced_mesh%lon_edge_left_area(buf_j)
      reduced_mesh%lon_edge_area(buf_j) = reduced_mesh%lon_edge_left_area(buf_j) + reduced_mesh%lon_edge_right_area(buf_j)
      ! Check reduced areas.
      if (abs(reduced_mesh%cell_area(buf_j) / raw_mesh%cell_area(j+buf_j) - reduced_mesh%reduce_factor) > 1.0d-12) then
        call log_error('Failed to reduce cell_area!', __FILE__, __LINE__)
      end if
      if (abs(reduced_mesh%subcell_area(1,buf_j) / raw_mesh%subcell_area(1,j+buf_j) - reduced_mesh%reduce_factor) > 1.0d-12) then
        call log_error('Failed to reduce subcell_area!', __FILE__, __LINE__)
      end if
      if (abs(reduced_mesh%subcell_area(2,buf_j) / raw_mesh%subcell_area(2,j+buf_j) - reduced_mesh%reduce_factor) > 1.0d-12) then
        call log_error('Failed to reduce subcell_area!', __FILE__, __LINE__)
      end if
      if (abs(reduced_mesh%lon_edge_area(buf_j) / raw_mesh%lon_edge_area(j+buf_j) - reduced_mesh%reduce_factor) > 1.0d-2) then
        call log_error('Failed to reduce lon_edge_area!', __FILE__, __LINE__)
      end if
    end do
    ! Vertex area
    do buf_j = -1, 1
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
      ! Check areas
      if (abs(reduced_mesh%vertex_area(buf_j) / raw_mesh%vertex_area(j+buf_j) - reduced_mesh%reduce_factor) > 1.0d-12) then
        call log_error('Failed to reduce vertex_area!', __FILE__, __LINE__)
      end if
      if (abs(reduced_mesh%lat_edge_area(buf_j) / raw_mesh%lat_edge_area(j+buf_j) - reduced_mesh%reduce_factor) > 1.0d-2) then
        call log_error('Failed to reduce lat_edge_area!', __FILE__, __LINE__)
      end if
    end do
#else
#endif

  end subroutine reduce_mesh

  subroutine allocate_reduced_static(reduced_mesh, reduced_static)

    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(inout) :: reduced_static

    allocate(reduced_static%ghs(reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:1,reduced_mesh%reduce_factor))

  end subroutine allocate_reduced_static

  subroutine allocate_reduced_state(reduced_mesh, reduced_state)

    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state

    allocate(reduced_state%u        (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-2:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%v        (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-2:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%gd       (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-2:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%pv       (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%pv_lon   (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%pv_lat   (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lon_t(reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lat_t(reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lon_n(reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lat_n(reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%m_vtx    (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%m_lon    (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%m_lat    (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%mf_lon_n (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%mf_lat_n (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%mf_lon_t (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%mf_lat_t (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-1:1,reduced_mesh%reduce_factor))

  end subroutine allocate_reduced_state

  subroutine allocate_reduced_tend(reduced_mesh, reduced_tend)

    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_tend_type), intent(inout) :: reduced_tend

    allocate(reduced_tend%mf_div_lon(reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub))
    allocate(reduced_tend%dpedlon   (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub))
    allocate(reduced_tend%dkedlon   (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub))

  end subroutine allocate_reduced_tend

  subroutine reduce_static(j, raw_mesh, raw_static, reduced_mesh, reduced_static)

    integer, intent(in) :: j
    type(mesh_type), intent(in) :: raw_mesh
    type(static_type), intent(in) :: raw_static
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(inout) :: reduced_static

    integer buf_j, move, raw_i, i

    if (reduced_mesh%reduce_factor == 0) return
    do buf_j = -1, 1
      do move = 1, reduced_mesh%reduce_factor
        call reduce_array(j, move, raw_mesh, raw_static%ghs(:,j+buf_j), reduced_mesh, reduced_static%ghs(:,buf_j,move))
      end do
    end do

  end subroutine reduce_static

  subroutine reduce_state(j, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer buf_j, move, raw_i, i

    if (reduced_mesh%reduce_factor == 0) return
    ! Reduce prognostic variables.
    do buf_j = -2, 2
      do move = 1, reduced_mesh%reduce_factor
        call reduce_array(j, move, raw_mesh, raw_state%u (:,j+buf_j), reduced_mesh, reduced_state%u (:,buf_j,move))
        call reduce_array(j, move, raw_mesh, raw_state%v (:,j+buf_j), reduced_mesh, reduced_state%v (:,buf_j,move))
        call reduce_array(j, move, raw_mesh, raw_state%gd(:,j+buf_j), reduced_mesh, reduced_state%gd(:,buf_j,move))
      end do
    end do
    ! Reduce diagnostic variables.
    do buf_j = -1, 1
      do move = 1, reduced_mesh%reduce_factor
        call diagnose_1(buf_j, move, reduced_mesh, reduced_state)
      end do
    end do

    do move = 1, reduced_mesh%reduce_factor
      call diagnose_2(0, move, reduced_mesh, reduced_state)
    end do

  end subroutine reduce_state

  subroutine reduce_array(j, move, raw_mesh, raw_array, reduced_mesh, reduced_array)

    integer, intent(in) :: j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    real(r8), intent(in) :: raw_array(raw_mesh%full_lon_lb:raw_mesh%full_lon_ub)
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    real(r8), intent(out) :: reduced_array(reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub)

    integer raw_i, i

    raw_i = move
    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
      reduced_array(i) = sum(raw_array(raw_i:raw_i+reduced_mesh%reduce_factor-1)) / reduced_mesh%reduce_factor
      raw_i = raw_i + reduced_mesh%reduce_factor
    end do

    call parallel_fill_halo(reduced_mesh%halo_width, reduced_array)

  end subroutine reduce_array

  subroutine diagnose_1(buf_j, move, reduced_mesh, reduced_state)

    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state

    integer i

    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
      ! Mass on vertex
#ifdef STAGGER_V_ON_POLE
      reduced_state%m_vtx(i,buf_j,move) = ((reduced_state%gd(i,buf_j-1,move) + reduced_state%gd(i+1,buf_j-1,move)) * reduced_mesh%subcell_area(2,buf_j-1) + &
                                           (reduced_state%gd(i,buf_j  ,move) + reduced_state%gd(i+1,buf_j  ,move)) * reduced_mesh%subcell_area(1,buf_j  )   &
                                          ) / reduced_mesh%vertex_area(buf_j) / g
#else
#endif
      ! Mass on zonal edge
      reduced_state%m_lon(i,buf_j,move) = (reduced_mesh%lon_edge_left_area (buf_j) * reduced_state%gd(i  ,buf_j,move) + &
                                           reduced_mesh%lon_edge_right_area(buf_j) * reduced_state%gd(i+1,buf_j,move)   &
                                          ) / reduced_mesh%lon_edge_area(buf_j) / g
      ! Mass on meridional edge
#ifdef STAGGER_V_ON_POLE
      reduced_state%m_lat(i,buf_j,move) = (reduced_mesh%lat_edge_up_area  (buf_j) * reduced_state%gd(i,buf_j  ,move) + &
                                           reduced_mesh%lat_edge_down_area(buf_j) * reduced_state%gd(i,buf_j-1,move)   &
                                          ) / reduced_mesh%lat_edge_area(buf_j) / g
#else
#endif
      ! Normal mass flux on zonal edge
      reduced_state%mf_lon_n(i,buf_j,move) = reduced_state%m_lon(i,buf_j,move) * reduced_state%u(i,buf_j,move)
      ! Normal mass flux on meridional edge
      reduced_state%mf_lat_n(i,buf_j,move) = reduced_state%m_lat(i,buf_j,move) * reduced_state%v(i,buf_j,move)
      ! PV on vertex
#ifdef STAGGER_V_ON_POLE
      reduced_state%pv(i,buf_j,move) = ((reduced_state%u(i  ,buf_j-1,move) * reduced_mesh%de_lon(buf_j-1) - &
                                         reduced_state%u(i  ,buf_j  ,move) * reduced_mesh%de_lon(buf_j  ) + &
                                         reduced_state%v(i+1,buf_j  ,move) * reduced_mesh%de_lat(buf_j  ) - &
                                         reduced_state%v(i  ,buf_j  ,move) * reduced_mesh%de_lat(buf_j  )   &
                                        ) / reduced_mesh%vertex_area(buf_j) + reduced_mesh%half_f(buf_j)        &
                                       ) / reduced_state%m_vtx(i,buf_j,move)
#else
#endif
    end do
    ! PV on edge
    ! TODO: Here we need to implement APVM.

    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%mf_lon_n(:,buf_j,move))
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%mf_lat_n(:,buf_j,move))

  end subroutine diagnose_1

  subroutine diagnose_2(buf_j, move, reduced_mesh, reduced_state)

    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state

    integer i

    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
      ! Tangent mass flux on zonal edge
#ifdef STAGGER_V_ON_POLE
      reduced_state%mf_lat_t(i,buf_j,move) = reduced_mesh%full_tangent_wgt(1,buf_j) * (reduced_state%mf_lat_n(i,buf_j  ,move) + reduced_state%mf_lat_n(i+1,buf_j  ,move)) + &
                                             reduced_mesh%full_tangent_wgt(2,buf_j) * (reduced_state%mf_lat_n(i,buf_j+1,move) + reduced_state%mf_lat_n(i+1,buf_j+1,move))
#else
#endif
      ! Tangent mass flux on meridional edge
#ifdef STAGGER_V_ON_POLE
      reduced_state%mf_lon_t(i,buf_j,move) = reduced_mesh%half_tangent_wgt(1,buf_j) * (reduced_state%mf_lon_n(i-1,buf_j-1,move) + reduced_state%mf_lon_n(i,buf_j-1,move)) + &
                                             reduced_mesh%half_tangent_wgt(2,buf_j) * (reduced_state%mf_lon_n(i-1,buf_j  ,move) + reduced_state%mf_lon_n(i,buf_j  ,move))
#else
#endif
    end do

    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%mf_lon_t(:,buf_j,move))
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%mf_lat_t(:,buf_j,move))

  end subroutine diagnose_2

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

    if (allocated(this%u       )) deallocate(this%u       )
    if (allocated(this%v       )) deallocate(this%v       )
    if (allocated(this%gd      )) deallocate(this%gd      )
    if (allocated(this%pv      )) deallocate(this%pv      )
    if (allocated(this%pv_lon  )) deallocate(this%pv_lon  )
    if (allocated(this%pv_lat  )) deallocate(this%pv_lat  )
    if (allocated(this%m_vtx   )) deallocate(this%m_vtx   )
    if (allocated(this%m_lon   )) deallocate(this%m_lon   )
    if (allocated(this%m_lat   )) deallocate(this%m_lat   )
    if (allocated(this%mf_lon_n)) deallocate(this%mf_lon_n)
    if (allocated(this%mf_lat_n)) deallocate(this%mf_lat_n)
    if (allocated(this%mf_lon_t)) deallocate(this%mf_lon_t)
    if (allocated(this%mf_lat_t)) deallocate(this%mf_lat_t)

  end subroutine reduced_state_final

  subroutine reduced_tend_final(this)

    type(reduced_tend_type), intent(inout) :: this

    if (allocated(this%mf_div_lon)) deallocate(this%mf_div_lon)
    if (allocated(this%dpedlon   )) deallocate(this%dpedlon   )
    if (allocated(this%dkedlon   )) deallocate(this%dkedlon   )

  end subroutine reduced_tend_final

end module reduce_mod
