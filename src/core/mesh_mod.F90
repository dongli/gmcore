module mesh_mod

  use flogger
  use const_mod
  use namelist_mod
  use sphere_geometry_mod

  implicit none

  private

  public mesh_type
  public mesh
  public create_meshes

  type mesh_type
    integer id
    integer halo_width
    integer num_full_lon
    integer num_half_lon
    integer num_full_lat
    integer num_half_lat
    integer full_lon_start_idx
    integer full_lon_end_idx
    integer full_lat_start_idx
    integer full_lat_end_idx
    integer full_lat_start_idx_no_pole
    integer full_lat_end_idx_no_pole
    integer half_lon_start_idx
    integer half_lon_end_idx
    integer half_lat_start_idx
    integer half_lat_end_idx
    integer half_lat_start_idx_no_pole
    integer half_lat_end_idx_no_pole
    integer full_lon_lb
    integer full_lon_ub
    integer half_lon_lb
    integer half_lon_ub
    integer full_lat_lb
    integer full_lat_ub
    integer half_lat_lb
    integer half_lat_ub
    real(r8) start_lon
    real(r8) end_lon
    real(r8) start_lat
    real(r8) end_lat
    real(r8) dlon
    real(r8) dlat
    real(r8) total_area
    real(r8), allocatable :: full_lon(:)
    real(r8), allocatable :: half_lon(:)
    real(r8), allocatable :: full_lat(:)
    real(r8), allocatable :: half_lat(:)
    real(r8), allocatable :: full_cos_lon(:)
    real(r8), allocatable :: half_cos_lon(:)
    real(r8), allocatable :: full_sin_lon(:)
    real(r8), allocatable :: half_sin_lon(:)
    real(r8), allocatable :: full_cos_lat(:)
    real(r8), allocatable :: half_cos_lat(:)
    real(r8), allocatable :: full_sin_lat(:)
    real(r8), allocatable :: half_sin_lat(:)
    ! For output
    real(r8), allocatable :: full_lon_deg(:)
    real(r8), allocatable :: half_lon_deg(:)
    real(r8), allocatable :: full_lat_deg(:)
    real(r8), allocatable :: half_lat_deg(:)
    ! Area for weighting
    real(r8), allocatable :: cell_area(:)
    real(r8), allocatable :: lon_edge_area(:)
    real(r8), allocatable :: lon_edge_left_area(:)
    real(r8), allocatable :: lon_edge_right_area(:)
    real(r8), allocatable :: lat_edge_area(:)
    real(r8), allocatable :: lat_edge_up_area(:)
    real(r8), allocatable :: lat_edge_down_area(:)
    real(r8), allocatable :: vertex_area(:)
    real(r8), allocatable :: subcell_area(:,:)
    ! Edge length
    real(r8), allocatable :: cell_lon_dist(:)
    real(r8), allocatable :: cell_lat_dist(:)
    real(r8), allocatable :: vertex_lon_dist(:)
    real(r8), allocatable :: vertex_lat_dist(:)
    ! Weight for constructing tangential wind
    real(r8), allocatable :: full_tangent_wgt(:,:)
    real(r8), allocatable :: half_tangent_wgt(:,:)
    ! Weight for dissipating potential enstrophy
    real(r8), allocatable :: full_upwind_beta(:)
    real(r8), allocatable :: half_upwind_beta(:)
    ! Coriolis parameters
    real(r8), allocatable :: full_f(:)
    real(r8), allocatable :: half_f(:)
  contains
    procedure :: init => mesh_init
    procedure :: has_south_pole => mesh_has_south_pole
    procedure :: has_north_pole => mesh_has_north_pole
    procedure :: is_south_pole => mesh_is_south_pole
    procedure :: is_north_pole => mesh_is_north_pole
    procedure :: is_pole => mesh_is_pole
    final :: mesh_final
  end type mesh_type

  type(mesh_type) mesh

contains

  subroutine mesh_init(this, num_lon, num_lat, id, halo_width, start_lon, end_lon, start_lat, end_lat)

    class(mesh_type), intent(inout)           :: this
    integer         , intent(in   )           :: num_lon
    integer         , intent(in   )           :: num_lat
    integer         , intent(in   ), optional :: id
    integer         , intent(in   ), optional :: halo_width
    real(r8)        , intent(in   ), optional :: start_lon
    real(r8)        , intent(in   ), optional :: end_lon
    real(r8)        , intent(in   ), optional :: start_lat
    real(r8)        , intent(in   ), optional :: end_lat

    real(r8) x(3), y(3), z(3), total_area
    integer i, j

    this%num_full_lon = num_lon
    this%num_half_lon = num_lon
#ifdef STAGGER_V_ON_POLE
    this%num_full_lat = num_lat - 1
    this%num_half_lat = num_lat
#else
    this%num_full_lat = num_lat
    this%num_half_lat = num_lat - 1
#endif

    this%full_lon_start_idx = 1
    this%full_lon_end_idx = this%num_full_lon
    this%full_lat_start_idx = 1
    this%full_lat_end_idx = this%num_full_lat
    this%half_lon_start_idx = 1
    this%half_lon_end_idx = this%num_half_lon
    this%half_lat_start_idx = 1
    this%half_lat_end_idx = this%num_half_lat

    this%id         = merge(id        ,  0         , present(id))
    this%halo_width = merge(halo_width,  1         , present(halo_width))
    this%start_lon  = merge(start_lon ,  0.0_r8     , present(start_lon))
    this%end_lon    = merge(end_lon   ,  2.0_r8 * pi, present(end_lon))
    this%start_lat  = merge(start_lat , -0.5_r8 * pi, present(start_lat))
    this%end_lat    = merge(end_lat   ,  0.5_r8 * pi, present(end_lat))
    this%total_area = radius**2 * (this%end_lon - this%start_lon) * (sin(this%end_lat) - sin(this%start_lat))

#ifdef STAGGER_V_ON_POLE
    this%full_lat_start_idx_no_pole = this%full_lat_start_idx
    this%full_lat_end_idx_no_pole   = this%full_lat_end_idx
    this%half_lat_start_idx_no_pole = merge(this%half_lat_start_idx + 1, this%half_lat_start_idx, this%has_south_pole())
    this%half_lat_end_idx_no_pole   = merge(this%half_lat_end_idx   - 1, this%half_lat_end_idx  , this%has_north_pole())
#else
    this%full_lat_start_idx_no_pole = merge(this%full_lat_start_idx + 1, this%full_lat_start_idx, this%has_south_pole())
    this%full_lat_end_idx_no_pole   = merge(this%full_lat_end_idx   - 1, this%full_lat_end_idx  , this%has_north_pole())
    this%half_lat_start_idx_no_pole = this%half_lat_start_idx
    this%half_lat_end_idx_no_pole   = this%half_lat_end_idx
#endif

    this%full_lon_lb = this%full_lon_start_idx - this%halo_width
    this%full_lon_ub = this%full_lon_end_idx   + this%halo_width
    this%full_lat_lb = this%full_lat_start_idx - 2
    this%full_lat_ub = this%full_lat_end_idx   + 2
    this%half_lon_lb = this%half_lon_start_idx - this%halo_width
    this%half_lon_ub = this%half_lon_end_idx   + this%halo_width
    this%half_lat_lb = this%half_lat_start_idx - 2
    this%half_lat_ub = this%half_lat_end_idx   + 2

    allocate(this%full_lon           (this%full_lon_lb:this%full_lon_ub)); this%full_lon            = inf
    allocate(this%half_lon           (this%half_lon_lb:this%half_lon_ub)); this%half_lon            = inf
    allocate(this%full_lat           (this%full_lat_lb:this%full_lat_ub)); this%full_lat            = inf
    allocate(this%half_lat           (this%half_lat_lb:this%half_lat_ub)); this%half_lat            = inf
    allocate(this%full_cos_lon       (this%full_lon_lb:this%full_lon_ub)); this%full_cos_lon        = inf
    allocate(this%half_cos_lon       (this%half_lon_lb:this%half_lon_ub)); this%half_cos_lon        = inf
    allocate(this%full_sin_lon       (this%full_lon_lb:this%full_lon_ub)); this%full_sin_lon        = inf
    allocate(this%half_sin_lon       (this%half_lon_lb:this%half_lon_ub)); this%half_sin_lon        = inf
    allocate(this%full_cos_lat       (this%full_lat_lb:this%full_lat_ub)); this%full_cos_lat        = inf
    allocate(this%half_cos_lat       (this%half_lat_lb:this%half_lat_ub)); this%half_cos_lat        = inf
    allocate(this%full_sin_lat       (this%full_lat_lb:this%full_lat_ub)); this%full_sin_lat        = inf
    allocate(this%half_sin_lat       (this%half_lat_lb:this%half_lat_ub)); this%half_sin_lat        = inf
    allocate(this%full_lon_deg       (this%full_lon_lb:this%full_lon_ub)); this%full_lon_deg        = inf
    allocate(this%half_lon_deg       (this%half_lon_lb:this%half_lon_ub)); this%half_lon_deg        = inf
    allocate(this%full_lat_deg       (this%full_lat_lb:this%full_lat_ub)); this%full_lat_deg        = inf
    allocate(this%half_lat_deg       (this%half_lat_lb:this%half_lat_ub)); this%half_lat_deg        = inf
    allocate(this%cell_area          (this%full_lat_lb:this%full_lat_ub)); this%cell_area           = 0.0_r8
    allocate(this%lon_edge_area      (this%full_lat_lb:this%full_lat_ub)); this%lon_edge_area       = 0.0_r8
    allocate(this%lon_edge_left_area (this%full_lat_lb:this%full_lat_ub)); this%lon_edge_left_area  = 0.0_r8
    allocate(this%lon_edge_right_area(this%full_lat_lb:this%full_lat_ub)); this%lon_edge_right_area = 0.0_r8
    allocate(this%lat_edge_area      (this%half_lat_lb:this%half_lat_ub)); this%lat_edge_area       = 0.0_r8
    allocate(this%lat_edge_up_area   (this%half_lat_lb:this%half_lat_ub)); this%lat_edge_up_area    = 0.0_r8
    allocate(this%lat_edge_down_area (this%half_lat_lb:this%half_lat_ub)); this%lat_edge_down_area  = 0.0_r8
    allocate(this%vertex_area        (this%half_lat_lb:this%half_lat_ub)); this%vertex_area         = 0.0_r8
    allocate(this%subcell_area     (2,this%full_lat_lb:this%full_lat_ub)); this%subcell_area        = 0.0_r8
    allocate(this%cell_lon_dist      (this%full_lat_lb:this%full_lat_ub)); this%cell_lon_dist       = inf
    allocate(this%cell_lat_dist      (this%half_lat_lb:this%half_lat_ub)); this%cell_lat_dist       = inf
    allocate(this%vertex_lon_dist    (this%half_lat_lb:this%half_lat_ub)); this%vertex_lon_dist     = inf
    allocate(this%vertex_lat_dist    (this%full_lat_lb:this%full_lat_ub)); this%vertex_lat_dist     = inf
    allocate(this%full_tangent_wgt (2,this%full_lat_lb:this%full_lat_ub)); this%full_tangent_wgt    = inf
    allocate(this%half_tangent_wgt (2,this%half_lat_lb:this%half_lat_ub)); this%half_tangent_wgt    = inf
    allocate(this%full_upwind_beta   (this%full_lat_lb:this%full_lat_ub)); this%full_upwind_beta    = inf
    allocate(this%half_upwind_beta   (this%half_lat_lb:this%half_lat_ub)); this%half_upwind_beta    = inf
    allocate(this%full_f             (this%full_lat_lb:this%full_lat_ub)); this%full_f              = inf
    allocate(this%half_f             (this%half_lat_lb:this%half_lat_ub)); this%half_f              = inf

    this%dlon = (this%end_lon - this%start_lon) / this%num_full_lon
    do i = this%full_lon_lb, this%full_lon_ub
      this%full_lon(i) = this%start_lon + (i - 1) * this%dlon
      this%half_lon(i) = this%full_lon(i) + 0.5_r8 * this%dlon
      this%full_lon_deg(i) = this%full_lon(i) * deg
      this%half_lon_deg(i) = this%half_lon(i) * deg
    end do

#ifdef STAGGER_V_ON_POLE
    this%dlat = (this%end_lat - this%start_lat) / this%num_full_lat
    do j = this%half_lat_lb, this%half_lat_ub
      this%half_lat(j) = this%start_lat + (j - 1) * this%dlat
      if (abs(this%half_lat(j)) < 1.0e-14) this%half_lat(j) = 0.0_r8
      this%half_lat_deg(j) = this%half_lat(j) * deg
      if (this%half_lat(j) < -pi05 .or. this%half_lat(j) > pi05) then
        this%half_lat(j) = inf
        this%half_lat_deg(j) = inf
      end if
    end do
    this%half_lat(this%num_half_lat) = this%end_lat
    this%half_lat_deg(this%num_half_lat) = this%end_lat * deg

    do j = this%full_lat_lb, this%full_lat_ub
      this%full_lat(j) = this%half_lat(j) + 0.5_r8 * this%dlat
      if (abs(this%full_lat(j)) < 1.0e-14) this%full_lat(j) = 0.0_r8
      this%full_lat_deg(j) = this%full_lat(j) * deg
      if (this%full_lat(j) < -pi05 .or. this%full_lat(j) > pi05) then
        this%full_lat(j) = inf
        this%full_lat_deg(j) = inf
      end if
    end do
#else
    this%dlat = (this%end_lat - this%start_lat) / this%num_half_lat
    do j = this%full_lat_lb, this%full_lat_ub
      this%full_lat(j) = this%start_lat + (j - 1) * this%dlat
      if (abs(this%full_lat(j)) < 1.0e-14) this%full_lat(j) = 0.0_r8
      this%full_lat_deg(j) = this%full_lat(j) * deg
      if (this%full_lat(j) < -pi05 .or. this%full_lat(j) > pi05) then
        this%full_lat(j) = inf
        this%full_lat_deg(j) = inf
      end if
    end do
    this%full_lat(this%num_full_lat) = this%end_lat
    this%full_lat_deg(this%num_full_lat) = this%end_lat * deg

    do j = this%half_lat_lb, this%half_lat_ub
      this%half_lat(j) = this%full_lat(j) + 0.5_r8 * this%dlat
      if (abs(this%half_lat(j)) < 1.0e-14) this%half_lat(j) = 0.0_r8
      this%half_lat_deg(j) = this%half_lat(j) * deg
      if (this%half_lat(j) < -pi05 .or. this%half_lat(j) > pi05) then
        this%half_lat(j) = inf
        this%half_lat_deg(j) = inf
      end if
    end do
#endif

    do i = this%full_lon_lb, this%full_lon_ub
      this%full_cos_lon(i) = cos(this%full_lon(i))
      this%full_sin_lon(i) = sin(this%full_lon(i))
    end do

    do i = this%half_lon_lb, this%half_lon_ub
      this%half_cos_lon(i) = cos(this%half_lon(i))
      this%half_sin_lon(i) = sin(this%half_lon(i))
    end do

    do j = this%half_lat_lb, this%half_lat_ub
      if (this%half_lat(j) >= -pi05 .and. this%half_lat(j) <= pi05) then
        this%half_cos_lat(j) = cos(this%half_lat(j))
        this%half_sin_lat(j) = sin(this%half_lat(j))
      end if
    end do

    do j = this%full_lat_lb, this%full_lat_ub
      if (this%full_lat(j) >= -pi05 .and. this%full_lat(j) <= pi05) then
        this%full_cos_lat(j) = cos(this%full_lat(j))
        this%full_sin_lat(j) = sin(this%full_lat(j))
      end if
    end do

    ! Ensure the values of cos_lat and sin_lat are expected at the Poles.
#ifdef STAGGER_V_ON_POLE
    if (this%has_south_pole()) then
      this%half_cos_lat(this%half_lat_start_idx) =  0.0_r8
      this%half_sin_lat(this%half_lat_start_idx) = -1.0_r8
    end if
    if (this%has_north_pole()) then
      this%half_cos_lat(this%half_lat_end_idx) = 0.0_r8
      this%half_sin_lat(this%half_lat_end_idx) = 1.0_r8
    end if
#else
    if (this%has_south_pole()) then
      this%full_cos_lat(this%full_lat_start_idx) =  0.0_r8
      this%full_sin_lat(this%full_lat_start_idx) = -1.0_r8
    end if
    if (this%has_north_pole()) then
      this%full_cos_lat(this%full_lat_end_idx) = 0.0_r8
      this%full_sin_lat(this%full_lat_end_idx) = 1.0_r8
    end if
#endif

#ifdef STAGGER_V_ON_POLE
    do j = this%full_lat_start_idx, this%full_lat_end_idx
      this%cell_area(j) = radius**2 * this%dlon * (this%half_sin_lat(j+1) - this%half_sin_lat(j))
      this%subcell_area(1,j) = radius**2 * 0.5d0 * this%dlon * (this%full_sin_lat(j) - this%half_sin_lat(j))
      this%subcell_area(2,j) = radius**2 * 0.5d0 * this%dlon * (this%half_sin_lat(j+1) - this%full_sin_lat(j))
      call cartesian_transform(this%full_lon(1), this%full_lat(j  ), x(1), y(1), z(1))
      call cartesian_transform(this%half_lon(1), this%half_lat(j  ), x(2), y(2), z(2))
      call cartesian_transform(this%half_lon(1), this%half_lat(j+1), x(3), y(3), z(3))
      this%lon_edge_left_area(j) = calc_area(x, y, z)
      this%lon_edge_right_area(j) = this%lon_edge_left_area(j)
      this%lon_edge_area(j) = this%lon_edge_left_area(j) + this%lon_edge_right_area(j)
    end do

    do j = this%half_lat_start_idx, this%half_lat_end_idx
      if (this%is_south_pole(j)) then
        this%vertex_area(j) = radius**2 * this%dlon * (this%full_sin_lat(j) + 1)
      else if (this%is_north_pole(j)) then
        this%vertex_area(j) = radius**2 * this%dlon * (1 - this%full_sin_lat(j-1))
      else
        this%vertex_area(j) = radius**2 * this%dlon * (this%full_sin_lat(j) - this%full_sin_lat(j-1))
        call cartesian_transform(this%full_lon(2), this%full_lat(j  ), x(1), y(1), z(1))
        call cartesian_transform(this%half_lon(1), this%half_lat(j  ), x(2), y(2), z(2))
        call cartesian_transform(this%half_lon(2), this%half_lat(j  ), x(3), y(3), z(3))
        this%lat_edge_up_area(j) = calc_area_with_last_small_arc(x, y, z)
        call cartesian_transform(this%full_lon(2), this%full_lat(j-1), x(1), y(1), z(1))
        call cartesian_transform(this%half_lon(2), this%half_lat(j  ), x(2), y(2), z(2))
        call cartesian_transform(this%half_lon(1), this%half_lat(j  ), x(3), y(3), z(3))
        this%lat_edge_down_area(j) = calc_area_with_last_small_arc(x, y, z)
        this%lat_edge_area(j) = this%lat_edge_up_area(j) + this%lat_edge_down_area(j)
      end if
    end do
#else
    do j = this%full_lat_start_idx, this%full_lat_end_idx
      if (this%is_south_pole(j)) then
        this%cell_area(j) = radius**2 * this%dlon * (this%half_sin_lat(j) + 1.0d0)
        this%subcell_area(2,j) = radius**2 * 0.5d0 * this%dlon * (this%half_sin_lat(j) + 1.0d0)
      else if (this%is_north_pole(j)) then
        this%cell_area(j) = radius**2 * this%dlon * (1.0 - this%half_sin_lat(j-1))
        this%subcell_area(1,j) = radius**2 * 0.5d0 * this%dlon * (1.0d0 - this%half_sin_lat(j-1))
      else
        this%cell_area(j) = radius**2 * this%dlon * (this%half_sin_lat(j) - this%half_sin_lat(j-1))
        this%subcell_area(1,j) = radius**2 * 0.5d0 * this%dlon * (this%full_sin_lat(j) - this%half_sin_lat(j-1))
        this%subcell_area(2,j) = radius**2 * 0.5d0 * this%dlon * (this%half_sin_lat(j) - this%full_sin_lat(j))
        call cartesian_transform(this%full_lon(1), this%full_lat(j  ), x(1), y(1), z(1))
        call cartesian_transform(this%half_lon(1), this%half_lat(j-1), x(2), y(2), z(2))
        call cartesian_transform(this%half_lon(1), this%half_lat(j  ), x(3), y(3), z(3))
        this%lon_edge_left_area(j) = calc_area(x, y, z)
        this%lon_edge_right_area(j) = this%lon_edge_left_area(j)
        this%lon_edge_area(j) = this%lon_edge_left_area(j) + this%lon_edge_right_area(j)
      end if
    end do

    do j = this%half_lat_start_idx, this%half_lat_end_idx
      this%vertex_area(j) = radius**2 * this%dlon * (this%full_sin_lat(j+1) - this%full_sin_lat(j))
      call cartesian_transform(this%full_lon(2), this%full_lat(j+1), x(1), y(1), z(1))
      call cartesian_transform(this%half_lon(1), this%half_lat(j  ), x(2), y(2), z(2))
      call cartesian_transform(this%half_lon(2), this%half_lat(j  ), x(3), y(3), z(3))
      this%lat_edge_up_area(j) = calc_area_with_last_small_arc(x, y, z)
      call cartesian_transform(this%full_lon(2), this%full_lat(j), x(1), y(1), z(1))
      call cartesian_transform(this%half_lon(2), this%half_lat(j), x(2), y(2), z(2))
      call cartesian_transform(this%half_lon(1), this%half_lat(j), x(3), y(3), z(3))
      this%lat_edge_down_area(j) = calc_area_with_last_small_arc(x, y, z)
      this%lat_edge_area(j) = this%lat_edge_up_area(j) + this%lat_edge_down_area(j)
    end do
#endif

    total_area = 0.0d0
    do j = this%full_lat_start_idx, this%full_lat_end_idx
      total_area = total_area + this%cell_area(j) * this%num_full_lon
    end do
    if (abs((this%total_area - total_area) / this%total_area) > 1.0d-12) then
      call log_error('Failed to calculate cell area!', __FILE__, __LINE__)
    end if

    total_area = 0.0d0
    do j = this%half_lat_start_idx, this%half_lat_end_idx
      total_area = total_area + this%vertex_area(j) * this%num_full_lon
    end do
    if (abs((this%total_area - total_area) / this%total_area) > 1.0d-12) then
      call log_error('Failed to calculate vertex area!', __FILE__, __LINE__)
    end if

    total_area = 0.0d0
    do j = this%full_lat_start_idx, this%full_lat_end_idx
      total_area = total_area + sum(this%subcell_area(:,j)) * this%num_full_lon * 2
    end do
    if (abs((this%total_area - total_area) / this%total_area) > 1.0d-12) then
      call log_error('Failed to calculate subcell area!', __FILE__, __LINE__)
    end if

    do j = this%full_lat_start_idx, this%full_lat_end_idx
      if (abs((this%cell_area(j) - 2.0d0 * sum(this%subcell_area(:,j))) / this%cell_area(j)) > 1.0d-12) then
        call log_error('Failed to calculate subcell area!', __FILE__, __LINE__)
      end if
    end do

#ifdef STAGGER_V_ON_POLE
    do j = this%half_lat_start_idx, this%half_lat_end_idx
      if (this%is_south_pole(j)) then
        if (abs((this%vertex_area(j) - 2.0d0 * this%subcell_area(1,j)) / this%vertex_area(j)) > 1.0d-12) then
          call log_error('Failed to calculate subcell area!', __FILE__, __LINE__)
        end if
      else if (this%is_north_pole(j)) then
        if (abs((this%vertex_area(j) - 2.0d0 * this%subcell_area(2,j-1)) / this%vertex_area(j)) > 1.0d-12) then
          call log_error('Failed to calculate subcell area!', __FILE__, __LINE__)
        end if
      else
        if (abs((this%vertex_area(j) - 2.0d0 * (this%subcell_area(2,j-1) + this%subcell_area(1,j))) / this%vertex_area(j)) > 1.0d-12) then
          call log_error('Failed to calculate subcell area!', __FILE__, __LINE__)
        end if
      end if
    end do
#else
    do j = this%half_lat_start_idx, this%half_lat_end_idx
      if (abs((this%vertex_area(j) - 2.0d0 * (this%subcell_area(2,j) + this%subcell_area(1,j+1))) / this%vertex_area(j)) > 1.0d-12) then
        call log_error('Failed to calculate subcell area!', __FILE__, __LINE__)
      end if
    end do
#endif

    total_area = 0.0d0
    do j = this%full_lat_start_idx_no_pole, this%full_lat_end_idx_no_pole
      total_area = total_area + this%lon_edge_area(j) * this%num_full_lon
    end do
    do j = this%half_lat_start_idx_no_pole, this%half_lat_end_idx_no_pole
      total_area = total_area + this%lat_edge_area(j) * this%num_full_lon
    end do
    if (abs((this%total_area - total_area) / this%total_area) > 1.0d-11) then
      call log_error('Failed to calculate edge area!', __FILE__, __LINE__)
    end if

    do j = this%full_lat_start_idx_no_pole, this%full_lat_end_idx_no_pole
      this%vertex_lat_dist(j) = this%dlat * radius
      this%cell_lon_dist(j) = 2.0d0 * this%lon_edge_area(j) / this%vertex_lat_dist(j)
    end do

    do j = this%half_lat_start_idx_no_pole, this%half_lat_end_idx_no_pole
      this%vertex_lon_dist(j) = radius * this%half_cos_lat(j) * mesh%dlon
      this%cell_lat_dist(j) = 2.0d0 * this%lat_edge_area(j) / this%vertex_lon_dist(j)
    end do

    !  ____________________                 ____________________                  ____________________                  ____________________                 
    ! |          |         |               |          |         |                |          |         |                |          |         |                
    ! |          |         |               |          |         |                |          |         |                |          |         |                
    ! |          |         |               |          |         |                |          |         |                |          |         |                
    ! |          |         |               |          |         |                |          |         |                |          |         |                
    ! |_____o____|____o____|   j           |_____o____|____*____|   j            |_____*____|____o____|   j            |_____o____|____o____|   j             
    ! |          |////|////|               |          |////|    |                |/////|    |         |                |     |    |         |                
    ! |          |/3//|/2//|               |          |////|    |                |/////|    |         |                |     |    |         |                
    ! |          x---------|   j           |          x---------|   j            |-----|----x         |   j            |-----|----x         |   j            
    ! |          |    |/1//|               |          |    |    |                |/////|////|         |                |     |////|         |                
    ! |_____o____|____*____|   j - 1       |_____o____|____o____|   j - 1        |_____o____|____o____|   j - 1        |_____*____|____o____|   j - 1        
    !       i    i   i+1                         i    i   i+1                          i    i   i+1                          i
    !
    !
    !       [ 1    As_1 + As_2 + As_3]
    ! w = - [--- - ------------------]
    !       [ 2        A_{i+1,j}     ]
    !
    !

    select case (tangent_wgt_scheme)
    case ('classic')
      do j = this%full_lat_start_idx_no_pole, this%full_lat_end_idx_no_pole
#ifdef STAGGER_V_ON_POLE
        this%full_tangent_wgt(1,j) = this%vertex_lon_dist(j  ) / this%cell_lon_dist(j) * 0.25d0
        this%full_tangent_wgt(2,j) = this%vertex_lon_dist(j+1) / this%cell_lon_dist(j) * 0.25d0
#else
        this%full_tangent_wgt(1,j) = this%vertex_lon_dist(j-1) / this%cell_lon_dist(j) * 0.25d0
        this%full_tangent_wgt(2,j) = this%vertex_lon_dist(j  ) / this%cell_lon_dist(j) * 0.25d0
#endif
      end do

      do j = this%half_lat_start_idx_no_pole, this%half_lat_end_idx_no_pole
#ifdef STAGGER_V_ON_POLE
        this%half_tangent_wgt(1,j) = this%vertex_lat_dist(j-1) / this%cell_lat_dist(j) * 0.25d0
        this%half_tangent_wgt(2,j) = this%vertex_lat_dist(j  ) / this%cell_lat_dist(j) * 0.25d0
#else
        this%half_tangent_wgt(1,j) = this%vertex_lat_dist(j  ) / this%cell_lat_dist(j) * 0.25d0
        this%half_tangent_wgt(2,j) = this%vertex_lat_dist(j+1) / this%cell_lat_dist(j) * 0.25d0
#endif
      end do
    case ('thuburn09')
      do j = this%full_lat_start_idx_no_pole, this%full_lat_end_idx_no_pole
#ifdef STAGGER_V_ON_POLE
        this%full_tangent_wgt(1,j) = this%vertex_lon_dist(j  ) / this%cell_lon_dist(j) * this%subcell_area(2,j  ) / this%cell_area(j  )
        this%full_tangent_wgt(2,j) = this%vertex_lon_dist(j+1) / this%cell_lon_dist(j) * this%subcell_area(1,j  ) / this%cell_area(j  )
#else
        this%full_tangent_wgt(1,j) = this%vertex_lon_dist(j-1) / this%cell_lon_dist(j) * this%subcell_area(2,j  ) / this%cell_area(j  )
        this%full_tangent_wgt(2,j) = this%vertex_lon_dist(j  ) / this%cell_lon_dist(j) * this%subcell_area(1,j  ) / this%cell_area(j  )
#endif
      end do

      do j = this%half_lat_start_idx_no_pole, this%half_lat_end_idx_no_pole
#ifdef STAGGER_V_ON_POLE
        this%half_tangent_wgt(1,j) = this%vertex_lat_dist(j-1) / this%cell_lat_dist(j) * this%subcell_area(1,j-1) / this%cell_area(j-1)
        this%half_tangent_wgt(2,j) = this%vertex_lat_dist(j  ) / this%cell_lat_dist(j) * this%subcell_area(2,j  ) / this%cell_area(j  )
#else
        this%half_tangent_wgt(1,j) = this%vertex_lat_dist(j  ) / this%cell_lat_dist(j) * this%subcell_area(1,j  ) / this%cell_area(j  )
        this%half_tangent_wgt(2,j) = this%vertex_lat_dist(j+1) / this%cell_lat_dist(j) * this%subcell_area(2,j+1) / this%cell_area(j+1)
#endif
      end do
    end select

    do j = this%full_lat_start_idx, this%full_lat_end_idx
      this%full_upwind_beta(j) = 4 / pi**2 * this%full_lat(j)**2
    end do
    do j = this%half_lat_start_idx, this%half_lat_end_idx
      this%half_upwind_beta(j) = 4 / pi**2 * this%half_lat(j)**2
    end do

    do j = this%full_lat_start_idx, this%full_lat_end_idx
      this%full_f(j) = 2.0_r8 * omega * this%full_sin_lat(j)
    end do
    do j = this%half_lat_start_idx, this%half_lat_end_idx
      this%half_f(j) = 2.0_r8 * omega * this%half_sin_lat(j)
    end do

    call log_notice('Mesh module is initialized.')

  end subroutine mesh_init

  logical function mesh_has_south_pole(this) result(res)

    class(mesh_type), intent(in) :: this

    res = this%start_lat == -0.5_r8 * pi

  end function mesh_has_south_pole

  logical function mesh_has_north_pole(this) result(res)

    class(mesh_type), intent(in) :: this

    res = this%end_lat == 0.5_r8 * pi

  end function mesh_has_north_pole

  logical function mesh_is_south_pole(this, j) result(res)

    class(mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = this%has_south_pole() .and. j == 1

  end function mesh_is_south_pole

  logical function mesh_is_north_pole(this, j) result(res)

    class(mesh_type), intent(in) :: this
    integer, intent(in) :: j

#ifdef STAGGER_V_ON_POLE
    res = this%has_north_pole() .and. j == this%num_half_lat
#else
    res = this%has_north_pole() .and. j == this%num_full_lat
#endif

  end function mesh_is_north_pole

  logical function mesh_is_pole(this, j) result(res)

    class(mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = this%is_south_pole(j) .or. this%is_north_pole(j)

  end function mesh_is_pole

  subroutine mesh_final(this)

    type(mesh_type), intent(inout) :: this

    if (allocated(this%full_lon           )) deallocate(this%full_lon           )
    if (allocated(this%full_lat           )) deallocate(this%full_lat           )
    if (allocated(this%half_lon           )) deallocate(this%half_lon           )
    if (allocated(this%half_lat           )) deallocate(this%half_lat           )
    if (allocated(this%full_cos_lon       )) deallocate(this%full_cos_lon       )
    if (allocated(this%half_cos_lon       )) deallocate(this%half_cos_lon       )
    if (allocated(this%full_sin_lon       )) deallocate(this%full_sin_lon       )
    if (allocated(this%half_sin_lon       )) deallocate(this%half_sin_lon       )
    if (allocated(this%full_cos_lat       )) deallocate(this%full_cos_lat       )
    if (allocated(this%half_cos_lat       )) deallocate(this%half_cos_lat       )
    if (allocated(this%full_sin_lat       )) deallocate(this%full_sin_lat       )
    if (allocated(this%half_sin_lat       )) deallocate(this%half_sin_lat       )
    if (allocated(this%full_lon_deg       )) deallocate(this%full_lon_deg       )
    if (allocated(this%half_lon_deg       )) deallocate(this%half_lon_deg       )
    if (allocated(this%full_lat_deg       )) deallocate(this%full_lat_deg       )
    if (allocated(this%half_lat_deg       )) deallocate(this%half_lat_deg       )
    if (allocated(this%cell_area          )) deallocate(this%cell_area          )
    if (allocated(this%lon_edge_area      )) deallocate(this%lon_edge_area      )
    if (allocated(this%lon_edge_left_area )) deallocate(this%lon_edge_left_area )
    if (allocated(this%lon_edge_right_area)) deallocate(this%lon_edge_right_area)
    if (allocated(this%lat_edge_area      )) deallocate(this%lat_edge_area      )
    if (allocated(this%lat_edge_up_area   )) deallocate(this%lat_edge_up_area   )
    if (allocated(this%lat_edge_down_area )) deallocate(this%lat_edge_down_area )
    if (allocated(this%vertex_area        )) deallocate(this%vertex_area        )
    if (allocated(this%subcell_area       )) deallocate(this%subcell_area       )
    if (allocated(this%cell_lon_dist      )) deallocate(this%cell_lon_dist      )
    if (allocated(this%cell_lat_dist      )) deallocate(this%cell_lat_dist      )
    if (allocated(this%vertex_lon_dist    )) deallocate(this%vertex_lon_dist    )
    if (allocated(this%vertex_lat_dist    )) deallocate(this%vertex_lat_dist    )
    if (allocated(this%full_tangent_wgt   )) deallocate(this%full_tangent_wgt   )
    if (allocated(this%half_tangent_wgt   )) deallocate(this%half_tangent_wgt   )
    if (allocated(this%full_upwind_beta   )) deallocate(this%full_upwind_beta   )
    if (allocated(this%half_upwind_beta   )) deallocate(this%half_upwind_beta   )
    if (allocated(this%full_f             )) deallocate(this%full_f             )
    if (allocated(this%half_f             )) deallocate(this%half_f             )

  end subroutine mesh_final

  subroutine create_meshes()

    call mesh%init(num_lon, num_lat, halo_width=4)

  end subroutine create_meshes

end module mesh_mod
