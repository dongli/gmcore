module state_mod

  use const_mod
  use mesh_mod
  use allocator_mod

  implicit none

  private

  public state_type
  public states
  public create_states

  type state_type
    type(mesh_type), pointer :: mesh => null()
    ! Prognostic variables
    real(r8), allocatable, dimension(:,:) :: u
    real(r8), allocatable, dimension(:,:) :: v
    real(r8), allocatable, dimension(:,:) :: gd
    ! Diagnostic variables
    real(r8), allocatable, dimension(:,:) :: pv
    real(r8), allocatable, dimension(:,:) :: m_vtx
    real(r8), allocatable, dimension(:,:) :: m_lon
    real(r8), allocatable, dimension(:,:) :: m_lat
    real(r8), allocatable, dimension(:,:) :: mf_lon_n
    real(r8), allocatable, dimension(:,:) :: mf_lat_n
    real(r8), allocatable, dimension(:,:) :: mf_lon_t
    real(r8), allocatable, dimension(:,:) :: mf_lat_t
    real(r8), allocatable, dimension(:,:) :: pv_lon
    real(r8), allocatable, dimension(:,:) :: pv_lat
    real(r8), allocatable, dimension(:,:) :: dpv_lon_t
    real(r8), allocatable, dimension(:,:) :: dpv_lat_n
    real(r8), allocatable, dimension(:,:) :: dpv_lat_t
    real(r8), allocatable, dimension(:,:) :: dpv_lon_n
    real(r8), allocatable, dimension(:,:) :: ke_cell
    real(r8), allocatable, dimension(:,:) :: pvc_lon ! Corrected PV on U grids
    real(r8), allocatable, dimension(:,:) :: pvc_lat ! Corrected PV on V grids
    real(r8) total_m
    real(r8) total_ke
    real(r8) total_e
    real(r8) total_av
    real(r8) total_pe
  contains
    procedure :: init => state_init
    procedure :: clear => state_clear
    final :: state_final
  end type state_type

  type(state_type), allocatable, target :: states(:)

contains

  subroutine create_states()

    integer i

    if (.not. allocated(states)) then
      allocate(states(0:2))
      do i = lbound(states, 1), ubound(states, 1)
        call states(i)%init(mesh)
      end do
    end if

  end subroutine create_states

  subroutine state_init(this, mesh)

    class(state_type), intent(inout)         :: this
    type(mesh_type  ), intent(in   ), target :: mesh

    call this%clear()

    this%mesh => mesh

    call allocate_array(mesh, this%u        , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%v        , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%gd       , full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%pv       , half_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%m_vtx    , half_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%m_lon    , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%m_lat    , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%mf_lon_n , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%mf_lat_n , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%mf_lon_t , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%mf_lat_t , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%pv_lon   , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%pv_lat   , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%dpv_lon_t, half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%dpv_lon_n, half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%dpv_lat_t, full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%dpv_lat_n, full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%ke_cell  , full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%pvc_lon  , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%pvc_lat  , full_lon=.true., half_lat=.true.)

  end subroutine state_init

  subroutine state_clear(this)

    class(state_type), intent(inout) :: this

    if (allocated(this%u        )) deallocate(this%u        )
    if (allocated(this%v        )) deallocate(this%v        )
    if (allocated(this%gd       )) deallocate(this%gd       )
    if (allocated(this%pv       )) deallocate(this%pv       )
    if (allocated(this%m_vtx    )) deallocate(this%m_vtx    )
    if (allocated(this%m_lon    )) deallocate(this%m_lon    )
    if (allocated(this%m_lat    )) deallocate(this%m_lat    )
    if (allocated(this%mf_lon_n )) deallocate(this%mf_lon_n )
    if (allocated(this%mf_lat_n )) deallocate(this%mf_lat_n )
    if (allocated(this%mf_lon_t )) deallocate(this%mf_lon_t )
    if (allocated(this%mf_lat_t )) deallocate(this%mf_lat_t )
    if (allocated(this%pv_lon   )) deallocate(this%pv_lon   )
    if (allocated(this%pv_lat   )) deallocate(this%pv_lat   )
    if (allocated(this%dpv_lon_t)) deallocate(this%dpv_lon_t)
    if (allocated(this%dpv_lat_n)) deallocate(this%dpv_lat_n)
    if (allocated(this%dpv_lat_t)) deallocate(this%dpv_lat_t)
    if (allocated(this%dpv_lon_n)) deallocate(this%dpv_lon_n)
    if (allocated(this%ke_cell  )) deallocate(this%ke_cell  )
    if (allocated(this%pvc_lon  )) deallocate(this%pvc_lon  )
    if (allocated(this%pvc_lat  )) deallocate(this%pvc_lat  )

  end subroutine state_clear

  subroutine state_final(this)

    type(state_type), intent(inout) :: this

    call this%clear()

  end subroutine state_final

end module state_mod
