module parallel_mod

  use mpi
  use mesh_mod

  implicit none

  private

  public parallel_init
  public parallel_final
  public parallel_fill_halo
  public parallel_zero_halo
  public parallel_overlay_inner_halo
  public parallel_zonal_sum

  type parallel_info_type
    integer comm
    integer num_proc
    integer rank
    character(30) proc_name
  end type parallel_info_type

  type(parallel_info_type) parallel_info

  interface parallel_fill_halo
    module procedure parallel_fill_halo_1d_r8_1
    module procedure parallel_fill_halo_1d_r8_2
    module procedure parallel_fill_halo_2d_r8
  end interface parallel_fill_halo

  interface parallel_zero_halo
    module procedure parallel_zero_halo_1d_r8
  end interface parallel_zero_halo

  interface parallel_zonal_sum
    module procedure paralle_zonal_sum_0d_r8
  end interface parallel_zonal_sum

contains

  subroutine parallel_init()

    integer ierr
    integer n

    parallel_info%comm = MPI_COMM_WORLD

    ! call MPI_INIT(ierr)
    ! call MPI_CART_CREATE(MPI_COMM_WORLD, 2, [2, 2], [.true., .false.], .true., parallel_info%comm, ierr)
    ! call MPI_COMM_SIZE(parallel_info%comm, parallel_info%num_proc, ierr)
    ! call MPI_COMM_RANK(parallel_info%comm, parallel_info%rank, ierr)
    ! call MPI_GET_PROCESSOR_NAME(parallel_info%proc_name, n, ierr)

  end subroutine parallel_init

  subroutine parallel_final()

    integer ierr

    ! call MPI_FINALIZE(ierr)

  end subroutine parallel_final

  subroutine parallel_fill_halo_1d_r8_1(halo_width, array, left_halo, right_halo)

    integer, intent(in   )           :: halo_width
    real(8), intent(inout)           :: array(:)
    logical, intent(in   ), optional :: left_halo
    logical, intent(in   ), optional :: right_halo

    integer i, m, n

    if (merge(left_halo, .true., present(left_halo))) then
      m = lbound(array, 1) - 1
      n = ubound(array, 1) - 2 * halo_width
      do i = 1, halo_width
        array(m+i) = array(n+i)
      end do
    end if

    if (merge(right_halo, .true., present(right_halo))) then
      m = ubound(array, 1) - halo_width
      n = lbound(array, 1) + halo_width - 1
      do i = 1, halo_width
        array(m+i) = array(n+i)
      end do
    end if

  end subroutine parallel_fill_halo_1d_r8_1

  subroutine parallel_fill_halo_1d_r8_2(mesh, array, left_halo, right_halo)

    type(mesh_type), intent(in   )           :: mesh
    real(8)        , intent(inout)           :: array(:)
    logical        , intent(in   ), optional :: left_halo
    logical        , intent(in   ), optional :: right_halo

    integer i, m, n

    if (merge(left_halo, .true., present(left_halo))) then
      m = lbound(array, 1) - 1
      n = ubound(array, 1) - 2 * mesh%halo_width
      do i = 1, mesh%halo_width
        array(m+i) = array(n+i)
      end do
    end if

    if (merge(right_halo, .true., present(right_halo))) then
      m = ubound(array, 1) - mesh%halo_width
      n = lbound(array, 1) + mesh%halo_width - 1
      do i = 1, mesh%halo_width
        array(m+i) = array(n+i)
      end do
    end if

  end subroutine parallel_fill_halo_1d_r8_2

  subroutine parallel_fill_halo_2d_r8(mesh, array, left_halo, right_halo, top_halo, bottom_halo)

    type(mesh_type), intent(in   )           :: mesh
    real(8)        , intent(inout)           :: array(:,:)
    logical        , intent(in   ), optional :: left_halo
    logical        , intent(in   ), optional :: right_halo
    logical        , intent(in   ), optional :: top_halo
    logical        , intent(in   ), optional :: bottom_halo

    integer i, j, m, n

    if (merge(left_halo, .true., present(left_halo))) then
      m = lbound(array, 1) - 1
      n = ubound(array, 1) - 2 * mesh%halo_width
      do j = lbound(array, 2), ubound(array, 2)
        do i = 1, mesh%halo_width
          array(m+i,j) = array(n+i,j)
        end do
      end do
    end if

    if (merge(right_halo, .true., present(right_halo))) then
      m = ubound(array, 1) - mesh%halo_width
      n = lbound(array, 1) + mesh%halo_width - 1
      do j = lbound(array, 2), ubound(array, 2)
        do i = 1, mesh%halo_width
          array(m+i,j) = array(n+i,j)
        end do
      end do
    end if

  end subroutine parallel_fill_halo_2d_r8

  subroutine parallel_zero_halo_1d_r8(mesh, array, left_halo, right_halo)

    type(mesh_type), intent(in   )           :: mesh
    real(8)        , intent(inout)           :: array(mesh%full_lon_lb:mesh%full_lon_ub)
    logical        , intent(in   ), optional :: left_halo
    logical        , intent(in   ), optional :: right_halo

    if (merge(left_halo, .false., present(left_halo))) then
      array(mesh%full_lon_lb:mesh%full_lon_start_idx-1) = 0.0d0
    end if

    if (merge(right_halo, .false., present(right_halo))) then
      array(mesh%full_lon_end_idx+1:mesh%full_lon_ub) = 0.0d0
    end if

  end subroutine parallel_zero_halo_1d_r8

  subroutine parallel_overlay_inner_halo(mesh, array, left_halo, right_halo)

    type(mesh_type), intent(in   )           :: mesh
    real(8)        , intent(inout)           :: array(mesh%full_lon_lb:mesh%full_lon_ub)
    logical        , intent(in   ), optional :: left_halo
    logical        , intent(in   ), optional :: right_halo

    integer i

    if (merge(left_halo, .false., present(left_halo))) then
      do i = mesh%full_lon_start_idx, mesh%full_lon_start_idx + mesh%halo_width - 2
        array(i) = array(i) + array(mesh%full_lon_end_idx+i-mesh%full_lon_start_idx+1)
      end do
    end if

  end subroutine parallel_overlay_inner_halo

  subroutine paralle_zonal_sum_0d_r8(value)

    real(8), intent(inout) :: value

  end subroutine paralle_zonal_sum_0d_r8

end module parallel_mod
