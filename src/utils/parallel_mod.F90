module parallel_mod

  use mpi
  use mesh_mod

  implicit none

  private

  public parallel_init
  public parallel_final
  public parallel_fill_halo
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

  subroutine parallel_fill_halo_1d_r8_1(halo_width, array, all_halo, left_halo, right_halo)

    integer, intent(in   )           :: halo_width
    real(8), intent(inout)           :: array(:)
    logical, intent(in   ), optional :: all_halo
    logical, intent(in   ), optional :: left_halo
    logical, intent(in   ), optional :: right_halo

    integer i, m, n

    if (merge(all_halo, .true., present(all_halo)) .or. merge(left_halo, .true., present(left_halo))) then
      m = lbound(array, 1) - 1
      n = ubound(array, 1) - 2 * halo_width
      do i = 1, halo_width
        array(m+i) = array(n+i)
      end do
    end if

    if (merge(all_halo, .true., present(all_halo)) .or. merge(right_halo, .true., present(right_halo))) then
      m = ubound(array, 1) - halo_width
      n = lbound(array, 1) + halo_width - 1
      do i = 1, halo_width
        array(m+i) = array(n+i)
      end do
    end if

  end subroutine parallel_fill_halo_1d_r8_1

  subroutine parallel_fill_halo_1d_r8_2(mesh, array, all_halo, left_halo, right_halo)

    type(mesh_type), intent(in   )           :: mesh
    real(8)        , intent(inout)           :: array(:)
    logical        , intent(in   ), optional :: all_halo
    logical        , intent(in   ), optional :: left_halo
    logical        , intent(in   ), optional :: right_halo

    integer i, m, n

    if (merge(all_halo, .true., present(all_halo)) .or. merge(left_halo, .true., present(left_halo))) then
      m = lbound(array, 1) - 1
      n = ubound(array, 1) - 2 * mesh%halo_width
      do i = 1, mesh%halo_width
        array(m+i) = array(n+i)
      end do
    end if

    if (merge(all_halo, .true., present(all_halo)) .or. merge(right_halo, .true., present(right_halo))) then
      m = ubound(array, 1) - mesh%halo_width
      n = lbound(array, 1) + mesh%halo_width - 1
      do i = 1, mesh%halo_width
        array(m+i) = array(n+i)
      end do
    end if

  end subroutine parallel_fill_halo_1d_r8_2

  subroutine parallel_fill_halo_2d_r8(mesh, array, all_halo, left_halo, right_halo, top_halo, bottom_halo)

    type(mesh_type), intent(in   )           :: mesh
    real(8)        , intent(inout)           :: array(:,:)
    logical        , intent(in   ), optional :: all_halo
    logical        , intent(in   ), optional :: left_halo
    logical        , intent(in   ), optional :: right_halo
    logical        , intent(in   ), optional :: top_halo
    logical        , intent(in   ), optional :: bottom_halo

    integer i, j, m, n

    if (merge(all_halo, .true., present(all_halo)) .or. merge(left_halo, .true., present(left_halo))) then
      m = lbound(array, 1) - 1
      n = ubound(array, 1) - 2 * mesh%halo_width
      do j = lbound(array, 2), ubound(array, 2)
        do i = 1, mesh%halo_width
          array(m+i,j) = array(n+i,j)
        end do
      end do
    end if

    if (merge(all_halo, .true., present(all_halo)) .or. merge(right_halo, .true., present(right_halo))) then
      m = ubound(array, 1) - mesh%halo_width
      n = lbound(array, 1) + mesh%halo_width - 1
      do j = lbound(array, 2), ubound(array, 2)
        do i = 1, mesh%halo_width
          array(m+i,j) = array(n+i,j)
        end do
      end do
    end if

  end subroutine parallel_fill_halo_2d_r8

  subroutine paralle_zonal_sum_0d_r8(value)

    real(8), intent(inout) :: value

  end subroutine paralle_zonal_sum_0d_r8

end module parallel_mod
