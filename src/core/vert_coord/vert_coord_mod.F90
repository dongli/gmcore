module vert_coord_mod

  use const_mod
  use namelist_mod
  use sigma_coord_mod
  use hybrid_coord_mod

  implicit none

  private

  public vert_coord_init
  public vert_coord_final
  public vert_coord_calc_ph
  public vert_coord_calc_ph_lev
  public vert_coord_calc_dphdt_lev

  interface
    pure real(r8) function vert_coord_calc_ph_interface(k, phs)
      import r8
      integer, intent(in) :: k
      real(r8), intent(in) :: phs
    end function vert_coord_calc_ph_interface

    pure real(r8) function vert_coord_calc_ph_lev_interface(k, phs)
      import r8
      integer, intent(in) :: k
      real(r8), intent(in) :: phs
    end function vert_coord_calc_ph_lev_interface

    pure real(r8) function vert_coord_calc_dphdt_lev_interface(k, dphsdt)
      import r8
      integer, intent(in) :: k
      real(r8), intent(in) :: dphsdt
    end function vert_coord_calc_dphdt_lev_interface
  end interface

  procedure(vert_coord_calc_ph_interface), pointer :: vert_coord_calc_ph
  procedure(vert_coord_calc_ph_lev_interface), pointer :: vert_coord_calc_ph_lev
  procedure(vert_coord_calc_dphdt_lev_interface), pointer :: vert_coord_calc_dphdt_lev

contains

  subroutine vert_coord_init(num_lev, namelist_file, template)

    integer, intent(in) :: num_lev
    character(*), intent(in), optional :: namelist_file
    character(*), intent(in), optional :: template

    select case (vert_coord_scheme)
    case ('sigma')
      call sigma_coord_init(num_lev, namelist_file, template)
      vert_coord_calc_ph => sigma_coord_calc_ph
      vert_coord_calc_ph_lev => sigma_coord_calc_ph_lev
      vert_coord_calc_dphdt_lev => sigma_coord_calc_dphdt_lev
    case ('hybrid')
      call hybrid_coord_init(num_lev, namelist_file, template)
      vert_coord_calc_ph => hybrid_coord_calc_ph
      vert_coord_calc_ph_lev => hybrid_coord_calc_ph_lev
      vert_coord_calc_dphdt_lev => hybrid_coord_calc_dphdt_lev
    end select

  end subroutine vert_coord_init

  subroutine vert_coord_final()

    select case (vert_coord_scheme)
    case ('sigma')
      call sigma_coord_final()
    case ('hybrid')
      call hybrid_coord_final()
    end select

  end subroutine vert_coord_final

end module vert_coord_mod
