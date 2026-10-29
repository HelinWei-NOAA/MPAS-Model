module mpas_gfs_base_pressure_option
  implicit none
  private

  public :: gfs_fixed_base_mode
  public :: gfs_mynn_style_mode
  public :: gfs_local_dp_mode
  public :: gfs_pressure_mode

  logical, save :: initialized = .false.
  integer, save :: pressure_mode = 0

contains

  subroutine initialize_gfs_pressure_mode()
    character(len=32) :: value
    integer :: status

    if (initialized) return

    value = ''
    call get_environment_variable( &
         'MPAS_GFS_PRESSURE_MODE', value, status=status)

    if (status == 0) then
      select case (trim(adjustl(value)))
      case ('0', 'original')
        pressure_mode = 0
      case ('1', 'fixed_base')
        pressure_mode = 1
      case ('2', 'mynn_style')
        pressure_mode = 2
      case ('3', 'local_dp')
        pressure_mode = 3
      case default
        error stop 'MPAS_GFS_PRESSURE_MODE must be 0, 1, 2, or 3'
      end select
    endif

    initialized = .true.
  end subroutine initialize_gfs_pressure_mode


  integer function gfs_pressure_mode()
    call initialize_gfs_pressure_mode()
    gfs_pressure_mode = pressure_mode
  end function gfs_pressure_mode


  logical function gfs_fixed_base_mode()
    call initialize_gfs_pressure_mode()
    gfs_fixed_base_mode = (pressure_mode == 1)
  end function gfs_fixed_base_mode


  logical function gfs_mynn_style_mode()
    call initialize_gfs_pressure_mode()
    gfs_mynn_style_mode = (pressure_mode == 2)
  end function gfs_mynn_style_mode

  logical function gfs_local_dp_mode()
    call initialize_gfs_pressure_mode()
    gfs_local_dp_mode = (pressure_mode == 3)
  end function gfs_local_dp_mode

end module mpas_gfs_base_pressure_option
