module mpas_sfc_diff_wrapper_mod

  use mpas_kind_types, only: RKIND
  use sfc_diff, only: sfc_diff_run

  implicit none

  public :: mpas_call_sfc_diff
  public :: gfs_sfcl_available
  public :: gfs_sfcl_garea, gfs_sfcl_zvfun, gfs_sfcl_zorl
  public :: gfs_sfcl_rb_lnd, gfs_sfcl_fm_lnd, gfs_sfcl_fh_lnd

  logical, save :: gfs_sfcl_available = .false.
  real(kind=RKIND), allocatable, save :: gfs_sfcl_garea(:)
  real(kind=RKIND), allocatable, save :: gfs_sfcl_zvfun(:)
  real(kind=RKIND), allocatable, save :: gfs_sfcl_zorl(:)
  real(kind=RKIND), allocatable, save :: gfs_sfcl_rb_lnd(:)
  real(kind=RKIND), allocatable, save :: gfs_sfcl_fm_lnd(:)
  real(kind=RKIND), allocatable, save :: gfs_sfcl_fh_lnd(:)

contains

  subroutine ensure_gfs_sfcl_cache(im)
    integer, intent(in) :: im

    if (allocated(gfs_sfcl_garea)) then
       if (size(gfs_sfcl_garea) /= im) then
          deallocate(gfs_sfcl_garea, gfs_sfcl_zvfun, gfs_sfcl_zorl)
          deallocate(gfs_sfcl_rb_lnd, gfs_sfcl_fm_lnd, gfs_sfcl_fh_lnd)
       endif
    endif

    if (.not. allocated(gfs_sfcl_garea)) then
       allocate(gfs_sfcl_garea(im), gfs_sfcl_zvfun(im), gfs_sfcl_zorl(im))
       allocate(gfs_sfcl_rb_lnd(im), gfs_sfcl_fm_lnd(im), gfs_sfcl_fh_lnd(im))
    endif
  end subroutine ensure_gfs_sfcl_cache

  subroutine mpas_call_sfc_diff(num_cells, num_lev, dt, areaCell,        &
                                psfc, t1, q1, z1, wspd,                 &
                                prsl1, prslki, prsik1, prslk1,          &
                                sigmaf, vegtype, shdmax, ivegsrc,       &
                                xland, xice,                            &
                                znt, z0h, tsk, u10, v10,                &
                                fm, fh, cm, ch, stress, rb_lnd_out,     &
                                errmsg, errflg)

    integer, intent(in) :: num_cells, num_lev
    integer, intent(in) :: ivegsrc
    integer, intent(in) :: vegtype(num_cells)

    real(kind=RKIND), intent(in) :: dt
    real(kind=RKIND), intent(in) :: areaCell(num_cells)
    real(kind=RKIND), intent(in) :: psfc(num_cells)
    real(kind=RKIND), intent(in) :: t1(num_cells)
    real(kind=RKIND), intent(in) :: q1(num_cells)
    real(kind=RKIND), intent(in) :: z1(num_cells)
    real(kind=RKIND), intent(in) :: wspd(num_cells)
    real(kind=RKIND), intent(in) :: prsl1(num_cells)
    real(kind=RKIND), intent(in) :: prslki(num_cells)
    real(kind=RKIND), intent(in) :: prsik1(num_cells)
    real(kind=RKIND), intent(in) :: prslk1(num_cells)
    real(kind=RKIND), intent(in) :: sigmaf(num_cells)
    real(kind=RKIND), intent(in) :: shdmax(num_cells)
    real(kind=RKIND), intent(in) :: xland(num_cells)
    real(kind=RKIND), intent(in) :: xice(num_cells)

    real(kind=RKIND), intent(inout) :: znt(num_cells)
    real(kind=RKIND), intent(inout) :: z0h(num_cells)
    real(kind=RKIND), intent(in) :: tsk(num_cells)
    real(kind=RKIND), intent(in) :: u10(num_cells), v10(num_cells)
    real(kind=RKIND), intent(inout) :: fm(num_cells), fh(num_cells)
    real(kind=RKIND), intent(inout) :: cm(num_cells), ch(num_cells)
    real(kind=RKIND), intent(inout) :: stress(num_cells)
    real(kind=RKIND), intent(inout) :: rb_lnd_out(num_cells)
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    integer :: i
    real(kind=RKIND), parameter :: grav = 9.80665_RKIND
    real(kind=RKIND), parameter :: rd = 287.0_RKIND
    real(kind=RKIND), parameter :: rv = 461.5_RKIND
    real(kind=RKIND), parameter :: eps = rd / rv
    real(kind=RKIND), parameter :: epsm1 = eps - 1.0_RKIND
    real(kind=RKIND), parameter :: rvrdm1 = rv / rd - 1.0_RKIND

    logical :: redrag, use_oceanuv, thsfc_loc
    integer :: sfc_z0_type
    logical, allocatable :: flag_iter(:), flag_lakefreeze(:)
    logical, allocatable :: wet(:), dry(:), icy(:)
    real(kind=RKIND), allocatable :: z0pert(:), ztpert(:)
    real(kind=RKIND), allocatable :: lakefrac(:), fice(:)
    real(kind=RKIND), allocatable :: u1(:), v1(:), usfco(:), vsfco(:)
    real(kind=RKIND), allocatable :: tskin_wat(:), tskin_lnd(:), tskin_ice(:)
    real(kind=RKIND), allocatable :: tsurf_wat(:), tsurf_lnd(:), tsurf_ice(:)
    real(kind=RKIND), allocatable :: z0rl_wat(:), z0rl_lnd(:), z0rl_ice(:), z0rl_wav(:)
    real(kind=RKIND), allocatable :: ustar_wat(:), ustar_lnd(:), ustar_ice(:)
    real(kind=RKIND), allocatable :: cm_wat(:), cm_lnd(:), cm_ice(:)
    real(kind=RKIND), allocatable :: ch_wat(:), ch_lnd(:), ch_ice(:)
    real(kind=RKIND), allocatable :: rb_wat(:), rb_lnd(:), rb_ice(:)
    real(kind=RKIND), allocatable :: stress_wat(:), stress_lnd(:), stress_ice(:)
    real(kind=RKIND), allocatable :: fm_wat(:), fm_lnd(:), fm_ice(:)
    real(kind=RKIND), allocatable :: fh_wat(:), fh_lnd(:), fh_ice(:)
    real(kind=RKIND), allocatable :: fm10_wat(:), fm10_lnd(:), fm10_ice(:)
    real(kind=RKIND), allocatable :: fh2_wat(:), fh2_lnd(:), fh2_ice(:)
    real(kind=RKIND), allocatable :: ztmax_wat(:), ztmax_lnd(:), ztmax_ice(:)
    real(kind=RKIND), allocatable :: zvfun(:)

    errmsg = ''
    errflg = 0

    call ensure_gfs_sfcl_cache(num_cells)
    gfs_sfcl_available = .false.

    allocate(flag_iter(num_cells), flag_lakefreeze(num_cells))
    allocate(wet(num_cells), dry(num_cells), icy(num_cells))
    allocate(z0pert(num_cells), ztpert(num_cells))
    allocate(lakefrac(num_cells), fice(num_cells))
    allocate(u1(num_cells), v1(num_cells), usfco(num_cells), vsfco(num_cells))
    allocate(tskin_wat(num_cells), tskin_lnd(num_cells), tskin_ice(num_cells))
    allocate(tsurf_wat(num_cells), tsurf_lnd(num_cells), tsurf_ice(num_cells))
    allocate(z0rl_wat(num_cells), z0rl_lnd(num_cells), z0rl_ice(num_cells), z0rl_wav(num_cells))
    allocate(ustar_wat(num_cells), ustar_lnd(num_cells), ustar_ice(num_cells))
    allocate(cm_wat(num_cells), cm_lnd(num_cells), cm_ice(num_cells))
    allocate(ch_wat(num_cells), ch_lnd(num_cells), ch_ice(num_cells))
    allocate(rb_wat(num_cells), rb_lnd(num_cells), rb_ice(num_cells))
    allocate(stress_wat(num_cells), stress_lnd(num_cells), stress_ice(num_cells))
    allocate(fm_wat(num_cells), fm_lnd(num_cells), fm_ice(num_cells))
    allocate(fh_wat(num_cells), fh_lnd(num_cells), fh_ice(num_cells))
    allocate(fm10_wat(num_cells), fm10_lnd(num_cells), fm10_ice(num_cells))
    allocate(fh2_wat(num_cells), fh2_lnd(num_cells), fh2_ice(num_cells))
    allocate(ztmax_wat(num_cells), ztmax_lnd(num_cells), ztmax_ice(num_cells))
    allocate(zvfun(num_cells))

    redrag = .false.
    use_oceanuv = .false.
    thsfc_loc = .false.
    sfc_z0_type = 0

    flag_iter = .true.
    flag_lakefreeze = .false.
    z0pert = 0.0_RKIND
    ztpert = 0.0_RKIND
    lakefrac = 0.0_RKIND
    fice = max(0.0_RKIND, min(1.0_RKIND, xice))
    u1 = u10
    v1 = v10
    usfco = 0.0_RKIND
    vsfco = 0.0_RKIND

    tskin_wat = tsk
    tskin_lnd = tsk
    tskin_ice = tsk
    tsurf_wat = tsk
    tsurf_lnd = tsk
    tsurf_ice = tsk

    do i = 1, num_cells
       ! MPAS xland convention is normally 1=land, 2=water.
       dry(i) = (xland(i) < 1.5_RKIND)
       wet(i) = .not. dry(i)
       icy(i) = (fice(i) > 0.5_RKIND)

       ! GFS sfc_diff expects z0rl in cm and internally uses 0.01*z0rl.
       z0rl_lnd(i) = max(znt(i), 1.0e-6_RKIND) * 100.0_RKIND
       z0rl_wat(i) = max(znt(i), 1.0e-6_RKIND) * 100.0_RKIND
       z0rl_ice(i) = max(znt(i), 1.0e-6_RKIND) * 100.0_RKIND
       z0rl_wav(i) = max(znt(i), 1.0e-6_RKIND) * 100.0_RKIND

       ustar_lnd(i) = max(stress(i), 0.0_RKIND)
       ustar_wat(i) = ustar_lnd(i)
       ustar_ice(i) = ustar_lnd(i)
       cm_lnd(i) = max(cm(i), 1.0e-8_RKIND)
       cm_wat(i) = cm_lnd(i)
       cm_ice(i) = cm_lnd(i)
       ch_lnd(i) = max(ch(i), 1.0e-8_RKIND)
       ch_wat(i) = ch_lnd(i)
       ch_ice(i) = ch_lnd(i)
       rb_lnd(i) = rb_lnd_out(i)
       rb_wat(i) = rb_lnd(i)
       rb_ice(i) = rb_lnd(i)
       stress_lnd(i) = max(stress(i), 0.0_RKIND)
       stress_wat(i) = stress_lnd(i)
       stress_ice(i) = stress_lnd(i)
       fm_lnd(i) = max(fm(i), 1.0e-6_RKIND)
       fm_wat(i) = fm_lnd(i)
       fm_ice(i) = fm_lnd(i)
       fh_lnd(i) = max(fh(i), 1.0e-6_RKIND)
       fh_wat(i) = fh_lnd(i)
       fh_ice(i) = fh_lnd(i)
       fm10_lnd(i) = fm_lnd(i)
       fm10_wat(i) = fm_lnd(i)
       fm10_ice(i) = fm_lnd(i)
       fh2_lnd(i) = fh_lnd(i)
       fh2_wat(i) = fh_lnd(i)
       fh2_ice(i) = fh_lnd(i)
       ztmax_lnd(i) = max(z0h(i), 1.0e-6_RKIND)
       ztmax_wat(i) = ztmax_lnd(i)
       ztmax_ice(i) = ztmax_lnd(i)
       zvfun(i) = 0.0_RKIND
    enddo

    call sfc_diff_run(num_cells, rvrdm1, eps, epsm1, grav,              &
         psfc, t1, q1, z1, areaCell, wspd,                              &
         prsl1, prslki, prsik1, prslk1,                                 &
         sigmaf, vegtype, shdmax, ivegsrc,                              &
         z0pert, ztpert, flag_iter, redrag,                             &
         flag_lakefreeze, lakefrac, fice,                               &
         u10, v10, sfc_z0_type, u1, v1, usfco, vsfco, use_oceanuv,      &
         wet, dry, icy, thsfc_loc,                                      &
         tskin_wat, tskin_lnd, tskin_ice,                               &
         tsurf_wat, tsurf_lnd, tsurf_ice,                               &
         z0rl_wat, z0rl_lnd, z0rl_ice, z0rl_wav,                        &
         ustar_wat, ustar_lnd, ustar_ice,                               &
         cm_wat, cm_lnd, cm_ice,                                       &
         ch_wat, ch_lnd, ch_ice,                                       &
         rb_wat, rb_lnd, rb_ice,                                       &
         stress_wat, stress_lnd, stress_ice,                           &
         fm_wat, fm_lnd, fm_ice,                                       &
         fh_wat, fh_lnd, fh_ice,                                       &
         fm10_wat, fm10_lnd, fm10_ice,                                 &
         fh2_wat, fh2_lnd, fh2_ice,                                    &
         ztmax_wat, ztmax_lnd, ztmax_ice,                              &
         zvfun, errmsg, errflg)

    if (errflg /= 0) return

    do i = 1, num_cells
       if (dry(i)) then
          znt(i) = max(0.01_RKIND*z0rl_lnd(i), 1.0e-6_RKIND)
          z0h(i) = max(ztmax_lnd(i), 1.0e-6_RKIND)
          stress(i) = max(stress_lnd(i), 0.0_RKIND)
          cm(i) = max(cm_lnd(i), 1.0e-8_RKIND)
          ch(i) = max(ch_lnd(i), 1.0e-8_RKIND)
          rb_lnd_out(i) = rb_lnd(i)
          fm(i) = max(fm_lnd(i), 1.0e-6_RKIND)
          fh(i) = max(fh_lnd(i), 1.0e-6_RKIND)
       elseif (icy(i)) then
          znt(i) = max(0.01_RKIND*z0rl_ice(i), 1.0e-6_RKIND)
          z0h(i) = max(ztmax_ice(i), 1.0e-6_RKIND)
          stress(i) = max(stress_ice(i), 0.0_RKIND)
          cm(i) = max(cm_ice(i), 1.0e-8_RKIND)
          ch(i) = max(ch_ice(i), 1.0e-8_RKIND)
          rb_lnd_out(i) = rb_ice(i)
          fm(i) = max(fm_ice(i), 1.0e-6_RKIND)
          fh(i) = max(fh_ice(i), 1.0e-6_RKIND)
       else
          znt(i) = max(0.01_RKIND*z0rl_wat(i), 1.0e-6_RKIND)
          z0h(i) = max(ztmax_wat(i), 1.0e-6_RKIND)
          stress(i) = max(stress_wat(i), 0.0_RKIND)
          cm(i) = max(cm_wat(i), 1.0e-8_RKIND)
          ch(i) = max(ch_wat(i), 1.0e-8_RKIND)
          rb_lnd_out(i) = rb_wat(i)
          fm(i) = max(fm_wat(i), 1.0e-6_RKIND)
          fh(i) = max(fh_wat(i), 1.0e-6_RKIND)
       endif

       gfs_sfcl_garea(i) = areaCell(i)
       gfs_sfcl_zvfun(i) = zvfun(i)
       gfs_sfcl_zorl(i) = znt(i)
       gfs_sfcl_rb_lnd(i) = rb_lnd(i)
       gfs_sfcl_fm_lnd(i) = fm_lnd(i)
       gfs_sfcl_fh_lnd(i) = fh_lnd(i)
       print*,'inside sfcl wrapper',gfs_sfcl_garea(i),gfs_sfcl_zvfun(i),gfs_sfcl_zorl(i),gfs_sfcl_rb_lnd(i),gfs_sfcl_fm_lnd(i),gfs_sfcl_fh_lnd(i)
    enddo

    gfs_sfcl_available = .true.

  end subroutine mpas_call_sfc_diff

end module mpas_sfc_diff_wrapper_mod
