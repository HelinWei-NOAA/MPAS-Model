module mpas_satmedmfvdifq_wrapper_mod
! COMBINED DEBUG PACKAGE: qfx direct, z0 capped, pressure/vertical/unit checks.

  use mpas_kind_types, only: RKIND
  use satmedmfvdifq, only: satmedmfvdifq_run
  implicit none

  type mpas_satmedmfvdifq_config_type
    logical :: sa3dtke      = .false.
    logical :: tte_edmf     = .false.
    logical :: dspheat      = .true.
    logical :: use_oceanuv  = .false.
    logical :: do_canopy    = .false.
    logical :: cplaqm       = .false.
    logical :: gen_tend     = .false.
    logical :: ldiag3d      = .false.

    integer :: sfc_rlm = 0
    integer :: tc_pbl  = 0
    integer :: use_lpt = 0

    real(kind=RKIND) :: xkzm_m = 1.0_RKIND
    real(kind=RKIND) :: xkzm_h = 1.0_RKIND
    real(kind=RKIND) :: xkzm_s = 1.0_RKIND
    real(kind=RKIND) :: dspfac = 1.0_RKIND
    real(kind=RKIND) :: bl_upfr = 0.13_RKIND
    real(kind=RKIND) :: bl_dnfr = 0.1_RKIND
    real(kind=RKIND) :: rlmx = 300.0_RKIND
    real(kind=RKIND) :: elmx = 300.0_RKIND
  end type mpas_satmedmfvdifq_config_type

contains

  subroutine mpas_call_satmedmfvdifq(nCells, nVertLevels, ntrac, dt,      &
                                     z_mid, z_int, areaCell,             &
                                     u_mpas, v_mpas, t_mpas,             &
                                     qv_mpas, qc_mpas, qi_mpas, tke_mpas,&
                                     p_mid, p_int, exner_mid,            &
                                     sw_heat, lw_heat, coszen,           &
                                     skin_temp, shflx, qfx_mpas, stress_in, &
                                     z0_mpas, u10, v10, vegfra_in, rb_in, fm_in, fh_in,    &
                                     hpbl_out, kpbl_out,                 &
                                     ten_t_out, ten_u_out, ten_v_out,    &
                                     ten_qv_out, ten_qc_out, ten_qi_out, &
                                     cfg, errmsg, errflg)

    integer, intent(in) :: nCells, nVertLevels, ntrac
    real(kind=RKIND), intent(in) :: dt
    real(kind=RKIND), intent(in) :: z_mid(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: z_int(nCells,nVertLevels+1)
    real(kind=RKIND), intent(in) :: areaCell(nCells)
    real(kind=RKIND), intent(in) :: u_mpas(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: v_mpas(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: t_mpas(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: qv_mpas(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: qc_mpas(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: qi_mpas(nCells,nVertLevels)
    real(kind=RKIND), intent(inout) :: tke_mpas(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: p_mid(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: p_int(nCells,nVertLevels+1)
    real(kind=RKIND), intent(in) :: exner_mid(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: sw_heat(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: lw_heat(nCells,nVertLevels)
    real(kind=RKIND), intent(in) :: coszen(nCells)
    real(kind=RKIND), intent(in) :: skin_temp(nCells)
    real(kind=RKIND), intent(in) :: shflx(nCells)
    real(kind=RKIND), intent(in) :: qfx_mpas(nCells)
    real(kind=RKIND), intent(in) :: stress_in(nCells)
    real(kind=RKIND), intent(in) :: z0_mpas(nCells)
    real(kind=RKIND), intent(in) :: u10(nCells), v10(nCells)
    real(kind=RKIND), intent(in) :: vegfra_in(nCells),rb_in(nCells)
    real(kind=RKIND), intent(in) :: fm_in(nCells), fh_in(nCells)
    real(kind=RKIND), intent(out) :: hpbl_out(nCells)
    integer, intent(out) :: kpbl_out(nCells)
    real(kind=RKIND), intent(out) :: ten_t_out(nCells,nVertLevels)
    real(kind=RKIND), intent(out) :: ten_u_out(nCells,nVertLevels)
    real(kind=RKIND), intent(out) :: ten_v_out(nCells,nVertLevels)
    real(kind=RKIND), intent(out) :: ten_qv_out(nCells,nVertLevels)
    real(kind=RKIND), intent(out) :: ten_qc_out(nCells,nVertLevels)
    real(kind=RKIND), intent(out) :: ten_qi_out(nCells,nVertLevels)
    type(mpas_satmedmfvdifq_config_type), intent(in) :: cfg
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    integer :: im, km
    integer :: i, k
    integer :: kcheck, ibad, kbad
    integer :: nbad_tke_input, nbad_tke_output
    integer :: ntqv, ntcw, ntiw, ntrw, ntke
    integer :: index_of_temperature, index_of_x_wind
    integer :: index_of_y_wind, index_of_process_pbl

    real(kind=RKIND), parameter :: grav  = 9.80665_RKIND
    real(kind=RKIND), parameter :: pi    = 3.14159265358979323846_RKIND
    real(kind=RKIND), parameter :: rd    = 287.0_RKIND
    real(kind=RKIND), parameter :: cp    = 1004.0_RKIND
    real(kind=RKIND), parameter :: rv    = 461.5_RKIND
    real(kind=RKIND), parameter :: hvap  = 2.5e6_RKIND
    real(kind=RKIND), parameter :: hfus  = 3.3358e5_RKIND
    real(kind=RKIND), parameter :: fv    = rv/rd - 1.0_RKIND
    real(kind=RKIND), parameter :: eps   = rd/rv
    real(kind=RKIND), parameter :: epsm1 = eps - 1.0_RKIND
    real(kind=RKIND), parameter :: z0lo = 0.1_RKIND
    real(kind=RKIND), parameter :: z0up = 1.0_RKIND

    real(kind=RKIND), allocatable :: rtg(:,:,:), q1(:,:,:)
    real(kind=RKIND), allocatable :: u1(:,:), v1(:,:), t1(:,:)
    real(kind=RKIND), allocatable :: swh(:,:), hlw(:,:), xmu(:)
    real(kind=RKIND), allocatable :: garea(:), zvfun(:), sigmaf(:)
    real(kind=RKIND), allocatable :: psk(:), rbsoil(:), zorl(:)
    real(kind=RKIND), allocatable :: u10m(:), v10m(:), fm(:), fh(:)
    real(kind=RKIND), allocatable :: tsea(:), heat(:), evap(:)
    real(kind=RKIND), allocatable :: stress(:), spd1(:)
    real(kind=RKIND), allocatable :: prsi(:,:), del(:,:), prsl(:,:), prslk(:,:)
    real(kind=RKIND), allocatable :: phii(:,:), phil(:,:)
    real(kind=RKIND), allocatable :: dusfc(:), dvsfc(:), dtsfc(:), dqsfc(:)
    real(kind=RKIND), allocatable :: hpbl(:), dkt(:,:), dku(:,:), tkeh(:,:)
    real(kind=RKIND), allocatable :: def_1(:,:), def_2(:,:), def_3(:,:)
    real(kind=RKIND), allocatable :: dku3d_h(:,:), dku3d_e(:,:)
    real(kind=RKIND), allocatable :: ten_t(:,:), ten_u(:,:), ten_v(:,:)
    real(kind=RKIND), allocatable :: dv(:,:), du(:,:), tdt(:,:)
    real(kind=RKIND), allocatable :: dtend(:,:,:)
    real(kind=RKIND), allocatable :: claie(:), cfch(:), cfrt(:), cclu(:), cpopu(:)
    real(kind=RKIND) :: rho1, tem1, tem2
    integer, allocatable :: kpbl(:), kinver(:), dtidx(:,:)

    im = nCells
    km = nVertLevels

    errmsg = ''
    errflg = 0

    ! Tracer layout for the UFS routine.
    ntqv = 1
    ntcw = 2
    ntiw = 3
    ntrw = 0
    ntke = ntrac

    if (ntrac < 4) then
      errmsg = 'mpas_call_satmedmfvdifq: ntrac must be at least 4: qv,qc,qi,tke'
      errflg = 1
      return
    endif

    allocate(rtg(im,km,ntrac), q1(im,km,ntrac))
    allocate(u1(im,km), v1(im,km), t1(im,km))
    allocate(swh(im,km), hlw(im,km), xmu(im))
    allocate(garea(im), zvfun(im), sigmaf(im))
    allocate(psk(im), rbsoil(im), zorl(im))
    allocate(u10m(im), v10m(im), fm(im), fh(im))
    allocate(tsea(im), heat(im), evap(im), stress(im), spd1(im))
    allocate(prsi(im,km+1), del(im,km), prsl(im,km), prslk(im,km))
    allocate(phii(im,km+1), phil(im,km))
    allocate(dusfc(im), dvsfc(im), dtsfc(im), dqsfc(im))
    allocate(hpbl(im), dkt(im,km), dku(im,km), tkeh(im,km))
    allocate(def_1(im,km), def_2(im,km), def_3(im,km))
    allocate(dku3d_h(im,km), dku3d_e(im,km))
    allocate(ten_t(im,km), ten_u(im,km), ten_v(im,km))
    allocate(dv(im,km), du(im,km), tdt(im,km))
    allocate(kpbl(im), kinver(im))
    allocate(dtidx(3,1), dtend(im,km,3))
    allocate(claie(im), cfch(im), cfrt(im), cclu(im), cpopu(im))

    rtg = 0.0_RKIND
    q1 = 0.0_RKIND
    dkt = 0.0_RKIND
    dku = 0.0_RKIND
    tkeh = 0.0_RKIND
    def_1 = 0.0_RKIND
    def_2 = 0.0_RKIND
    def_3 = 0.0_RKIND
    dku3d_h = 0.0_RKIND
    dku3d_e = 0.0_RKIND
    ten_t = 0.0_RKIND
    ten_u = 0.0_RKIND
    ten_v = 0.0_RKIND
    dv = 0.0_RKIND
    du = 0.0_RKIND
    tdt = 0.0_RKIND
    dtend = 0.0_RKIND
    dtidx = 0
    claie = 0.0_RKIND
    cfch  = 0.0_RKIND
    cfrt  = 0.0_RKIND
    cclu  = 0.0_RKIND
    cpopu = 0.0_RKIND

    do k = 1, km
      do i = 1, im
        u1(i,k) = u_mpas(i,k)
        v1(i,k) = v_mpas(i,k)
        t1(i,k) = t_mpas(i,k)

        q1(i,k,ntqv) = max(qv_mpas(i,k), 1.0e-12_RKIND)
        q1(i,k,ntcw) = max(qc_mpas(i,k), 0.0_RKIND)
        q1(i,k,ntiw) = max(qi_mpas(i,k), 0.0_RKIND)
        q1(i,k,ntke) = max(tke_mpas(i,k), 1.0e-9_RKIND)

        prsl(i,k)  = p_mid(i,k)
        prslk(i,k) = exner_mid(i,k)

        ! UFS routine expects geopotential, not geometric height.
        phil(i,k) = grav * z_mid(i,k)

        swh(i,k) = sw_heat(i,k)
        hlw(i,k) = lw_heat(i,k)
      enddo
    enddo

    do k = 1, km+1
      do i = 1, im
        prsi(i,k) = p_int(i,k)

        ! UFS routine internally does zi=phii/grav, so pass geopotential.
        phii(i,k) = grav * z_int(i,k)
      enddo
    enddo

    do k = 1, km
      do i = 1, im
        del(i,k) = abs(prsi(i,k) - prsi(i,k+1))
      enddo
    enddo

    do i = 1, im
      garea(i) = areaCell(i)

      xmu(i) = max(coszen(i), 0.0_RKIND)

      ! UFS code uses z0 = 0.01*zorl, so zorl is in cm.
      zorl(i) = min(max(z0_mpas(i), 1.0e-4_RKIND), 1.0e-3_RKIND) * 100.0_RKIND

      rho1 = prsl(i,1) / (rd * max(t1(i,1), 180.0_RKIND))
      rho1 = min(max(rho1, 0.2_RKIND), 1.5_RKIND)

      tsea(i) = min(max(skin_temp(i), 180.0_RKIND), 340.0_RKIND)

      ! TKE-EDMF expects kinematic sensible/moisture fluxes.
      ! shflx: W m-2 -> K m s-1
      heat(i) = shflx(i)/(rho1*cp)
      ! qfx_mpas: kg m-2 s-1 -> kg kg-1 m s-1
      evap(i) = qfx_mpas(i)/rho1

      ! Temporary safety caps for debugging.
      heat(i) = min(max(heat(i), -1.0_RKIND), 1.0_RKIND)
      evap(i) = min(max(evap(i), -1.0e-3_RKIND), 1.0e-3_RKIND)

      stress(i) = min(max(stress_in(i), 1.0e-6_RKIND), 5.0_RKIND)
      spd1(i) = max(sqrt(u1(i,1)**2 + v1(i,1)**2), 1.0_RKIND)

      u10m(i) = u10(i)
      v10m(i) = v10(i)
      rbsoil(i) = min(max(rb_in(i), -10.0_RKIND), 10.0_RKIND)
      fm(i)     = min(max(fm_in(i), 1.0e-6_RKIND), 10.0_RKIND)
      fh(i)     = min(max(fh_in(i), 1.0e-6_RKIND), 10.0_RKIND)

      ! If MPAS does not have rbsoil, start neutral.
!     rbsoil(i) = 0.0_RKIND

      ! psk is surface Exner. Use lowest layer as fallback.
      psk(i) = prslk(i,1)

      ! Need real mappings later from vegetation/roughness data.
      sigmaf(i) = min(max(vegfra_in(i), 0.0_RKIND), 1.0_RKIND)
!     zvfun(i) = 1.0_RKIND
      !-------------------------------------------------------
! Compute zvfun: function of surface roughness and vegetation
!-------------------------------------------------------

! z0 from MPAS is in meters; UFS uses zorl in cm
! Here we stay consistent with MPAS units (meters)
      tem1 = (z0_mpas(i) - z0lo) / (z0up - z0lo)
      tem1 = min(max(tem1, 0.0_RKIND), 1.0_RKIND)
      tem2 = max(sigmaf(i), 0.1_RKIND)
      zvfun(i) = sqrt(tem1 * tem2)

      ! No inversion limiter initially.
      kinver(i) = km

      kpbl(i) = 1
      hpbl(i) = 0.0_RKIND
    enddo

    index_of_temperature = 1
    index_of_x_wind = 2
    index_of_y_wind = 3
    index_of_process_pbl = 1


    !write(0,*) 'TKE-EDMF input min/max: rb=', minval(rbsoil), maxval(rbsoil), &
    !    ' fm=', minval(fm), maxval(fm), ' fh=', minval(fh), maxval(fh)
    !write(0,*) 'TKE-EDMF input min/max: zorl=', minval(zorl), maxval(zorl), &
    !    ' stress=', minval(stress), maxval(stress), ' spd1=', minval(spd1), maxval(spd1)
    !write(0,*) 'TKE-EDMF input min/max: heat=', minval(heat), maxval(heat), &
    !    ' evap=', minval(evap), maxval(evap), ' sigmaf=', minval(sigmaf), maxval(sigmaf)
    !write(0,*) 'TKE-EDMF input min/max: p=', minval(prsl), maxval(prsl), &
    !    ' t=', minval(t1), maxval(t1), ' qv=', minval(q1(:,:,ntqv)), maxval(q1(:,:,ntqv))
   !call flush(0)


    ! Wrapper-side vertical orientation check.
    ! If k=1 is lowest layer, z_mid(:,1) should be small and prsl(:,1) should be large.
    !write(0,*) 'WRAPPER VERT CHECK: z_mid k1/km=', minval(z_mid(:,1)), maxval(z_mid(:,1)), &
    !    minval(z_mid(:,km)), maxval(z_mid(:,km))
    !write(0,*) 'WRAPPER VERT CHECK: prsl k1/km=', minval(prsl(:,1)), maxval(prsl(:,1)), &
    !    minval(prsl(:,km)), maxval(prsl(:,km))
   !call flush(0)

    if (minval(prsl) < 0.1_RKIND) then
   !   !write(0,*) 'INFO full-column pressure includes upper-atmosphere levels: prsl=', minval(prsl), maxval(prsl), &
   !              ' t1=', minval(t1), maxval(t1)
   !   call flush(0)
    endif

    ! Lower-column sanity check. Avoid flagging the whole model top;
    ! TKE-EDMF should mainly depend on the lower PBL-relevant levels.
    kcheck = min(km, 30)
    ibad = 0
    kbad = 0

    do k = 1, kcheck
      do i = 1, im
        if (prsl(i,k) /= prsl(i,k) .or. t1(i,k) /= t1(i,k) .or. &
            exner_mid(i,k) /= exner_mid(i,k) .or. q1(i,k,ntqv) /= q1(i,k,ntqv)) then
           ibad = i
           kbad = k
        endif

        if (prsl(i,k) < 1000.0_RKIND .or. exner_mid(i,k) <= 1.0e-6_RKIND .or. &
            t1(i,k) < 120.0_RKIND .or. t1(i,k) > 360.0_RKIND .or. &
            q1(i,k,ntqv) < 0.0_RKIND) then
           ibad = i
           kbad = k
        endif
      enddo
    enddo

    if (ibad > 0) then
       !write(0,*) 'BAD lower-column TKE-EDMF input before satmedmfvdifq_run'
       !write(0,*) '  local i,k=', ibad, kbad
       !write(0,*) '  prsl=', prsl(ibad,kbad), ' t1=', t1(ibad,kbad), &
       !          ' qv=', q1(ibad,kbad,ntqv), ' exner=', exner_mid(ibad,kbad)
       !write(0,*) '  lower-column ranges: prsl=', minval(prsl(:,1:kcheck)), maxval(prsl(:,1:kcheck)), &
       !          ' t1=', minval(t1(:,1:kcheck)), maxval(t1(:,1:kcheck)), &
       !          ' qv=', minval(q1(:,1:kcheck,ntqv)), maxval(q1(:,1:kcheck,ntqv)), &
       !          ' exner=', minval(exner_mid(:,1:kcheck)), maxval(exner_mid(:,1:kcheck))
      !call flush(0)

       ten_t_out(:,:)  = 0.0_RKIND
       ten_u_out(:,:)  = 0.0_RKIND
       ten_v_out(:,:)  = 0.0_RKIND
       ten_qv_out(:,:) = 0.0_RKIND
       ten_qc_out(:,:) = 0.0_RKIND
       ten_qi_out(:,:) = 0.0_RKIND
       hpbl_out(:) = 0.0_RKIND
       kpbl_out(:) = 1
       errflg = 0
       errmsg = ''
       return
    endif

    if (any(rbsoil /= rbsoil) .or. any(fm /= fm) .or. any(fh /= fh) .or. &
        any(zorl /= zorl) .or. any(heat /= heat) .or. any(evap /= evap)) then
       !write(0,*) 'BAD TKE-EDMF input: NaN detected before satmedmfvdifq_run'
      !call flush(0)
       ten_t_out(:,:)  = 0.0_RKIND
       ten_u_out(:,:)  = 0.0_RKIND
       ten_v_out(:,:)  = 0.0_RKIND
       ten_qv_out(:,:) = 0.0_RKIND
       ten_qc_out(:,:) = 0.0_RKIND
       ten_qi_out(:,:) = 0.0_RKIND
       hpbl_out(:) = 0.0_RKIND
       kpbl_out(:) = 1
       errflg = 0
       errmsg = ''
       return
    endif


    ! Unit/range checks for roughness and surface flux conversion.
    if (any(z0_mpas < 0.0_RKIND) .or. any(z0_mpas > 1.0_RKIND)) then
       !write(0,*) 'WARNING z0_mpas suspicious [m]: ', minval(z0_mpas), maxval(z0_mpas)
      !call flush(0)
    endif

    if (any(abs(shflx) > 2000.0_RKIND)) then
       !write(0,*) 'WARNING shflx suspicious [W/m2]: ', minval(shflx), maxval(shflx)
      !call flush(0)
    endif

    if (any(abs(qfx_mpas) > 1.0e-3_RKIND)) then
       !write(0,*) 'WARNING qfx suspicious [kg/m2/s]: ', minval(qfx_mpas), maxval(qfx_mpas)
      !call flush(0)
    endif

    !write(0,*) 'UNIT CHECK TKE wrapper: z0_m=', minval(z0_mpas), maxval(z0_mpas), &
      !  ' zorl_cm=', minval(zorl), maxval(zorl), &
      !  ' hfx_Wm2=', minval(shflx), maxval(shflx), &
      !  ' heat_Kms=', minval(heat), maxval(heat), &
      !  ' qfx_kgm2s=', minval(qfx_mpas), maxval(qfx_mpas), &
      !  ' evap_ms=', minval(evap), maxval(evap)
   !call flush(0)

    call satmedmfvdifq_run(im, km, ntrac, ntcw, ntrw, ntiw, ntke,       &
         grav, pi, rd, cp, rv, hvap, hfus, fv, eps, epsm1,             &
         def_1, def_2, def_3, cfg%sa3dtke, dku3d_h, dku3d_e,           &
         dv, du, tdt, rtg, u1, v1, t1, q1,                             &
         swh, hlw, xmu, garea, zvfun, sigmaf,                          &
         psk, rbsoil, zorl, u10m, v10m, fm, fh,                        &
         tsea, heat, evap, stress, spd1, kpbl,                         &
         prsi, del, prsl, prslk, phii, phil, dt, cfg%tte_edmf,         &
         cfg%dspheat, dusfc, dvsfc, dtsfc, dqsfc, hpbl, dkt, dku, tkeh,&
         kinver, cfg%xkzm_m, cfg%xkzm_h, cfg%xkzm_s, cfg%dspfac,       &
         cfg%bl_upfr, cfg%bl_dnfr, cfg%rlmx, cfg%elmx,                 &
         cfg%sfc_rlm, cfg%tc_pbl, cfg%use_lpt,                         &
         cfg%do_canopy, cfg%cplaqm, claie, cfch, cfrt, cclu, cpopu,    &
         ntqv, dtend, dtidx, index_of_temperature,                     &
         index_of_x_wind, index_of_y_wind, index_of_process_pbl,       &
         cfg%gen_tend, cfg%ldiag3d, errmsg, errflg)

    if (errflg /= 0) then
       !write(0,*) 'TKE-EDMF returned errflg=', errflg, ' errmsg=', trim(errmsg)
      !call flush(0)
       return
    endif

    !write(0,*) 'TKE-EDMF raw output min/max: tdt=', minval(tdt), maxval(tdt), &
     !   ' du=', minval(du), maxval(du), ' dv=', minval(dv), maxval(dv)
    !write(0,*) 'TKE-EDMF raw output min/max: rtg_qv=', minval(rtg(:,:,ntqv)), maxval(rtg(:,:,ntqv)), &
     !   ' rtg_qc=', minval(rtg(:,:,ntcw)), maxval(rtg(:,:,ntcw)), &
     !   ' rtg_qi=', minval(rtg(:,:,ntiw)), maxval(rtg(:,:,ntiw))
    !write(0,*) 'TKE-EDMF raw output min/max: hpbl=', minval(hpbl), maxval(hpbl), &
     !   ' kpbl=', minval(kpbl), maxval(kpbl), ' tke=', minval(q1(:,:,ntke)), maxval(q1(:,:,ntke))
   !call flush(0)

    do k = 1, km
      do i = 1, im
        if (exner_mid(i,k) > 1.0e-8_RKIND) then
          ten_t_out(i,k) = tdt(i,k) / exner_mid(i,k)
        else
          ten_t_out(i,k) = 0.0_RKIND
        endif

        ten_u_out(i,k) = du(i,k)
        ten_v_out(i,k) = dv(i,k)

        ten_qv_out(i,k) = rtg(i,k,ntqv)
        ten_qc_out(i,k) = rtg(i,k,ntcw)
        ten_qi_out(i,k) = rtg(i,k,ntiw)

        ! Temporary safety caps for debugging coupled gfs_sfcl + TKE-EDMF.
        ten_t_out(i,k)  = min(max(ten_t_out(i,k),  -2.0e-4_RKIND), 2.0e-4_RKIND)
        ten_qv_out(i,k) = min(max(ten_qv_out(i,k), -1.0e-6_RKIND), 1.0e-6_RKIND)
        ten_qc_out(i,k) = min(max(ten_qc_out(i,k), -1.0e-7_RKIND), 1.0e-7_RKIND)
        ten_qi_out(i,k) = min(max(ten_qi_out(i,k), -1.0e-8_RKIND), 1.0e-8_RKIND)

        tke_mpas(i,k) = min(max(q1(i,k,ntke), 1.0e-9_RKIND), 10.0_RKIND)

      enddo
    enddo

    !write(0,*) 'TKE-EDMF final tendencies: rth=', minval(ten_t_out), maxval(ten_t_out), &
     !   ' ru=', minval(ten_u_out), maxval(ten_u_out), ' rv=', minval(ten_v_out), maxval(ten_v_out)
    !write(0,*) 'TKE-EDMF final q/tke: qv=', minval(ten_qv_out), maxval(ten_qv_out), &
     !   ' qc=', minval(ten_qc_out), maxval(ten_qc_out), ' qi=', minval(ten_qi_out), maxval(ten_qi_out), &
     !   ' tke=', minval(tke_mpas), maxval(tke_mpas)
  ! call flush(0)

    if (any(ten_t_out /= ten_t_out) .or. any(ten_u_out /= ten_u_out) .or. &
        any(ten_v_out /= ten_v_out) .or. any(ten_qv_out /= ten_qv_out) .or. &
        any(tke_mpas /= tke_mpas)) then
       !write(0,*) 'BAD TKE-EDMF output: NaN detected after satmedmfvdifq_run'
  !    call flush(0)
       errflg = 1
       errmsg = 'NaN in TKE-EDMF wrapper output'
       return
    endif

    do i = 1, im
      hpbl_out(i) = hpbl(i)
      kpbl_out(i) = kpbl(i)
    enddo

  end subroutine mpas_call_satmedmfvdifq

end module mpas_satmedmfvdifq_wrapper_mod
