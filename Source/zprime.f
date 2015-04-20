program zprime

  ! calculates the cross section and generates distributions for
  ! pp -> tt,
  ! pp -> tt -> bw^+bbarw^- -> bbbare^+nue^-nubarc
  ! (future:) pp -> bw^+bbarw^- -> b bbar e^+ nu qqbar'
  ! uses adapted madgraph functions.
  ! uses cteq6 and mrs99 pdf subroutines.
  ! authors: declan millar, stefano moretti.

  use mathematics, only: pi
  use configuration
  use modelling
  use kinematics
  use integration
  use distributions  

  implicit none

  external differential_cross_section

  real :: poltot(-1:1, -1:1), polchi(-1:1, -1:1)
  real :: spattot(n_fb_asymmetries, -1:1), spatchi(n_fb_asymmetries, -1:1)

  real :: avgi, chi2a, cross, error, sd, stantot, alfas, qcdl4
  integer :: iab, ndimensions, iasy, icteq, iphel,ispat,jphel

  ! branching ratio for t->bev=bmv=btv (with qcd corrections?)
  real :: brtbln = 0.10779733d0

  ! branching ratio for t->bev=bmv=btv=1/9 (tree level)
  ! real brtbln/0.11111111d0/

  ! branching ratio for t->beq=bmq=btq=6/9 (tree level)
  ! real :: brtbeq=0.66666666d0

  integer :: i, j, k
  integer :: today(3), now(3)
  double precision :: start_time, finish_time

  call cpu_time(start_time)

  call read_config

  call modify_config

  open(unit = 10, file = 'Output/'//output_file, status = "replace", action = "write")
  close(10)

  call create_distributions
  call initialise_distributions
  call disable_distributions

  s = collider_energy*collider_energy

  ! (pdfs are intrinsically linked to the value of lamda_qcd; alpha_qcd)
  if (structure_function == 1) qcdl4 = 0.326d0
  if (structure_function == 2) qcdl4 = 0.326d0
  if (structure_function == 3) qcdl4 = 0.326d0
  if (structure_function == 4) qcdl4 = 0.215d0
  if (structure_function == 5) qcdl4 = 0.300d0
  if (structure_function == 6) qcdl4 = 0.300d0
  if (structure_function == 7) qcdl4 = 0.300d0
  if (structure_function == 8) qcdl4 = 0.229d0
  if (structure_function == 9) qcdl4 = 0.383d0

  ! qcdl4 is qcd lambda4 (to match pdf fits).
  rlambdaqcd4 = qcdl4

  ! initialise cteq grids.
  if (structure_function <= 4) then
    icteq = structure_function
    call setctq6(icteq)
  end if

  ! use appropriately evolved alphas.
  if (structure_function == 1) nloops = 2
  if (structure_function == 2) nloops = 2
  if (structure_function == 3) nloops = 1
  if (structure_function == 4) nloops = 1
  if (structure_function == 5) nloops = 1
  if (structure_function == 6) nloops = 1
  if (structure_function == 7) nloops = 1
  if (structure_function == 8) nloops = 1
  if (structure_function == 9) nloops = 1

  ! initialise madgraph - masses and coupling constants of particles
  call initialise_standard_model
  call initialise_zprimes

  ! dimensions
  if (final_state == 0) then
    ndimensions = 3
  else if (final_state > 0) then
    ndimensions = 15
  end if

  ! if nprn<0 no print-out
  nprn = 0

  if (final_state == 0) then
    m3 = tmass
    m4 = tmass
    m5 = 0.d0
    m6 = 0.d0
    m7 = 0.d0
    m8 = 0.d0

    ! integrates on:
    ! x(3) = (x1 - tau)/(1 - tau),
    ! x(2) = (ecm - m3 - m4)/(ecm_max - m3 - m4),
    ! x(1) = cos(theta3_cm)

    ! limits:
    do i = 3, 2, -1
      xl(i) = 0.d0
      xu(i) = 1.d0
    end do
    do i = 1, 1
      xl(i) = -1.d0
      xu(i) = 1.d0
    end do

  else if (final_state >= 1) then
    m3 = bmass
    m4 = bmass
    m5 = emass
    m6 = nuemass
    m7 = emass
    m8 = nuemass

    ! integrates on:         
    ! x(15) = (x1 - tau)/(1 - tau),
    ! x(14) = (ecm - m3 - m4 - m5 - m6 - m7 - m8)
    !        /(ecm_max - m3 - m4 - m5 - m6 - m7 - m8),
    ! x(13) = (xx356 - xx356min)/(xx356max - xx356min),
    ! where xx356 = arctg((m356**2 - m3**2)/m3/gamt),
    ! x(12) = (xx478 - xx478min)/(xx478max - xx478min),
    ! where xx478 = arctg((m478**2 - m3**2)/m3/gamt),
    ! x(11) = (xx56 - xx56min)/(xx56max - xx56min),
    ! where xx56 = arctg((m56**2 - rm_w**2)/rm_w/gamw),
    ! x(10) = (xx78 - xx78min)/(xx78max - xx78min),
    ! where xx78 = arctg((m78**2 - rm_w**2)/rm_w/gamw),
    ! x(9) = cos(theta_cm_356) = -cos(theta_cm_478)
    ! x(8) = cos(theta56_cm_356),
    ! x(7) = cos(theta78_cm_478),
    ! x(6) = cos(theta5_cm_56),
    ! x(5) = cos(theta7_cm_78),
    ! x(4) = fi56_cm_356,
    ! x(3) = fi78_cm_478,
    ! x(2) = fi5_cm_56,
    ! x(1) = fi8_cm_78;

    ! limits:
    do i = 15, 14, -1
      xl(i) = 0.d0
      xu(i) = 1.d0
    end do
    do i = 13, 10, -1
      xl(i) = 0.d0
      xu(i) = 1.d0
    end do
    do i = 9, 5, -1
      xl(i) = -1.d0
      xu(i) = 1.d0
    end do
    do i = 4, 1, -1
      xl(i) = 0.d0
      xu(i) = 2.d0*pi
    end do
  end if

  call generate_bins

  ! output information before integration
  open(unit = 10, file = 'Output/'//output_file, status = "old", action = "write", position="append")
  write(10,*) '====================================================='
  call idate(today)     ! today(1) = day, (2) = month, (3) = year
  call itime(now)       ! now(1) = hour, (2) = minute, (3) = second
  write(10,*) 'date ', today(3), today(2), today(1)
  write(10,*) 'time ', now(1), now(2), now(3)
  write(10,*) '-----------------------------------------------------'
  write(10,*) 'process'
  if (initial_state == 0) then
    if (final_state == 0) then
      write(10,*) 'pp #rightarrow t#bar{t}', &
               ' #times br(t#rightarrow bl#nu)^{2}'
    else if (final_state == 1) then
      write(10,*) 'pp #rightarrow t#bar{t}', &
               '#rightarrow b#bar{b} w^{+}w^{-}', &
               '#rightarrow b#bar{b} l^{+}l^{-} #nu#bar{#nu}'
    else if (final_state == 2) then
      write(10,*) 'pp #rightarrow t#bar{t}', &
               '#rightarrow b#bar{b} w^{+}w^{-}', &
               '#rightarrow b#bar{b} q#bar{q} l #nu'
    else if (final_state == 3) then
      write(10,*) 'pp #rightarrow t#bar{t}', &
               '#rightarrow b#bar{b} w^{+}w^{-}', &
               "#rightarrow b#bar{b} q#bar{q}q'#bar{q}'"
    end if
  else if (initial_state == 1) then
    if (final_state == 0) then
      write(10,*) 'p#bar{p} #rightarrow t#bar{t}', &
               ' #times br(t#rightarrow bl#nu)^{2}'
    else if (final_state == 1) then
      write(10,*) 'p#bar{p} #rightarrow t#bar{t}', &
               '#rightarrow b#bar{b} w^{+}w^{-}', &
               '#rightarrow b#bar{b} l^{+}l^{-} #nu#bar{#nu}'
    else if (final_state == 2) then
      write(10,*) 'p#bar{p} #rightarrow t#bar{t}', &
               '#rightarrow b#bar{b} w^{+}w^{-}', &
               '#rightarrow b#bar{b} q#bar{q} l #nu'
    else if (final_state == 3) then
      write(10,*) 'p#bar{p} #rightarrow t#bar{t}', &
               '#rightarrow b#bar{b} w^{+}w^{-}', &
               "#rightarrow b#bar{b} q#bar{q}q'#bar{q}'"
    end if
  end if
  write(10,*) '-----------------------------------------------------'
  write(10,*) 'notes'
  write(10,*) 'units: gev'
  write(10,*) 'quarks: all massless except t, b.'
  if (structure_function == 1) write(10,*) 'pdfs: cteq6m.'
  if (structure_function == 2) write(10,*) 'pdfs: cteq6d.'
  if (structure_function == 3) write(10,*) 'pdfs: cteq6l.'
  if (structure_function == 4) write(10,*) 'pdfs: cteq6l1.'
  if (structure_function == 5) write(10,*) 'pdfs: mrs99 (cor01).'
  if (structure_function == 6) write(10,*) 'pdfs: mrs99 (cor02).'
  if (structure_function == 7) write(10,*) 'pdfs: mrs99 (cor03).'
  if (structure_function == 8) write(10,*) 'pdfs: mrs99 (cor04).'
  if (structure_function == 9) write(10,*) 'pdfs: mrs99 (cor05).'
  if ((final_state >= 1) .and. (use_nwa == 0)) write(10,*) 'tops: off-shell.'
  if ((final_state >= 1) .and. (use_nwa == 1)) write(10,*) 'tops: nwa.'
  write(10,*) 'bsm model_name: ', model_name
  if (include_qcd == 1) write(10,*) 'qcd: on '
  if (include_qcd == 0) write(10,*) 'qcd: off'
  if (include_ew == 1) write(10,*) 'ew:  on '
  if (include_ew == 0) write(10,*) 'ew:  off'
  if (include_bsm == 1) write(10,*) 'bsm: on '
  if (include_bsm == 0) write(10,*) 'bsm: off'
  if (interference == 0) write(10,*) 'interference: none'
  if (interference == 1) write(10,*) 'interference: sm'
  if (interference == 2) write(10,*) 'interference: full'
  if (interference == 3) write(10,*) 'interference: no square terms.'
  if (phase_space_only == 1) write(10,*) 'phase space only'
  if (symmetrise_x1x2 == 1) write(10,*) 'symmetrical integration over x1<->x2'
  if (symmetrise_costheta_t == 1) write(10,*) 'symmetrical integration over cost'
  if (include_errors == 1) write(10,*) 'Distribution errors included.'
  if (include_errors == 0) write(10,*) 'Distribution errors excluded.'
  write(10,*) 'seed: ', seed
  write(10,*) '-----------------------------------------------------'
  write(10,*) 'parameters'
  write(10,*) '#sqrt{s}              ', collider_energy
  write(10,*) 'at |y| <              ', abs(ytmax)
  write(10,*) 'loops a_s evaluated at', nloops
  write(10,*) 'a_{s}(m_{z})          ', alfas(rm_z, rlambdaqcd4, nloops)
  write(10,*) '#lambda_{qcd}(4)      ', qcdl4
  write(10,*) 'm_{b}                 ', fmass(12)
  write(10,*) '#gamma_{b}            ', fwidth(12)
  write(10,*) 'm_{t}                 ', fmass(11)
  write(10,*) '#gamma_{t}            ', fwidth(11)
  write(10,*) 'm_{z}                 ', rm_z
  write(10,*) '#gamma_{z}            ', gamma_z
  write(10,*) 'm_{w}                 ', rm_w
  write(10,*) '#gamma_{w}            ', gamma_w
  write(10,*) 'm_{h}                 ', rm_h
  write(10,*) '#gamma_{h}            ', gamma_h
  write(10,*) '-----------------------------------------------------'
  write(10,*) 'zprime parameters'
  do i = 1, 5
    if (rmzp(i) > 0) then
      write(10,*) 'z#prime               ', i
      write(10,*) 'm_{z#prime}           ', rmzp(i)
      write(10,*) '#gamma_{z#prime}      ', gamzp(i)
      write(10,*) 'g_{p}                 ', gp(i)
      write(10,*) 'g_{v}^{u}             ', gv_u(i)
      write(10,*) 'g_{a}^{u}             ', ga_u(i)
      write(10,*) 'g_{v}^{d}             ', gv_d(i)
      write(10,*) 'g_{a}^{d}             ', ga_d(i)
      write(10,*)
    end if
  end do
  write(10,*) '-----------------------------------------------------'
  write(10,*) 'cuts'
  write(10,*) '-----------------------------------------------------'
  close(10)

  ! reset counter
  npoints = 0

  ! reset 
  if (final_state == 0) then
    do i = 1, 20
      resl(i) = 0.d0
      standdevl(i) = 0.d0
      cnorm(i) = 0.d0
      do iphel = -1, +1, 2
        do jphel = -1, +1, 2
          xsec_polar(i, iphel, jphel) = 0.d0
          error_polar(i, iphel, jphel) = 0.d0
        end do
      end do
      do ispat = 1, n_fb_asymmetries
        do iasy = -1, +1, 2
          xsec_fb(ispat, i, iasy) = 0.d0
          error_fb(ispat, i, iasy) = 0.d0
        end do
      end do
    end do
  end if

  ! integrate 
  print *, 'Starting integration'
  it = 0 
  call vegas(ndimensions, differential_cross_section, avgi, sd, chi2a)

  print *, 'Finished integration'  

  if (final_state == 0 .and. use_branching_ratio == 1) then
    ! multiply by branching ratios
    avgi = avgi*brtbln*brtbln
    sd = sd*brtbln*brtbln
  end if

  if (final_state == 2) then
    ! multiply by dilepton to semi-hadronic conversion
    avgi = avgi*12
    sd = sd*12
  end if

  if (final_state == 3) then
    ! multiply by dilepton to fully hadronic conversion
    avgi = avgi*36
    sd = sd*36
  end if

  ! collect total cross-section
  cross = avgi
  sigma = avgi
  error = sd

  ! print integrated cross section
  open(unit = 10, file = 'Output/'//output_file, status = "old", action = "write", position="append")
  write(10,*) 'integrated cross section'
  if (cross == 0.d0) then
    write(10,*) 'sigma = 0  ! check permitted gauge sectors.'
    stop
  else
    write(10,*) 'sigma (pb)', 'error (same units)'
    write(10,*) cross, error
    write(10,*) '(using ', npoints, ' points)'
  end if

  ! re-weight distributions for different iterations
  stantot = 0.d0
  do i = 1, it
    stantot = stantot + 1.d0/standdevl(i)/standdevl(i)
  end do
  do i = 1, it
    standdevl(i) = standdevl(i)*standdevl(i)*stantot
  end do
  do i = 1, it
    cnorm(i) = resl(i)*standdevl(i)
  end do

  ! total asymmetries
  ! collect polarised cross sections.
  if (include_asymmetries == 1) then
    if (final_state == 0) then
      do iphel = -1, +1, 2
        do jphel = -1, +1, 2
          do i = 1, it
            xsec_polar(i, iphel, jphel) = xsec_polar(i, iphel, jphel) &
                                          *avgi/cnorm(i)
            error_polar(i, iphel, jphel) = xsec_polar(i, iphel, jphel) &
                                           *sd/cnorm(i)
          end do
          poltot(iphel, jphel) = 0.d0
          polchi(iphel, jphel) = 0.d0
          do i = 1, it
            poltot(iphel, jphel) = poltot(iphel, jphel) &
                                   + xsec_polar(i, iphel, jphel)
            polchi(iphel, jphel) = polchi(iphel, jphel) &
                                   + error_polar(i, iphel, jphel)
          end do
          polchi(iphel, jphel) = polchi(iphel, jphel) &
          /poltot(iphel, jphel)
          !        polchi(iphel,jphel)=
          ! & sqrt(abs(polchi(iphel,jphel)
          ! &         -poltot(iphel,jphel)**2*dfloat(ncall)))
          ! & /dfloat(ncall)
        end do
      end do
    end if

    ! collect unpolarised spatial asymmetry
!     do ispat = 1, n_fb_asymmetries
!       if (o_asym(ispat + 3) == 0) then
!         continue
!       else
!         do iab = -1,+1, 2
!           do i = 1, it
!             xsec_fb(ispat, i, iab) = xsec_fb(ispat, i, iab) &
!                                      *avgi/cnorm(i)
!             error_fb(ispat, i, iab) = xsec_fb(ispat, i, iab) &
!                                       *sd/cnorm(i)
!           end do
!           spattot(ispat, iab) = 0.d0
!           spatchi(ispat, iab) = 0.d0
!           do i = 1, it
!             spattot(ispat, iab) = spattot(ispat, iab) &
!                                   + xsec_fb(ispat, i, iab)
!             spatchi(ispat, iab) = spatchi(ispat, iab) &
!                                   + error_fb(ispat, i, iab)
!           end do
!           spatchi(ispat, iab) = spatchi(ispat, iab) &
!                                 /spattot(ispat, iab)
!           !         spatchi(iasy)=
!           !    & sqrt(abs(spatchi(iasy)
!           !    &         -spattot(iasy)**2*dfloat(ncall)))
!           !    & /dfloat(ncall)
!         end do
!       end if
!     end do

    ! define asymmetries
    if (final_state == 0) then
      ! all
      atot(1) = (poltot(+1, +1) - poltot(+1, -1) &
                 - poltot(-1, +1) + poltot(-1, -1)) &
                /cross
      atoterr(1) = (polchi(+1, +1) + polchi(+1, -1) &
                    +polchi(-1, +1) + polchi(-1, -1)) &
                   /4.d0*atot(1)
      ! al
      atot(2) = (poltot(-1, -1) - poltot(+1, -1) &
                 + poltot(-1, +1) - poltot(+1, +1)) &
                /cross
      atoterr(2) = (polchi(-1, -1) + polchi(+1, -1) &
                    +polchi(-1, +1) + polchi(+1, +1)) &
                   /4.d0*atot(2)
      ! apv
      atot(3) = (poltot(-1, -1) - poltot(+1, +1)) &
                /cross/2.d0
      atoterr(3) = (polchi(-1, -1) + polchi(+1, +1)) &
                   /2.d0*atot(3)
    end if

!     do iasy = 4, n_asymmetries
!       ispat = iasy - 3
!       if (o_asym(iasy) > 0) then
!         atot(iasy) = (spattot(ispat, +1) - spattot(ispat, -1))/cross
!         atoterr(iasy) = sd/avgi*atot(iasy)
!       end if
!     end do

    ! print asymmetries
    write(10,*) 'total asymmetries'
!     if (final_state == 0) then
!       write(10,*) 'ALL:                    uncertainty (same units):'
!       write(10,*) atot(1), atoterr(1)
!       write(10,*) 'AL:                     uncertainty (same units):'
!       write(10,*) atot(2), atoterr(2)
!       write(10,*) 'APV:                    uncertainty (same units):'
!       write(10,*) atot(3), atoterr(3)
!       write(10,*) 'AFB:                    uncertainty (same units):'
!       write(10,*) atot(4), atoterr(4)
!       write(10,*) 'AFB*:                   uncertainty (same units):'
!       write(10,*) atot(5), atoterr(5)
!       write(10,*) 'AtRFB:                  uncertainty (same units):'
!       write(10,*) atot(6), atoterr(6)
!       write(10,*) "AttbRFB:                uncertainty (same units):"
!       write(10,*) atot(7), atoterr(7)
!       write(10,*) "ARFB/:                  uncertainty (same units):"
!       write(10,*) atot(8), atoterr(8)
!     else if (final_state > 0) then
!       write(10,*) 'A_l:                    uncertainty (same units):'
!       write(10,*) atot(9), atoterr(9)
!     end if
  end if

  call finalise_distributions

  close(10)
  call cpu_time(finish_time)
  print*, 'Program complete'
  print '("Time = ",f6.3," seconds.")', finish_time-start_time
  stop
end program zprime
