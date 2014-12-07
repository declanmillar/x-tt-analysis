! ======================================================================
      program ppbbllnn_BSM
! ----------------------------------------------------------------------
! Header
  ! Authors: Declan Millar, Stefano Moretti.

  ! Calculates the cross section and generates distributions for
  !   pp -> tt,
  !   pp -> tt -> bW^+bbarW^- -> bbbare^+nue^-nubarc
  ! (Future:) pp -> bW^+bbarW^- -> b bbar e^+ nu qqbar'
  ! Uses adapted Madgraph functions.
  ! Uses cteq6 and mrs99 PDF subroutines.
! ----------------------------------------------------------------------
! Declarations
  ! implicit
      implicit real*8 (a-h,o-z)

  !  LHE format
    ! Initialization information
      integer maxpup
      parameter (maxpup=100)
      integer idbmup,pdfgup,pdfsup,idwtup,nprup,lprup 
      double precision ebmup,xsecup,xerrup,xmaxup 
      common/heprup/idbmup(2),ebmup(2),pdfgup(2),pdfsup(2),
     & idwtup,nprup,xsecup(maxpup),xerrup(maxpup),
     & xmaxup(maxpup),lprup(maxpup)
    ! Information on each separate event
      integer maxnup
      parameter (maxnup=500)
      integer nup,idprup,idup,istup,mothup,icolup
      double precision xwgtup,scalup,aqedup,aqcdup,pup,vtimup,spinup
      common/hepeup/nup,idprup,xwgtup,scalup,aqedup,aqcdup, 
     & idup(maxnup),istup(maxnup),mothup(2,maxnup),icolup(2,maxnup),
     & pup(5,maxnup),vtimup(maxnup),spinup(maxnup)

  ! Global variables     
  !   vegas
      common/bveg1/ncall,itmx,nprn,ndev,xl(100),xu(100),acc
      common/bveg2/it,ndo,si,swgt,schi,xi(50,100)
      common/rndm/iseed
      common/reslocal/resl(20),standdevl(20)
  !   Kinematics
      common/par/rm3,rm4,rm5,rm6,rm7,rm8,s
      common/limfac/fac
      common/EW/a_em,s2w
      common/final/ifinal,ipmax
      common/top/rmt,gamt
      common/W/rmW,gamW
      common/Z/rmZ,gamZ
      common/H/rmH,gamH      
      common/stat/npoints
      common/coll/ecm_coll
      common/cuts/ytcut,yttcut
  !   Polarised/Spatial cross sections
      common/polarised/polcross(20,-1:1,-1:1),polerror(20,-1:1,-1:1)
      common/spatial/asycross(5,20,-1:1),asyerror(5,20,-1:1)  !nasy -3
  !   Permitted gauge sectors
      common/igauge/iQCD,iEW,iBSM
  !   Interference       
      common/interference/iint      
  !   Z' masses and VA/LR couplings
      common/Zp/rmZp(5),gamZp(5)
      common/Zpparam/paramZp(5)
      common/ZpAVcoup/gp(5),gV_d(5),gA_d(5),gV_u(5),gA_u(5)
      common/ZpLRcoup/gZpd(2,5),gZpu(2,5)
  !   Narrow width approximation (NWA)
      common/NWA/gNWA
  !   Structure functions
      common/partdist/istructure
      common/ALFASTRONG/rlambdaQCD4,nloops
      common/collider/icoll
  !   Switch for all distributions
      common/distros/idist
  !   Distributions in pTs of external particles
      common/ext_pT/pTmax(8),pTmin(8),pTw(8)
      common/dist_pT/xpT(8,500),fxpT(8,500,20),fxpTtot(8,500)
      common/inp_pT/m_pT(8)
      common/div_pT/ndiv_pT(8)
  !   Distributions in etas of external particles
      common/ext_eta/etamax(8),etamin(8),etaw(8)
      common/dist_eta/xeta(8,500),fxeta(8,500,20),fxetatot(8,500)
      common/inp_eta/m_eta(8)
      common/div_eta/ndiv_eta(8)
  !   Distributions in phis of external particles
      common/ext_phi/phimax(8),phimin(8),phiw(8)
      common/dist_phi/xphi(8,500),fxphi(8,500,20),fxphitot(8,500)
      common/inp_phi/m_phi(8)
      common/div_phi/ndiv_phi(8)
  !   Distribution in pT of the top
      common/ext_pT356/pT356max,pT356min,pT356w
      common/dist_pT356/xpT356(500),fxpT356(500,20),fxpT356tot(500)
      common/inp_pT356/m_pT356
      common/div_pT356/ndiv_pT356
  !   Distribution in pT of the anti-top
      common/ext_pT478/pT478max,pT478min,pT478w
      common/dist_pT478/xpT478(500),fxpT478(500,20),fxpT478tot(500)
      common/inp_pT478/m_pT478
      common/div_pT478/ndiv_pT478
  !   Distribution in ETmiss
      common/ext_ETmiss/ETmissmax,ETmissmin,ETmissw
      common/dist_ETmiss/xETmiss(500),fxETmiss(500,20),fxETmisstot(500)
      common/inp_ETmiss/m_ETmiss
      common/div_ETmiss/ndiv_ETmiss
  !   Distribution in eta of the top
      common/ext_eta356/eta356max,eta356min,eta356w
      common/dist_eta356/xeta356(500),fxeta356(500,20),fxeta356tot(500)
      common/inp_eta356/m_eta356
      common/div_eta356/ndiv_eta356
  !   Distribution in eta of the anti-top
      common/ext_eta478/eta478max,eta478min,eta478w
      common/dist_eta478/xeta478(500),fxeta478(500,20),fxeta478tot(500)
      common/inp_eta478/m_eta478
      common/div_eta478/ndiv_eta478
  !   Distribution in invarient mass of the top pair
      common/ext_rmass/rmassmax,rmassmin,rmassw
      common/dist_rmass/xrmass(500),fxrmass(500,20),fxrmasstot(500)
      common/inp_rmass/m_rmass
      common/div_rmass/ndiv_rmass
  !   Distribution in boost of top pair centre of mass frame
      common/ext_beta/betamax,betamin,betaw
      common/dist_beta/xbeta(500),fxbeta(500,20),fxbetatot(500)
      common/inp_beta/m_beta
      common/div_beta/ndiv_beta
  !   Distribution in cos(theta_t)
      common/ext_cost/costmax,costmin,costw
      common/dist_cost/xcost(500),fxcost(500,20),fxcosttot(500)
      common/inp_cost/m_cost
      common/div_cost/ndiv_cost
  !   Distribution in top energy
      common/ext_Et/Etmax,Etmin,Etw
      common/dist_Et/xEt(500),fxEt(500,20),fxEttot(500)
      common/inp_Et/m_Et
      common/div_Et/ndiv_Et
  !   Distributions in transverse variables
      integer ntrans
      parameter (ntrans=10)
      common/ext_trans/transmax(ntrans),transmin(ntrans),transw(ntrans)
      common/dist_trans/xtrans(ntrans,500),fxtrans(ntrans,500,20)
     &                                           ,fxtranstot(ntrans,500)
      common/inp_trans/m_trans(ntrans)
      common/div_trans/ndiv_trans(ntrans)
      dimension sfxtranstot(ntrans)
  !   Distribution in phi_l (lepton azimuthal angle)
      common/ext_fl/flmax,flmin,flw
      common/dist_fl/xfl(500),fxfl(500,20),fxfltot(500)
      common/inp_fl/m_fl
      common/div_fl/ndiv_fl
  !   Distribution in cos_phi_l
      common/ext_cosfl/cosflmax,cosflmin,cosflw
      common/dist_cosfl/xcosfl(500),fxcosfl(500,20),fxcosfltot(500)
      common/inp_cosfl/m_cosfl
      common/div_cosfl/ndiv_cosfl
  !   Distributions in asymmetries
      common/inp_asym/m_asy(8)  ! nasy      
  !   Distribution in sigp
      common/ext_sigp/sigpmax,sigpmin,sigpw
      common/dist_sigp/xsigp(1000),fxsigp(8,1000,20),fxsigptot(8,1000)
      common/inp_sigp/m_sigp
      common/div_sigp/ndiv_sigp
  !   Distribution in sigm
      common/ext_sigm/sigmmax,sigmmin,sigmw
      common/dist_sigm/xsigm(1000),fxsigm(8,1000,20),fxsigmtot(8,1000)
      common/inp_sigm/m_sigm
      common/div_sigm/ndiv_sigm

  ! Local variables
  !   Flag for Zp width specification
      dimension iwidth(5)
  !   Model name
      character*50 model
  !   Polarised/hemispherised cross sections
      dimension cnorm(20)
      dimension snorm(6)  !,ave(4) 
      dimension poltot(-1:1,-1:1),polchi(-1:1,-1:1)
      dimension asytot(5,-1:1),asychi(5,-1:1)
      dimension sfxpTtot(8),sfxetatot(8),sfxphitot(8)
      dimension sfxsigptot(8),sfxsigmtot(8)
      dimension asym_int(8)
      dimension Atot(8),Atoterr(8)

  ! Local constants
  !   pi
      parameter (pi=3.14159265358979323846d0)
  !   Unit conversion GeV -> nb
      parameter (conv=0.38937966d9)
  !   Date and time
      integer today(3), now(3)      
  !   Quark masses
      data rmu/0.00d0/,
     &     rmd/0.00d0/,
     &     rms/0.00d0/,
     &     rmc/0.00d0/,
     &     rmb/4.25d0/
  !   Higgs/gauge boson masses and widths
      data              gamW/2.08d0/
      data rmZ/91.19d0/,gamZ/2.50d0/
      data rmH/125.0d0/,gamH/0.31278d-2/
  !   Branching ratio for t->bln      
      real*8 BRtbln/0.10779733d0/
  !       real*8 BRtbln/1d0/

  ! External procedures
      external fxn
! ----------------------------------------------------------------------
! Read input files
  ! Read config file
  !   Permitted gauge sector flags
      read(5,*) iQCD
      read(5,*) iEW
      read(5,*) iBSM
  !   Interference flag
      read(5,*) iint
  !   Final state flag (ifinal=0: no top decay; ifinal=1: dileptonic top decay)      
      read(5,*) ifinal
  !   NWA flag (iNWA = 0: Actual top widths; iNWA = 1: tops in NWA)
      read(5,*) iNWA
  !   Name of model file
      read(5,*) model
  !   Collider energy     
      read(5,*) ecm_coll
  !   Collider flag (icoll=0: pp; icoll=1: ppbar)
      read(5,*) icoll
  !   PDFs
      read(5,*) istructure
  !   Cut on top rapidity 
      read(5,*) ytcut
  !   Cut on top pair boost
      read(5,*) yttcut
  !   Number of Vegas calls per iteration
      read(5,*) ncall
  !   Maximum number of Vegas iterations
      read(5,*) itmx
  !   Desired Accuracy (If negative, run maximum iterations.)
      read(5,*) acc
  !   Random number seed
      read(5,*) iseed
  !   Standard distributions flag
      read(5,*) idist
  !   Transverse mass distributions flag
      read(5,*) itdist
  !   Asymmetry distributions flag
      read(5,*) iadist
  !   Outout in lhe format
      read(5,*) ilhe
  !   !   Manually sum over costheta
  !       read(5,*) isycost

  ! Interpret config
    ! Number of external lines
      if(ifinal.eq.0)then 
        ipmax=4
      else if(ifinal.eq.1)then 
        ipmax=8
      else 
        write(*,*)'Invalid final state identifier!'
        stop
      end if
    ! NWA only for six-body final state
      if(ifinal.eq.0) iNWA=0
    ! itmx no more than 20.      
      if(itmx.gt.20)then
        write(*,*)'itmx does not have to exceed 20!'
        stop
      end if
    ! No event weighting for *true* event generation
      if(ilhe.eq.1)then
        itmx=1
      end if
    ! Extract model filename (Remove white space.)
      imodel = len(model)
      do while(model(imodel:imodel).eq.'') 
        imodel = imodel-1
      end do

  ! Read model file
      open(unit=42,file='Models/'//model(1:imodel)//'.mdl',status='old')
      read(42,*) rmZp
      read(42,*) gamZp      
      read(42,*) gp 
      read(42,*) paramZp
      read(42,*) gV_u
      read(42,*) gA_u
      read(42,*) gV_d
      read(42,*) gA_d      
  !   Check whether width has been specified
  !   (If gamZp is zero, the function widthZp is used instead.)
      do i=1,5
        if ((gamZp(i).eq.0d0).and.(rmZp(i).gt.0d0)) then
          iwidth(i) = 0
        else
          iwidth(i) = 1
        end if
      enddo
! ----------------------------------------------------------------------
! Distributions Setup
  ! (Set flags, binning range and divisions.)
  ! pT distributions
      do ip=1,ipmax
        m_pT(ip)=idist
        pTmax(ip)=7000.d0/(1+icoll*6)
        pTmin(ip)=0.d0
        ndiv_pT(ip)=175
  ! eta distributions
        m_eta(ip)=idist
        etamax(ip)=+10
        etamin(ip)=-10
        ndiv_eta(ip)=50
  ! phi distributions
        m_phi(ip)=idist
        phimax(ip)=+pi
        phimin(ip)=-pi
        ndiv_phi(ip)=100
      end do
  !   missing transverse momentum
      m_ETmiss=idist
      ETmissmax=7000.d0/(1+icoll*6)
      ETmissmin=0.d0
      ndiv_ETmiss=175      
  !   top transverse momentum
      m_pT356=idist
      pT356max=7000.d0/(1+icoll*6)
      pT356min=0.d0
      ndiv_pT356=175
  !   anti-top transverse momentum
      m_pT478=idist
      pT478max=7000.d0/(1+icoll*6)
      pT478min=0.d0
      ndiv_pT478=175
  !   2to6 top pseudorapidity
      m_eta356=idist
      eta356max=+10
      eta356min=-10
      ndiv_eta356=50
  !   2to6 anti-top pseudorapidity
      m_eta478=idist
      eta478max=+10
      eta478min=-10
      ndiv_eta478=50
  !   invarient mass of tt pair (always on)
      m_rmass=1
      rmassmax=14000.d0/(1+icoll*6)
      rmassmin=0.d0
      ndiv_rmass=500
  !   boost of parton CoM
      m_beta=idist
      betamax=1000.d0
      betamin=0.d0
      ndiv_beta=100
  !   costheta
      m_cost=idist
      costmax=+1.d0
      costmin=-1.d0
      ndiv_cost=50
  !   top energy
      m_Et=idist
      Etmax=7000.d0/(1+icoll*6)
      Etmin=0.d0
      ndiv_Et=175

  !   transverse variables
      do itrans=1,8
        m_trans(itrans)=itdist
      end do
  !   sum of tranvserse energy
      transmax(1)=ecm_col/2/(1+icoll*6)
      transmin(1)=0.d0
      ndiv_trans(1)=175
  !   invarient mass of the visible decay products of the tt pair
      transmax(2)=ecm_coll/(1+icoll*6)
      transmin(2)=0.d0
      ndiv_trans(2)=500   
  !   transverse mass
      transmax(3)=ecm_coll/(1+icoll*6)
      transmin(3)=0.d0
      ndiv_trans(3)=175
  !   transverse mass
      transmax(4)=ecm_coll/(1+icoll*6)
      transmin(4)=0.d0
      ndiv_trans(4)=175
  !   transverse mass
      transmax(5)=ecm_coll/(1+icoll*6)
      transmin(5)=0.d0
      ndiv_trans(5)=175
  !   contransverse mass 1
      transmax(6)=ecm_coll/(1+icoll*6)
      transmin(6)=0.d0
      ndiv_trans(6)=175     
  !   contransverse mass 2
      transmax(7)=ecm_coll/(1+icoll*6)
      transmin(7)=0.d0
      ndiv_trans(7)=175    
  !   contransverse mass 3
      transmax(8)=ecm_coll/(1+icoll*6)
      transmin(8)=0.d0
      ndiv_trans(8)=175            
  !   lepton transverse mass 
      transmax(9)=500.d0/(1+icoll*6)
      transmin(9)=0.d0
      ndiv_trans(9)=175
  !   lepton contransverse mass 
      transmax(10)=500.d0/(1+icoll*6)
      transmin(10)=0.d0
      ndiv_trans(10)=175      

  !   phi_l
      m_fl=iadist
      flmax=+2*pi
      flmin=0
      ndiv_fl=100
  !   cosphi_l
      m_cosfl=iadist
      cosflmax=+1.d0
      cosflmin=-1.d0
      ndiv_cosfl=100
  !   sigp
      m_sigp=iadist
      sigpmax=rmassmax
      sigpmin=rmassmin
      ndiv_sigp=ndiv_rmass/10
  !   sigm
      m_sigm=iadist
      sigmmax=rmassmax
      sigmmin=rmassmin
      ndiv_sigm=ndiv_rmass/10      

  !   asymmetries
      do i_asym=1,8   ! N_asym
        m_asy(i_asym)=iadist
      end do

  !   Turn off 2->6 only distributions
      if (ifinal.eq.0)then
        do ip=5,8
          m_pT(i)   = 0
          m_eta(i)  = 0
          m_phi(i)  = 0
        end do
        m_pT356  = 0
        m_eta356 = 0
        m_pT478  = 0
        m_eta478 = 0
        m_ETmiss = 0
        m_HT     = 0
        m_rM_T   = 0
        m_rM_CT  = 0
        m_rMlCT  = 0
        m_rMvis  = 0
        m_fl     = 0
        m_cosfl  = 0
        m_asy(8) = 0  ! turn off A_l
      end if
  !   Turn off 2->2 only distributions
      if (ifinal.eq.1)then
        m_asy(1) = 0 ! turn off A_LL
        m_asy(2) = 0 ! turn off A_L
        m_asy(3) = 0 ! turn off A_PV
      end if    
! ----------------------------------------------------------------------
! Set-up physics

  ! Collider CM energy squared.      
      s=ecm_coll*ecm_coll

  ! Factor outside integration
  !   Conversion GeV^-2 -> pb
      fac=conv
  !   Azimuthal angle integrated out (No initial transverse polarisation.)
      fac=fac*2.d0*pi      

  ! QCDL4 is QCD LAMBDA4 (to match PDF fits).
  ! (PDFs are intrinsically linked to the value of lamda_QCD; alpha_QCD) 
      if(istructure.eq.1)qcdl4=0.326d0
      if(istructure.eq.2)qcdl4=0.326d0
      if(istructure.eq.3)qcdl4=0.326d0
      if(istructure.eq.4)qcdl4=0.215d0
      if(istructure.eq.5)qcdl4=0.300d0
      if(istructure.eq.6)qcdl4=0.300d0
      if(istructure.eq.7)qcdl4=0.300d0
      if(istructure.eq.8)qcdl4=0.229d0
      if(istructure.eq.9)qcdl4=0.383d0
      rlambdaQCD4=QCDL4

  ! Initialise CTEQ grids.
      if(istructure.le.4)then
        icteq=istructure
        call SetCtq6(ICTEQ)
      end if

  ! Use appropriately evolved alphas.
      if(istructure.eq.1)nloops=2
      if(istructure.eq.2)nloops=2
      if(istructure.eq.3)nloops=1
      if(istructure.eq.4)nloops=1
      if(istructure.eq.5)nloops=1
      if(istructure.eq.6)nloops=1
      if(istructure.eq.7)nloops=1
      if(istructure.eq.8)nloops=1
      if(istructure.eq.9)nloops=1

  ! Initialize MadGraph for MEs
      rmt=175.d0
      gamt=1.55d0
      gNWA=gamt
      if(iNWA.eq.1)gamt=1.d-5
      call initialize(rmt,gamt)

  ! some EW parameters.
      a_em=1.d0/128.d0
      s2w=.2320d0
      rmW=rmZ*sqrt(1.d0-s2w)
      sw=sqrt(s2w)
      c2w=1.d0-s2w
      cw=sqrt(c2w)

  ! Calculate Zp couplings
      call coupZp

  ! Calculate sequential Zp widths
      do i=1,5
        if (iwidth(i).eq.0) gamZp(i)=
     &               widthZp(rmW,rmZ,rmZp(i),a_em,s2w,rlambdaQCD4,nloop)
      end do
! ---------------------------------------------------------------------- 
! VEGAS parameters
  ! Dimensions of integration
      if(ifinal.eq.0)then
        ndim=3
      else if(ifinal.eq.1)then
        ndim=15
      end if
  !   (If nprn<0 no print-out.)
      nprn=0
      if(ifinal.eq.0)then
  !   Final state masses     
        rm3=rmt
        rm4=rmt
        rm5=0.d0
        rm6=0.d0
        rm7=0.d0
        rm8=0.d0
  !   Integrates on:
  !   x(3)=(x1-tau)/(1-tau),
  !   x(2)=(ecm-rm3-rm4)/(ecm_max-rm3-rm4),
  !   x(1)=cos(theta3_cm)
  !   Limits:
        do i=3,2,-1
          xl(i)=0.d0
          xu(i)=1.d0
        end do
  !         if(isycost.eq.1)then  ! might not work
  !           nctpoints = 200
  !           do i=1,1
  !             xl(i)=0.d0
  !             xu(i)=0.d0
  !           end do
  !         else
  !           nctpoints = 0
        do i=1,1
          xl(i)=-1.d0
          xu(i)=1.d0
        end do
  !         end if

      else if(ifinal.eq.1)then
  !   Final state masses
        rm3=rmb
        rm4=rmb
        rm5=0.d0
        rm6=0.d0
        rm7=0.d0
        rm8=0.d0
  !   Integrates on:
  !   x(9)=cos(theta_cm_356)=-cos(theta_cm_478),
  !   x(15)=(x1-tau)/(1-tau),
  !   x(14)=(ecm-rm3-rm4-rm5-rm6-rm7-rm8)
  !        /(ecm_max-rm3-rm4-rm5-rm6-rm7-rm8),
  !   x(13)=(XX356-XX356min)/(XX356max-XX356min),
  !   where XX356=arctg((rm356**2-rmt**2)/rmt/gamt),
  !   x(12)=(XX478-XX478min)/(XX478max-XX478min),
  !   where XX478=arctg((rm478**2-rmt**2)/rmt/gamt),
  !   x(11)=(XX56-XX56min)/(XX56max-XX56min),
  !   where XX56=arctg((rm56**2-rmW**2)/rmW/gamW),
  !   x(10)=(XX78-XX78min)/(XX78max-XX78min),
  !   where XX78=arctg((rm78**2-rmW**2)/rmW/gamW),
  !   
  !   x(8)=cos(theta56_cm_356),
  !   x(7)=cos(theta78_cm_478),
  !   x(6)=cos(theta5_cm_56),
  !   x(5)=cos(theta7_cm_78),
  !   x(4)=fi56_cm_356,
  !   x(3)=fi78_cm_478,
  !   x(2)=fi5_cm_56,
  !   x(1)=fi8_cm_78;
  !   Limits:
        do i=15,14,-1
          xl(i)=0.d0
          xu(i)=1.d0
        end do
        do i=13,10,-1
          xl(i)=0.d0
          xu(i)=1.d0
        end do
        do i=9,5,-1
          xl(i)=-1.d0
          xu(i)=1.d0
        end do
        do i=4,1,-1
          xl(i)=0.d0
          xu(i)=2.d0*pi
        end do
      end if
! ----------------------------------------------------------------------
! Generate bins
  ! (Finds bin width, finds midpoints.)

      do ip=3,ipmax
        if(m_pT(ip).eq.1)then
          pTw(ip)=(pTmax(ip)-pTmin(ip))/ndiv_pT(ip)
          do j=1,ndiv_pT(ip)
            xpT(ip,j)=pTmin(ip)+pTw(ip)*(j-1)+pTw(ip)/2.d0
          end do
        end if
        if(m_eta(ip).eq.1)then
          etaw(ip)=(etamax(ip)-etamin(ip))/ndiv_eta(ip)
          do j=1,ndiv_eta(ip)
            xeta(ip,j)=etamin(ip)+etaw(ip)*(j-1)+etaw(ip)/2.d0
          end do
        end if
        if(m_phi(ip).eq.1)then
          phiw(ip)=(phimax(ip)-phimin(ip))/ndiv_phi(ip)
          do j=1,ndiv_phi(ip)
            xphi(ip,j)=phimin(ip)+phiw(ip)*(j-1)+phiw(ip)/2.d0
          end do
        end if
      end do

      if(m_ETmiss.eq.1)then
        ETmissw=(ETmissmax-ETmissmin)/ndiv_ETmiss
        do i=1,ndiv_ETmiss
          xETmiss(i)=ETmissmin+ETmissw*(i-1)+ETmissw/2.d0
        end do
      end if

      if(m_pT356.eq.1)then
        pT356w=(pT356max-pT356min)/ndiv_pT356
        do i=1,ndiv_pT356
          xpT356(i)=pT356min+pT356w*(i-1)+pT356w/2.d0
        end do
      end if

      if(m_pT478.eq.1)then
        pT478w=(pT478max-pT478min)/ndiv_pT478
        do i=1,ndiv_pT478
          xpT478(i)=pT478min+pT478w*(i-1)+pT478w/2.d0
        end do
      end if      

      if(m_eta356.eq.1)then
        eta356w=(eta356max-eta356min)/ndiv_eta356
        do i=1,ndiv_eta356
          xeta356(i)=eta356min+eta356w*(i-1)+eta356w/2.d0
        end do
      end if

      if(m_eta478.eq.1)then
        eta478w=(eta478max-eta478min)/ndiv_eta478
        do i=1,ndiv_eta478
          xeta478(i)=eta478min+eta478w*(i-1)+eta478w/2.d0
        end do
      end if

      if(m_rmass.eq.1)then
        rmassw=(rmassmax-rmassmin)/ndiv_rmass
        do i=1,ndiv_rmass
          xrmass(i)=rmassmin+rmassw*(i-1)+rmassw/2.d0
        end do
      end if

      if(m_beta.eq.1)then
        betaw=(betamax-betamin)/ndiv_beta
        do i=1,ndiv_beta
          xbeta(i)=betamin+betaw*(i-1)+betaw/2.d0
        end do
      end if
      if(m_cost.eq.1)then
        costw=(costmax-costmin)/ndiv_cost
        do i=1,ndiv_cost
          xcost(i)=costmin+costw*(i-1)+costw/2.d0
        end do
      end if

      if(m_Et.eq.1)then
        Etw=(Etmax-Etmin)/ndiv_Et
        do i=1,ndiv_Et
          xEt(i)=Etmin+Etw*(i-1)+Etw/2.d0
        end do
      end if

      do itrans=1,ntrans
        if(m_trans(itrans).eq.1)then
          transw(itrans)=(transmax(itrans)-transmin(itrans))
     &                                               /ndiv_trans(itrans)
          do i=1,ndiv_trans(itrans)
            xtrans(itrans,i)=transmin(itrans)+transw(itrans)*(i-1)
     &                                              +transw(itrans)/2.d0
          end do
        end if
      end do

      if(m_fl.eq.1)then
        flw=(flmax-flmin)/ndiv_fl
        do i=1,ndiv_fl
          xfl(i)=flmin+flw*(i-1)+flw/2.d0
        end do
      end if

      if(m_cosfl.eq.1)then
        cosflw=(cosflmax-cosflmin)/ndiv_cosfl
        do i=1,ndiv_cosfl
          xcosfl(i)=cosflmin+cosflw*(i-1)+cosflw/2.d0
        end do
      end if

      if(m_sigp.eq.1)then
        sigpw=(sigpmax-sigpmin)/ndiv_sigp
        do i=1,ndiv_sigp
          xsigp(i)=sigpmin+sigpw*(i-1)+sigpw/2.d0
        end do
      end if

      if(m_sigm.eq.1)then
        sigmw=(sigmmax-sigmmin)/ndiv_sigm
        do i=1,ndiv_sigm
          xsigm(i)=sigmmin+sigmw*(i-1)+sigmw/2.d0
        end do
      end if
! ----------------------------------------------------------------------
! Output information before integration
      write(*,*)'====================================================='
      call idate(today)   ! today(1)=day, (2)=month, (3)=year
      call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
      write(*,*)'DATE ',today(3),today(2),today(1)
      write(*,*)'TIME ',now(1),now(2),now(3)
      write(*,*)'-----------------------------------------------------'
      write(*,*)'PROCESS'
      if(icoll.eq.0)then
        if(ifinal.eq.0)
     &    write(*,*)'pp #rightarrow t#bar{t}',
     &               ' #times BR(t#rightarrow bl#nu)^{2}'
        if(ifinal.eq.1)
     &    write(*,*)'pp #rightarrow t#bar{t}',
     &               '#rightarrow b#bar{b} W^{+}W^{-}',
     &               '#rightarrow b#bar{b} l^{+}l^{-} #nu#bar{#nu}'
      else if(icoll.eq.1)then
        if(ifinal.eq.0)
     &    write(*,*)'p#bar{p} #rightarrow t#bar{t}',
     &               ' #times BR(t#rightarrow bl#nu)^{2}'
        if(ifinal.eq.1) 
     &    write(*,*)'p#bar{p} #rightarrow t#bar{t}',
     &               '#rightarrow b#bar{b} W^{+}W^{-}',
     &               '#rightarrow b#bar{b} l^{+}l^{-} #nu#bar{#nu}'
      end if
      write(*,*)'-----------------------------------------------------'
      write(*,*)'NOTES'            
      write(*,*)'Units: GeV'
      write(*,*)'Quarks: all massless except t, b.'           
      if(istructure.eq.1)write(*,*)'PDFs: cteq6m.'
      if(istructure.eq.2)write(*,*)'PDFs: cteq6d.'
      if(istructure.eq.3)write(*,*)'PDFs: cteq6l.'
      if(istructure.eq.4)write(*,*)'PDFs: cteq6l1.'
      if(istructure.eq.5)write(*,*)'PDFs: mrs99 (cor01).'
      if(istructure.eq.6)write(*,*)'PDFs: mrs99 (cor02).'
      if(istructure.eq.7)write(*,*)'PDFs: mrs99 (cor03).'
      if(istructure.eq.8)write(*,*)'PDFs: mrs99 (cor04).'
      if(istructure.eq.9)write(*,*)'PDFs: mrs99 (cor05).'
      if((ifinal.eq.1).and.(iNWA.eq.0))write(*,*)'Tops: off-shell.'
      if((ifinal.eq.1).and.(iNWA.eq.1))write(*,*)'Tops: NWA.'
      write(*,*)'BSM model: ',model
      if(iQCD.eq.1)write(*,*)'QCD: On '
      if(iQCD.eq.0)write(*,*)'QCD: Off'
      if(iEW.eq.1) write(*,*)'EW:  On '
      if(iEW.eq.0) write(*,*)'EW:  Off'
      if(iBSM.eq.1)write(*,*)'BSM: On '
      if(iBSM.eq.0)write(*,*)'BSM: Off'   
      if(iint.eq.0)write(*,*)'Interference: None'
      if(iint.eq.1)write(*,*)'Interference: SM'
      if(iint.eq.2)write(*,*)'Interference: Full'
      if(iint.eq.3)write(*,*)'Interference: No square terms.'      
      write(*,*)'-----------------------------------------------------'      
      write(*,*)'PARAMETERS'
      write(*,*)'#sqrt{s}              ',ecm_coll
      write(*,*)'at |y| <              ',abs(ytcut)
      write(*,*)'Loops a_s evaluated at',nloops
      write(*,*)'a_{s}(M_{Z})          ',alfas(rmZ,rLambdaQCD4,nloops)
      write(*,*)'#Lambda_{QCD}(4)      ',QCDL4
      write(*,*)'m_{b}                 ',rmb    
      write(*,*)'#Gamma_{b}            ',0.d0   
      write(*,*)'m_{t}                 ',rmt    
      write(*,*)'#Gamma_{t}            ',gamt   
      write(*,*)'m_{Z}                 ',rmZ    
      write(*,*)'#Gamma_{Z}            ',gamZ   
      write(*,*)'m_{W}                 ',rmW    
      write(*,*)'#Gamma_{W}            ',gamW   
      write(*,*)'m_{H}                 ',rmH    
      write(*,*)'#Gamma_{H}            ',gamH
      write(*,*)'-----------------------------------------------------' 
      write(*,*)'ZPRIME PARAMETERS'
      do i=1,5
        if(rmZp(i).gt.0)then
          write(*,*)'Z#prime',i 
          write(*,*)'m_{Z#prime}           ',rmZp(i) 
          write(*,*)'iwidth:               ',iwidth(i)   
          write(*,*)'#Gamma_{Z#prime}      ',gamZp(i)
          write(*,*)'g_{p}                 ',gp(i)
          write(*,*)'g_{V}^{u}             ',gV_u(i)
          write(*,*)'g_{A}^{u}             ',gA_u(i)
          write(*,*)'g_{V}^{d}             ',gV_d(i)
          write(*,*)'g_{A}^{d}             ',gA_d(i)
          write(*,*)         
        end if
      end do
! ----------------------------------------------------------------------
! Integration
  ! Section header
      write(*,*)'-----------------------------------------------------'
      write(*,*)'INTEGRATION'
  !   Reset counter
      npoints=0  
  !   Reset various iterative quantities
      if(ifinal.eq.0)then
        do i=1,20
          resl(i)=0.d0
          standdevl(i)=0.d0
          cnorm(i)=0.d0
          do iphel=-1,+1,2
            do jphel=-1,+1,2
              polcross(i,iphel,jphel)=0.d0
              polerror(i,iphel,jphel)=0.d0
            end do
          end do
          do jasy=1,5   ! nasy-3
            do iasy=-1,+1,2
              asycross(jasy,i,iasy)=0.d0
              asyerror(jasy,i,iasy)=0.d0
            end do
          end do 
        end do
      end if
  !   Integrate
      it=0
      call vegas(ndim,fxn,avgi,sd,chi2a)
      if(ifinal.eq.0)then
  !   Multiply by branching ratios (if ifinal = 0)      
        avgi=avgi*(BRtbln)**2
        sd=sd*(BRtbln)**2
      end if
  !  Collect total cross-section
      cross=avgi
      error=sd


    ! Print integrated cross section
      write(*,*)'-----------------------------------------------------'
      write(*,*)'INTEGRATED CROSS-SECTION'
      if(cross.eq.0d0)then
        write(*,*)'sigma = 0! Check permitted gauge sectors.'
        stop
      else      
        write(*,*)'sigma (pb)','error (same units)'
        write(*,*)cross,error
        write(*,*)'(using ',npoints,' points)'
      end if
  ! Re-weight distributions for different iterations     
      stantot=0.d0
        do i=1,it
          stantot=stantot+1.d0/standdevl(i)/standdevl(i)
        end do
        do i=1,it
          standdevl(i)=standdevl(i)*standdevl(i)*stantot
        end do
        do i=1,it
          cnorm(i)=resl(i)*standdevl(i)
        end do
! ----------------------------------------------------------------------
! Total asymmetries
  ! Collect polarised cross sections.
      if(iadist.eq.1)then
        if(ifinal.eq.0)then  
    
          do iphel=-1,+1,2
            do jphel=-1,+1,2
              do i=1,it
                polcross(i,iphel,jphel)=polcross(i,iphel,jphel)
     &                               *avgi/cnorm(i)
                polerror(i,iphel,jphel)=polcross(i,iphel,jphel)
     &                               *sd/cnorm(i)
              end do
              poltot(iphel,jphel)=0.d0
              polchi(iphel,jphel)=0.d0
              do i=1,it
                poltot(iphel,jphel)=poltot(iphel,jphel)
     &                       +polcross(i,iphel,jphel)
                polchi(iphel,jphel)=polchi(iphel,jphel)
     &                       +polerror(i,iphel,jphel)
              end do
              polchi(iphel,jphel)=polchi(iphel,jphel)
     &                         /poltot(iphel,jphel)
    !          polchi(iphel,jphel)=
    !   & sqrt(abs(polchi(iphel,jphel)
    !   &         -poltot(iphel,jphel)**2*dfloat(ncall)))
    !   & /dfloat(ncall)
            end do
          end do
        end if

    ! Collect unpolarised spatial asymmetry
        do iasy=1,5   ! nasy-3
    !         write(*,*)iasy,m_asy(iasy+3)
          if(m_asy(iasy+3).eq.0)then
            continue
          else
            do iAB=-1,+1,2
              do i=1,it
                asycross(iasy,i,iAB)=asycross(iasy,i,iAB)
     &                            *avgi/cnorm(i)
                asyerror(iasy,i,iAB)=asycross(iasy,i,iAB)
     &                            *sd/cnorm(i)
                write(*,*)asycross(iasy,i,iAB)
              end do        
              asytot(iasy,iAB)=0.d0
              asychi(iasy,iAB)=0.d0
              do i=1,it   ! add up each iteration
                asytot(iasy,iAB)=asytot(iasy,iAB)
     &                      +asycross(iasy,i,iAB)
                asychi(iasy,iAB)=asychi(iasy,iAB)
     &                      +asyerror(iasy,i,iAB)
              end do
    !             write(*,*)'bork1011'
              asychi(iasy,iAB)=asychi(iasy,iAB)
     &                        /asytot(iasy,iAB)
    !             write(*,*)'bork1014'
    !           asychi(iasy)=
    !      & sqrt(abs(asychi(iasy)
    !      &         -asytot(iasy)**2*dfloat(ncall)))
    !      & /dfloat(ncall)
            end do
          end if
        end do

    ! Define asymmetries
        if(ifinal.eq.0)then
    ! ALL
          Atot(1)=
     &          +(poltot(+1,+1)-poltot(+1,-1)
     &           -poltot(-1,+1)+poltot(-1,-1))
     &          /cross
          Atoterr(1)=
     &             +(polchi(+1,+1)+polchi(+1,-1)
     &              +polchi(-1,+1)+polchi(-1,-1))
     &             /4.d0*Atot(1)
    ! AL
          Atot(2)=
     &          +(poltot(-1,-1)-poltot(+1,-1) 
     &           +poltot(-1,+1)-poltot(+1,+1))
     &          /cross
          Atoterr(2)=
     &            +(polchi(-1,-1)+polchi(+1,-1) 
     &             +polchi(-1,+1)+polchi(+1,+1))
     &            /4.d0*Atot(2)
    ! APV
          Atot(3)=
     &          +(poltot(-1,-1)-poltot(+1,+1))
     &          /cross/2.d0
          Atoterr(3)=
     &             +(polchi(-1,-1)+polchi(+1,+1))
     &             /2.d0*Atot(3)
        end if
    ! AFBcm
        Atot(4)=
     &          +(asytot(1,+1)-asytot(1,-1))
     &          /cross
          Atoterr(4)=
     &            +sd/avgi*Atot(4)
    ! AtFB
        Atot(5)=
     &          +(asytot(2,+1)-asytot(2,-1))
     &          /cross
        Atoterr(5)=
     &            +sd/avgi*Atot(5)

        if(m_asy(6).gt.0)then
    ! A
          Atot(6)=
     &          +(asytot(3,+1)-asytot(3,-1))
     &          /cross
          Atoterr(6)=
     &            +sd/avgi*Atot(6)
        end if

        if(m_asy(7).gt.0)then
    ! A'
          Atot(7)=
     &          +(asytot(4,+1)-asytot(4,-1))
     &          /cross
          Atoterr(7)=
     &            +sd/avgi*Atot(7)
        end if

        
        if(m_asy(8).gt.0)then
    ! A_l
          Atot(8)=
     &          +(asytot(5,+1)-asytot(5,-1))
     &          /cross
          Atoterr(8)=
     &            +sd/avgi*Atot(8)
        end if
      

  ! Print Asymmetries     
        write(*,*)'TOTAL ASYMMETRIES'
        if(ifinal.eq.0)then
          write(*,*)'ALL:                  uncertainty (same units):'
          write(*,*)Atot(1),Atoterr(1) 
          write(*,*)'AL:                   uncertainty (same units):'
          write(*,*)Atot(2),Atoterr(2) 
          write(*,*)'APV:                  uncertainty (same units):'
          write(*,*)Atot(3),Atoterr(3) 
          write(*,*)'AFB*:                 uncertainty (same units):'
          write(*,*)Atot(4),Atoterr(4)
          write(*,*)'AtFB:                    uncertainty (same units):'
          write(*,*)Atot(5),Atoterr(5)
          write(*,*)'A:                    uncertainty (same units):'
          write(*,*)Atot(6),Atoterr(6)
          write(*,*)"A':                    uncertainty (same units):"
          write(*,*)Atot(7),Atoterr(7)
        else if(ifinal.gt.0)then
          write(*,*)'A_l:                  uncertainty (same units):'
          write(*,*)Atot(6),Atoterr(6) 
        end if
      end if
! ----------------------------------------------------------------------
! Plot Distributions
  ! Section header
      write(*,*)'-----------------------------------------------------'
      write(*,*)'HISTOGRAMS'
      do ip=3,8
  ! Plot distributions in pT
        if(m_pT(ip).eq.1)then        
          sfxpTtot(ip)=0d0
          do j=1,ndiv_pT(ip)
            fxpTtot(ip,j)=0.d0
            do i=1,it
              fxpT(ip,j,i)=fxpT(ip,j,i)*avgi/cnorm(i)/pTw(ip)
              fxpTtot(ip,j)=fxpTtot(ip,j)+fxpT(ip,j,i)
            end do
            sfxpTtot(ip)=sfxpTtot(ip)+fxpTtot(ip,j)*pTw(ip)
          end do          
          write(*,*)'HISTOGRAM'
          write(*,'(A,I1)')'pT',ip
          write(*,'(A,I1,A)')'d#sigma-/dp_{T}(',ip,')--[pb/GeV]'
          write(*,'(A,I1,A)')'p_{T}(',ip,')--[GeV]'
          do i=1,ndiv_pT(ip)
            write(*,*)xpT(ip,i),fxpTtot(ip,i)
          end do
          write(*,*)'END'
        end if
  ! Plot distributions in eta
        if(m_eta(ip).eq.1)then
          sfxetatot(ip)=0d0
          do j=1,ndiv_eta(ip)
            fxetatot(ip,j)=0.d0
            do i=1,it
              fxeta(ip,j,i)=fxeta(ip,j,i)*avgi/cnorm(i)/etaw(ip)
              fxetatot(ip,j)=fxetatot(ip,j)+fxeta(ip,j,i)            
            end do
            sfxetatot(ip)=sfxetatot(ip)+fxetatot(ip,j)*etaw(ip)
          end do
          write(*,*)'HISTOGRAM'
          write(*,'(A,I1)')'eta',ip
          write(*,'(A,I1,A)')'d#sigma-/d#eta(',ip,')--[pb]'
          write(*,'(A,I1,A)')'#eta(',ip,')'          
          do i=1,ndiv_eta(ip)
            write(*,*)xeta(ip,i),fxetatot(ip,i)
          end do
          write(*,*)'END'
        end if
  ! Plot distributions in phi        
        if(m_phi(ip).eq.1)then
          sfxphitot(ip)=0d0
          do j=1,ndiv_phi(ip)
            fxphitot(ip,j)=0.d0
            do i=1,it
              fxphi(ip,j,i)=fxphi(ip,j,i)*avgi/cnorm(i)/phiw(ip)
              fxphitot(ip,j)=fxphitot(ip,j)+fxphi(ip,j,i)            
            end do
            sfxphitot(ip)=sfxphitot(ip)+fxphitot(ip,j)*phiw(ip)
          end do
          write(*,*)'HISTOGRAM'
          write(*,'(A,I1)')'phi',ip
          write(*,'(A,I1,A)')'d#sigma-/d#phi(',ip,')--[pb/rad]'
          write(*,'(A,I1,A)')'#phi(',ip,')--[rad]'          
          do i=1,ndiv_phi(ip)
            write(*,*)xphi(ip,i),fxphitot(ip,i)
          end do
          write(*,*)'END'
        end if
      end do
  
  ! Plot distribution in ETmiss
      if(m_ETmiss.eq.1)then
        sfxETmisstot=0d0
        do j=1,ndiv_ETmiss
          fxETmisstot(j)=0.d0
          do i=1,it
            fxETmiss(j,i)=fxETmiss(j,i)*avgi/cnorm(i)/ETmissw
            fxETmisstot(j)=fxETmisstot(j)+fxETmiss(j,i)            
          end do
          sfxETmisstot=sfxETmisstot+fxETmisstot(j)*ETmissw
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'ETmiss'
        write(*,*)'d#sigma-/dp_{Tmiss}--[pb/GeV]'
        write(*,*)'p_{T}(miss)--[GeV]'
        do i=1,ndiv_ETmiss
          write(*,*)xETmiss(i),fxETmisstot(i)
        end do
        write(*,*)'END'
      end if

  ! Plot distribution in pT356
      if(m_pT356.eq.1)then
        sfxpT356tot=0d0
        do j=1,ndiv_pT356
          fxpT356tot(j)=0.d0
          do i=1,it
            fxpT356(j,i)=fxpT356(j,i)*avgi/cnorm(i)/pT356w
            fxpT356tot(j)=fxpT356tot(j)+fxpT356(j,i)            
          end do
          sfxpT356tot=sfxpT356tot+fxpT356tot(j)*pT356w
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'pT356'
        write(*,*)'d#sigma-/dp_{T}--[pb/GeV]'
        write(*,*)'p_{T}(t)--[GeV]'
        do i=1,ndiv_pT356
          write(*,*)xpT356(i),fxpT356tot(i)
        end do
        write(*,*)'END'
      end if
  ! Plot distribution in pT478
      if(m_pT478.eq.1)then  
        sfxpT478tot=0d0
        do j=1,ndiv_pT478
          fxpT478tot(j)=0.d0
          do i=1,it
            fxpT478(j,i)=fxpT478(j,i)*avgi/cnorm(i)/pT478w
            fxpT478tot(j)=fxpT478tot(j)+fxpT478(j,i)            
          end do
          sfxpT478tot=sfxpT478tot+fxpT478tot(j)*pT478w
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'pT478'
        write(*,*)'d#sigma-/dp_{T}--[pb/GeV]'
        write(*,*)'p_{T}(#bar{t})--[GeV]'
        do i=1,ndiv_pT478
          write(*,*)xpT478(i),fxpT478tot(i)
        end do
        write(*,*)'END'
      end if      
  
  ! Plot distribution in eta356
      if(m_eta356.eq.1)then

        do j=1,ndiv_eta356
          do i=1,it
            fxeta356(j,i)=fxeta356(j,i)*avgi/cnorm(i)/eta356w
          end do
        end do
        do j=1,ndiv_eta356
          fxeta356tot(j)=0.d0 
          do i=1,it
            fxeta356tot(j)=fxeta356tot(j)+fxeta356(j,i)
          end do
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'eta356'
        write(*,*)'d#sigma-/d#eta--[pb]'
        write(*,*)'#eta(t)'
        do i=1,ndiv_eta356
          write(*,*)xeta356(i),fxeta356tot(i)
        end do
        write(*,*)'END'
      end if
  ! Plot distribution in eta478
      if(m_eta478.eq.1)then

        sfxeta478tot=0d0
        do j=1,ndiv_eta478
          fxeta478tot(j)=0.d0
          do i=1,it
            fxeta478(j,i)=fxeta478(j,i)*avgi/cnorm(i)/eta478w
            fxeta478tot(j)=fxeta478tot(j)+fxeta478(j,i)            
          end do
          sfxeta478tot=sfxeta478tot+fxeta478tot(j)*eta478w
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'eta478'
        write(*,*)'d#sigma-/d#eta--[pb]'
        write(*,*)'#eta(#bar{t})'
        do i=1,ndiv_eta478
          write(*,*)xeta478(i),fxeta478tot(i)
        end do
        write(*,*)'END'
      end if      
  ! Plot distribution in Mtt
      if(m_rmass.eq.1)then

        sfxrmasstot=0d0
        do j=1,ndiv_rmass
          fxrmasstot(j)=0.d0
          do i=1,it
            fxrmass(j,i)=fxrmass(j,i)*avgi/cnorm(i)/rmassw
            fxrmasstot(j)=fxrmasstot(j)+fxrmass(j,i)            
          end do
          sfxrmasstot=sfxrmasstot+fxrmasstot(j)*rmassw
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'Mtt'
        write(*,*)'d#sigma-/dM_{tt}--[pb/GeV]'
        write(*,*)'M_{tt}--[GeV]'
        do i=1,ndiv_rmass
          write(*,*)xrmass(i),fxrmasstot(i)
        end do
        write(*,*)'END'
      end if
  ! Plot distribution in beta.
      if(m_beta.eq.1)then

        sfxbetatot=0d0
        do j=1,ndiv_beta
          fxbetatot(j)=0.d0
          do i=1,it
            fxbeta(j,i)=fxbeta(j,i)*avgi/cnorm(i)/betaw
            fxbetatot(j)=fxbetatot(j)+fxbeta(j,i)            
          end do
          sfxbetatot=sfxbetatot+fxbetatot(j)*betaw
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'Beta'
        write(*,*)'d#sigma-/d#Beta_{t}--[pb]'
        write(*,*)'#Beta_{t}'
        do i=1,ndiv_beta
          write(*,*)xbeta(i),fxbetatot(i)
        end do
        write(*,*)'END'
      end if  
  ! Plot distribution in cost
      if(m_cost.eq.1)then

        sfxcosttot=0d0
        do j=1,ndiv_cost
          fxcosttot(j)=0.d0
          do i=1,it
            fxcost(j,i)=fxcost(j,i)*avgi/cnorm(i)/costw
            fxcosttot(j)=fxcosttot(j)+fxcost(j,i)            
          end do
          sfxcosttot=sfxcosttot+fxcosttot(j)*costw
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'cost'
        write(*,*)'d#sigma-/dcos#theta--[pb]'
        write(*,*)'cos#theta'
        do i=1,ndiv_cost
          write(*,*)xcost(i),fxcosttot(i)
        end do
        write(*,*)'END'
      end if
  ! Plot distribution in fl
      if(m_fl.eq.1)then

        sfxfltot=0d0
        do j=1,ndiv_fl
          fxfltot(j)=0.d0
          do i=1,it
            fxfl(j,i)=fxfl(j,i)*avgi/cnorm(i)/flw
            fxfltot(j)=fxfltot(j)+fxfl(j,i)            
          end do
          sfxfltot=sfxfltot+fxfltot(j)*flw
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'fl'
        write(*,*)'d#sigma-/d#phi_{l}--[pb]'
        write(*,*)'#phi_{l}--[-rad-]'
        do i=1,ndiv_fl
          write(*,*)xfl(i),fxfltot(i)
        end do
        write(*,*)'END'
      end if 
  ! Plot distribution in cosfl
      if(m_cosfl.eq.1)then

        sfxcosfltot=0d0
        do j=1,ndiv_cosfl
          fxcosfltot(j)=0.d0
          do i=1,it
            fxcosfl(j,i)=fxcosfl(j,i)*avgi/cnorm(i)/cosflw
            fxcosfltot(j)=fxcosfltot(j)+fxcosfl(j,i)            
          end do
          sfxcosfltot=sfxcosfltot+fxcosfltot(j)*cosflw
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'cosfl'
        write(*,*)'d#sigma-/dcos#phi_{l}--[pb]'
        write(*,*)'cos#phi_{l}'
        do i=1,ndiv_cosfl
          write(*,*)xcosfl(i),fxcosfltot(i)
        end do
        write(*,*)'END'
      end if
  ! Plot distribution in Et
      if(m_Et.eq.1)then

        sfxEttot=0d0
        do j=1,ndiv_Et
          fxEttot(j)=0.d0
          do i=1,it
            fxEt(j,i)=fxEt(j,i)*avgi/cnorm(i)/Etw
            fxEttot(j)=fxEttot(j)+fxEt(j,i)            
          end do
          sfxEttot=sfxEttot+fxEttot(j)*Etw
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'Et'
        write(*,*)'d#sigma-/dE_{t}--[pb/GeV]'
        write(*,*)'E_{t}--[GeV]'
        do i=1,ndiv_Et
          write(*,*)xEt(i),fxEttot(i)
        end do
        write(*,*)'END'
      end if
  ! Plot distributions in all transverse variables
      do itrans=1,ntrans
        if(m_trans(itrans).eq.1)then
          sfxtranstot(itrans)=0d0
          do j=1,ndiv_trans(itrans)
            fxtranstot(itrans,j)=0.d0
            do i=1,it
              fxtrans(itrans,j,i)=fxtrans(itrans,j,i)
     &                           *avgi/cnorm(i)/transw(itrans)
              fxtranstot(itrans,j)=fxtranstot(itrans,j)
     &                            +fxtrans(itrans,j,i)            
            end do
            sfxtranstot(itrans)=sfxtranstot(itrans)+
     &                               fxtranstot(itrans,j)*transw(itrans)
          end do
          write(*,*)'HISTOGRAM'
          if (itrans.eq.1)then
            write(*,*)'Mvis'
            write(*,*)'d#sigma-/dM_{vis}--[pb/GeV]'
            write(*,*)'M_{T}--[GeV]'
          else if (itrans.eq.2)then
            write(*,*)'HT'
            write(*,*)'d#sigma-/dH_{T}--[pb/GeV]'
            write(*,*)'H_{T}--[GeV]'          
          else if (itrans.eq.3)then
            write(*,*)'M_T1'
            write(*,*)'d#sigma-/dM_{T1}--[pb/GeV]'
            write(*,*)'M_{T}--[GeV]'
          else if (itrans.eq.4)then
            write(*,*)'M_T2'
            write(*,*)'d#sigma-/dM_{T2}--[pb/GeV]'
            write(*,*)'M_{T}--[GeV]'
          else if (itrans.eq.5)then
            write(*,*)'M_T3'
            write(*,*)'d#sigma-/dM_{T3}--[pb/GeV]'
            write(*,*)'M_{T}--[GeV]'
          else if (itrans.eq.6)then
            write(*,*)'MlT'
            write(*,*)'d#sigma-/dM^{l}_{T}--[pb/GeV]'
            write(*,*)'M_{T}--[GeV]'
          else if (itrans.eq.7)then
            write(*,*)'M_CT1'
            write(*,*)'d#sigma-/dM_{T1}--[pb/GeV]'
            write(*,*)'M_{T}--[GeV]'
          else if (itrans.eq.8)then
            write(*,*)'M_CT2'
            write(*,*)'d#sigma-/dM_{T2}--[pb/GeV]'
            write(*,*)'M_{T}--[GeV]'
          else if (itrans.eq.9)then
            write(*,*)'M_CT3'
            write(*,*)'d#sigma-/dM_{T3}--[pb/GeV]'
            write(*,*)'M_{T}--[GeV]'
          else if (itrans.eq.10)then
            write(*,*)'MlCT'
            write(*,*)'d#sigma-/dM^{l}_{CT}--[pb/GeV]'
            write(*,*)'M^{l}_{CT}--[GeV]'
          else
            continue
          end if
          do i=1,ndiv_trans(itrans)
            write(*,*)xtrans(itrans,i),fxtranstot(itrans,i)
          end do
          write(*,*)'END'
        end if
      end do
  ! Plot distributions in all asymmetries
      if((m_sigp.eq.1).and.(m_sigm.eq.1))then
        do jasy=1,8
          if(m_asy(jasy).eq.0)then
            continue
          else
            ! snorm(jasy)=0.d0
            sfxsigptot(jasy)=0d0
            do j=1,ndiv_sigp
              fxsigptot(jasy,j)=0.d0
              do i=1,it
                fxsigp(jasy,j,i)=fxsigp(jasy,j,i)*avgi/cnorm(i)/sigpw
                fxsigptot(jasy,j)=fxsigptot(jasy,j)+fxsigp(jasy,j,i)
              end do
              sfxsigptot(jasy)=sfxsigptot(jasy)+fxsigptot(jasy,j)*sigpw
            end do
            sfxsigmtot(jasy)=0d0           
            do j=1,ndiv_sigm
              do i=1,it
                fxsigm(jasy,j,i)=fxsigm(jasy,j,i)*avgi/cnorm(i)/sigmw
                fxsigmtot(jasy,j)=fxsigmtot(jasy,j)+fxsigm(jasy,j,i)
              end do
              sfxsigmtot(jasy)=sfxsigmtot(jasy)+fxsigmtot(jasy,j)*sigmw
            end do
            write(*,*)'HISTOGRAM'
            if(jasy.eq.1)then
                write(*,*)'ALL'
                write(*,*)'A_{LL}'
            else if(jasy.eq.2)then
                write(*,*)'AL'
                write(*,*)'A_{L}'
            else if(jasy.eq.3)then
                write(*,*)'APV'
                write(*,*)'A_{PV}'
            else if(jasy.eq.4)then
                write(*,*)'AFBcm'
                write(*,*)'A_{FB^{*}}'
            else if(jasy.eq.5)then
                write(*,*)'AtFB'
                write(*,*)'A^{t}_{FB}'
            else if(jasy.eq.6)then
                write(*,*)'A'
                write(*,*)'A'
            else if(jasy.eq.7)then
                write(*,*)'Ap'
                write(*,*)"A'"
            else if(jasy.eq.8)then
                write(*,*)'A_l'
                write(*,*)'A_{l}'
            end if            
            write(*,*)'M_{tt}'
            ndiv_sig=(ndiv_sigm+ndiv_sigp)/2
            do i=1,ndiv_sig
              if(fxsigptot(jasy,i)+fxsigmtot(jasy,i).eq.0.d0)then
                write(*,*)(xsigm(i)+xsigp(i))/2.d0,0.d0
      !             snorm(jasy)=snorm(jasy)+0.d0
              else  
                write(*,*)(xsigm(i)+xsigp(i))/2.d0,
     &                   (fxsigptot(jasy,i)-fxsigmtot(jasy,i))/
     &                   (fxsigptot(jasy,i)+fxsigmtot(jasy,i))
      !               snorm(jasy)=snorm(jasy)+
      !      &               (fxsigptot(jasy,i)-fxsigmtot(jasy,i))/
      !      &               (fxsigptot(jasy,i)+fxsigmtot(jasy,i))
      !      &               *fxrmasstot(i)*rmassw/avgi      
              end if
            end do
            asym_int(jasy)=(sfxsigptot(jasy)-sfxsigmtot(jasy))/
     &                     (sfxsigptot(jasy)+sfxsigmtot(jasy))          
            write(*,*)'END'
      !           write(*,*)'(Total Asymmetry:',asym_int(jasy),')'
      !           write(*,*)'(Integrated Asymmetry:',snorm(jasy),' )'
          end if
        end do
      end if
! ----------------------------------------------------------------------
! Check distributions
      diff_max=1E-12
      n_error=0
      do ip=3,ipmax
        if(m_pT(ip).eq.1)then
          if(abs(cross-sfxpTtot(ip))>diff_max)then
            write(*,*)'pT',ip,' Error:',sfxpTtot(ip)
            n_error=n_error+1
          end if
        end if
        if(m_eta(ip).eq.1)then
          if(abs(cross-sfxetatot(ip))>diff_max)then
            write(*,*)'eta',ip,' Error:',sfxetatot(ip)
            n_error=n_error+1
          end if
        end if
        if(m_phi(ip).eq.1)then
          if(abs(cross-sfxphitot(ip))>diff_max)then
            write(*,*)'phi',ip,' Error:',sfxphitot(ip)
            n_error=n_error+1
          end if
        end if
      end do
      if(m_pT356.eq.1)then
        if(abs(cross-sfxpT356tot)>diff_max)then
          write(*,*)'pT356 Error:',sfxpT356tot
          n_error=n_error+1
        end if
      end if
      if(m_pT478.eq.1)then
        if(abs(cross-sfxpT478tot)>diff_max)then
          write(*,*)'pT478 Error:',sfxpT478tot
          n_error=n_error+1
        end if
      end if
      if(m_ETmiss.eq.1)then
        if(abs(cross-sfxETmisstot)>diff_max)then
          write(*,*)'ETmiss Error:',sfxETmisstot
          n_error=n_error+1
        end if
      end if
      if(m_rmass.eq.1)then
        if(abs(cross-sfxrmasstot)>diff_max)then
          write(*,*)'rmass Error:',sfxrmasstot
          n_error=n_error+1
        end if
      end if
      if(m_beta.eq.1)then
        if(abs(cross-sfxbetatot)>diff_max*10)then
          write(*,*)'beta Error:',sfxbetatot
          n_error=n_error+1
        end if
      end if
      if(m_cost.eq.1)then
        if(abs(cross-sfxcosttot)>diff_max)then
          write(*,*)'cost Error:',sfxcosttot
          n_error=n_error+1
        end if
      end if
      if(m_Et.eq.1)then
        if(abs(cross-sfxEttot)>diff_max)then
          write(*,*)'Et Error:',sfxEttot
          n_error=n_error+1
        end if
      end if
      if(m_fl.eq.1)then
        if(abs(cross-sfxfltot)>diff_max)then
          write(*,*)'fl Error:',sfxfltot
          n_error=n_error+1
        end if
      end if
      if(m_cosfl.eq.1)then
        if(abs(cross-sfxcosfltot)>diff_max)then
          write(*,*)'cosfl Error:',sfxcosfltot
          n_error=n_error+1
        end if
      end if
      do iasy=1,8 !nasy
        if(m_asy(iasy).eq.0)then
          continue
        else
          if(abs(Atot(iasy)-asym_int(iasy))>diff_max)then
            write(*,*)'A Error:',iasy,asym_int(iasy)
            n_error=n_error+1
          end if
        end if
      end do
      write(*,*)'INTEGRATION ERRORS:',n_error
! ----------------------------------------------------------------------
! End program
      write(*,*)'CLOSE'
      write(*,*)'======================================================'
      stop
      end
! ======================================================================