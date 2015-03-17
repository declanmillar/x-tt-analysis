function sggff_qcd(p1, p2, p3, p4, lam3, lam4)

! !function generated by madgraph
! returns amplitude squared summed/avg over colors
! and helicities
! for the point in phase space p1,p2,p3,p4,lam3,lam4

! for process : g g  -> t t~

  implicit none

! constants
  integer :: iq ! incoming fermion species (see initialise for IDs)
  integer ::    nexternal,   ncomb
  parameter (nexternal=4, ncomb= 16)

! arguments

  integer :: lam3,lam4
  real :: p1(0:3),p2(0:3),p3(0:3),p4(0:3)

! local variables

  integer :: nhel(nexternal,ncomb),ntry
  real :: t
  real :: sggff_qcd
  real :: ggff_qcd
  integer :: ihel
  logical :: goodhel(ncomb)
  data goodhel/ncomb* .FALSE. /
  data ntry/0/
  data (nhel(ihel,  1),ihel=1,4) / -1, -1, -1, -1/
  data (nhel(ihel,  2),ihel=1,4) / -1, -1, -1,  1/
  data (nhel(ihel,  3),ihel=1,4) / -1, -1,  1, -1/
  data (nhel(ihel,  4),ihel=1,4) / -1, -1,  1,  1/
  data (nhel(ihel,  5),ihel=1,4) / -1,  1, -1, -1/
  data (nhel(ihel,  6),ihel=1,4) / -1,  1, -1,  1/
  data (nhel(ihel,  7),ihel=1,4) / -1,  1,  1, -1/
  data (nhel(ihel,  8),ihel=1,4) / -1,  1,  1,  1/
  data (nhel(ihel,  9),ihel=1,4) /  1, -1, -1, -1/
  data (nhel(ihel, 10),ihel=1,4) /  1, -1, -1,  1/
  data (nhel(ihel, 11),ihel=1,4) /  1, -1,  1, -1/
  data (nhel(ihel, 12),ihel=1,4) /  1, -1,  1,  1/
  data (nhel(ihel, 13),ihel=1,4) /  1,  1, -1, -1/
  data (nhel(ihel, 14),ihel=1,4) /  1,  1, -1,  1/
  data (nhel(ihel, 15),ihel=1,4) /  1,  1,  1, -1/
  data (nhel(ihel, 16),ihel=1,4) /  1,  1,  1,  1/
! ----------
! begin code
! ----------
  sggff_qcd = 0d0
  ntry=ntry+1
  do ihel=1,ncomb
  !          if (goodhel(ihel) .or. ntry .lt. 10) then
    t=ggff_qcd(iq,p1, p2, p3, p4,lam3,lam4,nhel(1,ihel))
    sggff_qcd = sggff_qcd + t
  !              if (t .gt. 0d0 .and. .not. goodhel(ihel)) then
  !                  goodhel(ihel)=.true.
  !                  write(*,*) ihel,t
  !              endif
  !          endif
  enddo
  sggff_qcd = sggff_qcd /  4d0
end function sggff_qcd
      
      
function ggff_qcd(iq, p1, p2, p3, p4,lam3,lam4,nhel)

! function generated by madgraph
! returns amplitude squared summed/avg over colors
! for the point in phase space p1,p2,p3,p4,...
! and helicity nhel(1),nhel(2),....

! for process : g g  -> t t~

  implicit none

! constants
  real :: ggff_qcd

  integer :: iq ! incoming fermion species (see initialise for IDs)
  integer ::    ngraphs,    neigen,    nexternal
  parameter (ngraphs=  3,neigen=  2,nexternal=4)
  real ::     zero
  parameter (zero=0d0)

! arguments

  real :: p1(0:3),p2(0:3),p3(0:3),p4(0:3)
  integer :: lam3,lam4
  integer :: nhel(nexternal)


! local variables

  integer :: i,j
  real :: eigen_val(neigen), eigen_vec(ngraphs,neigen)
  complex*16 ztemp
  complex*16 amp(ngraphs)
  complex*16 w1(6)  , w2(6)  , w3(6)  , w4(6)  , w5(6)
  complex*16 w6(6)  , w7(6)

! global variables

  real ::           gg(2), g
  common /coupqcd/ gg,    g
  real ::            fmass(12), fwidth(12)
  common /fermions/ fmass,     fwidth

! color data

  data eigen_val(1  )/       7.2916666666666588e-02 /
  data eigen_vec(1  ,1  )/   7.0710678118654768e-01 /
  data eigen_vec(2  ,1  )/   7.0710678118654735e-01 /
  data eigen_vec(3  ,1  )/   0.0000000000000000e+00 /
  data eigen_val(2  )/       2.8125000000000000e-01 /
  data eigen_vec(1  ,2  )/  -4.0824829046386313e-01 /
  data eigen_vec(2  ,2  )/   4.0824829046386285e-01 /
  data eigen_vec(3  ,2  )/   8.1649658092772603e-01 /
! ----------
! begin code
! ----------
  if((nhel(3) == lam3) .AND. (nhel(4) == lam4))then
    continue
  else
    ggff_qcd = 0.d0
    return
  end if

  call vxxxxx(p1  , zero,nhel(1  ),-1,w1  )
  call vxxxxx(p2  , zero,nhel(2  ),-1,w2  )
  call oxxxxx(p3  ,fmass(11 ),nhel(3  ), 1,w3  )
  call ixxxxx(p4  ,fmass(11 ),nhel(4  ),-1,w4  )
  call fvoxxx(w3  ,w2  ,gg,fmass(11 ),fwidth(11 ),w5  )
  call iovxxx(w4  ,w5  ,w1  ,gg,amp(1  ))
  call fvoxxx(w3  ,w1  ,gg,fmass(11 ),fwidth(11 ),w6  )
  call iovxxx(w4  ,w6  ,w2  ,gg,amp(2  ))
  call jggxxx(w1  ,w2  ,g,w7  )
  call iovxxx(w4  ,w3  ,w7  ,gg,amp(3  ))
  ggff_qcd = 0.d0
  do i = 1, neigen
    ztemp = (0.d0,0.d0)
    do j = 1, ngraphs
      ztemp = ztemp + eigen_vec(j,i)*amp(j)
    enddo
    ggff_qcd =ggff_qcd+ztemp*eigen_val(i)*conjg(ztemp)
  enddo
!      call gaugecheck(amp,ztemp,eigen_vec,eigen_val,ngraphs,neigen)
end function ggff_qcd