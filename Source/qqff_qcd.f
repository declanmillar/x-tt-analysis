function sqqff_qcd(iq,p1,p2,p3,p4,lam3,lam4)

  ! function generated by madgraph
  ! returns amplitude squared summed/avg over colors
  ! and helicities
  ! for the point in phase space p1,p2,p3,p4,...
  
  ! for process : u u~  -> t t~

  implicit none

  ! constants
  real :: sqqff_qcd
  integer :: iq ! incoming fermion species (see initialise for IDs)
  integer ::    nexternal,   ncomb
  parameter (nexternal=4, ncomb= 16)

  ! arguments
  integer :: lam3,lam4
  real :: p1(0:3),p2(0:3),p3(0:3),p4(0:3)

  ! local variables
  integer :: nhel(nexternal,ncomb),ntry
  real :: t
  real :: qqff_qcd
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

  sqqff_qcd = 0d0
  ntry=ntry+1
  do ihel=1,ncomb
  !          if (goodhel(ihel) .or. ntry .lt. 10) then
    t=qqff_qcd(iq,p1, p2, p3, p4,lam3,lam4,nhel(1,ihel))
    sqqff_qcd = sqqff_qcd + t
  !              if (t .gt. 0d0 .and. .not. goodhel(ihel)) then
  !                  goodhel(ihel)=.true.
  !                  write(*,*) ihel,t
  !              endif
  !          endif
  enddo
  sqqff_qcd = sqqff_qcd /  4d0
end function sqqff_qcd
      
      
function qqff_qcd(iq,p1, p2, p3, p4,lam3,lam4,nhel)

  ! function generated by madgraph
  ! returns amplitude squared summed/avg over colors
  ! for the point in phase space p1,p2,p3,p4,...
  ! and helicity nhel(1),nhel(2),....

  ! for process : u u~  -> t t~

  use quantum_field_theory, only: g, gg, fmass, fwidth
  
  implicit none

  ! constants
  real :: qqff_qcd
  integer :: iq ! incoming fermion species (see initialise for IDs)
  integer ::    ngraphs,    neigen,    nexternal
  parameter (ngraphs=  1,neigen=  1,nexternal=4)
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

! color data

  data eigen_val(1  )/       2.2222222222222221e-01 /
  data eigen_vec(1  ,1  )/  -1.0000000000000000e+00 /
! ----------
! begin code
! ----------
  if((nhel(3) == lam3) .AND. (nhel(4) == lam4))then
    continue
  else
    qqff_qcd = 0.d0
    return
  end if

  call ixxxxx(p1  ,fmass(iq ),nhel(1  ), 1,w1  )
  call oxxxxx(p2  ,fmass(iq ),nhel(2  ),-1,w2  )
  call oxxxxx(p3  ,fmass(11 ),nhel(3  ), 1,w3  )
  call ixxxxx(p4  ,fmass(11 ),nhel(4  ),-1,w4  )
  call jioxxx(w1  ,w2  ,gg,zero,zero,w5  )
  call iovxxx(w4  ,w3  ,w5  ,gg,amp(1  ))
  qqff_qcd = 0.d0
  do i = 1, neigen
    ztemp = (0.d0,0.d0)
    do j = 1, ngraphs
      ztemp = ztemp + eigen_vec(j,i)*amp(j)
    enddo
    qqff_qcd =qqff_qcd+ztemp*eigen_val(i)*conjg(ztemp)
  enddo
!      call gaugecheck(amp,ztemp,eigen_vec,eigen_val,ngraphs,neigen)
end function qqff_qcd
