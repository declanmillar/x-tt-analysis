module class_Histogram
  use Configuration, only: include_errors
  use kinematics, only: sigma
  implicit none
  private

  integer :: include_errors = 1
  integer :: it = 20
  real :: sigma = 1
  real :: cnorm(20)
  real :: diff_max = 1e-12

  type, public :: Histogram
    private
    character (len = 50), public :: name
    character (len = 50), public :: ytitle
    character (len = 50), public :: xtitle
    real, public :: low
    real, public :: up
    integer, public :: ndiv

    real, private :: binw = 0
    real :: x(500)=0,fx(500,20)=0,fxtot(500)=0
    real :: sfxtot=0
    real :: sumw2(500,20)=0, sumw2tot(500)=0

    contains
      procedure :: bin_width => histogram_bin_width 
      procedure :: midpoints => histogram_midpoints
      procedure :: initialise => histogram_initialise
      procedure :: finalise => histogram_finalise
      procedure :: fill => histogram_fill
      procedure :: collate => histogram_collate
      procedure :: check => histogram_check
      procedure :: print => histogram_print
  end type Histogram

  interface Histogram
    module procedure new_Histogram
  end interface Histogram
contains

  function new_Histogram(this)
    class(Histogram), intent(inout) :: this
    type(Histogram) new_Histogram
  end function new_Histogram

  subroutine histogram_initialise(this)
    class(Histogram), intent(inout) :: this
    call this%bin_width
    call this%midpoints
  end subroutine histogram_initialise

  subroutine histogram_finalise(this)
    class(Histogram), intent(inout) :: this
    call this%collate
    call this%check
    call this%print
  end subroutine histogram_finalise

  subroutine histogram_bin_width(this)
    class(Histogram), intent(inout) :: this
    this%binw = (this%up-this%low)/this%ndiv
  end subroutine histogram_bin_width

  subroutine histogram_midpoints(this)
    class(Histogram), intent(inout) :: this
    integer :: i
    do i = 1, this%ndiv
      this%x(i) = this%low + this%binw * (i - 1) + this%binw / 2.d0
    end do
  end subroutine histogram_midpoints

  subroutine histogram_fill(this, value, weight)
    class(Histogram), intent(inout) :: this
    real :: value, weight
    integer :: nbin
    integer :: i

    nbin = int((value - this%low)/this%binw) + 1
    if (nbin >= (this%ndiv + 1)) then
      continue
    else if (nbin < 1) then
      continue
    else
      this%fx(nbin,it) = this%fx(nbin,it) + weight
      if (include_errors == 1) this%sumw2(nbin, it) = this%sumw2(nbin, it) + weight*weight
    end if
  end subroutine histogram_fill

  subroutine histogram_collate(this)
    class(Histogram), intent(inout) :: this
    integer :: i, j

    do i = 1, 20
      cnorm(i) = 1
    end do

    do j = 1, this%ndiv
      this%fxtot(j) = 0.d0
      do i = 1, it
        this%fx(j,i) = this%fx(j,i)*sigma/cnorm(i)/this%binw
        this%fxtot(j) = this%fxtot(j)+this%fx(j,i)
        if (include_errors == 1) then 
          this%sumw2(j,i) = this%sumw2(j,i)*sigma/cnorm(i)/this%binw*sigma/cnorm(i)/this%binw
          this%sumw2tot(j) = this%sumw2tot(j)+this%sumw2(j,i)
        else
          this%sumw2tot(j) = 0
        end if
      end do
      this%sfxtot = this%sfxtot + this%fxtot(j)*this%binw
    end do
  end subroutine histogram_collate

  subroutine histogram_check(this)
    class(Histogram), intent(inout) :: this
      if(abs(sigma-this%sfxtot) > diff_max) then
        write(*,*) this%name, ' error: ', this%sfxtot
      end if
  end subroutine histogram_check

  subroutine histogram_print(this)
    class(Histogram), intent(inout) :: this
    integer :: i
    write(*,*)'DISTRIBUTION'
    write(*,*) this%name
    write(*,*) this%ytitle
    write(*,*) this%xtitle
    write(*,*) this%low
    write(*,*) this%up
    write(*,*) this%binw
    do i = 1, 10
      write(*,*) this%x(i), this%fxtot(i)
    end do
    write(*,*)'END'
  end subroutine histogram_print

end module class_Histogram

program histogram_test
  use class_Histogram
  implicit none

  integer i

  type(Histogram) :: hist1
  type(Histogram) :: hist2

  hist1 = Histogram("hist1","y","x",0,10,10)
  hist2 = Histogram("hist2","y","x",0,20,10)

  call hist1%initialise()
  call hist2%initialise()

  call hist1%fill(5.0, 1.0)
  call hist2%fill(3.3, 1.0)

  call hist1%finalise()
  call hist2%finalise()       
end program histogram_test