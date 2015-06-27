module configuration

  use integration, only: seed, itmx, ncall, acc

  implicit none

  ! parameters in config file
  real    :: collider_energy
  integer :: initial_state ! 0 = pp, 1 = ppbar
  integer :: final_state  ! 0 = no decay, 1 = dilepton, 2 = semi-had, 3 = full-had
  
  integer :: structure_function
  character(50) :: model_name
  character(100) :: ntuple_file
  character(50) :: output_file
  integer :: include_qcd
  integer :: include_ew
  integer :: include_bsm
  integer :: include_gg
  integer :: include_qq
  integer :: interference  
  integer :: use_nwa
  integer :: use_branching_ratio ! 0 = no, 1 = dilepton, 2 = semi-had, 3 = full-had
  integer :: additional_kinematics
  integer :: symmetrise_x1x2
  integer :: symmetrise_costheta_t
  integer :: symmetrise_costheta_5
  integer :: symmetrise_costheta_7
  integer :: print_distributions
  integer :: include_errors
  integer :: phase_space_only
  integer :: use_rambo
  integer :: map_phase_space
  integer :: verbose
  integer :: nfb
  real :: lambdaqcd4
  integer :: nloops
  real :: ecm_low, ecm_up

  ! derived parameters
  integer :: n_final
  integer :: tops_decay
  integer :: ffinal

  ! distributions in asymmetries
  integer, parameter :: n_asymmetries = 12
  integer :: o_asym(n_asymmetries)

  integer :: ixmax, jxmax, i5max, i7max
  
  public :: read_config
  public :: modify_config

  contains

    subroutine read_config

      ! read config file
      print*, "Reading config file..."

      read(5,"(a)") ntuple_file

      read(5,*) output_file
      
      read(5,*) initial_state ! 0 = pp, 1 = ppbar

      read(5,*) final_state ! 1 = no decay, 1 = dilepton, 2 = semilepton, 4 = full hadron
      
      read(5,*) model_name
      
      read(5,*) structure_function

      read(5,*) include_qcd

      read(5,*) include_ew

      read(5,*) include_bsm

      read(5,*) include_gg

      read(5,*) include_qq

      read(5,*) phase_space_only

      read(5,*) interference

      read(5,*) use_nwa ! 0:actual top widths,1: tops in nwa

      read(5,*) collider_energy

      read(5,*) seed

      read(5,*) itmx

      read(5,*) ncall

      read(5,*) acc

      read(5,*) use_rambo

      read(5,*) map_phase_space

      read(5,*) symmetrise_x1x2

      read(5,*) symmetrise_costheta_t

      read(5,*) symmetrise_costheta_5

      read(5,*) symmetrise_costheta_7

      read(5,*) verbose

      read(5,*) ecm_low

      read(5,*) ecm_up

      print*, "...complete."
      
    end subroutine read_config

    subroutine modify_config

      print*, "Modifying config file..."

      if (final_state <= 0) then
        symmetrise_costheta_5 = 0
        symmetrise_costheta_7 = 0
      end if

      if (final_state <= 0) then
        n_final = 4
      else
        n_final = 8
      end if

      if (final_state == -1) then
        ffinal = 1
        include_qcd = 0
      else
        ffinal = 11
      end if

      if (use_rambo == 1) then
        print*, "Setting additional_kinematics to zero."
        map_phase_space = 0
        additional_kinematics = 0
      end if

      if (final_state <= 0) use_nwa = 0
      if (final_state > 0) use_branching_ratio = 0

      if (itmx > 20 )then
        write(*,*) 'itmx does not have to exceed 20!'
        stop
      end if

      if (symmetrise_x1x2 == 1)then
        ixmax = 2
      else
        ixmax = 1
      end if

      if (symmetrise_costheta_t == 1)then
        jxmax = 2
      else
        jxmax = 1
      end if

      if (symmetrise_costheta_5 == 1)then
        i5max = 2
      else
        i5max = 1
      end if

      if (symmetrise_costheta_7 == 1)then
        i7max = 2
      else
        i7max = 1
      end if

      if (final_state <= 0)then
        tops_decay = 0
      else
        tops_decay = 1
      end if

      if (phase_space_only == 1)then
        include_qcd = 0
        include_ew = 0
        include_bsm = 0
      end if
      
      print*, "...complete."
    end subroutine modify_config

    subroutine debug(message)
      character(*) :: message
      if (verbose == 1) print*, message
    end subroutine debug

end module configuration
