!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

MODULE THEORYWAVES

!BOP
!\newpage
! !MODULE: theorywaves
!
! !AUTHOR:
!  
!
! !DESCRIPTION:
!\\
!\\
!  References:\\
!\\
!\\

! !USES:
!  CONSTANTS
!EOP
USE CONSTANTS, ONLY: GRAV, PI

  implicit none
!  private
!  save

!BOP

! !DEFINED PARAMETERS:

  ! Kind Types:
  ! Use double precision for floating point computations.
!  integer, parameter :: tw_r8       = selected_real_kind(15, 307)

  ! Global parameters:
  ! The constant 1 is used repeatedly. 
  ! The value for pi is needed.
!  real(tw_r8), parameter :: tw_zero = real(0,tw_r8),         &
!                            tw_one  = real(1,tw_r8)
!  real(tw_r8), parameter :: PI      = &
!                               3.14159265358979323846_tw_r8
!  real(tw_r8), parameter :: Gravity = &
!                               9.80616_tw_r8

END MODULE THEORYWAVES
