!> @file
!> @brief Contains MODULE THEORYWAVES.
!>
!> @author Paul Hall @date 6-May-2024
!>
#include "w3macros.h"
!/ ------------------------------------------------------------------- /
!>
!> @brief Contains wave model subroutine, twmodel.
!>
!> @author Paul Hall @date 6-May-2024 
!>

MODULE THEORYWAVES
  !/
  !/                  +-----------------------------------+
  !/                  | THEORYWAVES                       |
  !/                  |           Paul Hall               |
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         06-May-2024 |
  !/                  +-----------------------------------+
  !
  !  1. Purpose :
  !
  !  2. Variables and types :
  !
  !  3. Subroutines and functions :
  !
  !      Name      Type  Scope    Description
  !     ----------------------------------------------------------------
  !      W3WAVE    Subr. Public   Actual wave model.
  !     ----------------------------------------------------------------
  !
  !  4. Subroutines and functions used :
  !
  !      Name      Type  Module   Description
  !     ----------------------------------------------------------------
  !      W3SETx    Subr. W3xDATMD Point to data structure.
  !     ----------------------------------------------------------------
  !
  !  5. Remarks : Call to W3NMIN hidden behind W3_DEBUGRUN. This call
  !               currently only serves to warn when one or more procs
  !               have no active seapoints. It has been hid as this
  !               dramatically increases runtime performance.
  !
  !  6. Switches :
  !
  !       !/SHRD  Switch for shared / distributed memory architecture.
  !       !/DIST  Id.
  !       !/MPI   Id.
  !       !/OMPG  Id.
  !
  !       !/PR1   First order propagation schemes.
  !       !/PR2   ULTIMATE QUICKEST scheme.
  !       !/PR3   Averaged ULTIMATE QUICKEST scheme.
  !       !/SMC   UNO2 scheme on SMC grid.
  !
  !       !/S     Enable subroutine tracing.
  !       !/T     Test output.
  !       !/MPIT  Test output for MPI specific code.
  !
  !  7. Source code :
  !
  !/ ------------------------------------------------------------------- /
  use w3parall     , only : init_get_isea
  !module default
  implicit none
  !
  PUBLIC
  !/
CONTAINS
  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Run WAVEWATCH III for a given time interval.
  !>
  !> @details Currents are updated before winds as currents are used in wind
  !> and USTAR processing.
  !>
  !> Ice and water levels can be updated only once per call.
  !>
  !> If ice or water level time are undefined, the update
  !> takes place asap, otherwise around the "half-way point"
  !> between the old and new times.
  !>
  !> To increase accuracy, the calculation of the intra-spectral
  !> propagation is performed in two parts around the spatial propagation.
  !>
  !> @param[in] IMOD      Model number.
  !> @param[in] TEND      Ending time of integration.
  !> @param[in] STAMP     Write time stamp (optional, defaults to T).
  !> @param[in] NO_OUT    Skip output (optional, defaults to F).
  !> @param[in] ODAT
  !> @param[in] ID_LCOMM  Present only when using W3_OASIS.
  !> @param[in] TIMEN     Present only when using W3_OASIS.
  !>
  !> @author H. L. Tolman  @date 22-Mar-2021
  !>
  SUBROUTINE TWMODEL ( IMOD )
    !/
    !/                  +-----------------------------------+
    !/                  | THEORYWAVES                       |
    !/                  |           Paul Hall               |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         06-May-2024 |
    !/                  +-----------------------------------+
    !/
    !/    17-Mar-1999 : Distributed FORTRAN 77 version.     ( version 1.18 )
    !  1. Purpose :
    !
    !     Run WAVEWATCH III for a given time interval.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !       IMOD    Int.   I   Model number.
    !     ----------------------------------------------------------------
    !
    !     Local parameters : Flags
    !     ----------------------------------------------------------------
    !       FLOUTG  Log.  Flag for running W3OUTG.
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !     See module documentation.
    !
    !  5. Called by :
    !
    !     Any program shell or integrated model which uses WAVEWATCH III.
    !
    !  6. Error messages :
    !
    !  7. Remarks :
    !     - Currents are updated before winds as currents are used in wind
    !       and USTAR processing.
    !     - Ice and water levels can be updated only once per call.
    !     - If ice or water level time are undefined, the update
    !       takes place asap, otherwise around the "half-way point"
    !       betweem the old and new times.
    !     - To increase accuracy, the calculation of the intra-spectral
    !       propagation is performed in two parts around the spatial propagation.
    !
    !  8. Structure :
    !
    !     -----------------------------------------------------------
    !       0.  Initializations
    !         a Point to data structures
    !       1.  Check the consistency of the input.
    !         a Ending time versus initial time.
    !       2.  Determine next time from ending and output
    !           time and get corresponding time step.
    !       3.  Loop over time steps (see below).
    !       4.  Perform output to file if requested.
    !         a Check if time is output time.
    !        -------------- loop over output types ------------------
    !         d Perform output.                           ( W3IOxx )
    !         e Update next output time.
    !        -------------------- end loop --------------------------
    !       5.  Update log file.
    !       6.  If time is not ending time, branch back to 2.
    !     ----------------------------------------------------------- 
    use CONSTANTS     , only : GRAV, PI
    use w3gdatmd      , only : nseal, mapsf, MAPSTA, USSPF, NK, w3setg
    use w3idatmd      , only : HSL, UWX, UWY
    use w3idatmd      , only : WX0, WY0
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    INTEGER, INTENT(IN)           :: IMOD
    !/
    !/ ------------------------------------------------------------------- /
    !/ Local parameters :
    !/

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
    real, parameter         :: ZERO = 0.
    real, parameter         :: ONE = 1.
    integer                 :: n, jsea, isea, ix, iy, ib
    real                    :: sww, langmt, lasl, laslpj, alphal
    real                    :: u10, ustar

!      do jsea=1, nseal_cpl
      do jsea=1, nseal
        call init_get_isea(isea, jsea)
        ix  = mapsf(isea,1)
        iy  = mapsf(isea,2)
        if (mapsta(iy,ix) == 1) then 
           
!           u10 = SQRT((WX0(ix,iy)**2)+(WY0(ix,iy)**2))
!           ustar = SQRT((UWX(ix,iy)**2)+(UWY(ix,iy)**2))
!           sw_lamult(jsea) = EFactor_model(u10,ustar,HSL(ix,iy))
           LAMULT(jsea) = 4.2
        else 
           LAMULT(jsea)  = 1. 
        endif
      enddo
    end if

!EOC

!  function EFactor_model(u10, ustar, hbl)
!
!! This function returns the enhancement factor, given the 10-meter
!! wind (m/s), friction velocity (m/s) and the boundary layer depth (m).
!!
!! Qing Li, 160606
!
!! Input
!!    real(tw_r8), intent(in) :: &
!!        ! 10 meter wind (m/s)
!!        u10, &
!!        ! water-side surface friction velocity (m/s)
!!        ustar, &
!!        ! boundary layer depth (m)
!!        hbl
!    REAL, intent(in) :: u10, ustar, hbl
!
!! Local variables
!!    real(tw_r8) :: us_sl, lasl_sqr_i
!!    real(tw_r8) :: EFactor_model
!    REAL :: us_sl, lasl_sqr_i
!    REAL :: EFactor_model
!
!    if (u10 .gt. ZERO .and. ustar .gt. ZERO) then
!      ! surface layer averaged Stokes drift
!      us_sl = ustokes_SL_model(u10, hbl)
!      !
!      ! LaSL^{-2}
!      lasl_sqr_i = us_sl/ustar
!      !
!      ! enhancement factor (Li et al., 2016)
!      EFactor_model = sqrt(ONE &
!                 +ONE/1.5**2*lasl_sqr_i &
!                 +ONE/5.4**4*lasl_sqr_i**2)
!    else
!      ! otherwise set to one
!      EFactor_model = ONE
!    endif
!
!  end function EFactor_model
!  function ustokes_SL_model(u10, hbl)
!
!! This function returns the surface layer averaged Stokes drift, given
!! the 10-meter wind (m/s) and the boundary layer depth (m).
!!
!! Qing Li, 20180130
!
!! Input
!!    real(tw_r8), intent(in) :: &
!!        ! 10 meter wind (m/s)
!!        u10, &
!!        ! boundary layer depth (m)
!!        hbl
!    REAL :: u10, hbl
!! Local variables
!    ! parameters
!    REAL, parameter :: &
!        ! ratio of U19.5 to U10 (Holthuijsen, 2007)
!        u19p5_to_u10 = 1.075, &
!        ! ratio of mean frequency to peak frequency for
!        ! Pierson-Moskowitz spectrum (Webb, 2011)
!        fm_to_fp = 1.296, &
!        ! ratio of surface Stokes drift to U10
!        us_to_u10 = 0.0162, &
!        ! loss ratio of Stokes transport
!        r_loss = 0.667
!
!!    real(tw_r8) :: us, hm0, fm, fp, vstokes, kphil, kstar
!!    real(tw_r8) :: z0, z0i, r1, r2, r3, r4, tmp
!!    real(tw_r8) :: ustokes_SL_model
!    REAL :: us, hm0, fm, fp, vstokes, kphil, kstar
!    REAL :: z0, z0i, r1, r2, r3, r4, tmp
!    REAL :: ustokes_SL_model
!
!    if (u10 .gt. ZERO) then
!      ! surface Stokes drift
!      us = us_to_u10*u10
!      !
!      ! significant wave height from Pierson-Moskowitz
!      ! spectrum (Bouws, 1998)
!      hm0 = 0.0246*u10**2
!      !
!      ! peak frequency (PM, Bouws, 1998)
!      tmp = 2.0*PI*u19p5_to_u10*u10
!      fp = 0.877*GRAV/tmp
!      !
!      ! mean frequency
!      fm = fm_to_fp*fp
!      !
!      ! total Stokes transport (a factor r_loss is applied to account
!      !  for the effect of directional spreading, multidirectional waves
!      !  and the use of PM peak frequency and PM significant wave height
!      !  on estimating the Stokes transport)
!      vstokes = 0.125*PI*r_loss*fm*hm0**2
!      !
!      ! the general peak wavenumber for Phillips' spectrum
!      ! (Breivik et al., 2016) with correction of directional spreading
!      kphil = 0.176*us/vstokes
!      !
!      ! surface layer averaged Stokes dirft with Stokes drift profile
!      ! estimated from Phillips' spectrum (Breivik et al., 2016)
!      ! the directional spreading effect from Webb and Fox-Kemper, 2015
!      ! is also included
!      kstar = kphil*2.56
!      ! surface layer
!      z0 = 0.2*abs(hbl)
!      z0i = ONE/z0
!      ! term 1 to 4
!      r1 = (0.151/kphil*z0i-0.84) &
!            *(ONE-exp(-2.0*kphil*z0))
!      r2 = -(0.84+0.0591/kphil*z0i) &
!             *sqrt(2.0*PI*kphil*z0) &
!             *erfc(sqrt(2.0*kphil*z0))
!      r3 = (0.0632/kstar*z0i+0.125) &
!            *(ONE-exp(-2.0*kstar*z0))
!      r4 = (0.125+0.0946/kstar*z0i) &
!             *sqrt(2.0*PI*kstar*z0) &
!             *erfc(sqrt(2.0*kstar*z0))
!      ustokes_SL_model = us*(0.715+r1+r2+r3+r4)
!    else
!      ustokes_SL_model = ZERO
!    endif
!
!    end function ustokes_SL_model
!
END MODULE THEORYWAVES
