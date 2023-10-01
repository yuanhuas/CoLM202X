#include <define.h>

MODULE MOD_Urban_Irrigation
! -----------------------------------------------------------------------
! !DESCRIPTION:
! Calculate irrigation water use in the urban vegetated fraction
!
! Created by Zhuo Liu， 09/2023
!
! !REVISIONS:
!
! -----------------------------------------------------------------------
! !USE
 USE MOD_Precision
 USE MOD_TimeManager
 USE MOD_Namelist
 USE MOD_Const_Physical
 USE MOD_Vars_Global, only:PI

 IMPLICIT NONE
 SAVE

 PUBLIC :: UrbanIrrigation

CONTAINS

! Identify urban irrigation triggers and irrigation period
 SUBROUTINE UrbanIrrigationCheck( &
        idate          ,fveg           ,deltim         ,lai         ,&
        patchlonr      ,urb_irrig_trigger              ,urb_irrig_check)

   integer, intent(in) ::     &
        idate(3)               ! calendar (year, julian day, seconds)

   real(r8), intent(in) ::    &
        fveg,                 &! fraction of vegetation cover
        deltim,               &!
        lai,                  &!
        patchlonr              !

   logical, intent(out) ::    &
        urb_irrig_trigger,    &!
        urb_irrig_check        !

!-----------------------Local Variables------------------------------
   REAL(r8) :: &
        ldate(3),             &! patitial volume of ice and water of surface layer
        seconds_start_time,   &!
        londeg                 !

   INTEGER  :: &
        urb_irrig_start_month,&!
        urb_irrig_end_month,  &!
        month,                &!
        urb_irrig_start_time, &!
        urb_irrig_time,       &!
        day                    !

!Months when urban irrigation is turned on and off
        urb_irrig_start_month = 1
        urb_irrig_end_month   = 12
        urb_irrig_start_time  = 1

!Irrigation hours
        urb_irrig_time        = 4

!Identify fveg > 0 and lai > 0

   !TODO: reason?
   IF (fveg > 1.e-5 .and. lai > 1.e-5) THEN
      CALL julian2monthday(idate(1), idate(2), month, day)
!month
      IF (month >= urb_irrig_start_month .and. month <= urb_irrig_end_month) THEN
!hour
         IF (DEF_simulation_time%greenwich) THEN
            ! convert GMT time to local time
            londeg = patchlonr*180/PI
            CALL gmt2local(idate, londeg, ldate)
         ENDIF

         IF (DEF_simulation_time%greenwich) THEN
            seconds_start_time = ldate(3) - deltim - urb_irrig_start_time*3600.
         ELSE
            seconds_start_time = idate(3) - deltim - urb_irrig_start_time*3600.
         ENDIF
       ! identify trigger
         IF ((seconds_start_time >= 0._r8) .and. (seconds_start_time < deltim)) THEN
            urb_irrig_trigger = .true.
         ELSE
            urb_irrig_trigger = .false.
         ENDIF
       ! identify irrigation period
         IF ((seconds_start_time >= 0._r8) .and. (seconds_start_time < (deltim + urb_irrig_time*3600.))) THEN
            urb_irrig_check = .true.
         ELSE
            urb_irrig_check = .false.
         ENDIF

      ELSE
         urb_irrig_check = .false.
      ENDIF
   ELSE
      urb_irrig_check = .false.
   ENDIF

 END SUBROUTINE UrbanIrrigationCheck


! the amount of total irrigation was accounted for by calculating the soil moisture deficit
 SUBROUTINE UrbanIrrigationRequire(lbp  , nl_soil  ,  dz_gpersno   ,z_gpersno   ,t_gpersno,    porsl,&
                                   psi0 , bsw      ,  wliq_gpersno ,fveg  ,&
                                   urb_irrig_deficit)
   IMPLICIT NONE

   integer, intent(in) :: &
        lbp                     ,    &!
        nl_soil                       !

   real(r8), intent(in) :: &
        dz_gpersno   (lbp:nl_soil),  &! layer thickness [m]
        t_gpersno    (lbp:nl_soil),  &! soil/snow skin temperature [K]
        z_gpersno    (lbp:nl_soil),  &!
        porsl        (1:nl_soil),    &! soil porosity [-]
        psi0         (1:nl_soil),    &! saturated soil suction [mm] (NEGATIVE)
        bsw          (1:nl_soil),    &! clapp and hornbereger "b" parameter [-]
        wliq_gpersno (lbp:nl_soil),  &! liquid water [kg/m2]
        fveg                          !

   real(r8), intent(out) :: &
        urb_irrig_deficit             ! [mm]

!-----------------------Local Variables------------------------------
   INTEGER  :: &
        i                                  !
   REAL(r8) :: &
        urb_irrig_wliq_tot               ,&! [kg/m2]
        urb_irrig_wliq_target_tot        ,&! [kg/m2]
        urb_irrig_wliq_wilting_point_tot ,&! [kg/m2]
        urb_irrig_wliq_field_point_tot   ,&! [kg/m2]
        wliq_wilting_point(1:nl_soil)    ,&! [kg/m2]
        wliq_field_capacity(1:nl_soil)   ,&! [kg/m2]
        dzmm(lbp:nl_soil)                ,&! [mm]
        urb_irrig_threshold              ,&! [kg/m2]
        urb_irrig_threshold_fraction     ,&! [-]
        urb_irrig_max_depth                ! [m]

!  initialize local variables
        wliq_wilting_point(1:nl_soil)        = 0._r8
        wliq_field_capacity(1:nl_soil)       = 0._r8
        urb_irrig_wliq_tot                   = 0._r8
        urb_irrig_wliq_target_tot            = 0._r8
        urb_irrig_wliq_wilting_point_tot     = 0._r8
        urb_irrig_threshold                  = 0._r8
        urb_irrig_deficit                    = 0._r8

!  set irrigation depth and threshold fraction
        urb_irrig_max_depth                  = 0.6
        urb_irrig_threshold_fraction         = 1._r8

!  calculate wilting point and field capacity
   DO i = 1, nl_soil
      !write(*,*) porsl(i)
      IF (t_gpersno(i)>tfrz .and. porsl(i)>=1.e-6) THEN
         !dzmm(i)                = dz_gpersno(i) * 1000.
         wliq_wilting_point(i)  = dz_gpersno(i)*porsl(i)*((-1.5e5/psi0(i))**(-1./bsw(i)))*denh2o
         wliq_field_capacity(i) = dz_gpersno(i)*porsl(i)*((-3399./psi0(i))**(-1./bsw(i)))*denh2o
      ENDIF
   ENDDO

   ! calculate total irrigation needed in all soil layers
   DO i = 1, nl_soil
      IF (z_gpersno(i)<urb_irrig_max_depth .and. t_gpersno(i)>tfrz) THEN
         urb_irrig_wliq_tot               = urb_irrig_wliq_tot + wliq_gpersno(i)
         urb_irrig_wliq_wilting_point_tot = urb_irrig_wliq_wilting_point_tot + wliq_wilting_point(i)
         urb_irrig_wliq_target_tot        = urb_irrig_wliq_target_tot + wliq_field_capacity(i)
      ENDIF
   ENDDO

   !TODO: need only one time?
   urb_irrig_threshold = urb_irrig_wliq_wilting_point_tot + urb_irrig_threshold_fraction * &
                         (urb_irrig_wliq_target_tot - urb_irrig_wliq_wilting_point_tot)

   ! calculate the amount of water needed
   IF (urb_irrig_wliq_tot < urb_irrig_threshold) THEN
      urb_irrig_deficit = urb_irrig_threshold - urb_irrig_wliq_tot
   ELSE
      urb_irrig_deficit = 0
   ENDIF

 END SUBROUTINE UrbanIrrigationRequire

!convert the calculated irrigation water amount into irrigation rate mm--->mm/s
 SUBROUTINE UrbanIrrigation(lbp         ,nl_soil     ,idate    ,    deltim   ,    fveg ,      lai ,  patchlonr,&
                            dz_gpersno  ,z_gpersno   ,t_gpersno,    porsl    ,    psi0 ,      bsw ,  wliq_gpersno,&
                            urb_irrig_qflx)

   IMPLICIT NONE

   integer, intent(in) :: &
        lbp    ,                    &!
        nl_soil,                    &!
        idate(3)                     ! calendar (year, julian day, seconds)
   real(r8), intent(in) :: &
        deltim                     ,&!
        fveg                       ,&!
        lai                        ,&!
        patchlonr                  ,&!
        dz_gpersno  (lbp:nl_soil)  ,&! layer thickness [m]
        z_gpersno   (lbp:nl_soil)  ,&!
        t_gpersno   (lbp:nl_soil)  ,&! soil/snow skin temperature [K]
        porsl       (1  :nl_soil)  ,&! soil porosity [-]
        psi0        (1  :nl_soil)  ,&! saturated soil suction [mm] (NEGATIVE)
        bsw         (1  :nl_soil)  ,&! clapp and hornbereger "b" parameter [-]
        wliq_gpersno(lbp:nl_soil)    ! liquid water [kg/m2]

   real(r8), intent(out) :: &
        urb_irrig_qflx                  !

!-----------------------Local Variables------------------------------
   INTEGER  :: &
        urb_irrig_time,         &!
        urb_irrig_sum_count,    &!
        urb_irrig_nsteps_per_day !
   REAL(r8) :: &
        urb_irrig_rate,         &!
        urb_irrig_deficit        !
   REAL(r8) ,SAVE :: &
        urb_irrig_rate_const,   &!
        urb_irrig_deficit_left   !

   LOGICAL ::  &
        urb_irrig_trigger,      &!
        urb_irrig_check          !
!  irrigation hours
        urb_irrig_time      = 4
!  initialize variables
        urb_irrig_qflx      = 0._r8
        urb_irrig_rate      = 0._r8

   CALL UrbanIrrigationCheck(idate,fveg,deltim,lai,patchlonr,urb_irrig_trigger,urb_irrig_check)
   ! When irrigation is triggered at this moment, the irrigation amount and irrigation rate are calculated.
   ! At the next moment it is not calculated.

   IF (urb_irrig_trigger) THEN
      CALL UrbanIrrigationRequire(lbp  , nl_soil  ,  dz_gpersno   ,z_gpersno   ,t_gpersno,    porsl,&
                                  psi0 , bsw      ,  wliq_gpersno ,fveg,&
                                  urb_irrig_deficit)
      !urb_irrig_deficit = 10.
      IF (urb_irrig_deficit > 0) THEN
         !TODO: only one time needed.
         urb_irrig_nsteps_per_day = NINT(urb_irrig_time*3600./deltim)
         ! kg/m2 -->m -->mm -->mm/s
         urb_irrig_rate           = (urb_irrig_deficit*1000.)/(denh2o*deltim*urb_irrig_nsteps_per_day)
         urb_irrig_rate_const     = urb_irrig_rate
         urb_irrig_deficit_left   = urb_irrig_deficit
      ENDIF
   ENDIF

   IF ((urb_irrig_check) .and. urb_irrig_deficit_left > 0) THEN
      urb_irrig_qflx          = urb_irrig_rate_const
      urb_irrig_deficit_left  = urb_irrig_deficit_left - urb_irrig_rate_const*deltim
      IF (urb_irrig_deficit_left < 0.) THEN
         urb_irrig_deficit_left = 0
      ENDIF
   ELSE
      urb_irrig_rate_const = 0._r8
   ENDIF

 END SUBROUTINE UrbanIrrigation

END MODULE MOD_Urban_Irrigation
