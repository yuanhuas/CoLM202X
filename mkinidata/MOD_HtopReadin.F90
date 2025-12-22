#include <define.h>

MODULE MOD_HtopReadin

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: HTOP_readin

CONTAINS

   SUBROUTINE HTOP_readin (dir_landdata, lc_year)

! ===========================================================
! Read in the canopy tree top height
! ===========================================================

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Vars_Global
   USE MOD_Const_LC
   USE MOD_Const_PFT
   USE MOD_Vars_TimeInvariants
   USE MOD_LandPatch
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_LandPFT
   USE MOD_Vars_PFTimeInvariants
   USE MOD_Vars_PFTimeVariables
#endif
   USE MOD_NetCDFVector
#ifdef SinglePoint
   USE MOD_SingleSrfdata
#endif

   IMPLICIT NONE

   integer, intent(in) :: lc_year    ! which year of land cover data used
   character(len=256), intent(in) :: dir_landdata

   ! Local Variables
   character(len=256) :: c
   character(len=256) :: landdir, lndname, cyear
   integer :: i,j,t,p,ps,pe,m,n,npatch

   real(r8), allocatable :: htoplc  (:)
   real(r8), allocatable :: htoppft (:), hbotpft(:), cratio_pft(:)

      write(cyear,'(i4.4)') lc_year
      landdir = trim(dir_landdata) // '/htop/' // trim(cyear)

! ==============USGS scheme===========================
#ifdef LULC_USGS

      IF (p_is_worker) THEN
         DO npatch = 1, numpatch
            m = patchclass(npatch)

            htop(npatch) = htop0(m)
            hbot(npatch) = hbot0(m)

         ENDDO
      ENDIF

#endif
! ==============IGBP scheme===========================
#ifdef LULC_IGBP
#ifdef SinglePoint
      allocate (htoplc (numpatch))
      htoplc(:) = SITE_htop
#else
      lndname = trim(landdir)//'/htop_patches.nc'
      CALL ncio_read_vector (lndname, 'htop_patches', landpatch, htoplc)
#endif

      IF (p_is_worker) THEN
         DO npatch = 1, numpatch
            m = patchclass(npatch)

            htop(npatch) = htop0(m)
            hbot(npatch) = hbot0(m)

            ! trees or woody savannas
            IF ( m<6 .or. m==8 ) THEN
               ! 01/06/2020, yuan: adjust htop reading
               ! 11/15/2021, yuan: adjust htop setting
               htop(npatch) = max(2., htoplc(npatch))
               hbot(npatch) = htoplc(npatch)*hbot0(m)/htop0(m)
               hbot(npatch) = max(1., hbot(npatch))
            ENDIF

         ENDDO
      ENDIF

      IF (allocated(htoplc))   deallocate ( htoplc )
#endif

! ==============PFT/PC scheme===========================
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
#ifdef SinglePoint
! 单点模拟时，htoppft来自MOD_SingleSrfdata.F90中定义读取的SITE_htop_pfts变量
      IF (numpft > 0) THEN
         allocate(htoppft(numpft))
         htoppft = pack(SITE_htop_pfts, SITE_pctpfts > 0.)! only keep those with pct>0
      ENDIF
#else
      lndname = trim(landdir)//'/htop_pfts.nc'
      CALL ncio_read_vector (lndname, 'htop_pfts', landpft,   htoppft)
      lndname = trim(landdir)//'/hbot_pfts.nc'
      CALL ncio_read_vector (lndname, 'hbot_pfts', landpft,   hbotpft)
      lndname = trim(landdir)//'/cratio_pfts.nc'
      CALL ncio_read_vector (lndname, 'cratio_pfts', landpft,   cratio_pft)
#endif

      IF (p_is_worker) THEN
         DO npatch = 1, numpatch
            t = patchtype(npatch) ! land cover type, 0: soil, 1: urban, 2: wetland, 3: ice, 4: lake
            m = patchclass(npatch) ! index of land cover type

            IF (t == 0) THEN ! land cover type is soil, growing with vegetation
               ps = patch_pft_s(npatch)
               pe = patch_pft_e(npatch)

               DO p = ps, pe
                  n = pftclass(p) ! pft index

                  htop_p(p) = htop0_p(n) ! for non-tree, use the default height
                  hbot_p(p) = hbot0_p(n)
                  cratio_p(p) = 1.

                  ! for trees
                  ! 01/06/2020, yuan: adjust htop reading
                  ! 11/15/2021, yuan: adjust htop setting
                  ! 12/16/2025, xiangjy: adjust hbot reading
                  IF ( n>0 .and. n<9 ) THEN ! pft is tree
                     htop_p(p) = max(2., htoppft(p))
                     
                     !! ----case0 -----
                     ! hbot_p(p) = htoppft(p)*hbot0_p(n)/htop0_p(n) ! calcu with ratio
                     ! hbot_p(p) = max(1., hbot_p(p))
                     ! cratio_p(p) = 1.

                     !! ----case1 -----
                     ! hbot_p(p) = max(1., hbotpft(p))
                     ! cratio_p(p) = 1.

                     !! ----case2 -----
                     ! hbot_p(p) = htoppft(p)*hbot0_p(n)/htop0_p(n) ! calcu with ratio
                     ! hbot_p(p) = max(1., hbot_p(p))
                     ! cratio_p(p) = max(0., cratio_pft(p))
                     
                     !! ----case3 -----
                     hbot_p(p) = max(1., hbotpft(p))
                     cratio_p(p) = max(0., cratio_pft(p))

                  ENDIF
                  print*, '=====ipatch = ', npatch, '/', numpatch
                  print*, 'patch_pft_s(ipatch): ', ps
                  print*, 'patch_pft_e(ipatch): ', pe
                  print*, '-----pft index: ', n
                  print*, 'htop_p(p): ', htop_p(p)
                  print*, 'hbot_p(p): ', hbot_p(p)
                  print*, 'cratio_p(p): ', cratio_p(p)
               ENDDO

               htop(npatch) = sum(htop_p(ps:pe)*pftfrac(ps:pe)) ! weighted average
               hbot(npatch) = sum(hbot_p(ps:pe)*pftfrac(ps:pe)) ! but in fact, only htop_p & hbot_p used in 3DCanopyRadiation
               cratio(npatch) = sum(cratio_p(ps:pe)*pftfrac(ps:pe))
            ELSE ! if land cover type is not soil, use the default height
               htop(npatch) = htop0(m)
               hbot(npatch) = hbot0(m)
               cratio(npatch) = 1.
            ENDIF

         ENDDO
      ENDIF

      IF (allocated(htoppft)) deallocate(htoppft)
      IF (allocated(hbotpft)) deallocate(hbotpft)
      IF (allocated(cratio_pft)) deallocate(cratio_pft)
#endif

   END SUBROUTINE HTOP_readin

END MODULE MOD_HtopReadin
