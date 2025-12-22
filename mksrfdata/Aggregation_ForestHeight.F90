#include <define.h>

SUBROUTINE Aggregation_ForestHeight ( &
      gland, dir_rawdata, dir_model_landdata, lc_year)

!-----------------------------------------------------------------------
!  Global Forest Height
!     (http://lidarradar.jpl.nasa.gov/)
!      Simard, M., N. Pinto, J. B. Fisher, and A. Baccini, 2011: Mapping
!      forest canopy height globally with spaceborne lidar.
!      J. Geophys. Res., 116, G04021.
!
!  Created by Yongjiu Dai, 02/2014
!
! !REVISIONS:
!  Hua Yuan,      ?/2020 : for land cover land use classifications
!  Shupeng Zhang, 01/2022: porting codes to MPI parallel version
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_LandPatch
   USE MOD_Land2mWMO
   USE MOD_NetCDFVector
   USE MOD_NetCDFBlock
#ifdef RangeCheck
   USE MOD_RangeCheck
#endif
   USE MOD_AggregationRequestData
   USE MOD_Utils

   USE MOD_Const_LC
   USE MOD_5x5DataReadin
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_LandPFT
#endif

#ifdef SrfdataDiag
   USE MOD_SrfdataDiag
#endif

   IMPLICIT NONE

   ! arguments:
   integer, intent(in) :: lc_year
   type(grid_type),  intent(in) :: gland
   character(len=*), intent(in) :: dir_rawdata
   character(len=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ---------------------------------------------------------------
   character(len=256) :: landdir, lndname, cyear
   integer :: L, ipatch, p
   integer :: wmo_src

   type (block_data_real8_2d) :: tree_height
   real(r8), allocatable :: tree_height_patches(:), tree_height_one(:)

   ! for IGBP data
   character(len=256) :: dir_5x5, fname
   type (block_data_real8_2d) :: htop
   type (block_data_real8_3d) :: pftPCT,cdepth,hbot,cratio
   real(r8), allocatable :: htop_patches(:), htop_pfts(:), htop_pcs(:,:)
   real(r8), allocatable :: hbot_pfts(:) !, hbot_patches(:)
   real(r8), allocatable :: cratio_pfts(:) !, cratio_patches(:)
   real(r8), allocatable :: htop_one(:), area_one(:), pct_one(:,:)
   real(r8), allocatable :: cdepth_one(:,:), hbot_one(:,:), cratio_one(:,:)
   ! logical, allocatable :: pft_has_area(:)
   integer  :: ip, ipft
   real(r8) :: sumarea, patch_weight_sum, patch_bot_sum

#ifdef SrfdataDiag
   integer :: typpatch(N_land_classification+1), ityp
#ifndef CROP
   integer :: typpft  (N_PFT)
#else
   integer :: typpft  (N_PFT+N_CFT)
#endif
   integer :: typpc   (N_land_classification+1)
#endif

      write(cyear,'(i4.4)') lc_year
      landdir = trim(dir_model_landdata) // '/htop/' //trim(cyear) !构建srfdata输出目录
      ! print*,'landdir: ', trim(landdir)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(/, A)') 'Aggregate forest height ...'
         CALL system('mkdir -p ' // trim(adjustl(landdir)))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef LULC_USGS
      lndname = trim(dir_rawdata)//'/Forest_Height.nc'

      IF (p_is_io) THEN
         CALL allocate_block_data (gland, tree_height)
         CALL ncio_read_block (lndname, 'forest_height', gland, tree_height)

#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = tree_height)
#endif
      ENDIF

      IF (p_is_worker) THEN

         allocate (tree_height_patches (numpatch))

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)

            IF(L/=0 .and. L/=1 .and. L/=16 .and. L/=24)THEN
               ! NOT OCEAN(0)/URBAN and BUILT-UP(1)/WATER BODIES(16)/ICE(24)
               CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, &
                  data_r8_2d_in1 = tree_height, data_r8_2d_out1 = tree_height_one)
               tree_height_patches (ipatch) = median (tree_height_one, size(tree_height_one))
            ELSE
               tree_height_patches (ipatch) = -1.0e36_r8
            ENDIF
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('htop_patches ', tree_height_patches)
#endif

      lndname = trim(landdir)//'/htop_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'htop_patches', 'patch', landpatch, tree_height_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typpatch = (/(ityp, ityp = 0, N_land_classification)/)
      lndname  = trim(dir_model_landdata) // '/diag/htop_patch_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (tree_height_patches, landpatch%settyp, typpatch, m_patch2diag, &
         -1.0e36_r8, lndname, 'htop', compress = 6, write_mode = 'one', defval=0._r8, create_mode=.true.)
#endif

      IF (p_is_worker) THEN
         deallocate ( tree_height_patches )
      ENDIF
#endif


#ifdef LULC_IGBP
      IF (p_is_io) THEN
         CALL allocate_block_data (gland, htop)
      ENDIF

      IF (p_is_io) THEN

         dir_5x5 = trim(dir_rawdata) // trim(DEF_rawdata%htop%dir)
         fname = trim(DEF_rawdata%htop%fname)
         CALL read_5x5_data (dir_5x5, fname, gland, trim(DEF_rawdata%htop%vname), htop)
#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = htop)
#endif
      ENDIF

      IF (p_is_worker) THEN

         allocate (htop_patches (numpatch))

         DO ipatch = 1, numpatch

            IF (landpatch%settyp(ipatch) /= 0) THEN

               CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, &
                  area = area_one, data_r8_2d_in1 = htop, data_r8_2d_out1 = htop_one)

               where (htop_one < 0.) htop_one = 0.

               htop_patches(ipatch) = sum(htop_one * area_one) / sum(area_one)
            ENDIF

         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('HTOP_patches ', htop_patches)
#endif

      lndname = trim(landdir)//'/htop_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'htop_patches', 'patch', landpatch, htop_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typpatch = (/(ityp, ityp = 0, N_land_classification)/)
      lndname  = trim(dir_model_landdata) // '/diag/htop_patch_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (htop_patches, landpatch%settyp, typpatch, m_patch2diag, &
         -1.0e36_r8, lndname, 'htop', compress = 6, write_mode = 'one', defval=0._r8, create_mode=.true.)
#endif

      IF (p_is_worker) THEN
         IF (allocated(htop_patches)) deallocate (htop_patches)
         IF (allocated(htop_one    )) deallocate (htop_one    )
         IF (allocated(area_one    )) deallocate (area_one    )
      ENDIF
#endif

! ================PFT/PC htop
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF (p_is_io) THEN
         CALL allocate_block_data (grid_pft, pftPCT, N_PFT_modis, lb1 = 0) ! (grid, data, 第一维度数量, 第一维度下界默认为1)
         CALL allocate_block_data (grid_pft, cdepth, N_PFT_modis, lb1 = 0)
         CALL allocate_block_data (grid_pft, hbot  , N_PFT_modis, lb1 = 0)
         CALL allocate_block_data (grid_pft, cratio, N_PFT_modis, lb1 = 0)
         CALL allocate_block_data (gland   , htop)
      ENDIF

      IF (p_is_io) THEN

         dir_5x5 = trim(DEF_dir_rawdata) // trim(DEF_rawdata%htop%dir)
         fname = trim(DEF_rawdata%htop%fname)
         print*, 'Reading htop from dir: ', trim(DEF_rawdata%htop%dir)
         print*, 'Reading htop from vname: ', trim(DEF_rawdata%htop%vname)
         print*, 'Reading htop from fname: ', trim(DEF_rawdata%htop%fname)
         CALL read_5x5_data     (dir_5x5, fname, gland   , trim(DEF_rawdata%htop%vname), htop)

         ! 新增读取crown depth
         dir_5x5= trim(DEF_dir_rawdata) // trim(DEF_rawdata%cdepth%dir)
         fname = trim(DEF_rawdata%cdepth%fname)
         CALL read_5x5_data_pft (dir_5x5, fname, grid_pft, trim(DEF_rawdata%cdepth%vname), cdepth, lb=1, ub=8)
         print*, 'Reading cdepth from dir: ', trim(DEF_rawdata%cdepth%dir)
         print*, 'Reading cdepth from vname: ', trim(DEF_rawdata%cdepth%vname)

         ! 新增读取crown aspect ratio
         dir_5x5= trim(DEF_dir_rawdata) // trim(DEF_rawdata%cratio%dir)
         fname = trim(DEF_rawdata%cratio%fname)
         CALL read_5x5_data_pft (dir_5x5, fname, grid_pft, trim(DEF_rawdata%cratio%vname), cratio, lb=1, ub=8)
         print*, 'Reading cratio from dir: ', trim(DEF_rawdata%cratio%dir)
         print*, 'Reading cratio from vname: ', trim(DEF_rawdata%cratio%vname)

         dir_5x5= trim(DEF_dir_rawdata) // trim(DEF_rawdata%pft%dir)
         fname  = trim(DEF_rawdata%pft%fname) // '.' // trim(cyear)
         CALL read_5x5_data_pft (dir_5x5, fname, grid_pft, 'PCT_PFT', pftPCT)
         print*, 'Reading pftPCT from dir: ', trim(DEF_rawdata%pft%dir)
         print*, 'Reading pftPCT from fname: ', trim(DEF_rawdata%pft%fname) // '.' // trim(cyear)

#ifdef USEMPI
         CALL aggregation_data_daemon_multigrid (grid_in1 = grid_pft, data_r8_3d_in1 = pftPCT, n1_r8_3d_in1 = 16, &
            grid_in2 = gland, data_r8_2d_in2 = htop, grid_in3 = grid_pft, data_r8_3d_in3 = cdepth, n1_r8_3d_in3 = 16, &
            grid_in4 = grid_pft, data_r8_3d_in4 = cratio, n1_r8_3d_in4 = 16 )
#endif
      ENDIF

      IF (p_is_worker) THEN

         allocate (htop_patches   (numpatch)) ! numpatch=144, numpft=701
         allocate (htop_pfts      (numpft  ))
         ! allocate (hbot_patches   (numpatch))
         allocate (hbot_pfts      (numpft  ))
         ! allocate (cratio_patches (numpatch))
         allocate (cratio_pfts    (numpft  ))
         
         DO ipatch = 1, numpatch ! 外循环144个patch，内循环16个PFT

            IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
               wmo_src = wmo_source (landpatch%ielm(ipatch))

               ! set patch htop
               htop_patches(ipatch) = htop_patches(wmo_src)
               ! hbot_patches(ipatch) = hbot_patches(wmo_src)
               ! cratio_patches(ipatch) = cratio_patches(wmo_src)
               ! set pft htop at the same time
               ip = patch_pft_s(ipatch)
               htop_pfts(ip) = htop_patches(ipatch)
               ! hbot_pfts(ip) = hbot_patches(wmo_src)
               ! cratio_pfts(ip) = cratio_patches(wmo_src)

               CYCLE
            ENDIF
            ! 获取ipatch的返回数据，降维度，去掉了ipatch这个维度
            ! 从全局数据中提取当前进程负责的patch数据
            ! n_points是这个patch包含的网格点数，每个patch不同
            ! htop(nlat,nlon)——>htop_one(npoints)
            ! cdepth(n_pft=16,nlat,nlon)——>cdepth_one(n_pft,npoints)
            ! pftPCT(n_pft=16,nlat,nlon)——>pct_one(n_pft,npoints)
            CALL aggregation_request_data_multigrid(landpatch, ipatch, &
               grid_in1 = grid_pft, area = area_one, & 
               data_r8_3d_in1 = pftPCT, data_r8_3d_out1 = pct_one, n1_r8_3d_in1 = 16, lb1_r8_3d_in1 = 0, &
               grid_in2 = gland, data_r8_2d_in2 = htop, data_r8_2d_out2 = htop_one, &
               grid_in3 = grid_pft, data_r8_3d_in3 = cdepth, data_r8_3d_out3 = cdepth_one, n1_r8_3d_in3 = 16, lb1_r8_3d_in3 = 0, &
               grid_in4 = grid_pft, data_r8_3d_in4 = cratio, data_r8_3d_out4 = cratio_one, n1_r8_3d_in4 = 16, lb1_r8_3d_in4 = 0)
            where (htop_one < 0.) htop_one = 0.

            ! htop_patches (ipatch): 这个patch的平均冠层顶部高度
            ! htop_pfts (ip*ipatch): 这个PFT在它实际所占据面积上的平均冠层顶部高度
            ! htop_one (npoints): 这个patch中每个网格点的树高，维度就是这个ipatch中包含的网格点数
            htop_patches(ipatch) = sum(htop_one * area_one) / sum(area_one) !只是诊断变量
            pct_one = max(pct_one , 0.0)
            print*, '-----ipatch = ', ipatch, '/', numpatch
            print*, 'patch_pft_s(ipatch): ', patch_pft_s(ipatch)
            print*, 'patch_pft_e(ipatch): ', patch_pft_e(ipatch)
            ! ! 输出cdepth_one中小于等于0的个数
            ! print*, 'cdepth_one <= 0: ', count(cdepth_one <= 0.)
            ! print*, 'cdepth_one: ', cdepth_one
            hbot_one = spread(htop_one, dim=1, ncopies=16) - cdepth_one
            ! print*, 'hbot_one (no limited): ', hbot_one
            where (cdepth_one <= 0.) hbot_one = 0.
            where (cratio_one <= 0.) cratio_one = 0.
            print*, 'htop_one: ', htop_one
            ! print*, 'hbot_one: (cdepth_one <=0)', hbot_one
            ! print*, 'hbot_one(1,:): ', hbot_one(1,:)
            ! print*, 'hbot_one(2,:): ', hbot_one(2,:)
            ! print*, 'hbot_one(3,:): ', hbot_one(3,:)
            ! print*, 'cratio_one: (cratio_one <= 0) ', cratio_one
            ! print*, 'cratio_one(1,:): ', cratio_one(1,:)
            ! print*, 'cratio_one(2,:): ', cratio_one(2,:)
            ! print*, 'cratio_one(3,:): ', cratio_one(3,:)
            ! print*, 'pct_one: ', pct_one
            ! print*, 'pct_one >0: ', count(pct_one > 0.)
            ! print*, 'pct_one(1,:): ', pct_one(1,:)
            ! print*, 'pct_one(2,:): ', pct_one(2,:)
            ! print*, 'pct_one(3,:): ', pct_one(3,:)
            DO ip = patch_pft_s(ipatch), patch_pft_e(ipatch)
               p = landpft%settyp(ip)
               print*, '-----pft index: ', p
               print*, 'pct_pft: ', pct_one(p,:)
               print*, 'hbot: ', hbot_one(p+1,:)
               print*, 'cratio: ', cratio_one(p,:)
            ENDDO

#ifndef CROP
            IF (patchtypes(landpatch%settyp(ipatch)) == 0) THEN
#else
            IF (patchtypes(landpatch%settyp(ipatch)) == 0 .and. landpatch%settyp(ipatch)/=CROPLAND) THEN
#endif
               DO ip = patch_pft_s(ipatch), patch_pft_e(ipatch) ! 这个patch的第一个PFT索引~最后一个PFT索引，PFT在这个patch中的本地索引
                  p = landpft%settyp(ip) ! PFT索引ip对应的PFT类型编号，eg 02=BoNET
                  sumarea = sum(pct_one(p,:) * area_one) ! 这个PFT所占的有效面积
                  IF (sumarea > 0) THEN
                     htop_pfts(ip) = sum(htop_one * pct_one(p,:) * area_one) / sumarea 
                     hbot_pfts(ip) = sum(hbot_one(p+1,:) * pct_one(p,:) * area_one) / sumarea
                     cratio_pfts(ip) = sum(cratio_one(p,:) * pct_one(p,:) * area_one) / sumarea
                     ! hbot_one(1~16,npoints), pct_one(0~15,npoints)
                     ! 如果是作物或裸土，hbot_one本身就是0，hbot_pfts也是0
                     ! 没关系，mkinidata/MOD_HtopReadin.F90中，只使用了PFT=1~8tree的htop，其余PFT用不到，都是Const默认值
                  ELSE ! 暂时没有走这个选项的
                     ! htop_pfts(ip) = htop_patches(ipatch)
                     htop_pfts(ip) = htop_patches(ipatch)
                     hbot_pfts(ip) = 0.
                     cratio_pfts(ip) = 1.
                  ENDIF
               ENDDO
               ! IF (patch_weight_sum > 0) THEN
               !    hbot_patches(ipatch) = patch_bot_sum / patch_weight_sum ! 计算patch平均值
               ! ELSE
               !    hbot_patches(ipatch) = 0.
               ! ENDIF
               ! DO ip = patch_pft_s(ipatch), patch_pft_e(ipatch)
               !    p = landpft%settyp(ip)
               !    sumarea = sum(pct_one(p,:) * area_one)
               !    IF (sumarea <= 0.) THEN
               !       htop_pfts(ip) = htop_patches(ipatch)
               !       ! hbot_pfts(ip) = hbot_patches(ipatch) ! 为无面积的PFT赋值
               !       hbot_pfts(ip) = 0.
               !    ENDIF
               ! ENDDO
#ifdef CROP
            ELSEIF (landpatch%settyp(ipatch) == CROPLAND) THEN
               ip = patch_pft_s(ipatch)
               htop_pfts(ip) = htop_patches(ipatch)
               ! hbot_pfts(ip) = hbot_patches(ipatch)
               hbot_pfts(ip) = 0. ! 本来crop的默认hbot就是0
               cratio_pfts(ip) = 1.
#endif
            ENDIF
         ENDDO

#ifdef USEMPI
      CALL aggregation_worker_done_multigrid ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('HTOP_patches ', htop_patches)
      CALL check_vector_data ('HTOP_pfts    ', htop_pfts   )
      ! CALL check_vector_data ('HBOT_patches ', hbot_patches)
      CALL check_vector_data ('HBOT_pfts    ', hbot_pfts   )
      CALL check_vector_data ('CRATIO_pfts  ', cratio_pfts )
#endif

      lndname = trim(landdir)//'/htop_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'htop_patches', 'patch', landpatch, htop_patches, DEF_Srfdata_CompressLevel)

      ! lndname = trim(landdir)//'/hbot_patches.nc'
      ! CALL ncio_create_file_vector (lndname, landpatch)
      ! CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      ! CALL ncio_write_vector (lndname, 'hbot_patches', 'patch', landpatch, hbot_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typpatch = (/(ityp, ityp = 0, N_land_classification)/)
      lndname  = trim(dir_model_landdata) // '/diag/htop_patch_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (htop_patches, landpatch%settyp, typpatch, m_patch2diag, &
         -1.0e36_r8, lndname, 'htop', compress = 6, write_mode = 'one', defval=0._r8, create_mode=.true.)
      ! lndname  = trim(dir_model_landdata) // '/diag/hbot_patch_' // trim(cyear) // '.nc'
      ! CALL srfdata_map_and_write (hbot_patches, landpatch%settyp, typpatch, m_patch2diag, &
      !    -1.0e36_r8, lndname, 'hbot', compress = 6, write_mode = 'one', defval=0._r8, create_mode=.true.)
#endif

      lndname = trim(landdir)//'/htop_pfts.nc'
      CALL ncio_create_file_vector (lndname, landpft)
      CALL ncio_define_dimension_vector (lndname, landpft, 'pft')
      CALL ncio_write_vector (lndname, 'htop_pfts', 'pft', landpft, htop_pfts, DEF_Srfdata_CompressLevel)

      lndname = trim(landdir)//'/hbot_pfts.nc'
      CALL ncio_create_file_vector (lndname, landpft)
      CALL ncio_define_dimension_vector (lndname, landpft, 'pft')
      CALL ncio_write_vector (lndname, 'hbot_pfts', 'pft', landpft, hbot_pfts, DEF_Srfdata_CompressLevel)

      lndname = trim(landdir)//'/cratio_pfts.nc'
      CALL ncio_create_file_vector (lndname, landpft)
      CALL ncio_define_dimension_vector (lndname, landpft, 'pft')
      CALL ncio_write_vector (lndname, 'cratio_pfts', 'pft', landpft, cratio_pfts, DEF_Srfdata_CompressLevel)
#ifdef SrfdataDiag
#ifndef CROP
      typpft  = (/(ityp, ityp = 0, N_PFT-1)/)
#else
      typpft  = (/(ityp, ityp = 0, N_PFT+N_CFT-1)/)
#endif
      lndname = trim(dir_model_landdata) // '/diag/htop_pft_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (htop_pfts, landpft%settyp, typpft, m_pft2diag, &
         -1.0e36_r8, lndname, 'htop_pft', compress = 6, write_mode = 'one', defval=0._r8, create_mode=.true.)
      lndname = trim(dir_model_landdata) // '/diag/hbot_pft_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (hbot_pfts, landpft%settyp, typpft, m_pft2diag, &
         -1.0e36_r8, lndname, 'hbot_pft', compress = 6, write_mode = 'one', defval=0._r8, create_mode=.true.)
      lndname = trim(dir_model_landdata) // '/diag/cratio_pft_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (cratio_pfts, landpft%settyp, typpft, m_pft2diag, &
         -1.0e36_r8, lndname, 'cratio_pft', compress = 6, write_mode = 'one', defval=0._r8, create_mode=.true.)
#endif

      IF (p_is_worker) THEN
         IF (allocated(htop_patches)) deallocate (htop_patches)
         IF (allocated(htop_pfts   )) deallocate (htop_pfts   )
         IF (allocated(htop_one    )) deallocate (htop_one    )
         IF (allocated(pct_one     )) deallocate (pct_one     )
         IF (allocated(area_one    )) deallocate (area_one    )
         ! IF (allocated(hbot_patches)) deallocate (hbot_patches)
         IF (allocated(hbot_pfts   )) deallocate (hbot_pfts   )
         IF (allocated(cdepth_one  )) deallocate (cdepth_one  )
         IF (allocated(hbot_one    )) deallocate (hbot_one    )
         IF (allocated(cratio_pfts )) deallocate (cratio_pfts )
         IF (allocated(cratio_one  )) deallocate (cratio_one  )
      ENDIF
#endif

END SUBROUTINE Aggregation_ForestHeight
