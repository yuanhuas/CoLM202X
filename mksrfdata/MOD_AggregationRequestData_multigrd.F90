#include <define.h>

MODULE MOD_AggregationRequestData_multigrid

!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!    Aggregation Utilities.
!
!    On IO processes, a data daemon is running to provide data
!       at fine resolutions for worker processes.
!    On worker processes, request is sent to IO processes and
!       data is returned from IO processes.
!
!  Created by Shupeng Zhang, May 2023
!-----------------------------------------------------------------------

   IMPLICIT NONE

   PUBLIC :: aggregation_request_data_multigrid

#ifdef USEMPI
   PUBLIC :: aggregation_data_daemon_multigrid
   PUBLIC :: aggregation_worker_done_multigrid
#endif

! ---- subroutines ----
CONTAINS

#ifdef USEMPI
   SUBROUTINE aggregation_data_daemon_multigrid (grid_in1, data_r8_2d_in1, &
         data_r8_3d_in1, n1_r8_3d_in1  ,&
         grid_in2, data_r8_2d_in2, &
         data_r8_3d_in2, n1_r8_3d_in2)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_DataType

   IMPLICIT NONE

   type (grid_type), intent(in) :: grid_in1
   type (grid_type), intent(in) :: grid_in3

   ! 2D REAL data
   type (block_data_real8_2d), intent(in), optional :: data_r8_2d_in1
   type (block_data_real8_2d), intent(in), optional :: data_r8_2d_in2

   ! 3D REAL data
   integer, intent(in), optional :: n1_r8_3d_in1
   type (block_data_real8_3d), intent(in), optional :: data_r8_3d_in1

   integer, intent(in), optional :: n1_r8_3d_in2
   type (block_data_real8_3d), intent(in), optional :: data_r8_3d_in2

   ! Local Variables
   integer :: nreq, ireq, rmesg(3), isrc, idest, grdn
   integer :: xblk, yblk, xloc, yloc
   integer,  allocatable :: ylist(:), xlist(:)

   real(r8), allocatable :: sbuf_r8_1d(:), sbuf_r8_2d(:,:)
   integer , allocatable :: sbuf_i4_1d(:)

   logical,  allocatable :: worker_done (:)

      IF (p_is_io) THEN

         allocate (worker_done (0:p_np_worker-1))

         worker_done(:) = .false.
         DO WHILE (any(.not. worker_done))

            CALL mpi_recv (rmesg, 3, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc = rmesg(1)
            nreq = rmesg(2)
            grdn = rmesg(3)

            IF (nreq > 0) THEN

               allocate (xlist (nreq))
               allocate (ylist (nreq))

               CALL mpi_recv (xlist, nreq, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL mpi_recv (ylist, nreq, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               idest = isrc

               allocate (sbuf_r8_1d (nreq))

               IF (present(data_r8_2d_in1)) THEN
                  IF (grdn == 1) THEN
                     DO ireq = 1, nreq
                        xblk = grid_in1%xblk(xlist(ireq))
                        yblk = grid_in1%yblk(ylist(ireq))
                        xloc = grid_in1%xloc(xlist(ireq))
                        yloc = grid_in1%yloc(ylist(ireq))

                        sbuf_r8_1d(ireq) = data_r8_2d_in1%blk(xblk,yblk)%val(xloc,yloc)

                     ENDDO
                     CALL mpi_send (sbuf_r8_1d, nreq, MPI_REAL8, &
                        idest, mpi_tag_data, p_comm_glb, p_err)

                   ENDIF
               ENDIF

               IF (present(data_r8_2d_in2) .and. present(grid_in2)) THEN
                  IF (grdn == 2) THEN
                     DO ireq = 1, nreq
                        xblk = grid_in2%xblk(xlist(ireq))
                        yblk = grid_in2%yblk(ylist(ireq))
                        xloc = grid_in2%xloc(xlist(ireq))
                        yloc = grid_in2%yloc(ylist(ireq))

                        sbuf_r8_1d(ireq) = data_r8_2d_in2%blk(xblk,yblk)%val(xloc,yloc)

                     ENDDO

                     CALL mpi_send (sbuf_r8_1d, nreq, MPI_REAL8, &
                        idest, mpi_tag_data, p_comm_glb, p_err)

                  ENDIF
               ENDIF

               deallocate (sbuf_r8_1d)

               IF ( present(data_r8_3d_in1) .and. present(n1_r8_3d_in1) ) THEN

                  allocate (sbuf_r8_2d (n1_r8_3d_in1,nreq))
                  IF (grdn == 1) THEN
                     DO ireq = 1, nreq
                        xblk = grid_in1%xblk(xlist(ireq))
                        yblk = grid_in1%yblk(ylist(ireq))
                        xloc = grid_in1%xloc(xlist(ireq))
                        yloc = grid_in1%yloc(ylist(ireq))

                        sbuf_r8_2d(:,ireq) = data_r8_3d_in1%blk(xblk,yblk)%val(:,xloc,yloc)

                     ENDDO

                     CALL mpi_send (sbuf_r8_2d, n1_r8_3d_in1*nreq, MPI_REAL8, &
                        idest, mpi_tag_data, p_comm_glb, p_err)

                  ENDIF

                 deallocate (sbuf_r8_2d)
               ENDIF

               IF ( present(data_r8_3d_in2) .and. present(n1_r8_3d_in2) ) THEN

                  allocate (sbuf_r8_2d (n1_r8_3d_in2,nreq))
                  IF (grdn == 2) THEN
                     DO ireq = 1, nreq
                        xblk = grid_in2%xblk(xlist(ireq))
                        yblk = grid_in2%yblk(ylist(ireq))
                        xloc = grid_in2%xloc(xlist(ireq))
                        yloc = grid_in2%yloc(ylist(ireq))

                        sbuf_r8_2d(:,ireq) = data_r8_3d_in2%blk(xblk,yblk)%val(:,xloc,yloc)

                     ENDDO

                     CALL mpi_send (sbuf_r8_2d, n1_r8_3d_in2*nreq, MPI_REAL8, &
                        idest, mpi_tag_data, p_comm_glb, p_err)

                  ENDIF

                 deallocate (sbuf_r8_2d)
               ENDIF

               deallocate (ylist)
               deallocate (xlist)

            ELSE
               worker_done(p_itis_worker(isrc)) = .true.
            ENDIF

         ENDDO

         deallocate (worker_done)

      ENDIF

   END SUBROUTINE aggregation_data_daemon_multigrid


#endif

   !----------------------------------------------------
   SUBROUTINE aggregation_request_data_multigrid (     &
         pixelset, iset, grid_in1, area, &
         data_r8_2d_in1, data_r8_2d_out1, &
         data_r8_3d_in1, data_r8_3d_out1, n1_r8_3d_in1, lb1_r8_3d_in1, &
         grid_in2, &
         data_r8_2d_in2, data_r8_2d_out2, &
         data_r8_3d_in2, data_r8_3d_out2, n1_r8_3d_in2, lb1_r8_3d_in2, &
         filledvalue_i4)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Block
   USE MOD_Pixel
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_Mesh
   USE MOD_Pixelset
   USE MOD_Utils

   IMPLICIT NONE

   type (pixelset_type), intent(in) :: pixelset
   integer, intent(in)  :: iset

   type (grid_type), intent(in) :: grid_in1
   type (grid_type), intent(in) :: grid_in2

   real(r8), allocatable, intent(out), optional :: area(:)

   type (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in1
   real(r8), allocatable,      intent(out), optional :: data_r8_2d_out1 (:)

   type (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in2
   real(r8), allocatable,      intent(out), optional :: data_r8_2d_out2 (:)

   integer, intent(in), optional :: n1_r8_3d_in1, lb1_r8_3d_in1
   type (block_data_real8_3d), intent(in),  optional :: data_r8_3d_in1
   real(r8), allocatable,      intent(out), optional :: data_r8_3d_out1 (:,:)

   integer, intent(in), optional :: n1_r8_3d_in2, lb1_r8_3d_in2
   type (block_data_real8_3d), intent(in),  optional :: data_r8_3d_in2
   real(r8), allocatable,      intent(out), optional :: data_r8_3d_out2 (:,:)

   integer, intent(in), optional :: filledvalue_i4

   ! Local Variables
   integer :: totalreq, ireq, nreq1, nreq2, smesg(3), isrc, idest, iproc
   integer :: ilon, ilat, xblk, yblk, xloc, yloc, iloc, nx, ny, ix, iy, ig
   integer :: ie, ipxstt, ipxend, npxl, ipxl, lb1, xgrdthis, ygrdthis
   integer,  allocatable :: ylist1(:), xlist1(:), ipt1(:)
   integer,  allocatable :: ylist2(:), xlist2(:), ipt2(:)
   integer,  allocatable :: ibuf(:), rbuf_i4_1d(:)
   integer,  allocatable :: xsorted(:), ysorted(:), xy2d(:,:)
   real(r8), allocatable :: area2d(:,:), rbuf_r8_1d(:), rbuf_r8_2d(:,:)
   logical,  allocatable :: msk1(:), msk2(:)


      ie     = pixelset%ielm  (iset)
      ipxstt = pixelset%ipxstt(iset)
      ipxend = pixelset%ipxend(iset)
      npxl   = ipxend - ipxstt + 1

      allocate(xlist1 (npxl))
      allocate(ylist1 (npxl))

      allocate(xlist2 (npxl))
      allocate(ylist2 (npxl))

      IF (present(area)) allocate (area (npxl))

      totalreq = npxl
      DO ipxl = ipxstt, ipxend
         xlist1(ipxl-ipxstt+1) = grid_in1%xgrd(mesh(ie)%ilon(ipxl))
         ylist1(ipxl-ipxstt+1) = grid_in1%ygrd(mesh(ie)%ilat(ipxl))

         xlist2(ipxl-ipxstt+1) = grid_in2%xgrd(mesh(ie)%ilon(ipxl))
         ylist2(ipxl-ipxstt+1) = grid_in2%ygrd(mesh(ie)%ilat(ipxl))

         IF (present(area)) THEN
            area(ipxl-ipxstt+1) = areaquad (&
               pixel%lat_s(mesh(ie)%ilat(ipxl)), pixel%lat_n(mesh(ie)%ilat(ipxl)), &
               pixel%lon_w(mesh(ie)%ilon(ipxl)), pixel%lon_e(mesh(ie)%ilon(ipxl)) )
         ENDIF
      ENDDO


      IF (present(data_r8_2d_in1) .and. present(data_r8_2d_out1))  allocate (data_r8_2d_out1 (totalreq))
      IF (present(data_r8_2d_in2) .and. present(data_r8_2d_out2))  allocate (data_r8_2d_out2 (totalreq))

      IF (present(data_r8_3d_in1) .and. present(data_r8_3d_out1) .and. present(n1_r8_3d_in1)) THEN
         IF (present(lb1_r8_3d_in1)) THEN
            lb1 = lb1_r8_3d_in1
         ELSE
            lb1 = 1
         ENDIF
         allocate (data_r8_3d_out1 (lb1:lb1-1+n1_r8_3d_in1,totalreq))
      ENDIF

      IF (present(data_r8_3d_in2) .and. present(data_r8_3d_out2) .and. present(n1_r8_3d_in2)) THEN
         IF (present(lb1_r8_3d_in2)) THEN
            lb1 = lb1_r8_3d_in2
         ELSE
            lb1 = 1
         ENDIF
         allocate (data_r8_3d_out2 (lb1:lb1-1+n1_r8_3d_in2,totalreq))
      ENDIF

#ifdef USEMPI

      allocate (ipt1 (totalreq))
      allocate (msk1 (totalreq))

      allocate (ipt2 (totalreq))
      allocate (msk2 (totalreq))

      ipt1(:) = -1
      ipt2(:) = -1

      DO ireq = 1, totalreq
         xblk = grid_in1%xblk(xlist1(ireq))
         yblk = grid_in1%yblk(ylist1(ireq))
         ipt1(ireq) = gblock%pio(xblk,yblk)

         IF ( present(grid_in2) ) THEN
            xblk = grid_in2%xblk(xlist2(ireq))
            yblk = grid_in2%yblk(ylist2(ireq))
            ipt2(ireq) = gblock%pio(xblk,yblk)
         ENDIF
      ENDDO

      DO iproc = 0, p_np_io-1
         msk1 = (ipt1 == p_address_io(iproc))
         nreq1 = count(msk1)

         IF (nreq1 > 0) THEN

            smesg = (/p_iam_glb, nreq1, 1/)

            idest = p_address_io(iproc)
            CALL mpi_send (smesg, 3, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)

            allocate (ibuf (nreq1))

            ibuf = pack(xlist1(1:totalreq), msk1)
            CALL mpi_send (ibuf, nreq1, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

            ibuf = pack(ylist1(1:totalreq), msk1)
            CALL mpi_send (ibuf, nreq1, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

            isrc = idest

            allocate (rbuf_r8_1d (nreq1))

            IF (present(data_r8_2d_in1) .and. present(data_r8_2d_out1)) THEN
               CALL mpi_recv (rbuf_r8_1d, nreq1, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_1d, msk1, data_r8_2d_out1)
            ENDIF

            deallocate (rbuf_r8_1d)

            IF (present(data_r8_3d_in1) .and. present(data_r8_3d_out1) .and. present(n1_r8_3d_in1)) THEN
               allocate (rbuf_r8_2d (n1_r8_3d_in1,nreq1))
               CALL mpi_recv (rbuf_r8_2d, n1_r8_3d_in1*nreq1, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_2d, msk1, data_r8_3d_out1)
               deallocate (rbuf_r8_2d)
            ENDIF

            deallocate (ibuf)
         ENDIF

         IF ( present(grid_in2) ) THEN
            msk2 = (ipt2 == p_address_io(iproc))
            nreq2 = count(msk2)

            IF (nreq2 > 0) THEN
               smesg = (/p_iam_glb, nreq2, 2/)

               idest = p_address_io(iproc)
               CALL mpi_send (smesg, 3, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)

               allocate (ibuf (nreq2))

               ibuf = pack(xlist2(1:totalreq), msk2)
               CALL mpi_send (ibuf, nreq2, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

               ibuf = pack(ylist2(1:totalreq), msk2)
               CALL mpi_send (ibuf, nreq2, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

               isrc = idest

               allocate (rbuf_r8_1d (nreq2))

               IF (present(data_r8_2d_in2) .and. present(data_r8_2d_out2)) THEN
                  CALL mpi_recv (rbuf_r8_1d, nreq2, MPI_REAL8, &
                     isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
                  CALL unpack_inplace (rbuf_r8_1d, msk2, data_r8_2d_out2)
               ENDIF

               deallocate (rbuf_r8_1d)

               IF (present(data_r8_3d_in2) .and. present(data_r8_3d_out2)) THEN
                  allocate (rbuf_r8_2d (n1_r8_3d_in2,nreq2))
                  CALL mpi_recv (rbuf_r8_2d, n1_r8_3d_in2*nreq2, MPI_REAL8, &
                     isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
                  CALL unpack_inplace (rbuf_r8_2d, msk2, data_r8_3d_out2)
                  deallocate (rbuf_r8_2d)
               ENDIF

               deallocate (ibuf)
            ENDIF
         ENDIF
      ENDDO

      deallocate (xlist1)
      deallocate (ylist1)
      deallocate (ipt1  )
      deallocate (msk1  )

      IF ( present(grid_in2) ) THEN
         deallocate (xlist2)
         deallocate (ylist2)
         deallocate (ipt2  )
         deallocate (msk2  )
      ENDIF
#else

      DO ireq = 1, totalreq

         IF (present(data_r8_2d_in1) .and. present(data_r8_2d_out1)) THEN
            xblk = grid_in1%xblk(xlist1(ireq))
            yblk = grid_in1%yblk(ylist1(ireq))
            xloc = grid_in1%xloc(xlist1(ireq))
            yloc = grid_in1%yloc(ylist1(ireq))

            data_r8_2d_out1(ireq) = data_r8_2d_in1%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_2d_in2) .and. present(data_r8_2d_out2)) THEN
            xblk = grid_in2%xblk(xlist2(ireq))
            yblk = grid_in2%yblk(ylist2(ireq))
            xloc = grid_in2%xloc(xlist2(ireq))
            yloc = grid_in2%yloc(ylist2(ireq))

            data_r8_2d_out2(ireq) = data_r8_2d_in2%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_3d_in1) .and. present(data_r8_3d_out1) .and. present(n1_r8_3d_in1)) THEN
            xblk = grid_in1%xblk(xlist1(ireq))
            yblk = grid_in1%yblk(ylist1(ireq))
            xloc = grid_in1%xloc(xlist1(ireq))
            yloc = grid_in1%yloc(ylist1(ireq))

            data_r8_3d_out1(:,ireq) = data_r8_3d_in1%blk(xblk,yblk)%val(:,xloc,yloc)
         ENDIF

         IF (present(data_r8_3d_in2) .and. present(data_r8_3d_out2) .and. present(n1_r8_3d_in2)) THEN
            xblk = grid_in2%xblk(xlist2(ireq))
            yblk = grid_in2%yblk(ylist2(ireq))
            xloc = grid_in2%xloc(xlist2(ireq))
            yloc = grid_in2%yloc(ylist2(ireq))

            data_r8_3d_out2(:,ireq) = data_r8_3d_in2%blk(xblk,yblk)%val(:,xloc,yloc)
         ENDIF

      ENDDO

#endif

   END SUBROUTINE aggregation_request_data_multigrid

#ifdef USEMPI

   SUBROUTINE aggregation_worker_done_multigrid ()

   USE MOD_SPMD_Task

   IMPLICIT NONE

   integer :: smesg(3), iproc, idest

      IF (p_is_worker) THEN
         DO iproc = 0, p_np_io-1
            smesg = (/p_iam_glb, -1, 1/)
            idest = p_address_io(iproc)
            CALL mpi_send (smesg, 3, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)
          ENDDO
      ENDIF

   END SUBROUTINE aggregation_worker_done_multigrid

#endif

END MODULE MOD_AggregationRequestData_multigrid
