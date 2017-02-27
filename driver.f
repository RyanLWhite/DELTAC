c-----------------------------------------------------------------------
c     file driver.f.
c     tests deltac.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. driver_mod.
c     1. driver_run.
c     2. driver_dqscan.
c     3. driver_sxscan.
c     4. driver_conv.
c     5. driver_main.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subprogram 0. driver_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE driver_mod
      USE deltar_mod
      USE deltac_mod
      IMPLICIT NONE

      COMPLEX(r8), DIMENSION(2,2) :: fulldelta

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. driver_run.
c     controls program execution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE driver_run

      LOGICAL :: dqscan_flag=.FALSE.,sxscan_flag=.FALSE.,
     $     conv_flag=.FALSE.

      NAMELIST/run_list/dqscan_flag,sxscan_flag,conv_flag
c-----------------------------------------------------------------------
c     read input.
c-----------------------------------------------------------------------
      OPEN(UNIT=in_unit,FILE="deltac.in",STATUS="OLD")
      READ(in_unit,NML=run_list)
      CLOSE(UNIT=in_unit)
c-----------------------------------------------------------------------
c     run tests.
c-----------------------------------------------------------------------
      IF(dqscan_flag)CALL driver_dqscan
      IF(sxscan_flag)CALL driver_sxscan
      IF(conv_flag)CALL driver_conv
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE driver_run
c-----------------------------------------------------------------------
c     subprogram 2. driver_dqscan.
c     scans over DR and Q.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE driver_dqscan

      LOGICAL :: out,bin
      LOGICAL :: deltar_flag=.FALSE.,deltac_flag=.FALSE.
      INTEGER :: nslog,is,idr,ndr,count
      REAL(r8), PARAMETER :: eps=1e-10,epsd=1e-4,rmatch=3
      REAL(r8) :: e,f,g,h,k,m,taua,taur,slogmin,slogmax,slog,dslog,
     $     drmax,drmin,ddr,dr,dlim=2,v1=1
      COMPLEX(r8) :: s
      COMPLEX(r8), DIMENSION(2) :: deltar,deltac
      TYPE(resist_type) :: restype

      NAMELIST/dqscan_list/drmax,drmin,ndr,e,f,g,h,k,m,taua,taur,
     $     slogmin,slogmax,nslog,out,bin,rtol,atol,ntmax,fmax,
     $     deltar_flag,deltac_flag,dlim,v1
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(6x,"dr",9x,"e",10x,"f",10x,"g",10x,"h",10x,"k",10x,"m"
     $     //1p,7e11.3)
c-----------------------------------------------------------------------
c     read input.
c-----------------------------------------------------------------------
      OPEN(UNIT=in_unit,FILE="deltac.in",STATUS="OLD")
      READ(in_unit,NML=dqscan_list)
      CLOSE(UNIT=in_unit)
      WRITE(*,'(a,1p,e10.3)')" dlim = ",dlim
c-----------------------------------------------------------------------
c     open output files.
c-----------------------------------------------------------------------
      IF(deltar_flag)THEN
         OPEN(UNIT=binr_unit,FILE="deltar.bin",STATUS="REPLACE",
     $        FORM="UNFORMATTED")
      ENDIF
      IF(deltac_flag)THEN
         OPEN(UNIT=binc_unit,FILE="deltac.bin",STATUS="REPLACE",
     $        FORM="UNFORMATTED")
         CALL deltac_read_parameters("deltac.in")
      ENDIF
      OPEN(UNIT=x01_bin_unit,FILE="x01.bin",STATUS="REPLACE",
     $     FORM="UNFORMATTED")
c-----------------------------------------------------------------------
c     set restype.
c-----------------------------------------------------------------------
      CALL timer(0,out_unit)
      restype%f=f
      restype%g=g
      restype%h=h
      restype%k=k
      restype%m=m
      restype%v1=v1
      restype%taua=taua
      restype%taur=taur
c-----------------------------------------------------------------------
c     start loops over dr and s.
c-----------------------------------------------------------------------
      dslog=one/nslog
      ddr=(drmax-drmin)/ndr
      dr=drmin
      count=0
      DO idr=0,ndr
         restype%e=dr-f-h**2
         IF(ABS(dr) < eps)restype%e=restype%e+eps
         slog=slogmin
         DO is=0,nslog*(slogmax-slogmin)
            s=10**slog
c-----------------------------------------------------------------------
c     run deltar code and record output.
c-----------------------------------------------------------------------
            IF(deltar_flag)THEN
               CALL deltar_run(restype,s,deltar)
               IF(bin .AND. ABS(deltar(2)) < dlim)
     $              WRITE(binr_unit)REAL(slog,4),
     $              mylog(deltar(1)),REAL(deltar(2),4)
            ENDIF
c-----------------------------------------------------------------------
c     run deltac code and record output.
c-----------------------------------------------------------------------
            IF(deltac_flag)THEN
               CALL deltac_run(restype,s,deltac,fulldelta)
               IF(bin .AND. ABS(deltac(2)) < dlim)
     $              WRITE(binc_unit)REAL(slog,4),
     $              mylog(deltac(1)),REAL(deltac(2),4) !,
            ENDIF
c-----------------------------------------------------------------------
c     finish loops over dr and s.
c-----------------------------------------------------------------------
            slog=slog+dslog
            count=count+1
         ENDDO
         IF(deltar_flag)WRITE(binr_unit)
         IF(deltac_flag)WRITE(binc_unit)
         WRITE(x01_bin_unit)
         dr=dr+ddr
      ENDDO
      WRITE(*,'(a,i5)')" count = ",count
c-----------------------------------------------------------------------
c     close output files.
c-----------------------------------------------------------------------
      IF(deltar_flag)THEN
         CLOSE(UNIT=binr_unit)
      ENDIF
      IF(deltac_flag)THEN
         CLOSE(UNIT=binc_unit)
      ENDIF
      CLOSE(UNIT=x01_bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE driver_dqscan
c-----------------------------------------------------------------------
c     subprogram 3. driver_sxscan.
c     scans over Q and xmax.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE driver_sxscan

      LOGICAL :: out,bin
      INTEGER :: nslog,nxlog,is,ixmax,count
      REAL(r8), PARAMETER :: eps=1e-10
      REAL(r8) :: e,f,g,h,k,m,taua,taur,slogmin,slogmax,slog,dslog,
     $     xlogmax,xlogmin,dxlog,xlog,v1=1
      COMPLEX(r8) :: s
      COMPLEX(r8), DIMENSION(2) :: deltac
      TYPE(resist_type) :: restype

      NAMELIST/sxscan_list/e,f,g,h,k,m,taua,taur,
     $     xlogmax,xlogmin,nxlog,slogmin,slogmax,nslog,
     $     out,bin,rtol,atol,ntmax,fmax
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(6x,"dr",9x,"e",10x,"f",10x,"g",10x,"h",10x,"k",10x,"m"
     $     //1p,7e11.3)
c-----------------------------------------------------------------------
c     read input.
c-----------------------------------------------------------------------
      OPEN(UNIT=in_unit,FILE="deltac.in",STATUS="OLD")
      READ(in_unit,NML=sxscan_list)
      CLOSE(UNIT=in_unit)
      OPEN(UNIT=out_unit,FILE="sxscan.out",STATUS="REPLACE")
c-----------------------------------------------------------------------
c     open output files and read input.
c-----------------------------------------------------------------------
      OPEN(UNIT=binx_unit,FILE="deltax.bin",STATUS="REPLACE",
     $     FORM="UNFORMATTED")
      CALL deltac_read_parameters("deltac.in")
c-----------------------------------------------------------------------
c     set restype.
c-----------------------------------------------------------------------
      CALL timer(0,out_unit)
      restype%e=e
      restype%f=f
      restype%g=g
      restype%h=h
      restype%k=k
      restype%m=m
      restype%taua=taua
      restype%taur=taur
      restype%v1=v1
c-----------------------------------------------------------------------
c     start loops over q and xfac.
c-----------------------------------------------------------------------
      count=0
      dxlog=one/nxlog
      dslog=one/nslog
      slog=slogmin
      DO is=0,nslog*(slogmax-slogmin)
         s=10**slog
         xlog=xlogmin
         DO ixmax=0,nxlog*(xlogmax-xlogmin)
            xfac=10**xlog
c-----------------------------------------------------------------------
c     run deltac code and record output.
c-----------------------------------------------------------------------
            CALL deltac_run(restype,s,deltac,fulldelta)
            WRITE(binx_unit)
     $           REAL(LOG10(xfac),4),mylog(deltac(1)),mylog(deltac(2)),
     $           REAL(LOG10(xmax),4)
c-----------------------------------------------------------------------
c     finish loops over xlog and slog.
c-----------------------------------------------------------------------
            count=count+1
            xlog=xlog+dxlog
         ENDDO
         WRITE(binx_unit)
         slog=slog+dslog
      ENDDO
      WRITE(*,'(a,i5)')" count = ",count
      CLOSE(UNIT=binx_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE driver_sxscan
c-----------------------------------------------------------------------
c     subprogram 4. driver_conv.
c     graphic test for lsode convergence.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE driver_conv

      LOGICAL :: out,bin
      REAL(r8), PARAMETER :: eps=1e-10
      REAL(r8) :: e,f,g,h,k,m,taua,taur,xfac
      COMPLEX(r8) :: s
      COMPLEX(r8), DIMENSION(2) :: deltac
      TYPE(resist_type) :: restype

      NAMELIST/conv_list/e,f,g,h,k,m,taua,taur,xfac,s,out,bin,rtol,atol
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(6x,"dr",9x,"e",10x,"f",10x,"g",10x,"h",10x,"k",10x,"m"
     $     //1p,7e11.3)
c-----------------------------------------------------------------------
c     read input.
c-----------------------------------------------------------------------
      OPEN(UNIT=in_unit,FILE="deltac.in",STATUS="OLD")
      READ(in_unit,NML=conv_list)
      CLOSE(UNIT=in_unit)
c-----------------------------------------------------------------------
c     open output files and read input.
c-----------------------------------------------------------------------
      CALL deltac_read_parameters("deltac.in")
c-----------------------------------------------------------------------
c     set restype.
c-----------------------------------------------------------------------
      CALL timer(0,out_unit)
      restype%e=e
      restype%f=f
      restype%g=g
      restype%h=h
      restype%k=k
      restype%m=m
      restype%taua=taua
      restype%taur=taur
      restype%v1=1.0
c-----------------------------------------------------------------------
c     run deltac code and record output.
c-----------------------------------------------------------------------
c$$$      lsode_diagnose=.TRUE.
      OPEN(UNIT=conv_out_unit,FILE="conv.out",STATUS="REPLACE")
      OPEN(UNIT=conv_bin_unit,FILE="conv.bin",STATUS="REPLACE",
     $     FORM="UNFORMATTED")
      CALL deltac_run(restype,s,deltac,fulldelta)
      CLOSE(conv_out_unit)
      CLOSE(conv_bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE driver_conv
      END MODULE driver_mod
c-----------------------------------------------------------------------
c     subprogram 5. driver_main.
c     trivial main program.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      PROGRAM driver_main
      USE driver_mod
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      OPEN(UNIT=out_unit,FILE="driver.out",STATUS="REPLACE")
      CALL driver_run
      CALL program_stop("Normal termination.")
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      END PROGRAM driver_main
