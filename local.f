c-----------------------------------------------------------------------
c     file local.f
c     local defintions for most computers.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. local.
c     1. timer.
c     2. program_stop.
c     3. mylog.
c-----------------------------------------------------------------------
c     subprogram 0. local.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE local_mod
      USE io_mod
      IMPLICIT NONE

      CHARACTER(6) :: convert_type="none"
      LOGICAL, PARAMETER :: rewind_namel=.false.,single_pr=.false.
      INTEGER, PARAMETER ::
     $     r4=SELECTED_REAL_KIND(6,37),
     $     r8=SELECTED_REAL_KIND(13,307)
      REAL(r8), PARAMETER :: pi=3.1415926535897932385_r8,
     $     twopi=2*pi,pisq=pi*pi,mu0=4e-7_r8*pi,
     $     rtod=180/pi,dtor=pi/180,alog10=2.302585093_r8
      REAL(r8), PARAMETER :: 
     $     mp=1.672614e-27,me=9.1091e-31,
     $     c=2.997925e10,jzero=2.4048255577_r8
      REAL(r8), PARAMETER :: zero=0,one=1,two=2,half=.5
      COMPLEX(r8), PARAMETER :: ifac=(0,1)

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. timer.
c     handles machine-dependent timing statistics.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     delcarations.
c-----------------------------------------------------------------------
      SUBROUTINE timer(mode,unit)
       
      INTEGER, INTENT(IN) :: mode,unit

      CHARACTER(10) :: date,time,zone
      INTEGER, DIMENSION(8) :: values
      REAL(4), SAVE :: seconds
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(1x,a,1p,e10.3,a)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      IF(mode == 0)THEN
         CALL DATE_AND_TIME(date,time,zone,values)
         seconds=(values(5)*60+values(6))*60+values(7)+values(8)*1e-3
      ELSE
         CALL DATE_AND_TIME(date,time,zone,values)
         seconds=(values(5)*60+values(6))*60+values(7)+values(8)*1e-3
     $        -seconds
         WRITE(unit,10)"Total cpu time = ",seconds," seconds"
         WRITE(*,10)"Total cpu time = ",seconds," seconds"
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE timer
c-----------------------------------------------------------------------
c     subprogram 2. program_stop.
c     terminates program with message, calls timer, closes output file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE program_stop(message)

      CHARACTER(*), INTENT(IN) :: message
c-----------------------------------------------------------------------
c     write completion message.
c-----------------------------------------------------------------------
      CALL timer(1,out_unit)
      CLOSE(out_unit)
      WRITE(*,'(1x,2a)') 'PROGRAM_STOP => ', TRIM(message)
c-----------------------------------------------------------------------
c     write completion message.
c-----------------------------------------------------------------------
      STOP
      END SUBROUTINE program_stop
c-----------------------------------------------------------------------
c     subprogram 3. mylog.
c     returns bounded log of specific component.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION mylog(u) RESULT(ulog)

      COMPLEX(r8), INTENT(IN) :: u
      REAL(r4) :: ulog

      REAL, PARAMETER :: minlog=-15
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      IF(u == 0)THEN
         ulog=minlog
      ELSE
         ulog=LOG10(ABS(u))
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION mylog
      END MODULE local_mod
