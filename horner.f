c-----------------------------------------------------------------------
c     file horner.f.
c     Horner's method for polynomials and first two derivatives.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. horner_mod.
c     1. horner_1d.
c     2. horner_2d.
c     3. horner_3d.
c     4. horner_transform.
c     5. horner_pointeval.
c-----------------------------------------------------------------------
c     subprogram 0. horner_mod.
c-----------------------------------------------------------------------
      MODULE horner_mod
      USE local_mod
      USE jacobi_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. horner_1d.
c     1D polynomials.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE horner_1d(uu,order,x,u,du,ddu)

      REAL(8), DIMENSION(:,0:), INTENT(IN) :: uu
      INTEGER, INTENT(IN) :: order
      REAL(8), INTENT(IN) :: x
      REAL(8), DIMENSION(:), INTENT(OUT) :: u
      REAL(8), DIMENSION(:), INTENT(OUT) :: du
      REAL(8), DIMENSION(:), INTENT(OUT) :: ddu

      INTEGER :: ip,np
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      np=SIZE(uu,2)-1
      u=uu(:,np)
      IF(order > 0)du=np*uu(:,np)
      IF(order > 1)ddu=np*(np-1)*uu(:,np)
c-----------------------------------------------------------------------
c     iterate.
c-----------------------------------------------------------------------
      DO ip=np-1,0,-1

         u=u*x+uu(:,ip)

         IF(order < 1 .OR. ip < 1)CYCLE
         du=du*x+ip*uu(:,ip)

         IF(order < 2 .OR. ip < 2)CYCLE
         ddu=ddu*x+ip*(ip-1)*uu(:,ip)

      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE horner_1d
c-----------------------------------------------------------------------
c     subprogram 2. horner_2d.
c     2D polynomials.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE horner_2d(uu,order,x,y,u,du,ddu)

      REAL(8), DIMENSION(:,0:,0:), INTENT(IN) :: uu
      INTEGER, INTENT(IN) :: order
      REAL(8), INTENT(IN) :: x,y
      REAL(8), DIMENSION(:), INTENT(OUT) :: u
      REAL(8), DIMENSION(:,:), INTENT(OUT) :: du
      REAL(8), DIMENSION(:,:,:), INTENT(OUT) :: ddu

      INTEGER :: ip,np
      REAL(8), DIMENSION(SIZE(u,1)) :: term,dterm,ddterm
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      np=SIZE(uu,3)-1
      CALL horner_1d(uu(:,:,np),order,x,u,dterm,ddterm)
      IF(order > 0)THEN
         du(:,1)=dterm
         du(:,2)=np*u
      ENDIF
      IF(order > 1)THEN
         ddu(:,1,1)=ddterm
         ddu(:,1,2)=np*dterm
         ddu(:,2,2)=np*(np-1)*u
      ENDIF
c-----------------------------------------------------------------------
c     iterate.
c-----------------------------------------------------------------------
      DO ip=np-1,0,-1

         CALL horner_1d(uu(:,:,ip),order,x,term,dterm,ddterm)
         u=u*y+term

         IF(order < 1)CYCLE
         du(:,1)=du(:,1)*y+dterm
         IF(ip > 0)du(:,2)=du(:,2)*y+ip*term

         IF(order < 2)CYCLE
         ddu(:,1,1)=ddu(:,1,1)*y+ddterm
         IF(ip > 0)ddu(:,1,2)=ddu(:,1,2)*y+ip*dterm
         IF(ip > 1)ddu(:,2,2)=ddu(:,2,2)*y+ip*(ip-1)*term

      ENDDO
c-----------------------------------------------------------------------
c     symmetrize second derivatives.
c-----------------------------------------------------------------------
      ddu(:,2,1)=ddu(:,1,2)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE horner_2d
c-----------------------------------------------------------------------
c     subprogram 3. horner_3d.
c     3D polynomials.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE horner_3d(uu,order,x,y,z,u,du,ddu)

      REAL(8), DIMENSION(:,0:,0:,0:), INTENT(IN) :: uu
      INTEGER, INTENT(IN) :: order
      REAL(8), INTENT(IN) :: x,y,z
      REAL(8), DIMENSION(:), INTENT(OUT) :: u
      REAL(8), DIMENSION(:,:), INTENT(OUT) :: du
      REAL(8), DIMENSION(:,:,:), INTENT(OUT) :: ddu

      INTEGER :: ip,np
      REAL(8), DIMENSION(SIZE(u,1)) :: term
      REAL(8), DIMENSION(SIZE(u,1),3) :: dterm
      REAL(8), DIMENSION(SIZE(u,1),3,3) :: ddterm
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      np=SIZE(uu,4)-1
      CALL horner_2d(uu(:,:,:,np),order,x,y,u,dterm,ddterm)
      IF(order > 0)THEN
         du(:,1:2)=dterm(:,1:2)
         du(:,3)=np*u
      ENDIF
      IF(order > 1)THEN
         ddu(:,1:2,1:2)=ddterm(:,1:2,1:2)
         ddu(:,1:2,3)=np*dterm(:,1:2)
         ddu(:,3,3)=np*(np-1)*u
      ENDIF
c-----------------------------------------------------------------------
c     iterate.
c-----------------------------------------------------------------------
      DO ip=np-1,0,-1

         CALL horner_2d(uu(:,:,:,ip),order,x,y,term,dterm,ddterm)
         u=u*z+term

         IF(order < 1)CYCLE
         du(:,1:2)=du(:,1:2)*z+dterm(:,1:2)
         IF(ip > 0)du(:,3)=du(:,3)*z+ip*term

         IF(order < 2)CYCLE
         ddu(:,1:2,1:2)=ddu(:,1:2,1:2)*z+ddterm(:,1:2,1:2)
         IF(ip > 0)ddu(:,1:2,3)=ddu(:,1:2,3)*z+ip*dterm(:,1:2)
         IF(ip > 1)ddu(:,3,3)=ddu(:,3,3)*z+ip*(ip-1)*term

      ENDDO
c-----------------------------------------------------------------------
c     symmetrize second derivatives.
c-----------------------------------------------------------------------
      ddu(:,3,1:2)=ddu(:,1:2,3)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE horner_3d
c-----------------------------------------------------------------------
c     subprogram 4. horner_transform.
c     transforms basis amplitudes to polynomial coefficients.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE horner_transform(nodal,np,uu,name)

      LOGICAL, INTENT(IN) :: nodal
      INTEGER, INTENT(IN) :: np
      REAL(r8), DIMENSION(:,0:,0:,0:), INTENT(INOUT) :: uu
      CHARACTER(*), INTENT(IN) :: name

      INTEGER :: ip,jp,np1,info,iqty,nqty,ix,iy,iz,jx,jy,jz,nx,ny,nz
      INTEGER, DIMENSION(0:np) :: ipiv
      REAL(r8) :: z,z1
      REAL(r8), DIMENSION(0:np,0:np) :: pmat,tmat
      TYPE(jacobi_type) :: basis
c-----------------------------------------------------------------------
c     compute transformation matrix.
c-----------------------------------------------------------------------
      CALL jacobi_alloc(basis,np,nodal,.FALSE.,"gll")
      DO ip=0,np
         z=basis%qzero(ip)
         z1=(z+1)/2
         CALL jacobi_basis(z,basis)
         tmat(ip,:)=basis%pb
         pmat(ip,:)=(/(z1**jp,jp=0,np)/)
      ENDDO
      CALL jacobi_dealloc(basis)
      np1=np+1
      CALL dgetrf(np1,np1,pmat,np1,ipiv,info)
      CALL dgetrs('N',np1,np1,pmat,np1,ipiv,tmat,np1,info)
c-----------------------------------------------------------------------
c     compute sizes.
c-----------------------------------------------------------------------
      nqty=SIZE(uu,1)
      nx=SIZE(uu,2)/(np+1)
      ny=SIZE(uu,3)/(np+1)
      nz=SIZE(uu,4)/(np+1)
c-----------------------------------------------------------------------
c     start loops over quantities and cells.
c-----------------------------------------------------------------------
      DO iqty=1,nqty
         jx=0
         DO ix=1,nx
            jy=0
            DO iy=1,nx
               jz=0
               DO iz=1,nz
c-----------------------------------------------------------------------
c     transform in x.
c-----------------------------------------------------------------------
                  DO ip=0,np
                     DO jp=0,np
                        uu(iqty,jx:jx+np,jy+ip,jz+jp)
     $                       =MATMUL(tmat,uu(iqty,jx:jx+np,jy+ip,jz+jp))
                     ENDDO
                  ENDDO
c-----------------------------------------------------------------------
c     transform in y.
c-----------------------------------------------------------------------
                  DO ip=0,np
                     DO jp=0,np
                        uu(iqty,jx+ip,jy:jy+np,jz+jp)
     $                       =MATMUL(tmat,uu(iqty,jx+ip,jy:jy+np,jz+jp))
                     ENDDO
                  ENDDO
c-----------------------------------------------------------------------
c     transform in z.
c-----------------------------------------------------------------------
                  DO ip=0,np
                     DO jp=0,np
                        uu(iqty,jx+jp,jy+ip,jz:jz+np)
     $                       =MATMUL(tmat,uu(iqty,jx+jp,jy+ip,jz:jz+np))
                     ENDDO
                  ENDDO
c-----------------------------------------------------------------------
c     finish loops over cells.
c-----------------------------------------------------------------------
                  jz=jz+np1
               ENDDO
               jy=jy+np1
            ENDDO
            jx=jx+np1
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE horner_transform
c-----------------------------------------------------------------------
c     subprogram 6. horner_pointeval.
c     evaluate variable and its gradients at a given logical position.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE horner_pointeval(np,nxyz,kep,uu,uw,duw,dduw)

      INTEGER, INTENT(IN) :: np
      INTEGER, DIMENSION(3), INTENT(IN) :: nxyz
      REAL(r8), DIMENSION(ndim), INTENT(IN) :: kep
      REAL(r8), DIMENSION(:,0:,0:,0:), INTENT(IN) :: uu
      REAL(r8), DIMENSION(:), INTENT(OUT) :: uw
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: duw
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: dduw

      INTEGER, DIMENSION(3) :: ii,jj
      REAL(r8), DIMENSION(3) :: dx
c-----------------------------------------------------------------------
c     find grid cell, index, and offset.
c-----------------------------------------------------------------------
      ii=MIN(MAX(INT(kep),0),nxyz-1)
      jj=ii*(np+1)
      dx=kep-ii
c-----------------------------------------------------------------------
c     evaluate polynomials.
c-----------------------------------------------------------------------
      CALL horner_3d(uu(:,jj(1):jj(1)+np,jj(2):jj(2)+np,jj(3):jj(3)+np),
     $     2,dx(1),dx(2),dx(3),uw,duw,dduw)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE horner_pointeval
      END MODULE horner_mod
