      implicit real*8(a-h,o-z)
C
C     Program to read molecular coordinates in arbitrary orientation,
C     transform them to center of mass, align them along principal
C     axes of the inertia tensor, and possibly rescale them.
C
      Parameter (MAXAT=100)
      Dimension c(3,MAXAT),cp(3,MAXAT),ian(MAXAT),wt(MAXAT),scale(3),
     $  atmass(18),t(6),e(3),v(3,3),cm(3),jat(20)
      Character*80 Title
      character*2 iat
      Logical DoScal
c     atomic masses mostly from g03 output
      Data AtMass
     * /1.0078250D0, 4.002602D0, 6.941D0, 9.012182D0, 10.811D0,
     *  12.0107D0, 14.0067D0, 15.9994D0, 18.9984032D0, 20.1797D0,
     *  22.989770D0, 24.3050D0, 26.981538D0, 28.0855D0,
     *  30.973761D0, 32.065D0, 35.453D0, 39.948D0/
C     * /1.0078250D0, 4*0.d0,
C     * 12.0D0, 14.0030740d0, 15.9949146d0, 18.9984033d0,
C     *  2.014102d0,  4*0.d0,
C     * 13.00335094D0, 2*0.d0, 18.9984046D0/
      Data Conv, Bohr, Zero, Half, One, Two
     $ /505379.006D0, 0.5291771D0, 0.d0, 0.5d0, 1.d0, 2.d0/
C
C     Read the coordinates
C
      read(5,'(a)') title
      write(6,'(a)') title
C      read(5,*) nat, niso
      read(5,*) nat
      do i = 1, nat
        read(5,*) iat, (c(k,i),k=1,3)
        call getan(iat, ian(i))
      enddo
      write(6,'(/" Input Coordinates")')
      call bndang(nat,ian,c)
      doscal = .false.
      read(5,*,end=10) xpa, xpb, xpc
      doscal = .true.
c
c     make the inertia tensor and diagonalize it
c
   10 do i = 1, nat
        wt(i) = atmass(ian(i))
      enddo
      call inertia(nat,wt,c,cm,t,e,v)
c     write(6,'(/" Center of mass")')
c     write(6,102) cm
c     write(6,'(/" Inertia tensor")')
c     write(6,101) t
c     write(6,'(/" Eigenvalues of the Inertia tensor")')
c     write(6,102) e
c     write(6,'(/" Eigenvectors of the Inertia tensor")')
c     write(6,102) v
C
C     Calculate rotation constants and planar moments
C
      RCA = Conv/E(1)
      RCB = Conv/E(2)
      RCC = Conv/E(3)
      Pa = Half*(E(3)+E(2)-E(1))
      Pb = Half*(E(1)+E(3)-E(2))
      Pc = Half*(E(1)+E(2)-E(3))
      Write(6,104) RCA,E(1),PA,  RCB,E(2),PB,  RCC,E(3),PC
      print *, "moo"
      dtoau=0.529177249d0
      do i=1,3
          E(i) = E(i)/(dtoau*dtoau)
      enddo
      write(6,107) E(1), E(2), E(3)
C
C   Transform to principal coordinates
C
      write(6,'(/" Coordinates after principal axis transformation")')
      call trnsp(3,v)
      do i = 1, nat
c       Call DGEMV('t', 3, 3, one, v, 3, c(1,i), 1, zero, cp(1,i), 1)
        call mxv(3,3,v,c(1,i),cp(1,i))
        call getsy(ian(i),iat)
        write(6,106) iat, (cp(k,i),k=1,3)
      enddo
C
C     Use experimental moments to scale coordinates
C
      if(.not.doscal) goto 99
      write(6,'(/" Experimental second moments")')
      write(6,102) xpa, xpb, xpc
      scale(1) = sqrt(xpa/pa)
      scale(2) = sqrt(xpb/pb)
      scale(3) = sqrt(xpc/pc)
      rms = sqrt((scale(1)-1d0)**2+(scale(2)-1d0)**2
     $         + (scale(3)-1d0)**2)/sqrt(3d0)
      write(6,'(/" Scale factors and RMS deviation from unity")')
      write(6,102) scale, rms
      write(6,'(/" Coordinates after scaling")')
      do i = 1, nat
        do k = 1, 3
          cp(k,i) = cp(k,i)*scale(k)
        enddo
c       write(6,100) i, ian(i), (cp(k,i),k=1,3)
      enddo
      call bndang(nat,ian,cp)
      write(6,'(/" Moments after scaling")')
      call inertia(nat,wt,cp,cm,t,e,v)
      RCA = Conv/E(1)
      RCB = Conv/E(2)
      RCC = Conv/E(3)
      Pa = Half*(E(3)+E(2)-E(1))
      Pb = Half*(E(1)+E(3)-E(2))
      Pc = Half*(E(1)+E(2)-E(3))
      Write(6,104) RCA,E(1),PA,  RCB,E(2),PB,  RCC,E(3),PC
c
c     read(5,*,end=99) niso, (jat(k), k=1, niso)
C      if(niso.eq.0) call exit
C      read(5,*) (jat(k), k=1, niso)
C      do j = 1, niso
C        do i = 1, nat
C          wt(i) = atmass(ian(i),1)
C        enddo
C        wt(jat(j)) = atmass(ian(jat(j)),2)
C        write(6,105) jat(j), wt(jat(j))
C        call inertia(nat,wt,cp,cm,t,e,v)
C        RCA = Conv/E(1)
C        RCB = Conv/E(2)
C        RCC = Conv/E(3)
C        Pa = Half*(E(3)+E(2)-E(1))
C        Pb = Half*(E(1)+E(3)-E(2))
C        Pc = Half*(E(1)+E(2)-E(3))
C        Write(6,104) RCA,E(1),PA,  RCB,E(2),PB,  RCC,E(3),PC
C      enddo
   99 stop
c
  100 format(2i5,3f12.8)
  101 format(f12.6/2f12.6/3f12.6)
  102 format(4f12.6)
  103 format(3i5,3f12.8)
  104 format(1X/' Rotational constants/MHz and moments of inertia/',
     * '(amu A**2):'//
     * '  A =',F13.5,'   Ia =',F13.6,'   Pa =',F13.6/
     * '  B =',F13.5,'   Ib =',F13.6,'   Pb =',F13.6/
     * '  C =',F13.5,'   Ic =',F13.6,'   Pc =',F13.6)
  105 format(/' New mass for atom no',i3,' is',f15.8)
  106 format(5x,a5,3f15.6)
  107 format(1X/' Moments of inertial (amu bohr**2)',//
     * 4X,'XX= ',F10.5,', YY= ',F10.5,', ZZ= ',F10.5)
C
      end
*deck inertia
      subroutine inertia(nat,atmass,c,cm,t,e,v)
      implicit real*8(a-h,o-z)
      dimension atmass(*),c(3,*),cm(*),e(*),v(3,*),t(*),w(9)
C
C    Transformation to center of mass coordinates
C
      save zero, one
      data zero/0.d0/, one/1.d0/
c
      wtot = zero
      do k = 1, 3
        cm(k) = zero
      enddo
      do i = 1, nat
        wtot = wtot + atmass(i)
        call daxpy(3,atmass(i),c(1,i),1,cm,1)
      enddo
      call dscal(3,one/wtot,cm,1)
      do i = 1, nat
        call daxpy(3,-one,cm,1,c(1,i),1)
      enddo
C
C   Form inertia tensor and diagonalize it
C
      do k = 1, 6
        t(k) = zero
      enddo
      do i = 1, nat
        t(1) = t(1) + atmass(i)*(c(2,i)**2 + c(3,i)**2)
        t(2) = t(2) - atmass(i)*c(1,i)*c(2,i)
        t(3) = t(3) + atmass(i)*(c(1,i)**2 + c(3,i)**2)
        t(4) = t(4) - atmass(i)*c(1,i)*c(3,i)
        t(5) = t(5) - atmass(i)*c(2,i)*c(3,i)
        t(6) = t(6) + atmass(i)*(c(1,i)**2 + c(2,i)**2)
      enddo
      call tri2sq(3,3,v,t)
c     Call DSYEV( 'V', 'L', 3, V, 3, E, W, 9, INFO )
      call dsyevj3(v,w,e)
      call dcopy(9,w,1,v,1)
      return
      end
*deck tri2sq
      subroutine tri2sq(nd,n,a,b)
      implicit real*8(a-h,o-z)
      dimension b(*), a(nd,*)
c
c     convert triangular matrix to square
c
      ij = 0
      do i=1,n
        do j=1,i
          ij = ij + 1
          a(i,j) = b(ij)
          a(j,i) = b(ij)
        enddo
      enddo
      return
      end
*deck getsy
      subroutine getsy(ian,iat)
      implicit real*8(a-h,o-z)
      character*2 el(79),iat
c
c     get atomic symbol from atomic number
c
      data el/
     * 'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     * 'Na','Mg','Al','Si','P ','S ','Cl','Ar',
     * 60*'  ','Au'/
      iat = el(ian)
      return
      end
*deck getan
      subroutine getan(iat,ian)
      implicit real*8(a-h,o-z)
      character*2 el(79),iat
c
c     get atomic number from atomic symbol
c
      data el/
     * 'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     * 'Na','Mg','Al','Si','P ','S ','Cl','Ar',
     * 60*'  ','Au'/
      do i = 1,18
        if(iat.eq.el(i)) then
          ian = i
          return
        endif
      enddo
      ian = 0
      return
      end
*deck mxv
      subroutine mxv(n,ndim,c,v,vp)
      implicit real*8(a-h,o-z)
c
c     matrix vector product vp = c * v
c
      dimension c(ndim,*), v(*), vp(*)
c
      do i=1,n
        vp(i) = 0.d0
        do j=1,n
          vp(i) = vp(i) + c(i,j)*v(j)
        enddo
      enddo
      return
      end
*deck trnsp
      subroutine trnsp(n,a)
      implicit real*8(a-h,o-z)
      dimension a(n,*)
c
c     transpose square matrix
c
      do i=2,n
        do j=1,(i-1)
          t = a(i,j)
          a(i,j) = a(j,i)
          a(j,i) = t
        enddo
      enddo
      return
      end
c
*deck bndang
      subroutine bndang(natom,ian,c)
      implicit real*8(a-h,o-z)
      parameter (MAXATM=100) 
      character*2 iat
      dimension ian(*),c(3,*),d(MAXATM,MAXATM)
      dimension ib(MAXATM,5),nb(MAXATM)
c
c     input Cartesian coordinates and find bond lengths and angles
c
  100 format(5x,a5,3f15.6)
  101 format(18f8.3)
  104 format(5x,2i3,f10.6)
  105 format(5x,3i3,2x,f8.4)
  106 format(5x,4i3,2x,f9.4)
c
      if(natom.gt.MAXATM) then
        write(6,*)' too many atoms!'
        stop
      endif
c
c     read coordinates in Gaussian format
c
      write(6,'(/"                  Atomic Coordinates")')
      do k=1,natom
        call getsy(ian(k),iat)
        write(6,100) iat, (c(i,k),i=1,3)
      enddo
c
c     compute distance matrix
c
      do k=2,natom
        do j=1,k-1
          d(j,k) = sqrt((c(1,j)-c(1,k))**2 + (c(2,j)-c(2,k))**2 +
     $                  (c(3,j)-c(3,k))**2)
          d(k,j) = d(j,k)
        enddo
      enddo
c     write(6,*) '                  Distance Matrix'
c     do k=2,natom
c       write(6,101) (d(i,k),i=1,k-1)
c     enddo
c
c     find bonded atoms
c
      call bonds(natom,MAXATM,ian,nb,ib,d)
c     write(6,*) 'number of bonds'
c     write(6,*) (nb(k),k=1,natom)
c     do k=1,natom
c       write(6,*) k, ':',(ib(k,j),j=1,nb(k))
c     enddo
      write(6,'(/"         Bonds")')
      do k=1,natom
        do j=1,nb(k)
          if(ib(k,j).gt.k) write(6,104) k,ib(k,j),d(k,ib(k,j))
        enddo
      enddo
c
c     compute bond angles
c
      write(6,'(/"         Bond Angles")')
      do i=1,natom
        if(nb(i).gt.1) then
          do ja=1,nb(i)-1
            j = ib(i,ja)
            do ka=ja+1,nb(i)
              k = ib(i,ka)
              bonda = ang(j,i,k,c)
              write(6,105) j,i,k,bonda
            enddo
          enddo
        endif
      enddo
c
c     compute dihedral angles
c
      write(6,'(/"         Dihedral Angles")')
      do i=1,natom
        if(nb(i).ge.1) then
          do ja=1,nb(i)
            j = ib(i,ja)
            if(nb(j).ge.1) then
              do ka=1,nb(j)
                k = ib(j,ka)
                if(nb(k).ge.1.and.k.ne.i) then
                  do la=1,nb(k)
                    l = ib(k,la)
                    if(l.gt.i.and.l.ne.j) then
                      dih = dihed(i,j,k,l,c)
                      write(6,106) i,j,k,l,dih
                    endif
                  enddo
                endif
              enddo
            endif
          enddo
        endif
      enddo
      end
c
      subroutine bonds(natom,ndim,ian,nb,ib,d)
      implicit real*8(a-h,o-z)
      dimension ian(*),nb(*),ib(ndim,*),d(ndim,*)
c
c     find bonds
c
      do k=1,natom
        nb(k) = 0
        do j=1,natom
          djk = 1.6d0
          if(j.ne.k) then
            if(ian(j).eq.1.or.ian(k).eq.1) djk = 1.2d0
            if(d(j,k).lt.djk) then
              nb(k) = nb(k) + 1
              ib(k,nb(k)) = j
            endif
          endif
        enddo
      enddo
      return
      end
c
*deck dihed
      real*8 function dihed(i,j,k,l,c)
      implicit real*8(a-h,o-z)
c
c     compute i-j-k-l dihedral angle (in degrees)
c
      dimension c(3,*), b1(3), b2(3), b3(3), b4(3), b12(3), b34(3)
      save pi
      data pi/3.141592653589793d0/
c
      call vsub(b1,c(1,i),c(1,j))
      call vsub(b2,c(1,k),c(1,j))
      call vsub(b3,c(1,j),c(1,k))
      call vsub(b4,c(1,l),c(1,k))
      r2 = sqrt(dot(b2,b2))
      call cross(b12,b1,b2)
      call cross(b34,b3,b4)
      call cross(b1,b12,b34)
      dihed = atan2(dot(b1,b2), r2*dot(b12,b34))*180.d0/pi
      return
      end
c
*deck ang
      real*8 function ang(i,j,k,c)
      implicit real*8(a-h,o-z)
c
c     compute i-j-k bond angle (in degrees)
c
      dimension c(3,*),a(3),b(3)
      save pi
      data pi/3.141592653589793d0/
c
      call vsub(a,c(1,i),c(1,j))
      call vsub(b,c(1,k),c(1,j))
      ra = sqrt(dot(a,a))
      rb = sqrt(dot(b,b))
      costh = dot(a,b)/(ra*rb)
      ang = acos(costh)*180.d0/pi
      return
      end
c
*deck dot
      real*8 function dot(a,b)
      implicit real*8(a-h,o-z)
c
c     vector dot product
c
      dimension a(3),b(3)
      dot = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
      return
      end
c
*deck cross
      subroutine cross(c,a,b)
      implicit real*8(a-h,o-z)
c
c     vector cross product
c
      dimension a(3),b(3),c(3)
      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
      return
      end
c
*deck vsub
      subroutine vsub(c,a,b)
      implicit real*8(a-h,o-z)
c
c     c = a - b
c
      dimension a(3),b(3),c(3)
      do k=1,3
        c(k) = a(k) - b(k)
      enddo
      return
      end
*deck dsyevj3
* ----------------------------------------------------------------------------
* Numerical diagonalization of 3x3 matrcies
* Copyright (C) 2006  Joachim Kopp
* ----------------------------------------------------------------------------
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
* ----------------------------------------------------------------------------


* ----------------------------------------------------------------------------
      SUBROUTINE DSYEVJ3(A, Q, W)
* ----------------------------------------------------------------------------
* Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
* matrix A using the Jacobi algorithm.
* The upper triangular part of A is destroyed during the calculation,
* the diagonal elements are read but not destroyed, and the lower
* triangular elements are not referenced at all.
* ----------------------------------------------------------------------------
* Parameters:
*   A: The symmetric input matrix
*   Q: Storage buffer for eigenvectors
*   W: Storage buffer for eigenvalues
* ----------------------------------------------------------------------------
*     .. Arguments ..
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION Q(3,3)
      DOUBLE PRECISION W(3)

*     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 3 )
    
*     .. Local Variables ..
      DOUBLE PRECISION SD, SO
      DOUBLE PRECISION S, C, T
      DOUBLE PRECISION G, H, Z, THETA
      DOUBLE PRECISION THRESH
      INTEGER          I, X, Y, R

*     Initialize Q to the identitity matrix
*     --- This loop can be omitted if only the eigenvalues are desired ---
      DO 10 X = 1, N
        Q(X,X) = 1.0D0
        DO 11, Y = 1, X-1
          Q(X, Y) = 0.0D0
          Q(Y, X) = 0.0D0
   11   CONTINUE
   10 CONTINUE

*     Initialize W to diag(A)
      DO 20 X = 1, N
        W(X) = A(X, X)
   20 CONTINUE

*     Calculate SQR(tr(A))  
      SD = 0.0D0
      DO 30 X = 1, N
        SD = SD + ABS(W(X))
   30 CONTINUE
      SD = SD**2
 
*     Main iteration loop
      DO 40 I = 1, 50
*       Test for convergence
        SO = 0.0D0
        DO 50 X = 1, N
          DO 51 Y = X+1, N
            SO = SO + ABS(A(X, Y))
   51     CONTINUE
   50   CONTINUE
        IF (SO .EQ. 0.0D0) THEN
          RETURN
        END IF

        IF (I .LT. 4) THEN
          THRESH = 0.2D0 * SO / N**2
        ELSE
          THRESH = 0.0D0
        END IF

*       Do sweep
        DO 60 X = 1, N
          DO 61 Y = X+1, N
            G = 100.0D0 * ( ABS(A(X, Y)) )
            IF ( I .GT. 4 .AND. ABS(W(X)) + G .EQ. ABS(W(X))
     $                    .AND. ABS(W(Y)) + G .EQ. ABS(W(Y)) ) THEN
              A(X, Y) = 0.0D0
            ELSE IF (ABS(A(X, Y)) .GT. THRESH) THEN
*             Calculate Jacobi transformation
              H = W(Y) - W(X)
              IF ( ABS(H) + G .EQ. ABS(H) ) THEN
                T = A(X, Y) / H
              ELSE
                THETA = 0.5D0 * H / A(X, Y)
                IF (THETA .LT. 0.0D0) THEN
                  T = -1.0D0 / (SQRT(1.0D0 + THETA**2) - THETA)
                ELSE
                  T = 1.0D0 / (SQRT(1.0D0 + THETA**2) + THETA)
                END IF
              END IF

              C = 1.0D0 / SQRT( 1.0D0 + T**2 )
              S = T * C
              Z = T * A(X, Y)
              
*             Apply Jacobi transformation
              A(X, Y) = 0.0D0
              W(X)    = W(X) - Z
              W(Y)    = W(Y) + Z
              DO 70 R = 1, X-1
                T       = A(R, X)
                A(R, X) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   70         CONTINUE
              DO 80, R = X+1, Y-1
                T       = A(X, R)
                A(X, R) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   80         CONTINUE
              DO 90, R = Y+1, N
                T       = A(X, R)
                A(X, R) = C * T - S * A(Y, R)
                A(Y, R) = S * T + C * A(Y, R)
   90         CONTINUE

*             Update eigenvectors
*             --- This loop can be omitted if only the eigenvalues are desired ---
              DO 100, R = 1, N
                T       = Q(R, X)
                Q(R, X) = C * T - S * Q(R, Y)
                Q(R, Y) = S * T + C * Q(R, Y)
  100         CONTINUE
            END IF
   61     CONTINUE
   60   CONTINUE
   40 CONTINUE

      PRINT *, "DSYEVJ3: No convergence."
            
      END SUBROUTINE
* End of subroutine DSYEVJ3
*deck drot
      SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
*     .. Scalar Arguments ..
      DOUBLE PRECISION C,S
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     DROT applies a plane rotation.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*       code for both increments equal to 1
*
         DO I = 1,N
            DTEMP = C*DX(I) + S*DY(I)
            DY(I) = C*DY(I) - S*DX(I)
            DX(I) = DTEMP
         END DO
      ELSE
*
*       code for unequal increments or equal increments not equal
*         to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = C*DX(IX) + S*DY(IY)
            DY(IY) = C*DY(IY) - S*DX(IX)
            DX(IX) = DTEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
*deck dscal
      SUBROUTINE DSCAL(N,DA,DX,INCX)
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
*  Purpose
*  =======
*
*     DSCAL scales a vector by a constant.
*     uses unrolled loops for increment equal to one.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
*
*        code for increment equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               DX(I) = DA*DX(I)
            END DO
            IF (N.LT.5) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
            DX(I) = DA*DX(I)
            DX(I+1) = DA*DX(I+1)
            DX(I+2) = DA*DX(I+2)
            DX(I+3) = DA*DX(I+3)
            DX(I+4) = DA*DX(I+4)
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            DX(I) = DA*DX(I)
         END DO
      END IF
      RETURN
      END
*deck dcopy
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     DCOPY copies a vector, x, to a vector, y.
*     uses unrolled loops for increments equal to one.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,7)
         IF (M.NE.0) THEN
            DO I = 1,M
               DY(I) = DX(I)
            END DO
            IF (N.LT.7) RETURN
         END IF   
         MP1 = M + 1
         DO I = MP1,N,7
            DY(I) = DX(I)
            DY(I+1) = DX(I+1)
            DY(I+2) = DX(I+2)
            DY(I+3) = DX(I+3)
            DY(I+4) = DX(I+4)
            DY(I+5) = DX(I+5)
            DY(I+6) = DX(I+6)
         END DO
      ELSE      
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DY(IY) = DX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
*deck daxpy
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     DAXPY constant times a vector plus a vector.
*     uses unrolled loops for increments equal to one.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (DA.EQ.0.0d0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,4)
         IF (M.NE.0) THEN
            DO I = 1,M
               DY(I) = DY(I) + DA*DX(I)
            END DO
         END IF
         IF (N.LT.4) RETURN
         MP1 = M + 1
         DO I = MP1,N,4
            DY(I) = DY(I) + DA*DX(I)
            DY(I+1) = DY(I+1) + DA*DX(I+1)
            DY(I+2) = DY(I+2) + DA*DX(I+2)
            DY(I+3) = DY(I+3) + DA*DX(I+3)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
          DY(IY) = DY(IY) + DA*DX(IX)
          IX = IX + INCX
          IY = IY + INCY
         END DO
      END IF
      RETURN
      END
