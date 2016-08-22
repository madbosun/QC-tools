      implicit real*8(a-h,o-z)
C
C     Program to read molecular coordinates in arbitrary orientation,
C     transform them to center of mass, align them along principal
C     axes of the inertia tensor, and calculate distance from the CM.
C     Find Abelian point group and unique coordinates
C
      Parameter (MAXAT=1000)
      Dimension c(3,MAXAT), cp(3,MAXAT), ian(MAXAT), jan(MAXAT)
      dimension iua(MAXAT), indx(MAXAT), wt(MAXAT), d(MAXAT)
      dimension ip(MAXAT,7), atmass(100)
      dimension iop(7), t(6), e(3), v(3,3), cm(3), isym(3,7)
      dimension u(3,3)
      Character*80 Title
      character*32 arg
      character*8 lab(7)
      character*3 group
      character*2 ich(MAXAT), jch(MAXAT)
      logical ctest, dbg
      data isym/
     *   -1,1,1, 1,-1,1, 1,1,-1, -1,-1,1, -1,1,-1, 1,-1,-1, -1,-1,-1/
      data lab/
     *  'sig(yz) ','sig(xz) ','sig(xy) ','C2(z)   ',
     *  'C2(y)   ','C2(x)   ','i       '/
c     atomic masses from gamess
      Data AtMass/
     *   1.007825D+00,4.0026D+00,7.01600D+00,9.01218D+00,11.00931D+00,
     *   12.0D+00,14.00307D+00,15.99491D+00,18.99840D+00,19.99244D+00,
     *   22.9898D+00,23.98504D+00,26.98153D+00,27.97693D+00,
     *   30.97376D+00,31.97207D+00,34.96885D+00,39.948D+00,
     *   38.96371D+00,39.96259D+00,44.95592D+00,47.90D+00,50.9440D+00,
     *   51.9405D+00,54.9381D+00,55.9349D+00,58.9332D+00,57.9353D+00,
     *   62.9298D+00,63.9291D+00,68.9257D+00,73.9219D+00,74.9216D+00,
     *   79.9165D+00,78.9183D+00,83.9115D+00,
     *   84.9117D+00,87.9056D+00,89.9054D+00,89.9043D+00,92.9060D+00,
     *   97.9055D+00,97.0D+00,101.9037D+00,102.9048D+00,105.9032D+00,
     *   106.9041D+00,113.9036D+00,114.9041D+00,119.9022D+00,
     *   120.9038D+00,129.9067D+00,126.9044D+00,131.9042D+00,
     *   24*0.d0, 196.9666d0, 21*0.d0/
      Data Zero, One/0.d0, 1.d0/
      data tol/1.d-5/
C
C     Process command line
C
      dbg = .false.
      if(iargc().gt.0) then
        do i = 1, iargc()
          call getarg(i, arg)
          if (arg(1:1).eq.'d') then
            dbg = .true.
          else
            read(arg,'(f15.5)') tolx
            tol = tolx
          endif
        enddo
      endif
C      write(6,*) 'tol =', tol
C
C     Read the coordinates
C
      read(5,'(a)') title
C      write(6,'(a)') title
      read(5,*) nat
      do i = 1, nat
        read(5,*) ich(i), rtemp, (c(k,i),k=1,3)
        ian(i) = int(rtemp)
      enddo
C      write(6,*) 'Number of atoms', nat
      do i = 1, nat
C        write(6,100) i, ian(i), (c(k,i),k=1,3), atmass(ian(i))
c       write(6,111) i, ian(i), ity, (c(k,i),k=1,3)
      enddo
C
C     Make the inertia tensor and diagonalize it
C
      do i = 1, nat
        wt(i) = atmass(ian(i))
      enddo
      call inertia(nat,wt,c,cm,t,e,v)
C      write(6,'(/" Center of mass")')
C      write(6,102) cm
C      write(6,'(/" Inertia tensor")')
C      write(6,101) t
C      write(6,'(/" Eigenvalues of the Inertia tensor")')
C      write(6,102) e
C      write(6,'(/" Eigenvectors of the Inertia tensor")')
C      write(6,102) v
c debug
C      do i = 1, 3
C        do j = 1, 3
C          u(i,j) = zero
C          do k = 1, 3
C            u(i,j) = u(i,j) + v(k,i)*v(k,j)
C          enddo
C        enddo
C      enddo
C      write(6,'(/" U(dag)*U")')
C      write(6,102) u
c debug
C
C    Transform to principal axis coordinates
C
   10 continue
C      write(6,'(/" Coordinates after principal axis transformation")')
      do i = 1, nat
        call dgemv('t',3,3,one,v,3,c(1,i),1,zero,cp(1,i),1)
c       do k = 1, 3
c         cp(k,i) = ddot(3,v(1,k),1,c(1,i),1)
c       enddo
C        write(6,100) i, ian(i), (cp(k,i),k=1,3)
      enddo
c debug
      call inertia(nat,wt,cp,cm,t,e,v)
C      write(6,'(/" Inertia tensor")')
C      write(6,101) t
C      write(6,'(/" Eigenvalues of the Inertia tensor")')
C      write(6,102) e
C      write(6,'(/" Eigenvectors of the Inertia tensor")')
C      write(6,102) v
      vmax = 0.d0
      do j = 1, 3
        do k = 1, 3
          u(j,k) = abs(v(j,k))
        enddo
      enddo
      vmax = max(u(2,1),u(3,1),u(1,2),u(3,2),u(1,3),u(2,3))
C      write(6,*) 'vmax =', vmax
      if (vmax.gt.1.d-18) then
        call dcopy(3*nat,cp,1,c,1)
        goto 10
      endif
c debug
C
C     Calculate distance from center of mass
C
C      write(6,*) 'Distance from center of mass, in input order'
      do i = 1, nat
        d(i) = sqrt(cp(1,i)*cp(1,i) + cp(2,i)*cp(2,i) + cp(3,i)*cp(3,i))
C        write(6,110) i, d(i)
      enddo
C
C     Sort the distances and reorder the atoms
C
      call sortrx(nat,d,indx)
C      write(6,*) 'Coordinates ordered by center of mass distance'
      do j = 1, nat
        i = indx(j)
C        write(6,100) j, ian(i), (cp(k,i),k=1,3)
        jan(j) = ian(i)
        jch(j) = ich(i)
        do k = 1, 3
          c(k,j) = cp(k,i)
        enddo
      enddo
C
C     Find Abelian symmetry operations
C
      nop = 0
      do i = 1, 7
        iop(i) = 0
        call symop(nat,c,cp,isym(1,i),isym(2,i),isym(3,i))
C        if (dbg) write(6,*) 'after op ',i ,isym(1,i),isym(2,i),isym(3,i)
        do j = 1, nat
c         write(6,'(9f12.8)')(cp(k,j),k=1,3),(c(k,j)-cp(k,j),k=1,3),
c    $     (c(k,j)+cp(k,j),k=1,3)
        enddo
        if (ctest(nat,cp,c,ip(1,nop+1),tol,dbg)) then
          iop(i) = 1
          nop = nop + 1
C          write(6,'(1x,a8,''found'')') lab(i)
        endif
      enddo
C      write(6,*) nop, ' symmetry operations found'
C
C     Find symmetry unique atoms
C
      if(nop.gt.0) then
        do i = 1, nat
          iua(i) = 99999
          do k = 1, nop
            if(ip(i,k).lt.iua(i)) iua(i) = ip(i,k)
          enddo
        enddo
C        if (dbg) write(6,105)(iua(k),k=1,nat)
        nq = 1
        do i = 2, nat
          if(iua(i).gt.iua(nq)) then
            nq = nq + 1
            iua(nq) = iua(i)
          endif
        enddo
C        if (dbg) write(6,105)(iua(k),k=1,nq)
      else
        nq = nat
        do i = 1, nat
          iua(i) = i
        enddo
      endif
C
C     Determine point group
C
      group = 'C1 '
      if(nop.eq.7) then
        group = 'D2h'
      else if(nop.eq.3) then
        if(iop(4).eq.1.and.iop(5).eq.1.and.iop(6).eq.1) then
          group = 'D2 '
        else if(iop(7).eq.1) then
          group = 'C2h'
        else
          group = 'C2v'
        endif
      else if(nop.eq.1) then
        if(iop(7).eq.1) then
          group = 'Ci '
        else if(iop(4).eq.1.or.iop(5).eq.1.or.iop(6).eq.1) then
          group = 'C2 '
        else
          group = 'Cs '
        endif
      endif
C
C     Reorient if necessary
C
      if(group.eq.'C2h'.or.group.eq.'C2v'.or.group.eq.'C2 ') then
        if(iop(5).eq.1) call cswap(c,nat,1,3,2)
        if(iop(6).eq.1) call cswap(c,nat,2,3,1)
      else if(group.eq.'Cs ') then
        if(iop(1).eq.1) call cswap(c,nat,2,3,1)
        if(iop(2).eq.1) call cswap(c,nat,1,3,2)
      endif
C      write(6,*) 'Point group is ', group
C      write(6,*) 'Coordinates of unique atoms'
      write(6,*) '$data'
      write(6,*) trim(title)
      write(6,*) trim(group)
      if(group.ne.'C1') write(6,*)
      do i = 1, nq
        j = iua(i)
        write(6,106) jch(j), jan(j), (c(k,j),k=1,3)
      enddo
      write(6,*) '$end'
c
  100 format(2i5,3f12.8,f12.4)
c 101 format(f15.6/2f15.6/3f15.6)
  101 format(f25.14/2f25.14/3f25.14)
c 102 format(3f15.6)
c 102 format(3f25.14)
  102 format(3e15.6)
  105 format(18i3)
  106 format(1x,a2,i3,3x,3f14.8)
  110 format(' atom',i4,' distance ',f8.2)
C
      end
*deck ctest
      logical function ctest(nat,c,cp,ip,tol,dbg)
      implicit real*8(a-h,o-z)
      logical dbg
      dimension c(3,*),cp(3,*),ip(*)
C
C    Test for symmetry operation
C
      save zero
      data zero/0.d0/
c
      do i = 1, nat
        ip(i) = 0
      enddo
      do i = 1, nat
        if(ip(i).eq.0) then
          dmin = 100.d0
          do j = 1, nat
            d = zero
            do k = 1, 3
              d = d + (c(k,i) - cp(k,j))**2
            enddo
            d = sqrt(d)
            if (d.lt.dmin) then
              dmin = d
              jmin = j
            endif
          enddo
C          write(6,*) i, jmin, dmin, tol
          if(dmin.gt.tol) then
            ctest = .false.
C            if (dbg) write(6,'(''ctest: no match'')')
            return
          endif
c
C          if (dbg) write(6,'(''ctest: match'')')
C          if (dbg) write(6,*) i, jmin, dmin, tol
c         if (dbg) write(6,'(3f12.8)')(c(k,i),k=1,3)
c         if (dbg) write(6,'(3f12.8)')(cp(k,jmin),k=1,3)
          ip(i) = jmin
          ip(jmin) = i
        endif
      enddo
      ctest = .true.
C      if (dbg) write(6,'(''ctest:'',18i3)') (ip(k),k=1,nat)
      return
      end
*deck cswap
      subroutine cswap(c,nat,ix,iy,iz)
      implicit real*8(a-h,o-z)
      dimension c(3,*),t(3)
C
C     Interchange coordinates
C
C      write(6,*) 'cswap:',ix,iy,iz
      do i = 1, nat
        do k = 1, 3
          t(k) = c(k,i)
        enddo
        c(1,i) = t(ix)
        c(2,i) = t(iy)
        c(3,i) = t(iz)
      enddo
      return
      end
*deck symop
      subroutine symop(nat,c,cp,ix,iy,iz)
      implicit real*8(a-h,o-z)
      dimension c(3,*),cp(3,*)
C
C     Apply symmetry operation
C
      do i = 1, nat
        cp(1,i) = c(1,i)*ix
        cp(2,i) = c(2,i)*iy
        cp(3,i) = c(3,i)*iz
      enddo
      return
      end
*deck inertia
      subroutine inertia(nat,atmass,c,cm,t,e,v)
      implicit real*8(a-h,o-z)
      dimension atmass(*),c(3,*),cm(*),e(*),v(3,*),t(*),w(9),z(18)
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
c     lapack
c     call tri2sq(3,3,v,t)
c     Call DSYEV( 'V', 'L', 3, V, 3, E, W, 9, INFO )
c     jacobi
      call tri2sq(3,3,w,t)
c     call dsyevj3(w,v,e)
      call djac(3,3,w,v,e,z)
c     dum = 0.d0
c     idum = 1
c     lwrk = 18
c     call dgesvd('a','n',3,3,w,3,e,v,3,dum,idum,z,lwrk,info)
*     CALL DGESVD('A','N',N,N,A,N,D,U,N,DUM,IDUM,WRK,LWRK,INFO)
      return
      end
*deck tri2sq
      subroutine tri2sq(nd,n,a,b)
      implicit real*8(a-h,o-z)
      dimension a(nd,*), b(*)
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
*deck sortrx
      SUBROUTINE SORTRX(N,DATA,INDEX)
C===================================================================
C
C     SORTRX -- SORT, Real input, indeX output
C
C
C     Input:  N     INTEGER
C             DATA  REAL
C
C     Output: INDEX INTEGER (DIMENSION N)
C
C This routine performs an in-memory sort of the first N elements of
C array DATA, returning into array INDEX the indices of elements of
C DATA arranged in ascending order.  Thus,
C
C    DATA(INDEX(1)) will be the smallest number in array DATA;
C    DATA(INDEX(N)) will be the largest number in DATA.
C
C The original data is not physically rearranged.  The original order
C of equal input values is not necessarily preserved.
C
C===================================================================
C
C SORTRX uses a hybrid QuickSort algorithm, based on several
C suggestions in Knuth, Volume 3, Section 5.2.2.  In particular, the
C "pivot key" [my term] for dividing each subsequence is chosen to be
C the median of the first, last, and middle values of the subsequence;
C and the QuickSort is cut off when a subsequence has 9 or fewer
C elements, and a straight insertion sort of the entire array is done
C at the end.  The result is comparable to a pure insertion sort for
C very short arrays, and very fast for very large arrays (of order 12
C micro-sec/element on the 3081K for arrays of 10K elements).  It is
C also not subject to the poor performance of the pure QuickSort on
C partially ordered data.
C
C Created:  15 Jul 1986  Len Moss
C
C===================================================================
 
      INTEGER   N,INDEX(N)
      REAL*8    DATA(N)
 
      INTEGER   LSTK(31),RSTK(31),ISTK
      INTEGER   L,R,I,J,P,INDEXP,INDEXT
      REAL*8    DATAP
 
C     QuickSort Cutoff
C
C     Quit QuickSort-ing when a subsequence contains M or fewer
C     elements and finish off at end with straight insertion sort.
C     According to Knuth, V.3, the optimum value of M is around 9.
 
      INTEGER   M
      PARAMETER (M=9)
 
C===================================================================
C
C     Make initial guess for INDEX
 
      DO 50 I=1,N
         INDEX(I)=I
   50    CONTINUE
 
C     If array is short, skip QuickSort and go directly to
C     the straight insertion sort.
 
      IF (N.LE.M) GOTO 900
 
C===================================================================
C
C     QuickSort
C
C     The "Qn:"s correspond roughly to steps in Algorithm Q,
C     Knuth, V.3, PP.116-117, modified to select the median
C     of the first, last, and middle elements as the "pivot
C     key" (in Knuth's notation, "K").  Also modified to leave
C     data in place and produce an INDEX array.  To simplify
C     comments, let DATA[I]=DATA(INDEX(I)).
 
C Q1: Initialize
      ISTK=0
      L=1
      R=N
 
  200 CONTINUE
 
C Q2: Sort the subsequence DATA[L]..DATA[R].
C
C     At this point, DATA[l] <= DATA[m] <= DATA[r] for all l < L,
C     r > R, and L <= m <= R.  (First time through, there is no
C     DATA for l < L or r > R.)
 
      I=L
      J=R
 
C Q2.5: Select pivot key
C
C     Let the pivot, P, be the midpoint of this subsequence,
C     P=(L+R)/2; then rearrange INDEX(L), INDEX(P), and INDEX(R)
C     so the corresponding DATA values are in increasing order.
C     The pivot key, DATAP, is then DATA[P].
 
      P=(L+R)/2
      INDEXP=INDEX(P)
      DATAP=DATA(INDEXP)
 
      IF (DATA(INDEX(L)) .GT. DATAP) THEN
         INDEX(P)=INDEX(L)
         INDEX(L)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF
 
      IF (DATAP .GT. DATA(INDEX(R))) THEN
         IF (DATA(INDEX(L)) .GT. DATA(INDEX(R))) THEN
            INDEX(P)=INDEX(L)
            INDEX(L)=INDEX(R)
         ELSE
            INDEX(P)=INDEX(R)
         ENDIF
         INDEX(R)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF
 
C     Now we swap values between the right and left sides and/or
C     move DATAP until all smaller values are on the left and all
C     larger values are on the right.  Neither the left or right
C     side will be internally ordered yet; however, DATAP will be
C     in its final position.
 
  300 CONTINUE
 
C Q3: Search for datum on left >= DATAP
C
C     At this point, DATA[L] <= DATAP.  We can therefore start scanning
C     up from L, looking for a value >= DATAP (this scan is guaranteed
C     to terminate since we initially placed DATAP near the middle of
C     the subsequence).
 
         I=I+1
         IF (DATA(INDEX(I)).LT.DATAP) GOTO 300
 
  400 CONTINUE
 
C Q4: Search for datum on right <= DATAP
C
C     At this point, DATA[R] >= DATAP.  We can therefore start scanning
C     down from R, looking for a value <= DATAP (this scan is guaranteed
C     to terminate since we initially placed DATAP near the middle of
C     the subsequence).
 
         J=J-1
         IF (DATA(INDEX(J)).GT.DATAP) GOTO 400
 
C Q5: Have the two scans collided?
 
      IF (I.LT.J) THEN
 
C Q6: No, interchange DATA[I] <--> DATA[J] and continue
 
         INDEXT=INDEX(I)
         INDEX(I)=INDEX(J)
         INDEX(J)=INDEXT
         GOTO 300
      ELSE
 
C Q7: Yes, select next subsequence to sort
C
C     At this point, I >= J and DATA[l] <= DATA[I] == DATAP <= DATA[r],
C     for all L <= l < I and J < r <= R.  If both subsequences are
C     more than M elements long, push the longer one on the stack and
C     go back to QuickSort the shorter; if only one is more than M
C     elements long, go back and QuickSort it; otherwise, pop a
C     subsequence off the stack and QuickSort it.
 
         IF (R-J .GE. I-L .AND. I-L .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=J+1
            RSTK(ISTK)=R
            R=I-1
         ELSE IF (I-L .GT. R-J .AND. R-J .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=L
            RSTK(ISTK)=I-1
            L=J+1
         ELSE IF (R-J .GT. M) THEN
            L=J+1
         ELSE IF (I-L .GT. M) THEN
            R=I-1
         ELSE
C Q8: Pop the stack, or terminate QuickSort if empty
            IF (ISTK.LT.1) GOTO 900
            L=LSTK(ISTK)
            R=RSTK(ISTK)
            ISTK=ISTK-1
         ENDIF
         GOTO 200
      ENDIF
 
  900 CONTINUE
 
C===================================================================
C
C Q9: Straight Insertion sort
 
      DO 950 I=2,N
         IF (DATA(INDEX(I-1)) .GT. DATA(INDEX(I))) THEN
            INDEXP=INDEX(I)
            DATAP=DATA(INDEXP)
            P=I-1
  920       CONTINUE
               INDEX(P+1) = INDEX(P)
               P=P-1
               IF (P.GT.0) THEN
                  IF (DATA(INDEX(P)).GT.DATAP) GOTO 920
               ENDIF
            INDEX(P+1) = INDEXP
         ENDIF
  950    CONTINUE
 
C===================================================================
C
C     All done
 
      END
*deck djac
      SUBROUTINE DJAC( N, LDA, H, EVEC, EV, W )
*
*     Ivan Slapnicar, Faculty of Electrical Engineering, Mechanical
*     Engineering and Naval Architecture, R. Boskovica b.b,
*     58000 Split, Croatia, e-mail islapnicar at uni-zg.ac.mail.yu
*     This version:
*     June 30, 1992
*
*     .. Scalar Arguments
      INTEGER              N, LDA
*     ..
*     .. Array Arguments
      DOUBLE PRECISION     H( LDA, * ), EVEC( LDA, * ), EV( * ), W( * )
*     ..
*
*  Purpose
*  =======
*  
*  DJAC  computes all eigenvalues and and eigenvectors of a real
*  symmetric  N x N  matrix  H  using the standard Jacobi method
*  with delayed updates of the diagonal. For more information see 
*  
*     H. Rutishauser:
*     The Jacobi Method for Real Symmetric Matrices, 
*     Numer. Math. 9, 1-10 (1966)
*  
*  Only upper triangle of  H  is referenced.
*
*  Treshold is relative. Let  
*
*             (  A    CC  )  
*             (  CC    B  )
*
*  be the pivot submatrix. We perform the rotation if 
*
*    SQRT( ABS( CC ) ) > EPS * SQRT( ( ABS( A ) ) * SQRT( ABS( B ) ) 
*
*  where  EPS  is the machine precision. The algorithm has converged
*  if  N * ( N - 1 ) / 2  succesive rotations were not performed. 
*
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          Dimension of the symmetric matrix H. N >= 0.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array H.  LDA >= max(1,N)
*
*  H       (input) DOUBLE PRECISION array, dimension (LDA,N).
*          N x N  real symmetric matrix.
*
*  EVEC    (output) DOUBLE PRECISION array, dimension (LDA,N).
*          Contains the eigenvectors.
*
*  EV      (output) DOUBLE PRECISION array, dimansion (N).
*          contains eigenvalues of  H.
*
*  W       (workspace) DOUBLE PRECISION array, dimansion (N)
* 
*
*  Further Details
*  ===============
*
*  DJAC  uses BLAS1 subroutines  DROT, DCOPY and DAXPY, 
*             BLASJ subroutine  DROTJJ, and 
*             LAPACK auxiliary function DLAMCH.    
*
*  BLASJ contains BLAS-type routines for use by J-orthogonal Jacobi
*  methods. DLAMCH determines the machine precision.
*
*  =============================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, ISWEEP, NTEMP, NMAX
      DOUBLE PRECISION   EPS, A, B, CC, C, S, T 
*     ..
*     .. External Subroutines ..
      EXTERNAL           DROTJJ, DROT, DCOPY, DAXPY
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DSQRT, DABS
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( N .LT. 0 ) THEN
         INFO = 1
      ELSE IF( LDA .LT. MAX( 1,N ) ) THEN
         INFO = 2
      END IF
      IF( INFO .NE. 0 ) THEN
         WRITE(*,*) ' On entry to  DJAC  parameter number ', INFO 
         WRITE(*,*) ' had an illegal value. '
         STOP
      END IF
*
*     Quick return, if possible.
*
      IF( N .EQ. 0 ) RETURN
      IF( N .EQ. 1 ) THEN
         EVEC( 1, 1 ) = ONE
         EV( 1 ) = H( 1, 1 )
         RETURN
      END IF
*
*     Calculate the machine precision EPS, set  NMAX to size of the
*     sweep, and set  NTEMP  to 0.
*
      EPS = DLAMCH('E')
      NMAX = N * ( N - 1 ) / 2
      NTEMP = 0
*
*     Set the array  EVEC  to an identity matrix.
*
      DO 30 I = 1, N
         EVEC( I, I ) = ONE
      DO 30 J = I + 1, N
         EVEC( I, J ) = ZERO
         EVEC( J, I ) = ZERO
30    CONTINUE
*
*     Set EV to diagonal of  H.
*
      CALL DCOPY( N, H( 1,1 ), LDA + 1, EV, 1 )
*
*     Start the iterations.
*
      DO 40 ISWEEP = 1, 30
*
*     Set  W  to zero.
*
      DO 45 I = 1, N
         W( I ) = ZERO
45    CONTINUE 
*
*     Begin the main loop.
*
      DO 50 I=1,N-1
      DO 50 J=I+1,N
*
*     Calculate the pivot submatrix.
*
      A = H( I, I )
      B = H( J, J )
      CC = H( I, J )
*
*     Test whether to perform this rotation. If the rotation is
*     performed, set  NTEMP  to 0, otherwise  NTEMP = NTEMP + 1.
*     IF  NTEMP = NMAX, the algorithm has converged.
*
      IF( DABS( CC ) .LE. ( EPS * DSQRT( DABS( A ) ) * 
     $ DSQRT( DABS( B ) ) ) )  THEN
         NTEMP = NTEMP + 1 
         IF( NTEMP .GT. NMAX ) GO TO 100
         GO TO 50
      END IF
      NTEMP = 0
*
*     Calculate the Jacobi rotation parameters.
*
      CALL DROTJJ( A, B, CC, C, S, T, - ONE )
*
*     Update the upper triangle of  H  except the elements of the
*     pivot submatrix.
*
      CALL DROT( I-1, H( 1,I ), 1, H( 1,J ), 1, C, - S )
      CALL DROT( J-I-1, H( I, I + 1 ), LDA, H( I + 1, J ), 1, C, - S ) 
      CALL DROT( N-J, H( I, J + 1 ), LDA, H( J, J + 1 ), LDA, C, - S )
*
*     Update the elements of the pivot submatrix, and add the updates
*     of the diagonal to elements  W( I )  and  W( J ).
*
      CC = CC * T
      H( I, I ) = A - CC
      H( J, J ) = B + CC
      H( I, J ) = ZERO
      W( I ) = W( I ) - CC
      W( J ) = W( J ) + CC
*
*     Update the eigenvalue matrix.
*
      CALL DROT( N, EVEC( 1,I ), 1, EVEC( 1,J ), 1, C, - S )
*
*     End of the main loop.
*
50    CONTINUE
*
*     Add the updates of the diagonal stored in  W  to the diagonal
*     which was stored at the end of the previous sweep EV.
*
      CALL DAXPY( N, ONE, W, 1, EV, 1 )
*
*     Set diagonal of  H  to  EV, and set  W  to zero.
* 
      CALL DCOPY( N, EV, 1, H( 1,1 ), LDA + 1)
      DO 60 I = 1, N
         W( I ) = ZERO
60    CONTINUE
*
*     End of the  ISWEEP  loop.
*
40    CONTINUE
*
*     If this point is reached, the algorithm did not converge.
*
      WRITE(*,*) ' DJAC  executed ', ISWEEP, ' sweeps without', 
     $   ' convergence. '
      RETURN
*
*     Add the updates stored in  W  to  EV. 
*
100   CONTINUE 
      CALL DAXPY( N, ONE, W, 1, EV, 1 )
      RETURN
*
*     End of  DJAC.
*
      END
*deck drotjj
      SUBROUTINE DROTJJ( A, B, CC, C, S, T, HYP)
*
*     Ivan Slapnicar, Faculty of Electrical Engineering, Mechanical
*     Engineering and Naval Architecture, R. Boskovica b.b,
*     58000 Split, Croatia, e-mail islapnicar at uni-zg.ac.mail.yu
*     This version:
*     June 21, 1992
*
*     .. Scalar Arguments
      DOUBLE PRECISION     A, B, CC, C, S, T, HYP
*
*  Purpose
*  =======
*  
*  DROTJJ  constructs trigonometric or hyperbolic Jacobi plane rotation
*  of the form
*
*               (      C         S  )
*           R = (                   )
*               (  HYP * S       C  )
*
*  such that the matrix   
*
*               (  A    CC  )
*          R' * (           ) * R
*               (  CC    B  )
*
*  is diagonal.  If  HYP = - 1, the rotation is trigonometric, that is,
* 
*                    C ** 2 + S ** 2 = 1 .
*
*  If  HYP = 1, the rotation is hyperbolic, that is
* 
*                    C ** 2 - S ** 2 = 1 .
*
*  Arguments
*  =========
*
*  A, B, CC  (input) DOUBLE PRECISION
*            define the  2 x 2  input matrix.
*
*  C, S      (output) DOUBLE PRECISION
*            C and S  define the plane rotation. 
*            If  C = - 1, DROTJJ  failed to construct the rotation.
*            This can only happen in the hyperbolic case. Maximal
*            C  that can be calculated is of order of the fourth root
*            of the machine precision. 
*
*  T         (output) DOUBLE PRECISION
*            T = S / C
*
*  HYP       (input) DOUBLE PRECISION
*            If  HYP = 1.0D0,  DROTJJ computes a hyperbolic rotation.
*            If  HYP = -1.0D0,  DROTJJ computes a trigonometric rotation.
*  
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION      ABSCC, TMP, TMPBA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DSQRT, DABS  
*     ..
*     .. Executable Statements ..
*
      ABSCC = DABS( CC )
      TMP = 100.0D0 * ABSCC
      TMPBA = B + HYP * A
*
*     If TMP is small compared to TMPBA, use the approximate formulae.
*
      IF( ( TMPBA + TMP ) .EQ. TMPBA ) THEN
        T = - HYP * CC / TMPBA
        C = ONE
        S = T
        RETURN
      END IF
*
*     Calculate the tangent. If the rotation is hyperbolic, we have
*     to check whether it can be calculated.
*
      TMP = - HYP * 0.5D0 * TMPBA / CC
      TMPBA = TMP * TMP - HYP
      IF( TMPBA .GT. ZERO ) THEN
*
*     Standard formula for tangent.
*
         T = ONE / ( DABS( TMP ) + DSQRT( TMPBA ) )
      ELSE
*
*     Hyperbolic rotation cannot be calculated. Set  C  to - 1
*     and return.
*   
         C = - ONE
         RETURN
      END IF
*
*     Set the signum of the tangent.
*
      IF( TMP .LT. ZERO ) T = - T
*
*     Calculate cosine and sine.
*
      TMP = DSQRT( ONE - HYP * T * T )
      C = ONE / TMP
      S = T / TMP
      RETURN
*
*     End of  DROTJJ.
*
      END
