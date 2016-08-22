      program convert
c
c     convert Gaussian basis to CFOUR
c
      implicit real*8(a-h,o-z)
      parameter (MaxPrm=100,MaxBF=100,MaxL=6)
      character cline*80, name*2, lval*1, spd*1, basnm*16
      dimension is(40), ie(40), nex(MaxL), nbf(MaxL), spd(MaxL)
      dimension ex(MaxPrm,10), cc(100,100,6)
      save spd
      data spd/'S','P','D','F','G','H'/
c
c     lmax = maximum angular momentum in the basis set
c     nex(l) = no. of exponents of angular momentum l
c     nbf(l) = no. of contracted functions of angular momentum l
c
      lmax = 0
      do k=1,MaxL
        nex(k) = 0
        nbf(k) = 0
      enddo
      do i=1,MaxPrm
        do j=1,MaxBF
          do k=1,MaxL
            cc(i,j,k) = 0.d0
          enddo
        enddo
      enddo
c
      basnm = ' '
      if(iargc().eq.1) call getarg(1,basnm)
      read(5,'(a2)',end=99) name
    5 read(5,'(a)') cline
      call scan(cline,nw,is,ie)
      if(cline(is(1):ie(1)).eq.'****') goto 10
      read(cline(is(1):ie(1)),'(a)') lval
      read(cline(is(2):ie(2)),'(i3)') npr
      read(cline(is(3):ie(3)),'(f17.10)') scal
c     write(6,*) lval, npr, scal
      do l=1,MaxL
        if(lval.eq.spd(l)) lx = l
      enddo
      if(lx.gt.lmax) lmax = lx
      nbf(lx) = nbf(lx) + 1
        do i=1, npr
          read(5,*) z, c
c         write(6,*) z, c
          nx = 0
          if(nex(lx).gt.0) then
            nx = match(nex(lx),z,ex(1,lx))
          endif
          if(nx.gt.0) then
            cc(nx,nbf(lx),lx) = c
          else
            nex(lx) = nex(lx) + 1
            ex(nex(lx),lx) = z
            cc(nex(lx),nbf(lx),lx) = c
          endif
        enddo
      goto 5
   10 continue
c
c     write basis in ACES format
c
      write(6,*) name, ':', basnm
      write(6,*)
      write(6,'(i3)') lmax
      write(6,'(6i5)') (l-1,l=1,lmax)
      write(6,'(6i5)') (nbf(l),l=1,lmax)
      write(6,'(6i5)') (nex(l),l=1,lmax)
      do l=1, lmax
        write(6,*)
        write(6,'(6F14.7)') (ex(k,l),k=1,nex(l))
        write(6,*)
        do i=1, nex(l)
          write(6,'(6F11.7)') (cc(i,k,l),k=1,nbf(l))
        enddo
      enddo
   99 continue
  100 format(5f10.4)
  101 format(5f12.6)
  102 format(5f14.7)
      end
*deck match
      integer function match(n,x,c)
      implicit real*8(a-h,o-z)
c
c     Find match in array
c
      dimension c(*)
      do i = 1, n
        if(x.eq.c(i)) then
          match = i
c         write(6,*) 'match',i,x,c(i)
          return
        endif
      enddo
      match = 0
      return
      end
*deck scan
      subroutine scan(line,nw,is,ie)
c
c     find non-blank words in line
c
      character line*80
      dimension is(*), ie(*)
      nw = 0
      l = lenstr(line)
      if(l.eq.0) goto 10
      line(l+1:l+1) = ' '
c     write(6,*) 'len =', l
      js = iskpb(line,1)
    5   j = index(line(js:), ' ')
c       write(6,*) 'range =', js, js+j-2, line(js:js+j-2)
        nw = nw + 1
        is(nw) = js
        ie(nw) = js+j-2
        if(js+j.gt.l) goto 10
        js = iskpb(line,js+j)
        goto 5
  10  continue
      return
      end
*deck lenstr
      integer function lenstr(str)
c
c     Returns length of string ignoring trailing blanks
c
      character*(*) str
      do 15, i = len(str), 1, -1
         if(str(i:i) .ne. ' ') goto 20
15    continue
20    lenstr = i
      return
      end
*deck iskpb
      integer function iskpb(str,is)
c
c     Returns next non-blank char
c
      character*(*) str
      do 15, j = 1, len(str)
         i = is + j - 1
         if(str(i:i) .ne. ' ') goto 20
15    continue
      i = -1
20    iskpb = i
      return
      end
