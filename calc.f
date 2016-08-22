      implicit real*8(a-h,o-z)
      character cline*80 , op*1
      dimension is(40), ie(40), st(100)
c
c     hp style calculator for quantum chemists
c     johnm 7/09
c
      pi=3.141592653589793D+00
c     copied from molpro manual
      toev=27.2113961d0
      tok=3.157733d5
      tocm=219474.63067d0
      tohz=6.5796838999d15
      tokj=2625.500d0
      tokcal=627.5096d0
      toa=0.529177249d0
      tod=2.54158d0
c
      ik = 1
      call help
    5 read(5,'(a)') cline
c     write(6,'(a)') cline
      call scan(cline,nw,is,ie)
      do 20 n = 1, nw
        if(is(n).eq.ie(n).and.ichar(cline(is(n):ie(n))).lt.ichar('0'))
     $    then
          read(cline(is(n):ie(n)),'(a)') op
          call ckstak(ik-2)
          if(op.eq.'+') then
            r = st(ik-2) + st(ik-1)
          else if(op.eq.'-') then
            r = st(ik-2) - st(ik-1)
          else if(op.eq.'*') then
            r = st(ik-2) * st(ik-1)
          else if(op.eq.'/') then
            r = st(ik-2) / st(ik-1)
          endif
          st(ik-2) = r
          ik = ik - 1
          call ckstak(ik)
          write(6,100) r
        else if(ichar(cline(is(n):ie(n))).gt.ichar('@')) then
          if(cline(is(n):ie(n)).eq.'q') stop
          if(cline(is(n):ie(n)).eq.'s') then
            call ckstak(ik-2)
            tmp = st(ik-2)
            st(ik-2) = st(ik-1)
            st(ik-1) = tmp
            goto 20
            endif
          if(cline(is(n):ie(n)).eq.'c') then
            ik = 1
            goto 20
            endif
          if(cline(is(n):ie(n)).eq.'p') then
            call dump(ik,st)
            goto 20
            endif
          if(cline(is(n):ie(n)).eq.'h') then
            call help
            goto 20
            endif
          if(cline(is(n):ie(n)).eq.'pop') then
            call ckstak(ik-1)
            ik = ik - 1
            goto 20
            endif
          if(cline(is(n):ie(n)).eq.'exp') then
            call ckstak(ik-1)
            r = exp(st(ik-1))
          else if(cline(is(n):ie(n)).eq.'ln') then
            call ckstak(ik-1)
            r = log(st(ik-1))
          else if(cline(is(n):ie(n)).eq.'log') then
            call ckstak(ik-1)
            r = log10(st(ik-1))
          else if(cline(is(n):ie(n)).eq.'sqrt') then
            call ckstak(ik-1)
            r = sqrt(st(ik-1))
          else if(cline(is(n):ie(n)).eq.'sin') then
            call ckstak(ik-1)
            r = sin(st(ik-1))
          else if(cline(is(n):ie(n)).eq.'cos') then
            call ckstak(ik-1)
            r = cos(st(ik-1))
          else if(cline(is(n):ie(n)).eq.'tan') then
            call ckstak(ik-1)
            r = tan(st(ik-1))
          else if(cline(is(n):ie(n)).eq.'asin') then
            call ckstak(ik-1)
            r = asin(st(ik-1))
          else if(cline(is(n):ie(n)).eq.'acos') then
            call ckstak(ik-1)
            r = acos(st(ik-1))
          else if(cline(is(n):ie(n)).eq.'atan') then
            call ckstak(ik-1)
            r = atan(st(ik-1))
          else if(cline(is(n):ie(n)).eq.'inv') then
            call ckstak(ik-1)
            r = 1.d0/st(ik-1)
          else if(cline(is(n):ie(n)).eq.'sq') then
            call ckstak(ik-1)
            r = st(ik-1)**2
          else if(cline(is(n):ie(n)).eq.'xtoy') then
            call ckstak(ik-2)
            r = st(ik-2)**st(ik-1)
          else
            goto 10
          endif
          st(ik-1) = r
          write(6,100) st(ik-1)
          goto 20
   10     if(cline(is(n):ie(n)).eq.'pi') then
            st(ik) = pi
          else if(cline(is(n):ie(n)).eq.'toev') then
            st(ik) = toev
          else if(cline(is(n):ie(n)).eq.'toa') then
            st(ik) = toa
          else if(cline(is(n):ie(n)).eq.'tod') then
            st(ik) = tod
          else if(cline(is(n):ie(n)).eq.'tok') then
            st(ik) = tok
          else if(cline(is(n):ie(n)).eq.'tokj') then
            st(ik) = tokj
          else if(cline(is(n):ie(n)).eq.'tocm') then
            st(ik) = tocm
          else if(cline(is(n):ie(n)).eq.'tohz') then
            st(ik) = tohz
          else if(cline(is(n):ie(n)).eq.'tokcal') then
            st(ik) = tokcal
          else
            write(6,*) '?'
            goto 20
          endif
          write(6,100) st(ik)
          ik = ik + 1
          call ckstak(ik)
        else
          read(cline(is(n):ie(n)),*) st(ik)
          ik = ik + 1
          call ckstak(ik)
        endif
   20   continue
      goto 5
  100 format(g18.9)
      end
*deck dump
      subroutine dump(ik,s)
      implicit real*8(a-h,o-z)
      dimension s(*)
c     write(6,*) 'ik=',ik
      do i = 1, ik-1
        write(6,100) i, s(i)
      enddo
  100 format(i4,' : ',g18.9)
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
      end 
*deck ckstak
      subroutine ckstak(ik)
      if(ik.gt.100) then
        write(6,*) 'stack overflow'
        stop
      endif
      if(ik.lt.1) then
        write(6,*) 'stack underflow'
        stop
      endif
      end 
*deck help
      subroutine help
      write(6,100)
  100 format('HP style calculator for quantum chemists'//
     $ 'Binary Operators:'/'+, -, *, /, s (swap)'//
     $ 'Unary Operators:'/
     $ 'exp, sqrt, ln, log, sin, cos, tan, inv, sq, pop'/
     $ 'p (print stack), c (clear stack), h (help), q (quit)'//
     $ 'Constants:'/'pi=3.141592653589793  toev=27.2113961'/
     $ 'tok=3.157733d5        tocm=219474.63067'/
     $ 'tohz=6.5796838999d15  tokj=2625.500'/
     $ 'tokcal=627.5096       toa=0.529177249'/
     $ 'tod=2.54158d0'/)
      return
      end
