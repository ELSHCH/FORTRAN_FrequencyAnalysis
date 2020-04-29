	subroutine interp(k,m,t,time,y,ynew)
      integer k,m,ll,j,i
      real*8 time(m),t(k),ynew(m),y(k)
c	 y input array of length k
c	 t initial nodes sequence of length k
c	 ynew output array of length m
c	 time interpolation nodes sequence of length m
      time(1)=t(1)
      time(m)=t(k)
      h=(time(m)-time(1))/float(m)
      do 56, j=2,m
      time(j)=time(j-1)+h
   56 continue
c	 open(unit=7,file='datainit.dat')
c	do 34, i=1,n
c34	 write(7,16) t(i),y(i)
c	 close(unit=7)
      ynew(1)=y(1)
      i=2
      j=1
  125 do while (i.le.m)
      do while (j.le.k-1)
      if (time(i).eq.t(j)) then
      ynew(i)=y(j)
      ll=1
      elseif (time(i).eq.t(j+1)) then
      ynew(i)=y(j+1)
      ll=1
      elseif (time(i).gt.t(j) .and. time(i).lt.t(j+1)) then
c	   if (y(j+1).gt.y(j)) then
      ynew(i)=y(j+1)+(y(j)-y(j+1))/(t(j+1)-t(j))
     :*(t(j+1)-time(i))
      ll=1
c	   else
c	    ynew(i)=y(i)+(y(j+1)-y(j))/(t(j+1)-t(j))*(time(i)-t(j))
c	   ll=1
c	   endif
      else
      ll=0
      endif    
      if (ll.eq.1) then
      i=i+1
      goto 125
      else
      j=j+1
      endif
      end do
      end do
c	 open(unit=4,file='datainter.dat')
c	 do 33, i=1,m
c33	 write(4,16) time(i),ynew(i)
c	 close(unit=4)
c16	 format(1x,e13.7,3x,e13.7/)
c       print*,"finish"
      end subroutine
c
 
