	subroutine mean(nn,xinit,xm)
      integer nn
      real xinit(nn),xm,sum
      sum=0.
      do 13, i=1,nn
      sum=sum+xinit(i)
   13 continue
      xm=sum/float(nn)
      end subroutine 
