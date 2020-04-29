      subroutine norm(xx,n,sum)
      integer n
      double precision xx(1:n),sum
c	
      sum=0.
      do 15, i=1,n
      sum=xx(i)**2+sum
   15 continue
      sum=dsqrt(sum) 
      end
