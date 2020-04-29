	subroutine derivs(t,yi,yprime)
      double precision B,F,E,alpha,yprime(6),
     :yi(6),t,c1,c2,c3,c4,c5,c6,u1,u2,u3,pi,g
	parameter (pi=3.141592653589795)
c  
      external impampl
c
      call impampl(F,B,E,alpha)
	c1=yi(1)**2+yi(2)**2   
	c2=yi(1)*yi(5)-yi(2)*yi(4)
	c3=yi(1)**2-yi(2)**2 
	c4=yi(1)*yi(2)
	g=dmod(yi(3),2*pi) 
	c5=dcos(yi(3))
	c6=dsin(yi(3)) 
	u1=1+0.5*B
	u2=0.5*u1*c1
	u3=(alpha-1)*F
      yprime(1)=yi(4)+u2*yi(2)
      yprime(2)=yi(5)-u2*yi(1)
      yprime(3)=c1
      yprime(4)=-(2*yi(6)*yi(1)-u1*yi(1)*
     :c2-u2*yi(5)+2*F*yi(1)**3+6./32.*B**2*yi(1)*c1**2
     :+u3*c6*(yi(1)*(c3*c6+2*c4*c5)+c1*(yi(1)*c6+yi(2)*c5)))
	yprime(5)=-(2*yi(2)*yi(6)-u1*yi(2)*c2
     :+u2*yi(4)-2.*F*yi(2)**3+
     :6./32.*B**2*yi(2)*c1**2
     :+u3*c6*(yi(2)*(c3*c6+2*c4*c5)+c1*(-yi(2)*c6+yi(1)*c5)))
      yprime(6)=-u3*c1*(0.5*c3*dsin(2*yi(3))+c4*dcos(2*yi(3)))   
      end subroutine	
	  
	
