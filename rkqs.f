	SUBROUTINE rkqs(y,dydx,nl,xx,htry,eps,yscal,hnext,derivs)
      INTEGER nl,NMAX
      REAL*8 eps,hnext,htry,xx,dydx(nl),y(nl),df
      EXTERNAL derivs,rkck
      PARAMETER (NMAX=50)
c    USES derivs,rkck
      INTEGER i
      REAL*8 errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,
     *PSHRNK,ERRCON, yscal(NMAX)
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)
      h=htry
1     call rkck(y,dydx,nl,xx,h,ytemp,yerr,derivs)
c	if (df.gt.1e+15) return
      errmax=0.
      do 11 i=1,nl
c        errmax=max(errmax,abs(yerr(i)/yscal(i)))
         errmax=max(errmax,abs(yerr(i)))  
11    continue
      errmax=errmax/eps
      if(errmax.gt.1.)then
c        htemp=SAFETY*h*(errmax**PSHRNK)
         htemp=SAFETY*h*(errmax**PGROW)
c        h=sign(max(abs(htemp),0.1*abs(h)),h)
         h=sign(abs(htemp),h)
        xnew=xx+h
        if(xnew.eq.x)pause 'stepsize underflow in rkqs'
        goto 1
      else
c        if(errmax.gt.ERRCON)then
c          hnext=SAFETY*h*(errmax**PGROW)
c        else
	     hnext=SAFETY*h*(errmax**PGROW)
c          hnext=5.*h
c        endif
c        hdid=h
        xx=xx+h
        do 12 i=1,nl
          y(i)=ytemp(i)
12      continue
        return
      endif
      END
