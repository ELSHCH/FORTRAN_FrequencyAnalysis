      subroutine odeint(ystart,nvar,x1,x2,eps,h1,hmin,kount,xp,yp,
	*derivs,rkqs)
c	call odeint(ystart,n,t0,tf,eps,h1,hmin,nk,ttout,yout,derivs,rkqs)

      INTEGER nvar,KMAXX,MAXSTP,NMAX
      REAL*8 eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL rkqs,derivs
      PARAMETER (MAXSTP=1000000,NMAX=50,KMAXX=200,TINY=1.e-30)
      INTEGER i,kount,nstp
      real*8 h,hnext,xx,dydx(NMAX),xp(MAXSTP),y(nvar),
     *yp(MAXSTP,nvar),yscal(nvar),df
      xx=x1
c	hmin=0.001*h1
      h=sign(h1,x2-x1)
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
c      if (kmax.gt.0) xsav=x-2.*dxsav	
      do 16 nstp=1,MAXSTP
        call derivs(xx,y,dydx)
        do 12 i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
12      continue
c        if(kmax.gt.0)then
c          if(abs(x-xsav).gt.abs(dxsav)) then
c            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=xx
              do 13 i=1,nvar
                yp(kount,i)=y(i)
13            continue
c              xsav=x
c            endif
c          endif
c        endif
        if((xx+h-x2)*(xx+h-x1).gt.0.) h=x2-xx
        call rkqs(y,dydx,nvar,xx,h,eps,yscal,hnext,derivs)
c	  if (df.gt.1e+15) return 
c        if(hdid.eq.h)then
c          nok=nok+1
c        else
c          nbad=nbad+1
c        endif
        if((xx-x2)*(x2-x1).ge.0.)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
c          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=xx
            do 15 i=1,nvar
              yp(kount,i)=y(i)
15          continue
c          endif
          return
        endif
        if(abs(hnext).lt.hmin) then	  
        print*,'stepsize smaller than minimum in odeint'
        return
	  endif
        h=hnext
16    continue
c      pause 'too many steps in odeint'
      return
      END
