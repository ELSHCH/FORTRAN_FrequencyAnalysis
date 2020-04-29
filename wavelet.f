	module wavelet
      double precision sigma,xi,pi
      parameter (N2=4096*4,nvoice=40,noctave=12,nscale=480,
     :sigma=10,xi=6.,pi=3.141592653589795,par=5)
      integer N2,nvoice,nscale,noctave  
      contains
      subroutine t_f(bigN,ttime,xx,fx,freqx,nscalx)
c	Time Frequency analysis of trajectories using wavelets
c	 bigN number of points, is a power of two 
c	 noctave=log2(N2)-2
c	 nscale=noctave*nscale
c	 Input
c      xin(num) input sequence
c	 tt(num) time sequence
c	 Output
c      txn time sequence of length nscale*N2
c	 fxn frequency sequence of length nscale*N2
c	 lenx the last non-zero element in txn and fxn
c	 dif diffusion coefficient
      integer bigN,ll,nf
c	ridgex(bigN,nscale)
      real sigma2,echex(nscale),
     :fx(bigN),ttime(bigN),
     :tend,step,dif,freqx(nscale),nscalx(nscale,bigN)
      complex xx(bigN), mmx(bigN,nscale)
c
c      do 97, i=1,bigN
c      ttime(i)=(i-1)*tt(num)/float(bigN)
c   97 continue 
c      call interpl(num,bigN,tt,ttime,xin,xx)
c	print*, xx(20)
      sigma2=sigma*sigma;
c     Computation of Wavelet Transform
      call awt(bigN,xx,mmx,echex,noctave,nvoice,sigma2)
      do 70, i=1,bigN
      do 71, j=1,nscale
      nscalx(j,i)=mmx(i,j)*conjg(mmx(i,j))
   71 continue
   70 continue
      do 72,i=1,nscale
      freqx(i)=bigN*xi/((ttime(bigN)-ttime(1))*echex(i))
   72 continue
c	Choosing maximum frequency among set of all frequencies
c      do 64, i=1,bigN*nscale
c      tx(i)=0.
c	txn(i)=0.
c      fx(i)=0.
c	fxn(i)=0.
c   64 continue
      call ridge(bigN,nscale,nscalx,ttime,freqx,fx,5,1.)
	 
	print*,'ku-ku'
c      tend=ttime(bigN-600)
c      step=ttime(2)-ttime(1)
c	Diffusion of frequencies 
c      call simpson(tend,fx1,step,bigN-600,dif)
  112	format(1x,e13.7/)
      end subroutine
c
      subroutine awt(n,x,mx,echelle,noctave,nvoice,sigma2)
      integer nvoice,nscale,noctave,n,kscale			 
      real omega(n),sigma2,
     :scale,s,freq(n),Psi(n),echelle(noctave*nvoice)
      complex x(n),fftx(n),iffty(n),mx(n,noctave*nvoice)  
c     AWT -- Analytical Wavelet Transform with Gabor window
c     Inputs
c     x	     signal
c     nvoice  number of voices
c            default = 12
c     noctave number of octave
c            default = log2(n)-2
c     sigma2    first parameter, for Gaussian window, it could be variance,
c 	     default = 1
c     xi    shifting frequency, for Gaussian
c            default = 5
c     Outputs
c     mx     Analytical wavelet transform of x

c     Mallat, "A Wavelet Tour of Signal Processing";
c             4.3.3 Analytical Wavelet Transform
c
c     Fourier Transform of signal
c      open(unit=7,file='fourtrans.dat')
      do 50, i=1,n
      fftx(i)=x(i)	
   50	continue
      call four1(fftx,n,1)
      do 14, i=1,n/2-1
      omega(i)=(i-1)*2*pi/float(n)
      omega(n/2+i+1)=(-n/2+i)*2*pi/float(n)
   14	continue
      omega(n/2)=(n/2-1)*2*pi/float(n)
      omega(n/2+1)=pi
      nscale  = nvoice*noctave
      kscale  = 1
      scale   = 2
      do 15, i = 1,noctave
      do 16, j = 1,nvoice
      s=scale*2**(j/float(nvoice))
      do 17, k=1,n
	freq(k)=s*omega(k)-xi
c     This definition of Psi is what is given originally in Wavelab
c	Psi = realpow(4.*pi.*sigma2,1/4)*exp(-sigma2/2*freq.*freq);
c      This definition of Psi gives what is important in the computation
      Psi(k)=exp(-sigma2/2.*freq(k)**2)
      iffty(k)=fftx(k)*Psi(k)
   17	continue
c     This definition of Psi is the correct one
c     Psi = sqrt(s)*realpow(4.*pi.*sigma2,1/4)*exp(-sigma2/2*freq.*freq);
      do 18, k=n/2+2,n
      Psi(k)=0.
      iffty(k)=0.
   18	continue
      call four1(iffty,n,-1) 
      do 19, k=1,n 
      iffty(k)=iffty(k)/float(n)
      mx(k,kscale)=iffty(k)
   19	continue
      echelle(kscale)=s
      kscale=kscale+1
   16	continue
      scale=scale*2
   15 continue
c      do 56, i=1,nscale
c	write(7,113) real(mx(1400,i))
c  56	continue
c 113  format(1x,2e13.7/)
c      close(unit=7)
c     the matrix mx is ordered from high to low frequencies
      end subroutine
c
      subroutine four1(ft,nn,isign)
      integer isign,nn,i,istep,j,m,mmax,nix
      real dat(2*nn),tempi,tempr
      complex ft(nn)
c	 Replaces dat(1:2*nn) by its discrete Fourier transform,
c      if isign is input as 1; or replaces dat(1:2*nn) by nn times
c	 its inverse discrete Fourier transform, if isign is input as
c	 -1. dat is a complex array of length nn, equivalently, a real
c	 array of length 2*nn. nn MUST be an integer power of 2.
      double precision theta,wi,wpi,wpr,wr,wtemp
      do 55, i=1,nn
      dat(2*i)=aimag(ft(i))
      dat(2*i-1)=real(ft(i))
   55	continue
      nix=2*nn
      j=1
      do 11, i=1,nix,2
      if (j.gt.i) then
      tempr=dat(j)
      tempi=dat(j+1)
      dat(j)=dat(i)
      dat(j+1)=dat(i+1)
      dat(i)=tempr
      dat(i+1)=tempi
      endif
      m=nn
    1	if ((m.ge.2).and.(j.gt.m)) then
      j=j-m
      m=m/2
      goto 1
      endif
      j=j+m
   11	continue 
      mmax=2
    2	if (nix.gt.mmax) then
      istep=2*mmax
      theta=6.28318530717959d0/(isign*mmax)
      wpr=-2.d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      wr=1.d0
      wi=0.d0
      do 13, m=1,mmax,2
      do 12, i=m,nix,istep
      j=i+mmax
      tempr=wr*dat(j)-wi*dat(j+1)
      tempi=wr*dat(j+1)+wi*dat(j)
      dat(j)=dat(i)-tempr
      dat(j+1)=dat(i+1)-tempi
      dat(i)=dat(i)+tempr
      dat(i+1)=dat(i+1)+tempi
   12	continue
      wtemp=wr
      wr=wr*wpr-wi*wpi+wr
      wi=wi*wpr+wtemp*wpi+wi
   13 continue
      mmax=istep
      goto 2
      endif
      do 56, i=1,nn
      ft(i)=cmplx(dat(2*i-1),dat(2*i))
   56	continue
      end subroutine
c

c
      subroutine mean(nn,xinit,xm)
      integer nn
      real xinit(nn),xm,sum
      sum=0.
      do 13, i=1,nn
      sum=sum+xinit(i)
   13 continue
      xm=sum/float(nn)
      end subroutine 
c   
      subroutine interpl(k,m,t,time,y,ynew)
      integer k,m,ll,j,i
      real time(m),t(k),ynew_r(m),ynew_i(m),
    	:y_r(k),y_i(k)
      complex ynew(m),y(k)
c	 y input array of length k
c	 t initial nodes sequence of length k
c	 ynew output array of length m
c	 time interpolation nodes sequence of length m
      time(1)=t(1)
      time(m)=t(k)
      h=(time(m)-time(1))/float(m)
      do 56, j=2,m
      time(j)=time(j-1)+h
   56	continue
c	 open(unit=7,file='datainit.dat')
c	do 34, i=1,n
c34	 write(7,16) t(i),y(i)
c	 close(unit=7)
      do 89, i=1,k
      y_r(i)=real(y(i))
      y_i(i)=aimag(y(i))
   89	continue
      ynew_r(1)=y_r(1)
      ynew_i(1)=y_i(1)
      i=2
      j=1
  125	do while (i.le.m)
      do while (j.le.k-1)
      if (time(i).eq.t(j)) then
      ynew_r(i)=y_r(j)
      ynew_i(i)=y_i(i)
      ll=1
      elseif (time(i).eq.t(j+1)) then
      ynew_r(i)=y_r(j+1)
      ynew_i(i)=y_i(j+1)
      ll=1
      elseif (time(i).gt.t(j) .and. time(i).lt.t(j+1)) then
c	   if (y(j+1).gt.y(j)) then
      ynew_r(i)=y_r(j+1)+(y_r(j)-y_r(j+1))/(t(j+1)-t(j))
     :*(t(j+1)-time(i))
      ynew_i(i)=y_i(j+1)+(y_i(j)-y_i(j+1))/(t(j+1)-t(j))
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
      do 78, i=1,m
      ynew(i)=cmplx(ynew_r(i),ynew_i(i))
   78	continue
c	 open(unit=4,file='datainter.dat')
c	 do 33, i=1,m
c33	 write(4,16) time(i),ynew(i)
c	 close(unit=4)
c16	 format(1x,e13.7,3x,e13.7/)
c       print*,"finish"
      end subroutine

	subroutine rearrange(xseq,tseq,nseq)
      integer nseq,index
      real xseq(nseq),tseq(nseq),tmin
      l=0
      do while (l.le.(nseq-2))
      tmin=tseq(l+1)
      xval=xseq(l+1)
      index=l+1
      do 17, i=l+1,nseq-1
      if (tseq(i+1).lt.tmin) then
      tmin=tseq(i+1)
      xval=xseq(i+1)
      index=i+1
      endif
   17 continue
      l=l+1
c	Rearranging
      do 18, i=index-1,l,-1
      tseq(i+1)=tseq(i)
      xseq(i+1)=xseq(i)
   18 continue
      tseq(l)=tmin 
      xseq(l)=xval
      enddo
      end subroutine
c
      subroutine ridge(ncols,nrows,mmm,tvrem,freqx,yr,par,tau)
      integer tplus(nrows),tint(nrows),
c	localmaxima(ncols,nrows),
	:tminus(nrows),par,nscale,ind(nrows),indmax,lmax
      real mmm(nrows,ncols),tau,thresh,xsub(nrows),tvrem(ncols),
     :yr(ncols),freqx(nrows),xmax,x1(nrows),kmax(nrows),
     :kmax1(nrows),x2(nrows)	 
c     ridge -- Ridges of an Analytic Wavelet Transform
c     Inputs
c     m input matrix
c     par    	  parameter, 2*par is how many neighbours to
c     compare
c     tau       threshold of determination of the ridges  (default 0.5)
c     Outputs
c     localmaxima  local maxima of every column, binary matrix same size
c		  as m.
c	open(unit=8,file='lookupmatrix.dat')
c	open(unit=7,file='fourtrans2.dat')
      nscale=ncols*nrows
c      do 14, j=1,ncols
c      do 13, i=1,nrows
c      localmaxima(j,i)=0
c   13	continue
c   14 continue
      do 15, i=1,nrows
      tint(i)=i
   15	continue
      tplus(1)=tint(nrows)
      tminus(nrows)=tint(1)
      do 60, i=1,nrows-1
      tplus(i+1)=tint(i)
      tminus(i)=tint(i+1)
   60 continue 
      l1=1
      do 16, i=1,ncols
      do 20, k=1,nrows
      xsub(k)=mmm(k,i)
   20 continue 
	do 18, j=1,par
	do 44, l=1,nrows
	xsub(l) = max(xsub(tint(l)),xsub(tplus(l)),xsub(tminus(l)))
  44	continue
  18	continue
      thresh=xsub(1)
      do 65, l=1,nrows-1
      if (xsub(l+1).ge.thresh) then
      thresh=xsub(l+1)
      endif
   65 continue
      do 17, k=1,nrows
      if ((mmm(k,i).ge.thresh).and.(mmm(k,i).ge.xsub(k))) then 
    	yr(l1)=freqx(k)
	l1=l1+1
      endif
   17	continue  
   16 continue
c      do 19, i=1,nrows
c	write(8,114) localmaxima(1400,i)
c	write(7,113) mmm(i,1400)
c  19  continue
c 114	format(1x,i3/)
c 113  format(1x,e13.7/)
c      close(unit=8) 
c	close(unit=7) 
      end subroutine
c
      subroutine simpson(tint,seq1,h,nseq,sum)
      integer nseq
      real tint,seq1(nseq),meanseq,sum,h,seq2(nseq)
      call mean(nseq,seq1,meanseq)
      do 11, i=1,nseq
     	seq2(i)=abs(seq1(i)-meanseq)
   11	continue
      sum=0.
      i=1
      do while ((i+2).le.nseq)
      sum=sum+1/3.*h*(seq2(i)+4.*seq2(i+1)+seq2(i+2))
      i=i+2
      enddo 
      if ((i+1).eq.nseq) then
      sum=sum+h*seq2(nseq)+0.5*h*(seq2(nseq-1)-seq2(nseq))
      endif
      sum=sum/tint
      end subroutine
c     
      subroutine imageridges(m,ncols,nrows,t,freq,xr,yr,kk)
c     Subroutine locates frequency and time sequence
c	correponding to selected ridges
c     Inputs
c     m input matrix of size (nrows x ncols)
c     Outputs
c     (xr,yr)  the time-frequency plane
      integer ncols,nrows,indexx(nrows*ncols),indexy(ncols*nrows),
    	:m(nrows,ncols),kk,i,j
      real t(nrows),freq(ncols),xr(nrows),
	:yr(nrows)
c
      kk=1
      do 87, j=1,ncols
      do 88, i=1,nrows
      if (m(i,j).eq.1) then
      indexx(kk)=i
      indexy(kk)=j
      kk=kk+1
      endif
   88 continue
   87 continue
      kk=kk-1
      i=2
      xr(1)=t(indexx(1))
      yr(1)=freq(indexy(1))
      do while (i.le.nrows)
c	if (indexx(i).le.indexx(i-1)) then
c	i=i+1
c	else
      xr(i)=t(indexx(i))
      yr(i)=freq(indexy(i))
      i=i+1
c	endif
      end do
c	kk=i-1
      end	subroutine 
      end module
