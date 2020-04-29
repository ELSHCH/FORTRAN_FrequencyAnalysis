	subroutine simpson(seq1,nseq,sum)
      integer nseq
      real seq(nseq),sum
c      call mean(nseq,seq1,meanseq)
c      do 11, i=1,nseq
c     	seq2(i)=abs(seq1(i)-meanseq)
c   11	continue
      sum=0.
      i=1
      do while ((i+2).le.nseq)
      sum=sum+1/3.*(seq(i)+4.*seq(i+1)+seq(i+2))
      i=i+2
      enddo 
      if ((i+1).eq.nseq) then
      sum=sum+seq(nseq)+0.5*(seq(nseq-1)-seq(nseq))
      endif
c      sum=sum/tint
      end subroutine
