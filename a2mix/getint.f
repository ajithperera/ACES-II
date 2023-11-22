      Subroutine getint(onehao,ovrlp,twoint,trnfor,tmpmat,pmat,
     &                  unitar,buf,ibuf,norbs)
      Implicit Double Precision (a-h, o-z)
      Integer fw
      Parameter(luint=10)
      Dimension onehao(norbs,norbs),ovrlp(norbs,norbs),
     &  trnfor(norbs,norbs), buf(600),ibuf(600),pmat(norbs,norbs),
     &  twoint(norbs*norbs*norbs*norbs/8),tmpmat(norbs,norbs),
     &  unitar(norbs,norbs)

      Common /machsp/ iintln,ifltln,iintfp,ialone,ibitwd
      Common /files/ luout,moints
      nnm1o2(ix)=(ix*(ix-1))/2
      iexti(ix)=1 +(-1+aint(dsqrt(8.d0*ix+0.999d0)))/2
      iextj(ix)=ix-nnm1o2(iexti(ix))
      ilnbuf=600
C Open IIII
       open(luint,file='IIII   ',form='UNFORMATTED',
     &        access='SEQUENTIAL')
      rewind luint
C******* Begin Get Onehamil ***********
      nut=ilnbuf
      Call locate(luint,'ONEHAMIL')
      Do While (nut.eq.ilnbuf)
         read(luint) buf, ibuf, nut
         Do 10 int = 1, nut
            indi=iexti(ibuf(int))
            indj=iextj(ibuf(int))
            onehao(indi,indj)=buf(int)
            onehao(indj,indi)=buf(int)
 10      Continue
      End Do
C******* End Get onehao ***************

C******* Begin Get the Overlap **********
      nut = ilnbuf
      Call locate(luint,'OVERLAP ')
      Do While (nut.eq.ilnbuf)
         read(luint) buf, ibuf, nut
         Do 20 int = 1, nut
            indi=iexti(ibuf(int))
            indj=iextj(ibuf(int))
            ovrlp(indi,indj)=buf(int)
            ovrlp(indj,indi)=buf(int)
 20      Continue
      End Do
C ********* End Get the Overlap *************

C ********* Begin Get the 2E integrals ************
      If (iintfp.eq.2) Then
         fw=8
      End If
      If (iintfp.eq.1) Then
         fw=15
      End If
      nut = ilnbuf
      Call locate(luint,'TWOELSUP')
      Do While (nut.eq.ilnbuf)
         read(luint) buf, ibuf, nut
         Do 30 int = 1, nut 
            ind1=and(ishft(ibuf(int),-3*fw),2**fw-1)
            ind2=and(ishft(ibuf(int),-2*fw),2**fw-1)
            ind3=and(ishft(ibuf(int),-1*fw),2**fw-1)
            ind4=and(ibuf(int),2**fw-1)
            Call mysort(ind1,ind2,ind3,ind4,indtot)
            twoint(indtot)=buf(int)
 30      Continue
      End Do
C *********End Get the 2E integrals ************

C Close IIII
      close(luint,status='keep')

C**********Begin Get the Transformation *************
      Do 40 itmp=1,norbs
         Do 40 jtmp=1,norbs
            tmpmat(itmp,jtmp)=ovrlp(itmp,jtmp) 
 40      Continue
C tmpmat=S
      Call eig(tmpmat,unitar,1,norbs,1)
C tmpmat=s unitar=U
      Do 50 itmp=1,norbs
            tmpmat(itmp,itmp)=1/sqrt(tmpmat(itmp,itmp)) 
 50      Continue
C tmpmat=s^-1/2
C      call zero(ovrlp,norbs*norbs*iintfp)
C X=U.s^-1/2.U^t
      Call xgemm('n','n',norbs,norbs,norbs,1.D0,unitar,norbs,tmpmat,
     &           norbs,0.D0,pmat,norbs)
      Call xgemm('n','t',norbs,norbs,norbs,1.D0,pmat,norbs,unitar,
     &           norbs,0.D0,trnfor,norbs)
      call zero(unitar,norbs*norbs*iintfp)
C********** End Get the Transformation *************


      Return
      End 


