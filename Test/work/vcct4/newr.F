      subroutine rdov4(io,no,nu,ti,v)
      implicit double precision(a-h,o-z)
      character*8 lb1,lb2
      dimension ti(1),v(nu,nu,nu,nu)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      integer a,b,dirprd
c      write(6,*)'rdov4b entered:,io,no,nu:',io,no,nu
      if(io.eq.1)then
         if(nirrep.gt.1)then 
            call zeroma(v,1,no4)
            lb1='SOAOB2  '
            i1=1+no2
            i2=i1+no2
            itot=i2+nou
          call grdsym(nu,nu,nu,nu,no2,no2,v,ti,ti(i1),ti(i2),lb1,lb1,13)
         else
            call getlst(v,1,no2,1,1,13)
         endif
      endif
      if(io.eq.0)then
         ij=0
         do 10 i=1,nu
            do 10 j=1,nu
               ij=ij+1
               call getlst(ti,ij,1,1,1,233)
               iab=0
               do 8 a=1,nu
                  do 8 b=1,nu
                     iab=iab+1
                     v(a,b,i,j)=ti(iab)
 8             continue
 10      continue
         call mtrans(v,nu,7)
         nu4=nu*nu*nu*nu
 23      format(4f15.10)
      endif
      return
      end
      subroutine rrsymv(nfile,irec,na,a)
      implicit double precision(a-h,o-z)
      dimension a(na)
      read(nfile,rec=irec)a
      return
      end
      subroutine rdvem4(io,no,nu,ti,v)
      implicit double precision(a-h,o-z)
      integer dirprd
      character*8 lb1,lb2
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      dimension ti(1),v(1)
      nou=no*nu
      if(io.eq.0)then
         if(nirrep.gt.1)then
            call zeroma(v,1,nou3)
            lb1='SVAVB2  '
            lb2='SVBOA2  '
            i1=1+nu2
            i2=i1+nu2
            itot=i2+nou
          call grdsym(nu,nu,nu,no,nu2,nou,v,ti,ti(i1),ti(i2),lb1,lb2,30)
            call tranmd(v,nu,nu,nu,no,23)
         else
            call getlst(v,1,nou,1,1,30)
            call tranmd(v,nu,nu,nu,no,23)
         endif
      endif
      if (io.eq.1)then
         if(nirrep.gt.1)then
            call zeroma(v,1,no3u)
            lb1='SOAOB2  '
            lb2='SOAVB2  '
            i1=1+no2
            i2=i1+no2
            itot=i2+nou
          call grdsym(nu,nu,nu,no,no2,nou,v,ti,ti(i1),ti(i2),lb1,lb2,10)
            call tranmd(v,nu,nu,nu,no,23)
         else
            call getlst(v,1,nou,1,1,10)
            call tranmd(v,nu,nu,nu,no,23)
         endif
      endif
      return
      end
      subroutine rdvem3(io,no,nu,ti,v)
      implicit double precision(a-h,o-z)
      character*8 lb1,lb2
      integer dirprd
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      dimension ti(1),v(1)
      nou=no*nu
      if(io.eq.0)then
         if(nirrep.gt.1)then
            call zeroma(v,1,nou3)
            lb1='SVAVB2  '
            lb2='SVBOA2  '
            i1=1+nu2
            i2=i1+nu2
            itot=i2+nou
         call grdsym3(nu,nu,no,nu,nu2,nou,v,ti,ti(i1),ti(i2),lb1,lb2,30)
            call tranmd(v,nu,nu,no,nu,24)
         else
            call rdve3(no,nu,v)
            call tranmd(v,nu,nu,no,nu,24)
         endif
      endif
      if (io.eq.1)then
         if(nirrep.gt.1)then
            call zeroma(v,1,no3u)
            lb1='SOAOB2  '
            lb2='SOAVB2  '
            i1=1+no2
            i2=i1+no2
            itot=i2+nou
         call grdsym3(nu,nu,no,nu,no2,nou,v,ti,ti(i1),ti(i2),lb1,lb2,10)
            call tranmd(v,nu,nu,no,nu,24)
         else
            call rdvm3(nu,no,v)
            call tranmd(v,nu,nu,no,nu,24)
         endif
      endif
      return
      end
      subroutine ro2hpp(io,no,nu,ti,t2)
      implicit double precision(a-h,o-z)
      character*8 lb1,lb2
      integer dirprd
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      dimension ti(1),t2(nu,nu,no,no)
      if (io.eq.0) then
         if(nirrep.gt.1)then
            call zeroma(t2,1,no2u2)
            lb1='SVAVB2  '
            lb2='SOAOB2  '
            i1=1+nu2
            i2=i1+nu2
            itot=i2+no2
         call grdsym(nu,nu,no,no,nu2,no2,t2,ti,ti(i1),ti(i2),lb1,lb2,46)
            call insitu(nu,nu,no,no,ti,t2,13)
            call tranmd(t2,no,nu,nu,no,23)
         else
            CALL GETLST(t2,1,NO2,1,1,63)
            call insitu(nu,nu,no,no,ti,t2,13)
            call tranmd(t2,no,nu,nu,no,23)
         endif
      endif
cvhpe
      if (io.eq.1) then
         if(nirrep.gt.1)then
            call zeroma(t2,1,no2u2)
            lb1='SVAVB2  '
            lb2='SOAOB2  '
            i1=1+nu2
            i2=i1+nu2
            itot=i2+no2
        call grdsym(nu,nu,no,no,nu2,no2,t2,ti,ti(i1),ti(i2),lb1,lb2,16)
            call insitu(nu,nu,no,no,ti,t2,13)
            call tranmd(t2,no,nu,nu,no,23)
         else
            CALL GETLST(T2,1,no2,1,1,16)
            call insitu(nu,nu,no,no,ti,t2,13)
            call tranmd(t2,no,nu,nu,no,23)
         endif
      endif
cvhp
      if (io.eq.2) then
      call zeroma(t2,1,no2u2)
         if(nirrep.gt.1)then
            call zeroma(t2,1,no2u2)
            lb1='SVBOA2  '
            i1=1+nou
            i2=i1+nou
            itot=i2+nou
         call grdsym(nu,no,nu,no,nou,nou,t2,ti,ti(i1),ti(i2),lb1,lb1,25)
            call insitu(nu,no,nu,no,ti,t2,12)
         else
            CALL GETLST(T2,1,nou,1,1,25)
            call insitu(nu,no,nu,no,ti,t2,12)
         endif
      endif
      return
      end
      subroutine engy1(no,nu,ti,t2,vo,AA)
      implicit double precision (a-h,o-z)
      character*8 AA
      integer a,b
      common/en51/ec1
      common/mbpt5en/ec
      common/dgr22/e22
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension t2(no,nu,nu,no),vo(no,nu,nu,no),ti(1)
      data zero/0.0d+0/,two/2.0d+0/
      call ro2hpp(0,no,nu,ti,vo)
      x=zero
      do 10 i=1,no
         do 10 j=1,no
            do 10 a=1,nu
               do 10 b=1,nu
                  x=x+t2(i,a,b,j)*(two*vo(i,a,b,j)-vo(i,b,a,j))
 10   continue
      ec1=ec1+x
      if(aa.eq.'noprint ')then
      e22=e22+x
      else
      if(aa.eq.'t2ppt2sq')x=x+e22
      write(6,34)aa,x,ec1
      endif
 34   format(A8,' contribution =',f15.10,' Total=',f15.10)
      return
      end
      subroutine engy(no,nu,ti,t2,vo,aa)
      implicit double precision (a-h,o-z)
#include "machsp.com"
      integer a,b
      character*8 aa
      common/en51/ec1
      common/coreng/ecor
      common/mbpt5en/ec
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension t2(no,nu,nu,no),vo(no,nu,nu,no),ti(1)
      data zero/0.0d+0/,two/2.0d+0/
      call ro2hpp(0,no,nu,ti,vo)
      ec=zero
      do 10 i=1,no
      do 10 j=1,no
      do 10 a=1,nu
      do 10 b=1,nu
      ec=ec+t2(i,a,b,j)*(two*vo(i,a,b,j)-vo(i,b,a,j))
 10   continue
      write(6,34)aa,ec,(ec1+ec)
 34   format(A8,' contribution =',f15.10,' Total=',f15.10)
      ec=ec+ec1
 96   format(/,77('-'),/)
      write(6,96)
      write(6,99)ec
 99   format('Total nonit. T4 energy contribution:',f18.12)
 98   format('Total correlation energy:           ',f18.12)
 97   format('Total   energy :                    ',f18.12)
      call getrec(20,'JOBARC','SCFENEG ',iintfp,ESCF)
      call getrec(20,'JOBARC','TOTENERG',iintfp,ENTOT)
      ecor=ENTOT-ESCF
      write(6,98)ecor+ec
      ETOTPLUS5=ENTOT+EC
      write(6,97)etotplus5
      call putrec(20,'JOBARC','TOTENERG',iintfp,ETOTPLUS5)
      return
      end
      subroutine rdve3(no,nu,ve)
      implicit double precision (a-h,o-z)
      integer a
      dimension ve(nu,nu,no,nu)
      kk=0
      do 10 i=1,no
         do 10 a=1,nu
          kk=kk+1
          call getlst(ve(1,1,i,a),kk,1,1,1,30)
 10   continue
      return
      end
      subroutine rdvm3(no,nu,ve)
      implicit double precision (a-h,o-z)
      integer a
      dimension ve(no,no,nu,no)
      kk=0
      do 10 a=1,nu
      do 10 i=1,no
          kk=kk+1
          call getlst(ve(1,1,a,i),kk,1,1,1,10)
 10   continue
      return
      end
      SUBROUTINE INSI12(N1,N2,N3,TI,O2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TI(N1,N2,N3),O2(N2,N1,N3)
      N=N1*N2*N3
      CALL VECCOP(N,TI,O2)
      DO 10 I=1,N1
      DO 10 J=1,N2
      DO 10 K=1,N3
      O2(J,I,K)=TI(I,J,K)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE INSI13(N1,N2,N3,TI,O2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TI(N1,N2,N3),O2(N3,N2,N1)
      N=N1*N2*N3
      CALL VECCOP(N,TI,O2)
      DO 10 I=1,N1
      DO 10 J=1,N2
      DO 10 K=1,N3
      O2(K,J,I)=TI(I,J,K)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE INSI23(N1,N2,N3,TI,O2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TI(N1,N2,N3),O2(N1,N3,N2)
      N=N1*N2*N3
      CALL VECCOP(N,TI,O2)
      DO 10 I=1,N1
      DO 10 J=1,N2
      DO 10 K=1,N3
      O2(I,K,J)=TI(I,J,K)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE INSITU(N1,N2,N3,N4,TI,O2,IC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TI(N1,N2,N3),O2(N1,N2,N3,N4)
      DO 10 I=1,N4
      IF (IC.EQ.12)GOTO 12
      IF (IC.EQ.23)GOTO 23
      GOTO 13
 12   CALL INSI12(N1,N2,N3,TI,O2(1,1,1,I))
      GOTO 10
 23   CALL INSI23(N1,N2,N3,TI,O2(1,1,1,I))
      GOTO 10
 13   CALL INSI13(N1,N2,N3,TI,O2(1,1,1,I))
 10   CONTINUE
      RETURN
      END
      SUBROUTINE TRANMD(A,N1,N2,N3,N4,IJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1)
c      call ienter(43)
      n12=n1*n2
      n123=n12*n3
      IF(IJ.EQ.12) GO TO 12
      IF(IJ.EQ.13) GO TO 13
      IF(IJ.EQ.14) GO TO 14
      IF(IJ.EQ.23) GO TO 23
      IF(IJ.EQ.24) GO TO 24
      IF(IJ.EQ.34) GO TO 34
      IF(IJ.EQ.231) GO TO 231
      IF(IJ.EQ.312) GO TO 312
      IF(IJ.EQ.341) GO TO 341
      IF(IJ.EQ.413) GO TO 413
      IF(IJ.EQ.1234) GO TO 1234
      GOTO 100
   12 CONTINUE
      DO 10 L=1,N4     
         n123l=n123*(l-1)
      DO 10 K=1,N3      
         n12k=(k-1)*n12+n123l
      DO 10 I=1,N1      
         n1i=(i-1)*n1
      DO 10 J=1,I      
         ijkl=n12k+(j-1)*n1+i
         jikl=n12k+n1i+j
      X=A(ijkl)
      A(IJKL)=A(JIKL)
      A(jIKL)=X
   10 CONTINUE
      GO TO 100
   13 CONTINUE
      DO 20 L=1,N4      
         n123l=n123*(l-1)
      DO 20 I=1,N1      
         n12i=n123l+n12*(i-1)
      DO 20 J=1,N2      
         n1j=n1*(j-1)
      DO 20 K=1,I      
         ijkl=n123l+(k-1)*n12+n1j+i
         kjil=n12i+n1j+k
      X=A(IJKL)
      A(IJKL)=A(KJIL)
      A(KJIL)=X
   20 CONTINUE
      GO TO 100
   14 CONTINUE
      DO 25 I=1,N1
         n123i=n123*(i-1)
      DO 25 J=1,N2
         n1j=n1*(j-1)
      DO 25 K=1,N3
         n12k=n12*(k-1)+n1j
      DO 25 L=1,I
         ijkl=(l-1)*n123+n12k+i
         ljki=n123i     +n12k+l
      X=A(IJKL)
      A(IJKL)=A(LJKI)
      A(LJKI)= X
   25 CONTINUE
      GO TO 100
   23 CONTINUE
      DO 30 J=1,N2  
         j1=j-1
         n1j=n1*j1
         n12j=n12*j1
      DO 30 K=1,J   
         k1=k-1
         n1k=n1*k1
         n12k=n12*k1
      DO 30 L=1,N4
         n123l=n123*(l-1)
      DO 30 I=1,N1      
         n123li=n123l+i
         ijkl=n123li+n12k+n1j
         ikjl=n123li+n12j+n1k
      X=A(IJKL)
      A(IJKL)=A(IKJL)
      A(IKJL)=X
   30 CONTINUE
      GO TO 100
   24 CONTINUE
      DO 40 J=1,N2      
         j1=j-1
         n123j=n123*j1
         n1j  =n1*j1
      DO 40 L=1,J
         l1=l-1
         n123l=n123*l1
         n1l=n1*l1
      DO 40 K=1,N3
         n12k=n12*(k-1)
      DO 40 I=1,N1     
         n12ki=n12k+i
         ijkl=n123l+n12ki+n1j
         ilkj=n123j+n12ki+n1l
      X=A(IJKL)
      A(IJKL)=A(ILKJ)
      A(ILKJ)=X
   40 CONTINUE
      GO TO 100
   34 CONTINUE
      DO 50 K=1,N3      
         k1=k-1
         n12k=n12*k1
         n123k=n123*k1
      DO 50 L=1,K      
         l1=l-1
         n12l=n12*l1
         n123l=n123*l1
      DO 50 J=1,N2
         j1=j-1
         n1j=n1*(j-1)
      DO 50 I=1,N1    
         n1ji=n1j+i
         ijkl=n123l+n12k+n1ji
         ijlk=n123k+n12l+n1ji
      X=A(IJKL)
      A(IJKL)=A(IJLK)
      A(IJLK)=X
   50 CONTINUE
      GO TO 100
 231  CONTINUE
      DO 60 L=1,N4
         n123l=n123*(l-1)
      DO 60 J=1,N1
         j1=j-1
         n12j=n12*j1
         n1j=n1*j1
      DO 60 K=1,J
         n12jk=n12j+k
         k1=k-1
         n12k=n12*k1
         n1k=n1*k1
      DO 60 I=1,K
         i1=i-1
         ijkl=n123l+n12k+n1j+i
         jkil=n123l+i1*n12+n1k+j
         kijl=n123l+n12jk +i1*n1
      X=A(IJKL)
      A(IJKL)=A(JKIL)
      A(JKIL)=A(KIJL)
      A(KIJL)=X
      IF(J.EQ.K.OR.K.EQ.I) GOTO 60
         jikl=n123l+n12k  +n1*i1+j
         ikjl=n123l+n12j  +n1k  +i
         kjil=n123l+n12*i1+n1j  +k
      X=A(JIKL)
      A(JIKL)=A(IKJL)
      A(IKJL)=A(KJIL)
      A(KJIL)=X
 60   CONTINUE
      GOTO 100
 312  continue
      DO 70 L=1,N4
         n123l=n123*(l-1)
      DO 70 I=1,N1
         i1=i-1
         n1i=n1*i1
         n12i=n12*i1
      DO 70 J=1,I
         j1=j-1
         n1j =n1*j1
         n12j=n12*j1
      DO 70 K=1,J
         k1=k-1
         n1k=n1*k1
         n12k=n12*k1
         ijkl=n123l+n12k+n1j+i
         jkil=n123l+n12i+n1k+j
         kijl=n123l+n12j+n1i+k
      X=A(JKIL)
      A(JKIL)=A(IJKL)
      A(IJKL)=A(KIJL)
      A(KIJL)=X
      IF (I.EQ.J.OR.J.EQ.K)GOTO 70
         ikjl=n123l+n12j+n1k+i
         jikl=n123l+n12k+n1i+j
         kjil=n123l+n12i+n1j+k
      X=A(IKJL)
      A(IKJL)=A(JIKL)
      A(JIKL)=A(KJIL)
      A(KJIL)=X
 70   continue
      GOTO 100
 341  CONTINUE
      DO 80 L=1,N2
      DO 80 J=1,N1
      DO 80 K=1,J
      DO 80 I=1,K
         iljk=(k-1)*n123+(j-1)*n12+(l-1)*n1+i
         jlki=(i-1)*n123+(k-1)*n12+(l-1)*n1+j
         klij=(j-1)*n123+(i-1)*n12+(l-1)*n1+k
      X=A(ILJK)
      A(ILJK)=A(JLKI)
      A(JLKI)=A(KLIJ)
      A(KLIJ)=X
      IF(J.EQ.K.OR.K.EQ.I) GOTO 80
         ilkj=(j-1)*n123+(k-1)*n12+(l-1)*n1+i
         jlik=(k-1)*n123+(i-1)*n12+(l-1)*n1+j
         klji=(i-1)*n123+(j-1)*n12+(l-1)*n1+k
      X=A(JLIK)
      A(JLIK)=A(ILKJ)
      A(ILKJ)=A(KLJI)
      A(KLJI)=X
 80   CONTINUE
      GOTO 100
 413  CONTINUE
      DO 90 L=1,N2
      DO 90 I=1,N1
      DO 90 J=1,I
      DO 90 K=1,J
         jlki=(i-1)*n123+(k-1)*n12+(l-1)*n1+j
         iljk=(k-1)*n123+(j-1)*n12+(l-1)*n1+i
         klij=(j-1)*n123+(i-1)*n12+(l-1)*n1+k
      X=A(JLKI)
      A(JLKI)=A(ILJK)
      A(ILJK)=A(KLIJ)
      A(KLIJ)=X
      IF (I.EQ.J.OR.J.EQ.K)GOTO 90
         ilkj=(j-1)*n123+(k-1)*n12+(l-1)*n1+i
         jlik=(k-1)*n123+(i-1)*n12+(l-1)*n1+j
         klji=(i-1)*n123+(j-1)*n12+(l-1)*n1+k
      X=A(ILKJ)
      A(ILKJ)=A(JLIK)
      A(JLIK)=A(KLJI)
      A(KLJI)=X
 90   continue
      GO TO 100
 1234 CONTINUE
C      write(6,76)a
      DO 95 I=1,N1
      DO 95 J=1,N2
      DO 95 K=1,J
      DO 95 L=1,I
         ijkl=(l-1)*n123+(k-1)*n12+(j-1)*n1+i
         lkji=(i-1)*n123+(j-1)*n12+(k-1)*n1+l
      X=A(IJKL)
      A(IJKL)=A(LKJI)
      A(LKJI)=X
      if (i.eq.l.or.k.eq.j) goto 95
         ljki=(i-1)*n123+(k-1)*n12+(j-1)*n1+l
         ikjl=(l-1)*n123+(j-1)*n12+(k-1)*n1+i
      X=A(LJKI)
      A(LJKI)=A(IKJL)
      A(IKJL)=x
 95   CONTINUE
 100  CONTINUE
 77   FORMAT('IJKL:',4I3)
 76   format(4f15.10)
      RETURN
      END
      SUBROUTINE MTRANS(V,NU,ID)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION V(NU,NU,NU,NU)
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,
     *22,23),ID
    1 CONTINUE
      DO 100 K=1,NU
      DO 100 L=1,K
      DO 100 J=1,NU
      DO 100 I=1,NU
      X=V(I,J,K,L)
      V(I,J,K,L)=V(I,J,L,K)
      V(I,J,L,K)=X
  100 CONTINUE
      GO TO 1000
    2 CONTINUE
      DO 200 B=1,NU
      DO 200 C=1,B
      DO 200 D=1,C
      DO 200 A=1,NU
      X=V(A,D,B,C)
      V(A,D,B,C)=V(A,B,C,D)
      V(A,B,C,D)=V(A,C,D,B)
      V(A,C,D,B)=X
      IF(B.EQ.C.OR.C.EQ.D)GO TO 200
      X=V(A,B,D,C)
      V(A,B,D,C)=V(A,D,C,B)
      V(A,D,C,B)=V(A,C,B,D)
      V(A,C,B,D)=X
  200 CONTINUE
      GO TO 1000
    3 CONTINUE
      DO 300 B=1,NU
      DO 300 C=1,B
      DO 300 D=1,C
      DO 300 A=1,NU
      X=V(A,C,D,B)
      V(A,C,D,B)=V(A,B,C,D)
      V(A,B,C,D)=V(A,D,B,C)
      V(A,D,B,C)=X
      IF (B.EQ.C.OR.C.EQ.D)GO TO 300
      X=V(A,B,D,C)
      V(A,B,D,C)=V(A,C,B,D)
      V(A,C,B,D)=V(A,D,C,B)
      V(A,D,C,B)=X
  300 CONTINUE
      GO TO 1000
    4 CONTINUE
	call mtranrec(v,nu,3)
	call mtranrec(v,nu,8)
	goto 1000
      DO 400 A=1,NU
      DO 400 B=1,A
      DO 400 C=1,A
      LC=C
      IF (A.EQ.C)LC=B
      DO 400 D=1,LC
      X=V(D,C,A,B)
      V(D,C,A,B)=V(A,B,C,D)
      V(A,B,C,D)=V(C,D,B,A)
      V(C,D,B,A)=V(B,A,D,C)
      V(B,A,D,C)=X
      IF (C.EQ.D.OR.A.EQ.C.AND.B.EQ.D.OR.A.EQ.B)GO TO 400
      X=V(C,D,A,B)
      V(C,D,A,B)=V(A,B,D,C)
      V(A,B,D,C)=V(D,C,B,A)
      V(D,C,B,A)=V(B,A,C,D)
      V(B,A,C,D)=X
  400 CONTINUE
      GO TO 1000
    5 CONTINUE
	call mtranrec(v,nu,8)
	call mtranrec(v,nu,1)
	goto 1000
      DO 500 A=1,NU
      DO 500 C=1,A
      DO 500 D=1,C
      DO 500 B=1,NU
      X=V(C,B,D,A)
      V(C,B,D,A)=V(A,B,C,D)
      V(A,B,C,D)=V(D,B,A,C)
      V(D,B,A,C)=X
      IF (A.EQ.C.OR.C.EQ.D)GO TO 500
      X=V(A,B,D,C)
      V(A,B,D,C)=V(C,B,A,D)
      V(C,B,A,D)=V(D,B,C,A)
      V(D,B,C,A)=X
  500 CONTINUE
      GO TO 1000
    6 CONTINUE
      DO 600 D=1,NU
      DO 600 C=1,NU
      DO 600 A=1,NU
      DO 600 B=1,A
      X=V(B,A,C,D)
      V(B,A,C,D)=V(A,B,C,D)
      V(A,B,C,D)=X
  600 CONTINUE
      GO TO 1000
    7 CONTINUE
      DO 700 D=1,NU
      DO 700 B=1,NU
      DO 700 C=1,B
      DO 700 A=1,NU
      X=V(A,B,C,D)
      V(A,B,C,D)=V(A,C,B,D)
      V(A,C,B,D)=X
  700 CONTINUE
      GO TO 1000
    8 CONTINUE
      DO 800 D=1,NU
      DO 800 B=1,NU
      DO 800 A=1,NU
      DO 800 C=1,A
      X=V(A,B,C,D)
      V(A,B,C,D)=V(C,B,A,D)
      V(C,B,A,D)=X
  800 CONTINUE
      GO TO 1000
    9 CONTINUE
	call mtra1(v,nu,10)
	call mtranrec(v,nu,8)
	goto 1000
      DO 900 A=1,NU
      DO 900 B=1,NU
      DO 900 C=1,A
      LIMD=NU
      IF (A.EQ.C)LIMD=B
      DO 900 D=1,LIMD
      X=V(A,B,C,D)
      V(A,B,C,D)=V(C,D,A,B)
      V(C,D,A,B)=X
  900 CONTINUE
      GO TO 1000
 10   CONTINUE
	call mtra1(v,nu,10)
	goto 1000
      DO 950 A=1,NU
      DO 950 B=1,NU
      DO 950 C=1,NU
      DO 950 D=1,B
      X=V(A,B,C,D)
      V(A,B,C,D)=V(A,D,C,B)
      V(A,D,C,B)=X
 950  CONTINUE
      GOTO 1000
 11   CONTINUE
	call mtranrec(v,nu,1)
	call mtranrec(v,nu,8)
	call mtranrec(v,nu,1)
	goto 1000
      DO 960 A=1,NU
      DO 960 B=1,NU
      DO 960 C=1,NU
      DO 960 D=1,A
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,B,C,A)
      V(D,B,C,A)=X
 960  CONTINUE
      GOTO 1000
 12   CONTINUE
	call mtranrec(v,nu,7)
	call mtranrec(v,nu,6)
	goto 1000
      call tranmd(v,nu,nu,nu,nu,231)      
      GOTO 1000
 13   CONTINUE
	call mtranrec(v,nu,6)
	call mtranrec(v,nu,7)
	goto 1000
      call tranmd(v,nu,nu,nu,nu,312)
      GOTO 1000
 14   CONTINUE
	call mtranrec(v,nu,1)
	call mtranrec(v,nu,8)
      GO TO 1000
      DO 970 A=1,NU
      DO 970 B=1,NU
      DO 970 C=1,A
      DO 970 D=1,C
      X=V(D,B,A,C)
      V(D,B,A,C)=V(A,B,C,D)
      V(A,B,C,D)=V(C,B,D,A)
      V(C,B,D,A)=X
      IF (A.EQ.C.OR.C.EQ.D)GO TO 970
      X=V(D,B,C,A)
      V(D,B,C,A)=V(C,B,A,D)
      V(C,B,A,D)=V(A,B,D,C)
      V(A,B,D,C)=X
  970 CONTINUE
      GO TO 1000
 15   CONTINUE
      call mtra1(v,nu,10)
      call mtranrec(v,nu,6)
      GO TO 1000
      DO 980 A=1,NU
      DO 980 B=1,A
      DO 980 C=1,NU
      DO 980 D=1,B
      X=V(D,A,C,B)
      V(D,A,C,B)=V(A,B,C,D)
      V(A,B,C,D)=V(B,D,C,A)
      V(B,D,C,A)=X
      IF (A.EQ.B.OR.B.EQ.D)GO TO 980
      X=V(D,B,C,A)
      V(D,B,C,A)=V(B,A,C,D)
      V(B,A,C,D)=V(A,D,C,B)
      V(A,D,C,B)=X
 980  CONTINUE
      GO TO 1000
 16   CONTINUE
      call mtranrec(v,nu,6)
      call mtra1(v,nu,10)
      GO TO 1000
      DO 990 A=1,NU
      DO 990 B=1,A
      DO 990 C=1,NU
      DO 990 D=1,B
      X=V(B,D,C,A)
      V(B,D,C,A)=V(A,B,C,D)
      V(A,B,C,D)=V(D,A,C,B)
      V(D,A,C,B)=X
      IF (A.EQ.B.OR.B.EQ.D)GO TO 990
      X=V(A,D,C,B)
      V(A,D,C,B)=V(B,A,C,D)
      V(B,A,C,D)=V(D,B,C,A)
      V(D,B,C,A)=X
  990 CONTINUE
      GO TO 1000
 17   CONTINUE
      call mtranrec(v,nu,8)
      call mtranrec(v,nu,2)
      GO TO 1000
      DO 991 A=1,NU
      DO 991 B=1,A
      DO 991 C=1,A
      LC=C
      IF (A.EQ.C)LC=B
      DO 991 D=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,C,A,B)
      V(D,C,A,B)=V(B,A,D,C)
      V(B,A,D,C)=V(C,D,B,A)
      V(C,D,B,A)=X
      IF (C.EQ.D.OR.A.EQ.C.AND.B.EQ.D.OR.A.EQ.B)GO TO 991
      X=V(A,B,D,C)
      V(A,B,D,C)=V(C,D,A,B)
      V(C,D,A,B)=V(B,A,C,D)
      V(B,A,C,D)=V(D,C,B,A)
      V(D,C,B,A)=X
 991  CONTINUE
      GO TO 1000
 18   CONTINUE
	call mtranrec(v,nu,1)
	call mtranrec(v,nu,8)
	call mtranrec(v,nu,2)
      GO TO 1000
      DO 992 A=1,NU
      DO 992 B=1,NU
      DO 992 D=1,A
      LC=NU
      IF (A.EQ.d)LC=B
      DO 992 c=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,C,B,A)
      V(D,C,B,A)=X
 992  CONTINUE
      GO TO 1000
 19   CONTINUE
	call mtra1(v,nu,19)
	goto 1000
      DO 993 A=1,NU
c      write(6,*)'a=',a
      DO 993 B=1,a
      DO 993 C=1,nu
      DO 993 D=1,c
      X=V(A,B,C,D)
      V(A,B,C,D)=V(B,A,D,C)
      V(B,A,D,C)=X
      if(a.eq.b.or.c.eq.d)goto 9993
      x=v(a,b,d,c)
      v(a,b,d,c)=v(b,a,c,d)
      v(b,a,c,d)=x
 9993 continue
 993  CONTINUE
 909  format('a,b,c,d,v:',4i3,i8,5x,4i3,i8)
      GO TO 1000
 20   CONTINUE
      call mtranrec(v,nu,2)
      call mtranrec(v,nu,6)
      GO TO 1000
      DO 994 A=1,NU
      DO 994 B=1,A
      DO 994 C=1,A
      LD=B
      IF (A.EQ.B)Ld=C
      DO 994 D=1,LD
      X=V(A,B,C,D)
      V(A,B,C,D)=v(B,C,D,A)
      V(B,C,D,A)=V(C,D,A,B)
      V(C,D,A,B)=V(D,A,B,C)
      V(D,A,B,C)=X
      IF (B.EQ.D.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.C)GO TO 994
      X=V(A,D,C,B)
      V(A,D,C,B)=v(D,C,B,A)
      V(D,C,B,A)=v(C,B,A,D)
      V(C,B,A,D)=v(B,A,D,C)
      V(B,A,D,C)=X
 994  CONTINUE
      GO TO 1000
 21   CONTINUE
      call mtranrec(v,nu,6)
      call mtranrec(v,nu,3)
      GO TO 1000
      DO 995 A=1,NU
      DO 995 B=1,A
      DO 995 C=1,A
      LD=B
      IF (A.EQ.B)LD=C
      DO 995 D=1,LD
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,A,B,C)
      V(D,A,B,C)=V(C,D,A,B)
      V(C,D,A,B)=V(B,C,D,A)
      V(B,C,D,A)=X
      IF (B.EQ.D.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.C)GO TO 995
      X=V(A,D,C,B)
      V(A,D,C,B)=V(B,A,D,C)
      V(B,A,D,C)=V(C,B,A,D)
      V(C,B,A,D)=V(D,C,B,A)
      V(D,C,B,A)=X
 995  CONTINUE
      GO TO 1000
 22   CONTINUE
      call mtranrec(v,nu,3)
      call mtranrec(v,nu,6)
      GO TO 1000
      DO 996 A=1,NU
      DO 996 B=1,A
      DO 996 D=1,A
      LC=B
      IF (A.EQ.B)LC=D
      DO 996 C=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(B,D,A,C)
      V(B,D,A,C)=V(D,C,B,A)
      V(D,C,B,A)=V(C,A,D,B)
      V(C,A,D,B)=X
      IF (B.EQ.C.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.D)GO TO 996
      X=V(A,C,B,D)
      V(A,C,B,D)=V(C,D,A,B)
      V(C,D,A,B)=V(D,B,C,A)
      V(D,B,C,A)=V(B,A,D,C)
      V(B,A,D,C)=X
 996  CONTINUE
      GO TO 1000
 23   CONTINUE
      call mtranrec(v,nu,2)
      call mtranrec(v,nu,8)
      GO TO 1000
      DO 997 A=1,NU
      DO 997 B=1,A
      DO 997 D=1,A
      LC=B
      IF (A.EQ.B)LC=D
      DO 997 C=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(C,A,D,B)
      V(C,A,D,B)=V(D,C,B,A)
      V(D,C,B,A)=V(B,D,A,C)
      V(B,D,A,C)=X
      IF (B.EQ.C.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.D)GO TO 997
      X=V(A,C,B,D)
      V(A,C,B,D)=V(B,A,D,C)
      V(B,A,D,C)=V(D,B,C,A)
      V(D,B,C,A)=V(C,D,A,B)
      V(C,D,A,B)=X
 997  CONTINUE
      GO TO 1000
 1000 continue
      RETURN
      END
      SUBROUTINE MTRANREC(V,NU,ID)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION V(NU,NU,NU,NU)
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,
     *22,23),ID
    1 CONTINUE
      DO 100 K=1,NU
      DO 100 L=1,K
      DO 100 J=1,NU
      DO 100 I=1,NU
      X=V(I,J,K,L)
      V(I,J,K,L)=V(I,J,L,K)
      V(I,J,L,K)=X
  100 CONTINUE
      GO TO 1000
    2 CONTINUE
      DO 200 B=1,NU
      DO 200 C=1,B
      DO 200 D=1,C
      DO 200 A=1,NU
      X=V(A,D,B,C)
      V(A,D,B,C)=V(A,B,C,D)
      V(A,B,C,D)=V(A,C,D,B)
      V(A,C,D,B)=X
      IF(B.EQ.C.OR.C.EQ.D)GO TO 200
      X=V(A,B,D,C)
      V(A,B,D,C)=V(A,D,C,B)
      V(A,D,C,B)=V(A,C,B,D)
      V(A,C,B,D)=X
  200 CONTINUE
      GO TO 1000
    3 CONTINUE
      DO 300 B=1,NU
      DO 300 C=1,B
      DO 300 D=1,C
      DO 300 A=1,NU
      X=V(A,C,D,B)
      V(A,C,D,B)=V(A,B,C,D)
      V(A,B,C,D)=V(A,D,B,C)
      V(A,D,B,C)=X
      IF (B.EQ.C.OR.C.EQ.D)GO TO 300
      X=V(A,B,D,C)
      V(A,B,D,C)=V(A,C,B,D)
      V(A,C,B,D)=V(A,D,C,B)
      V(A,D,C,B)=X
  300 CONTINUE
      GO TO 1000
    4 CONTINUE
	call mtransnew(v,nu,3)
	call mtransnew(v,nu,8)
	goto 1000
      DO 400 A=1,NU
      DO 400 B=1,A
      DO 400 C=1,A
      LC=C
      IF (A.EQ.C)LC=B
      DO 400 D=1,LC
      X=V(D,C,A,B)
      V(D,C,A,B)=V(A,B,C,D)
      V(A,B,C,D)=V(C,D,B,A)
      V(C,D,B,A)=V(B,A,D,C)
      V(B,A,D,C)=X
      IF (C.EQ.D.OR.A.EQ.C.AND.B.EQ.D.OR.A.EQ.B)GO TO 400
      X=V(C,D,A,B)
      V(C,D,A,B)=V(A,B,D,C)
      V(A,B,D,C)=V(D,C,B,A)
      V(D,C,B,A)=V(B,A,C,D)
      V(B,A,C,D)=X
  400 CONTINUE
      GO TO 1000
    5 CONTINUE
      DO 500 A=1,NU
      DO 500 C=1,A
      DO 500 D=1,C
      DO 500 B=1,NU
      X=V(C,B,D,A)
      V(C,B,D,A)=V(A,B,C,D)
      V(A,B,C,D)=V(D,B,A,C)
      V(D,B,A,C)=X
      IF (A.EQ.C.OR.C.EQ.D)GO TO 500
      X=V(A,B,D,C)
      V(A,B,D,C)=V(C,B,A,D)
      V(C,B,A,D)=V(D,B,C,A)
      V(D,B,C,A)=X
  500 CONTINUE
      GO TO 1000
    6 CONTINUE
      DO 600 D=1,NU
      DO 600 C=1,NU
      DO 600 A=1,NU
      DO 600 B=1,A
      X=V(B,A,C,D)
      V(B,A,C,D)=V(A,B,C,D)
      V(A,B,C,D)=X
  600 CONTINUE
      GO TO 1000
    7 CONTINUE
      DO 700 D=1,NU
      DO 700 B=1,NU
      DO 700 C=1,B
      DO 700 A=1,NU
      X=V(A,B,C,D)
      V(A,B,C,D)=V(A,C,B,D)
      V(A,C,B,D)=X
  700 CONTINUE
      GO TO 1000
    8 CONTINUE
      DO 800 D=1,NU
      DO 800 B=1,NU
      DO 800 A=1,NU
      DO 800 C=1,A
      X=V(A,B,C,D)
      V(A,B,C,D)=V(C,B,A,D)
      V(C,B,A,D)=X
  800 CONTINUE
      GO TO 1000
    9 CONTINUE
      nu2=nu*nu
      GO TO 1000
 10   CONTINUE
      DO 950 A=1,NU
      DO 950 B=1,NU
      DO 950 C=1,NU
      DO 950 D=1,B
      X=V(A,B,C,D)
      V(A,B,C,D)=V(A,D,C,B)
      V(A,D,C,B)=X
 950  CONTINUE
      GOTO 1000
 11   CONTINUE
      DO 960 A=1,NU
      DO 960 B=1,NU
      DO 960 C=1,NU
      DO 960 D=1,A
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,B,C,A)
      V(D,B,C,A)=X
 960  CONTINUE
      GOTO 1000
 12   CONTINUE
      call tranmd(v,nu,nu,nu,nu,231)      
      GOTO 1000
 13   CONTINUE
      call tranmd(v,nu,nu,nu,nu,312)
      GOTO 1000
 14   CONTINUE
      DO 970 A=1,NU
      DO 970 B=1,NU
      DO 970 C=1,A
      DO 970 D=1,C
      X=V(D,B,A,C)
      V(D,B,A,C)=V(A,B,C,D)
      V(A,B,C,D)=V(C,B,D,A)
      V(C,B,D,A)=X
      IF (A.EQ.C.OR.C.EQ.D)GO TO 970
      X=V(D,B,C,A)
      V(D,B,C,A)=V(C,B,A,D)
      V(C,B,A,D)=V(A,B,D,C)
      V(A,B,D,C)=X
  970 CONTINUE
      GO TO 1000
 15   CONTINUE
      DO 980 A=1,NU
      DO 980 B=1,A
      DO 980 C=1,NU
      DO 980 D=1,B
      X=V(D,A,C,B)
      V(D,A,C,B)=V(A,B,C,D)
      V(A,B,C,D)=V(B,D,C,A)
      V(B,D,C,A)=X
      IF (A.EQ.B.OR.B.EQ.D)GO TO 980
      X=V(D,B,C,A)
      V(D,B,C,A)=V(B,A,C,D)
      V(B,A,C,D)=V(A,D,C,B)
      V(A,D,C,B)=X
 980  CONTINUE
      GO TO 1000
 16   CONTINUE
      DO 990 A=1,NU
      DO 990 B=1,A
      DO 990 C=1,NU
      DO 990 D=1,B
      X=V(B,D,C,A)
      V(B,D,C,A)=V(A,B,C,D)
      V(A,B,C,D)=V(D,A,C,B)
      V(D,A,C,B)=X
      IF (A.EQ.B.OR.B.EQ.D)GO TO 990
      X=V(A,D,C,B)
      V(A,D,C,B)=V(B,A,C,D)
      V(B,A,C,D)=V(D,B,C,A)
      V(D,B,C,A)=X
  990 CONTINUE
      GO TO 1000
 17   CONTINUE
      DO 991 A=1,NU
      DO 991 B=1,A
      DO 991 C=1,A
      LC=C
      IF (A.EQ.C)LC=B
      DO 991 D=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,C,A,B)
      V(D,C,A,B)=V(B,A,D,C)
      V(B,A,D,C)=V(C,D,B,A)
      V(C,D,B,A)=X
      IF (C.EQ.D.OR.A.EQ.C.AND.B.EQ.D.OR.A.EQ.B)GO TO 991
      X=V(A,B,D,C)
      V(A,B,D,C)=V(C,D,A,B)
      V(C,D,A,B)=V(B,A,C,D)
      V(B,A,C,D)=V(D,C,B,A)
      V(D,C,B,A)=X
 991  CONTINUE
      GO TO 1000
 18   CONTINUE
      DO 992 A=1,NU
      DO 992 B=1,NU
      DO 992 D=1,A
      LC=NU
      IF (A.EQ.d)LC=B
      DO 992 c=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,C,B,A)
      V(D,C,B,A)=X
 992  CONTINUE
      GO TO 1000
 19   CONTINUE
      DO 993 A=1,NU
c      write(6,*)'a=',a
      DO 993 B=1,a
      DO 993 C=1,nu
      DO 993 D=1,c
      X=V(A,B,C,D)
      V(A,B,C,D)=V(B,A,D,C)
      V(B,A,D,C)=X
      if(a.eq.b.or.c.eq.d)goto 9993
      x=v(a,b,d,c)
      v(a,b,d,c)=v(b,a,c,d)
      v(b,a,c,d)=x
 9993 continue
 993  CONTINUE
 909  format('a,b,c,d,v:',4i3,i8,5x,4i3,i8)
      GO TO 1000
 20   CONTINUE
      DO 994 A=1,NU
      DO 994 B=1,A
      DO 994 C=1,A
      LD=B
      IF (A.EQ.B)Ld=C
      DO 994 D=1,LD
      X=V(A,B,C,D)
      V(A,B,C,D)=v(B,C,D,A)
      V(B,C,D,A)=V(C,D,A,B)
      V(C,D,A,B)=V(D,A,B,C)
      V(D,A,B,C)=X
      IF (B.EQ.D.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.C)GO TO 994
      X=V(A,D,C,B)
      V(A,D,C,B)=v(D,C,B,A)
      V(D,C,B,A)=v(C,B,A,D)
      V(C,B,A,D)=v(B,A,D,C)
      V(B,A,D,C)=X
 994  CONTINUE
      GO TO 1000
 21   CONTINUE
      DO 995 A=1,NU
      DO 995 B=1,A
      DO 995 C=1,A
      LD=B
      IF (A.EQ.B)LD=C
      DO 995 D=1,LD
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,A,B,C)
      V(D,A,B,C)=V(C,D,A,B)
      V(C,D,A,B)=V(B,C,D,A)
      V(B,C,D,A)=X
      IF (B.EQ.D.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.C)GO TO 995
      X=V(A,D,C,B)
      V(A,D,C,B)=V(B,A,D,C)
      V(B,A,D,C)=V(C,B,A,D)
      V(C,B,A,D)=V(D,C,B,A)
      V(D,C,B,A)=X
 995  CONTINUE
      GO TO 1000
 22   CONTINUE
      DO 996 A=1,NU
      DO 996 B=1,A
      DO 996 D=1,A
      LC=B
      IF (A.EQ.B)LC=D
      DO 996 C=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(B,D,A,C)
      V(B,D,A,C)=V(D,C,B,A)
      V(D,C,B,A)=V(C,A,D,B)
      V(C,A,D,B)=X
      IF (B.EQ.C.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.D)GO TO 996
      X=V(A,C,B,D)
      V(A,C,B,D)=V(C,D,A,B)
      V(C,D,A,B)=V(D,B,C,A)
      V(D,B,C,A)=V(B,A,D,C)
      V(B,A,D,C)=X
 996  CONTINUE
      GO TO 1000
 23   CONTINUE
      DO 997 A=1,NU
      DO 997 B=1,A
      DO 997 D=1,A
      LC=B
      IF (A.EQ.B)LC=D
      DO 997 C=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(C,A,D,B)
      V(C,A,D,B)=V(D,C,B,A)
      V(D,C,B,A)=V(B,D,A,C)
      V(B,D,A,C)=X
      IF (B.EQ.C.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.D)GO TO 997
      X=V(A,C,B,D)
      V(A,C,B,D)=V(B,A,D,C)
      V(B,A,D,C)=V(D,B,C,A)
      V(D,B,C,A)=V(C,D,A,B)
      V(C,D,A,B)=X
 997  CONTINUE
      GO TO 1000
 1000 continue
      RETURN
      END
      SUBROUTINE MTRANSnew(V,NU,ID)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION V(NU,NU,NU,NU)
c      write(6,*)'entering mtrans:',nu,id
c      call ienter(2)
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,
     *22,23),ID
    1 CONTINUE
      DO 100 K=1,NU
      DO 100 L=1,K
      DO 100 J=1,NU
      DO 100 I=1,NU
      X=V(I,J,K,L)
      V(I,J,K,L)=V(I,J,L,K)
      V(I,J,L,K)=X
  100 CONTINUE
      GO TO 1000
    2 CONTINUE
      DO 200 B=1,NU
      DO 200 C=1,B
      DO 200 D=1,C
      DO 200 A=1,NU
      X=V(A,D,B,C)
      V(A,D,B,C)=V(A,B,C,D)
      V(A,B,C,D)=V(A,C,D,B)
      V(A,C,D,B)=X
      IF(B.EQ.C.OR.C.EQ.D)GO TO 200
      X=V(A,B,D,C)
      V(A,B,D,C)=V(A,D,C,B)
      V(A,D,C,B)=V(A,C,B,D)
      V(A,C,B,D)=X
  200 CONTINUE
      GO TO 1000
    3 CONTINUE
      DO 300 B=1,NU
      DO 300 C=1,B
      DO 300 D=1,C
      DO 300 A=1,NU
      X=V(A,C,D,B)
      V(A,C,D,B)=V(A,B,C,D)
      V(A,B,C,D)=V(A,D,B,C)
      V(A,D,B,C)=X
      IF (B.EQ.C.OR.C.EQ.D)GO TO 300
      X=V(A,B,D,C)
      V(A,B,D,C)=V(A,C,B,D)
      V(A,C,B,D)=V(A,D,C,B)
      V(A,D,C,B)=X
  300 CONTINUE
      GO TO 1000
    4 CONTINUE
	call mtranrec(v,nu,3)
	call mtranrec(v,nu,8)
	goto 1000
      DO 400 A=1,NU
      DO 400 B=1,A
      DO 400 C=1,A
      LC=C
      IF (A.EQ.C)LC=B
      DO 400 D=1,LC
      X=V(D,C,A,B)
      V(D,C,A,B)=V(A,B,C,D)
      V(A,B,C,D)=V(C,D,B,A)
      V(C,D,B,A)=V(B,A,D,C)
      V(B,A,D,C)=X
      IF (C.EQ.D.OR.A.EQ.C.AND.B.EQ.D.OR.A.EQ.B)GO TO 400
      X=V(C,D,A,B)
      V(C,D,A,B)=V(A,B,D,C)
      V(A,B,D,C)=V(D,C,B,A)
      V(D,C,B,A)=V(B,A,C,D)
      V(B,A,C,D)=X
  400 CONTINUE
      GO TO 1000
    5 CONTINUE
	call mtranrec(v,nu,8)
	call mtranrec(v,nu,1)
	goto 1000
      DO 500 A=1,NU
      DO 500 C=1,A
      DO 500 D=1,C
      DO 500 B=1,NU
      X=V(C,B,D,A)
      V(C,B,D,A)=V(A,B,C,D)
      V(A,B,C,D)=V(D,B,A,C)
      V(D,B,A,C)=X
      IF (A.EQ.C.OR.C.EQ.D)GO TO 500
      X=V(A,B,D,C)
      V(A,B,D,C)=V(C,B,A,D)
      V(C,B,A,D)=V(D,B,C,A)
      V(D,B,C,A)=X
  500 CONTINUE
      GO TO 1000
    6 CONTINUE
      DO 600 D=1,NU
      DO 600 C=1,NU
      DO 600 A=1,NU
      DO 600 B=1,A
      X=V(B,A,C,D)
      V(B,A,C,D)=V(A,B,C,D)
      V(A,B,C,D)=X
  600 CONTINUE
c      do 600 d=1,nu
c      do 600 c=1,nu
c      call transq(v(1,1,c,d),nu)
c 600  continue
      GO TO 1000
    7 CONTINUE
      DO 700 D=1,NU
      DO 700 B=1,NU
      DO 700 C=1,B
      DO 700 A=1,NU
      X=V(A,B,C,D)
      V(A,B,C,D)=V(A,C,B,D)
      V(A,C,B,D)=X
  700 CONTINUE
      GO TO 1000
    8 CONTINUE
      DO 800 D=1,NU
      DO 800 B=1,NU
      DO 800 A=1,NU
      DO 800 C=1,A
      X=V(A,B,C,D)
      V(A,B,C,D)=V(C,B,A,D)
      V(C,B,A,D)=X
  800 CONTINUE
      GO TO 1000
    9 CONTINUE
	call mtra1(v,nu,10)
	call mtranrec(v,nu,8)
	goto 1000
      DO 900 A=1,NU
      DO 900 B=1,NU
      DO 900 C=1,A
      LIMD=NU
      IF (A.EQ.C)LIMD=B
      DO 900 D=1,LIMD
      X=V(A,B,C,D)
      V(A,B,C,D)=V(C,D,A,B)
      V(C,D,A,B)=X
  900 CONTINUE
c      nu2=nu*nu
c      call transq1(v,nu2)
      GO TO 1000
 10   CONTINUE
	call mtra1(v,nu,10)
	goto 1000
      DO 950 A=1,NU
      DO 950 B=1,NU
      DO 950 C=1,NU
      DO 950 D=1,B
      X=V(A,B,C,D)
      V(A,B,C,D)=V(A,D,C,B)
      V(A,D,C,B)=X
 950  CONTINUE
      GOTO 1000
 11   CONTINUE
	call mtranrec(v,nu,1)
	call mtranrec(v,nu,8)
	call mtranrec(v,nu,1)
	goto 1000
      DO 960 A=1,NU
      DO 960 B=1,NU
      DO 960 C=1,NU
      DO 960 D=1,A
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,B,C,A)
      V(D,B,C,A)=X
 960  CONTINUE
      GOTO 1000
 12   CONTINUE
	call mtranrec(v,nu,7)
	call mtranrec(v,nu,6)
	goto 1000
      call tranmd(v,nu,nu,nu,nu,231)      
      GOTO 1000
 13   CONTINUE
	call mtranrec(v,nu,6)
	call mtranrec(v,nu,7)
	goto 1000
      call tranmd(v,nu,nu,nu,nu,312)
      GOTO 1000
 14   CONTINUE
	call mtranrec(v,nu,1)
	call mtranrec(v,nu,8)
      GO TO 1000
      DO 970 A=1,NU
      DO 970 B=1,NU
      DO 970 C=1,A
      DO 970 D=1,C
      X=V(D,B,A,C)
      V(D,B,A,C)=V(A,B,C,D)
      V(A,B,C,D)=V(C,B,D,A)
      V(C,B,D,A)=X
      IF (A.EQ.C.OR.C.EQ.D)GO TO 970
      X=V(D,B,C,A)
      V(D,B,C,A)=V(C,B,A,D)
      V(C,B,A,D)=V(A,B,D,C)
      V(A,B,D,C)=X
  970 CONTINUE
      GO TO 1000
 15   CONTINUE
      call mtra1(v,nu,10)
      call mtranrec(v,nu,6)
      GO TO 1000
      DO 980 A=1,NU
      DO 980 B=1,A
      DO 980 C=1,NU
      DO 980 D=1,B
      X=V(D,A,C,B)
      V(D,A,C,B)=V(A,B,C,D)
      V(A,B,C,D)=V(B,D,C,A)
      V(B,D,C,A)=X
      IF (A.EQ.B.OR.B.EQ.D)GO TO 980
      X=V(D,B,C,A)
      V(D,B,C,A)=V(B,A,C,D)
      V(B,A,C,D)=V(A,D,C,B)
      V(A,D,C,B)=X
 980  CONTINUE
      GO TO 1000
 16   CONTINUE
      call mtranrec(v,nu,6)
      call mtra1(v,nu,10)
      GO TO 1000
      DO 990 A=1,NU
      DO 990 B=1,A
      DO 990 C=1,NU
      DO 990 D=1,B
      X=V(B,D,C,A)
      V(B,D,C,A)=V(A,B,C,D)
      V(A,B,C,D)=V(D,A,C,B)
      V(D,A,C,B)=X
      IF (A.EQ.B.OR.B.EQ.D)GO TO 990
      X=V(A,D,C,B)
      V(A,D,C,B)=V(B,A,C,D)
      V(B,A,C,D)=V(D,B,C,A)
      V(D,B,C,A)=X
  990 CONTINUE
      GO TO 1000
 17   CONTINUE
      call mtranrec(v,nu,8)
      call mtranrec(v,nu,2)
      GO TO 1000
      DO 991 A=1,NU
      DO 991 B=1,A
      DO 991 C=1,A
      LC=C
      IF (A.EQ.C)LC=B
      DO 991 D=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,C,A,B)
      V(D,C,A,B)=V(B,A,D,C)
      V(B,A,D,C)=V(C,D,B,A)
      V(C,D,B,A)=X
      IF (C.EQ.D.OR.A.EQ.C.AND.B.EQ.D.OR.A.EQ.B)GO TO 991
      X=V(A,B,D,C)
      V(A,B,D,C)=V(C,D,A,B)
      V(C,D,A,B)=V(B,A,C,D)
      V(B,A,C,D)=V(D,C,B,A)
      V(D,C,B,A)=X
 991  CONTINUE
      GO TO 1000
 18   CONTINUE
	call mtranrec(v,nu,1)
	call mtranrec(v,nu,8)
	call mtranrec(v,nu,2)
      GO TO 1000
      DO 992 A=1,NU
      DO 992 B=1,NU
      DO 992 D=1,A
      LC=NU
      IF (A.EQ.d)LC=B
      DO 992 c=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,C,B,A)
      V(D,C,B,A)=X
 992  CONTINUE
      GO TO 1000
 19   CONTINUE
	call mtra1(v,nu,19)
	goto 1000
      DO 993 A=1,NU
c      write(6,*)'a=',a
      DO 993 B=1,a
      DO 993 C=1,nu
      DO 993 D=1,c
      X=V(A,B,C,D)
      V(A,B,C,D)=V(B,A,D,C)
      V(B,A,D,C)=X
      if(a.eq.b.or.c.eq.d)goto 9993
      x=v(a,b,d,c)
      v(a,b,d,c)=v(b,a,c,d)
      v(b,a,c,d)=x
 9993 continue
c      write(6,909)a,b,c,d,v(a,b,c,d),b,a,d,c,v(b,a,d,c)
 993  CONTINUE
 909  format('a,b,c,d,v:',4i3,i8,5x,4i3,i8)
      GO TO 1000
 20   CONTINUE
      call mtranrec(v,nu,2)
      call mtranrec(v,nu,6)
      GO TO 1000
      DO 994 A=1,NU
      DO 994 B=1,A
      DO 994 C=1,A
      LD=B
      IF (A.EQ.B)Ld=C
      DO 994 D=1,LD
      X=V(A,B,C,D)
      V(A,B,C,D)=v(B,C,D,A)
      V(B,C,D,A)=V(C,D,A,B)
      V(C,D,A,B)=V(D,A,B,C)
      V(D,A,B,C)=X
      IF (B.EQ.D.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.C)GO TO 994
      X=V(A,D,C,B)
      V(A,D,C,B)=v(D,C,B,A)
      V(D,C,B,A)=v(C,B,A,D)
      V(C,B,A,D)=v(B,A,D,C)
      V(B,A,D,C)=X
 994  CONTINUE
      GO TO 1000
 21   CONTINUE
      call mtranrec(v,nu,6)
      call mtranrec(v,nu,3)
      GO TO 1000
      DO 995 A=1,NU
      DO 995 B=1,A
      DO 995 C=1,A
      LD=B
      IF (A.EQ.B)LD=C
      DO 995 D=1,LD
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,A,B,C)
      V(D,A,B,C)=V(C,D,A,B)
      V(C,D,A,B)=V(B,C,D,A)
      V(B,C,D,A)=X
      IF (B.EQ.D.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.C)GO TO 995
      X=V(A,D,C,B)
      V(A,D,C,B)=V(B,A,D,C)
      V(B,A,D,C)=V(C,B,A,D)
      V(C,B,A,D)=V(D,C,B,A)
      V(D,C,B,A)=X
 995  CONTINUE
      GO TO 1000
 22   CONTINUE
      call mtranrec(v,nu,3)
      call mtranrec(v,nu,6)
      GO TO 1000
      DO 996 A=1,NU
      DO 996 B=1,A
      DO 996 D=1,A
      LC=B
      IF (A.EQ.B)LC=D
      DO 996 C=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(B,D,A,C)
      V(B,D,A,C)=V(D,C,B,A)
      V(D,C,B,A)=V(C,A,D,B)
      V(C,A,D,B)=X
      IF (B.EQ.C.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.D)GO TO 996
      X=V(A,C,B,D)
      V(A,C,B,D)=V(C,D,A,B)
      V(C,D,A,B)=V(D,B,C,A)
      V(D,B,C,A)=V(B,A,D,C)
      V(B,A,D,C)=X
 996  CONTINUE
      GO TO 1000
 23   CONTINUE
      call mtranrec(v,nu,2)
      call mtranrec(v,nu,8)
      GO TO 1000
      DO 997 A=1,NU
      DO 997 B=1,A
      DO 997 D=1,A
      LC=B
      IF (A.EQ.B)LC=D
      DO 997 C=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(C,A,D,B)
      V(C,A,D,B)=V(D,C,B,A)
      V(D,C,B,A)=V(B,D,A,C)
      V(B,D,A,C)=X
      IF (B.EQ.C.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.D)GO TO 997
      X=V(A,C,B,D)
      V(A,C,B,D)=V(B,A,D,C)
      V(B,A,D,C)=V(D,B,C,A)
      V(D,B,C,A)=V(C,D,A,B)
      V(C,D,A,B)=X
 997  CONTINUE
      GO TO 1000
 1000 continue
      RETURN
      END
      SUBROUTINE MTRA1(V,NU,IDENT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      LOGICAL AEB,AEC,AED,BEC,BED,CED,ACD,ABD
      DIMENSION V(1)
      NU2=NU*NU
      NU3=NU2*NU
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,
     *22,23),IDENT
    1 CONTINUE
      DO 100 J=1,NU
         j1=(j-1)*nu
         DO 101 K=1,NU
            k0=k-1
            k2=k0*nu2
            k3=k0*nu3
            DO 102 L=1,K
               l0=l-1
               jkl=j1+k2+l0*nu3
               jlk=j1+k3+l0*nu2
               DO 103 I=1,NU
                  ijkl=jkl+i
                  ijlk=jlk+i
                  X=V(ijkl)
                  V(ijkl)=V(ijlk)
                  V(ijlk)=X
 103           CONTINUE
 102        CONTINUE
 101     CONTINUE
 100  CONTINUE
      GO TO 1000
 2    CONTINUE
      DO 120 B=1,NU
         ib=b-1
         ib1=ib*nu
         ib2=ib*nu2
         ib3=ib*nu3
         DO 121 C=1,B
            BEC=B.EQ.C
            ic=c-1
            ic1=ic*nu
            ic2=ic*nu2
            ic3=ic*nu3
            DO 122 D=1,C
               CED=C.EQ.D
               id=d-1
               id1=id*nu
               id2=id*nu2
               id3=id*nu3
               idbc=id1+ib2+ic3
               ibcd=ib1+ic2+id3
               icdb=ic1+id2+ib3
               if(bec.OR.ced)GOTO 124
               ibdc=ib1+id2+ic3
               idcb=id1+ic2+ib3
               icbd=ic1+ib2+id3
 124           CONTINUE
               DO 123 A=1,NU
                  iadbc=idbc+a
                  iabcd=ibcd+a
                  iacdb=icdb+a
                  X=V(iadbc)
                  V(iadbc)=V(iabcd)
                  V(iabcd)=V(iacdb)
                  V(iacdb)=X
                  IF(BEC.OR.CED)GO TO 123
                  iabdc=a+ibdc
                  iadcb=a+idcb
                  iacbd=a+icbd
                  X=V(iabdc)
                  V(iabdc)=V(iadcb)
                  V(iadcb)=V(iacbd)
                  V(iacbd)=X
 123           CONTINUE
 122        CONTINUE
 121     CONTINUE
 120  CONTINUE
      GO TO 1000
    3 CONTINUE
      DO 300 B=1,NU
         ib=b-1
         ib1=ib*nu
         ib2=ib*nu2
         ib3=ib*nu3
         DO 301 C=1,B
            ic=c-1
            ic1=ic*nu
            ic2=ic*nu2
            ic3=ic*nu3
            DO 302 D=1,C
               id=d-1
               id1=id*nu
               id2=id*nu2
               id3=id*nu3
               icdb=ic1+id2+ib3
               ibcd=ib1+ic2+id3
               idbc=id1+ib2+ic3
               if(b.ne.c.and.c.ne.d)then
                  ibdc=ib1+id2+ic3
                  icbd=ic1+ib2+id3
                  idcb=id1+ic2+ib3
               endif
               DO 303 A=1,NU
                  iacdb=a+icdb
                  iabcd=a+ibcd
                  iadbc=a+idbc
                  X=V(iacdb)
                  V(IACDB)=V(IABCD)
                  V(IABCD)=V(IADBC)
                  V(IADBC)=X
                  IF (B.EQ.C.OR.C.EQ.D)GO TO 303
                  iabdc=a+ibdc
                  iacbd=a+icbd
                  iadcb=a+idcb
                  X=V(IABDC)
                  V(IABDC)=V(IACBD)
                  V(IACBD)=V(IADCB)
                  V(IADCB)=X
 303           CONTINUE
 302        CONTINUE
 301     CONTINUE
 300  CONTINUE
      GO TO 1000
    4 CONTINUE
      DO 400 A=1,NU
         ia=a-1
         ia1=ia*nu
         ia2=ia*nu2
         ia3=ia*nu3
         DO 401 B=1,A
            aeb=a.eq.b
            ib=b-1
            ib1=ib*nu
            ib2=ib*nu2
            ib3=ib*nu3
            DO 402 C=1,A
               aec=a.eq.c
               ic=c-1
               ic1=ic*nu
               ic2=ic*nu2
               ic3=ic*nu3
               icab=ic1+ia2+ib3
               iabc=a  +ib1+ic2
               icba=c  +ib2+ia3
               ibac=b  +ia1+ic3
               LC=C
               IF (AEC)LC=B
               id1=-nu
               id2=-nu2
               id3=-nu3
               DO 403 D=1,LC
                  ced=c.eq.d
                  bed=b.eq.d
c                  id=d-1
c                  id1=id*nu
c                  id2=id*nu2
c                  id3=id*nu3
                  id1=id1+nu
                  id2=id2+nu2
                  id3=id3+nu3
                  idcab=d   +icab
                  iabcd=iabc+id3
                  icdba=icba+id1
                  ibadc=ibac+id2
                  X=V(IDCAB)
                  V(IDCAB)=V(IABCD)
                  V(IABCD)=V(ICDBA)
                  V(ICDBA)=V(IBADC)
                  V(IBADC)=X
                  IF (CED.OR.AEC.AND.BED.OR.AEB)GO TO 403
                  icdab=c+id1+ia2+ib3
                  iabdc=a+ib1+id2+ic3
                  idcba=d+ic1+ib2+ia3
                  ibacd=b+ia1+ic2+id3
                  X=V(ICDAB)
                  V(ICDAB)=V(IABDC)
                  V(IABDC)=V(IDCBA)
                  V(IDCBA)=V(IBACD)
                  V(IBACD)=X
 403           CONTINUE
 402        CONTINUE
 401     CONTINUE
 400  CONTINUE
      GO TO 1000
    5 CONTINUE
      DO 500 A=1,NU
         ia=a-1
         ia2=ia*nu2
         ia3=ia*nu3
         DO 501 C=1,A
            AEC=A.EQ.C
            ic=c-1
            ic2=ic*nu2
            ic3=ic*nu3
            DO 502 D=1,C
               ACD=AEC.OR.C.EQ.D
               id=d-1
               id2=id*nu2
               id3=id*nu3
               icda=c+id2+ia3
               iacd=a+ic2+id3
               idac=d+ia2+ic3
               IF(ACD)GOTO 504
               iadc=a+id2+ic3
               icad=c+ia2+id3
               idca=d+ic2+ia3
 504           CONTINUE
               ib1=-nu
              DO 503 B=1,NU
                  ib1=ib1+nu
c                  ib=b-1
c                  ib1=ib*nu
                  icbda=ib1+icda
                  iabcd=ib1+iacd
                  idbac=ib1+idac
                  X=V(ICBDA)
                  V(ICBDA)=V(IABCD)
                  V(IABCD)=V(IDBAC)
                  V(IDBAC)=X
                  IF (ACD)GO TO 503
                  iabdc=ib1+iadc
                  icbad=ib1+icad
                  idbca=ib1+idca
                  X=V(IABDC)
                  V(IABDC)=V(ICBAD)
                  V(ICBAD)=V(IDBCA)
                  V(IDBCA)=X
 503           CONTINUE
 502        CONTINUE
 501     CONTINUE
 500  CONTINUE
      GO TO 1000
 6    CONTINUE
      DO 600 D=1,NU
         id=d-1
         id3=id*nu3
         DO 601 C=1,NU
            ic=c-1
            ic2=ic*nu2
            DO 602 A=1,NU
               ia=a-1
               ia1=ia*nu
               iacd =ia1+ic2+id3
               iacd_=a+ic2+id3
               ib1=-nu
               ib2=-nu2
               ib3=-nu3
               DO 603 B=1,A
                  ib1=ib1+nu
                  ib2=ib2+nu2
                  ib3=ib3+nu3
c                  ib=b-1
c                  ib1=ib*nu
c                  ib2=ib*nu2
c                  ib3=ib*nu3
                  ibacd=b  +iacd
                  iabcd=ib1+iacd_
                  X=V(IBACD)
                  V(IBACD)=V(IABCD)
                  V(IABCD)=X
 603           CONTINUE
 602        CONTINUE
 601     CONTINUE
 600  CONTINUE
      GO TO 1000
    7 CONTINUE
      DO 700 D=1,NU
         id=d-1
         id3=id*nu3
         DO 701 B=1,NU
            ib=b-1
            ib1=ib*nu
            ib2=ib*nu2
            DO 702 C=1,B
               ic=c-1
               ic1=ic*nu
               ic2=ic*nu2
               ibcd=ib1+ic2+id3
               icbd=ic1+ib2+id3
               DO 703 A=1,NU
                  iabcd=a+ibcd
                  iacbd=a+icbd
                  X=V(IABCD)
                  V(IABCD)=V(IACBD)
                  V(IACBD)=X
 703           CONTINUE
 702        CONTINUE
 701     CONTINUE
 700  CONTINUE
      GO TO 1000
 8    CONTINUE
      DO 800 D=1,NU
         id=d-1
         id3=id*nu3
         DO 801 B=1,NU
            ib=b-1
            ib1=ib*nu
            DO 802 A=1,NU
               ia=a-1
               ia2=ia*nu2
               iabd=a  +ib1+id3
               ibad=ib1+ia2+id3
               ic2=-nu2
               DO 803 C=1,A
                  ic2=ic2+nu2
                  iabcd=iabd+ic2
                  icbad=c   +ibad
                  X=V(IABCD)
                  V(IABCD)=V(ICBAD)
                  V(ICBAD)=X
 803           CONTINUE
 802        CONTINUE
 801     CONTINUE
 800  CONTINUE
      GO TO 1000
    9 CONTINUE
      DO 900 A=1,NU2
         ia=a-1
         ia2=ia*nu2
         ib2=-nu2
         DO 901 B=1,A
            ib2=ib2+nu2
            iab=a+ib2
            iba=b+ia2
            X=V(IAB)
            V(IAB)=V(IBA)
            V(IBA)=X
 901     CONTINUE
 900  CONTINUE
      GO TO 1000
 10   CONTINUE
      DO 950 B=1,NU
         ib=b-1
         ib1=ib*nu
         ib3=ib*nu3
         DO 951 D=1,B
            id=d-1
            id1=id*nu
            id3=id*nu3
            DO 952 C=1,NU
               ic=c-1
               ic2=ic*nu2
               ibcd=ib1+ic2+id3
               idcb=id1+ic2+ib3
               DO 953 A=1,NU
                  iabcd=a+ibcd
                  iadcb=a+idcb
                  X=V(IABCD)
                  V(IABCD)=V(IADCB)
                  V(IADCB)=X
 953           CONTINUE
 952        CONTINUE
 951     CONTINUE
 950  CONTINUE
      GOTO 1000
 11   CONTINUE
      DO 960 A=1,NU
         ia=a-1
         ia3=ia*nu3
         DO 961 D=1,A
            id=d-1
            id3=id*nu3
            ic2=-nu2
            iad=a+id3
            ida=d+ia3
            DO 962 C=1,NU
c               ic2=n2(c)
               ic2=ic2+nu2
c               ic=c-1
c               ic2=ic*nu2
               iacd=ic2+iad
               idca=ic2+ida
               ib1=-nu
               DO 963 B=1,NU
                  ib1=ib1+nu
c                  ib1=n1(b)
                  iabcd=ib1+iacd
                  idbca=ib1+idca
                  X=V(IABCD)
                  V(IABCD)=V(IDBCA)
                  V(IDBCA)=X
 963           CONTINUE
 962        CONTINUE
 961     CONTINUE
 960  CONTINUE
      GOTO 1000
 12   CONTINUE
      call tranmd(v,nu,nu,nu,nu,231)      
      GOTO 1000
 13   CONTINUE
      call tranmd(v,nu,nu,nu,nu,312)
      GOTO 1000
 14   CONTINUE
      DO 970 A=1,NU
         ia=a-1
         ia2=ia*nu2
         ia3=ia*nu3
         DO 971 C=1,A
            AEC=A.EQ.C
            ic=c-1
            ic2=ic*nu2
            ic3=ic*nu3
            DO 972 D=1,C
               id=d-1
               id2=id*nu2
               id3=id*nu3
               ACD=AEC.OR.C.EQ.D
               idac=d+ia2+ic3
               iacd=a+ic2+id3
               icda=c+id2+ia3
               IF (ACD)GOTO 974
               idca=d+ic2+ia3
               icad=c+ia2+id3
               iadc=a+id2+ic3
 974           CONTINUE
               ib1=-nu
               DO 973 B=1,NU
                  ib1=ib1+nu
c                  ib=b-1
c                  ib1=ib*nu
                  idbac=ib1+idac
                  iabcd=ib1+iacd
                  icbda=ib1+icda
                  X=V(IDBAC)
                  V(IDBAC)=V(IABCD)
                  V(IABCD)=V(ICBDA)
                  V(ICBDA)=X
                  IF (ACD)GO TO 973
                  idbca=ib1+idca
                  icbad=ib1+icad
                  iabdc=ib1+iadc
                  X=V(IDBCA)
                  V(IDBCA)=V(ICBAD)
                  V(ICBAD)=V(IABDC)
                  V(IABDC)=X
 973           CONTINUE
 972        CONTINUE
 971     CONTINUE
 970  CONTINUE
      GO TO 1000
 15   CONTINUE
      DO 980 A=1,NU
         ia=a-1
         ia1=ia*nu
         ia2=ia*nu2
         ia3=ia*nu3
         DO 981 B=1,A
            ib=b-1
            ib1=ib*nu
            ib2=ib*nu2
            ib3=ib*nu3
            AEB=A.EQ.B
            DO 982 D=1,B
               id=d-1
               id1=id*nu
               id2=id*nu2
               id3=id*nu3
               idab=d+ia1+ib3
               iabd=a+ib1+id3
               ibda=b+id1+ia3
               ABD=AEB.OR.B.EQ.D
               IF(ABD)GOTO 984
               idba=d+ib1+ia3
               ibad=b+ia1+id3
               iadb=a+id1+ib3
 984           CONTINUE
               ic2=-nu2
               DO 983 C=1,NU
                  ic2=ic2+nu2
c                  ic=c-1
c                  ic2=ic*nu2
                  idacb=ic2+idab
                  iabcd=ic2+iabd
                  ibdca=ic2+ibda
                  X=V(IDACB)
                  V(IDACB)=V(IABCD)
                  V(IABCD)=V(IBDCA)
                  V(IBDCA)=X
                  IF (ABD)GO TO 983
                  idbca=ic2+idba
                  ibacd=ic2+ibad
                  iadcb=ic2+iadb
                  X=V(IDBCA)
                  V(IDBCA)=V(IBACD)
                  V(IBACD)=V(IADCB)
                  V(IADCB)=X
 983           CONTINUE
 982        CONTINUE
 981     CONTINUE
 980  CONTINUE
      GO TO 1000
 16   CONTINUE
      DO 990 A=1,NU
         ia=a-1
         ia1=ia*nu
         ia2=ia*nu2
         ia3=ia*nu3
         DO 991 B=1,A
            AEB=A.EQ.B
            ib=b-1
            ib1=ib*nu
            ib2=ib*nu2
            ib3=ib*nu3
            DO 992 D=1,B
               ABD=AEB.OR.B.EQ.D
               id=d-1
               id1=id*nu
               id2=id*nu2
               id3=id*nu3
               ibda=b+id1+ia3
               iabd=a+ib1+id3
               idab=d+ia1+ib3
               IF (ABD)GO TO 994
               iadb=a+id1+ib3
               ibad=b+ia1+id3
               idba=d+ib1+ia3
 994           CONTINUE
               ic2=-nu2
               DO 993 C=1,NU
                  ic2=ic2+nu2
                  ibdca=ic2+ibda
                  iabcd=ic2+iabd
                  idacb=ic2+idab
                  X=V(IBDCA)
                  V(IBDCA)=V(IABCD)
                  V(IABCD)=V(IDACB)
                  V(IDACB)=X
                  IF (ABD)GO TO 993
                  iadcb=ic2+iadb
                  ibacd=ic2+ibad
                  idbca=ic2+idba
                  X=V(IADCB)
                  V(IADCB)=V(IBACD)
                  V(IBACD)=V(IDBCA)
                  V(IDBCA)=X
 993           CONTINUE
 992        CONTINUE
 991     CONTINUE
 990  CONTINUE
      GO TO 1000
 17   CONTINUE
      DO 170 A=1,NU
         ia=a-1
         ia1=ia*nu
         ia2=ia*nu2
         ia3=ia*nu3
         DO 171 B=1,A
            AEB=A.EQ.B
            ib=b-1
            ib1=ib*nu
            ib2=ib*nu2
            ib3=ib*nu3
            DO 172 C=1,A
               AEC=A.EQ.C
               ic=c-1
               ic1=ic*nu
               ic2=ic*nu2
               ic3=ic*nu3
               LD=C
               IF (AEC)LD=B
               iabc=a  +ib1+ic2
               icab=ic1+ia2+ib3
               ibac=b  +ia1+ic3
               icba=c  +ib2+ia3
               IF (AEB)GO TO 174
               iabc_=a  +ib1+ic3
               icab_=c  +ia2+ib3
               ibac_=b  +ia1+ic2
               icba_=ic1+ib2+ia3
 174           CONTINUE
               id1=-nu
               id2=-nu2
               id3=-nu3
               DO 173 D=1,LD
                  CED=C.EQ.D
                  BED=B.EQ.D
                  id1=id1+nu
                  id2=id2+nu2
                  id3=id3+nu3
c                  id=d-1
c                  id1=id*nu
c                  id2=id*nu2
c                  id3=id*nu3
                  iabcd=id3+iabc
                  idcab=d  +icab
                  ibadc=id2+ibac
                  icdba=id1+icba
                  X=V(IABCD)
                  V(IABCD)=V(IDCAB)
                  V(IDCAB)=V(IBADC)
                  V(IBADC)=V(ICDBA)
                  V(ICDBA)=X
                  IF (CED.OR.AEC.AND.BED.OR.AEB)GO TO 173
                  iabdc=id2+iabc_
                  icdab=id1+icab_
                  ibacd=id3+ibac_
                  idcba=d  +icba_
                  X=V(IABDC)
                  V(IABDC)=V(ICDAB)
                  V(ICDAB)=V(IBACD)
                  V(IBACD)=V(IDCBA)
                  V(IDCBA)=X
 173           CONTINUE
 172        CONTINUE
 171     CONTINUE
 170  CONTINUE
      GO TO 1000
 18   CONTINUE
      DO 180 A=1,NU
         ia=a-1
         ia3=ia*nu3
         DO 181 B=1,NU
            ib=b-1
            ib1=ib*nu
            ib2=ib*nu2
            DO 182 D=1,A
               id=d-1
               id3=id*nu3
               iabd=a+ib1+id3
               idba=d+ib2+ia3
               LC=NU
               IF (A.EQ.D)LC=B
               ic1=-nu
               ic2=-nu2
               DO 183 c=1,LC
                  ic1=ic1+nu
                  ic2=ic2+nu2
c                  ic=c-1
c                  ic1=ic*nu
c                  ic2=ic*nu2
                  iabcd=ic2+iabd
                  idcba=ic1+idba
                  X=V(IABCD)
                  V(IABCD)=V(IDCBA)
                  V(IDCBA)=X
 183           CONTINUE
 182        CONTINUE
 181     CONTINUE
 180  CONTINUE
      GO TO 1000
 19   CONTINUE
      DO 190 C=1,nu
         ic=c-1
         ic2=ic*nu2
         ic3=ic*nu3
         DO 191 D=1,c
            CED=C.EQ.D
            id=d-1
            id2=id*nu2
            id3=id*nu3
            DO 192 A=1,NU
               ia=a-1
               ia1=ia*nu
               iacd=a  +ic2+id3
               iadc=ia1+id2+ic3
               if(CED)goto 194
               iadc_=a  +id2+ic3
               iacd_=ia1+ic2+id3
 194           CONTINUE
               ib1=-nu
               DO 193 B=1,a
                  AEB=A.EQ.B
c                  ib=b-1
                  ib1=ib1+nu
c                  ib1=ib*nu
                  iabcd=ib1+iacd
                  ibadc=b  +iadc
                  X=V(IABCD)
                  V(IABCD)=V(IBADC)
                  V(IBADC)=X
                  if(AEB.OR.CED)goto 193
                  iabdc=ib1+iadc_
                  ibacd=b  +iacd_
                  x=v(iabdc)
                  v(iabdc)=v(ibacd)
                  v(ibacd)=x
 193           CONTINUE
 192        CONTINUE
 191     CONTINUE
 190  CONTINUE
      GO TO 1000
 20   CONTINUE
      DO 200 A=1,NU
         ia=a-1
         ia1=ia*nu
         ia2=ia*nu2
         ia3=ia*nu3
         DO 201 B=1,A
            AEB=A.EQ.B
            ib=b-1
            ib1=ib*nu
            ib2=ib*nu2
            ib3=ib*nu3
            DO 202 C=1,A
               AEC=A.EQ.C
               ic=c-1
               ic1=ic*nu
               ic2=ic*nu2
               ic3=ic*nu3
               iabc =a  +ib1+ic2
               ibca =b  +ic1+ia3
               icab =c  +ia2+ib3
               iabc_=ia1+ib2+ic3
               LD=B
               IF (A.EQ.B)Ld=C
               IF(AEC)GOTO 204
               iacb=a  +ic2+ib3
               icba=ic1+ib2+ia3
               icba_=c +ib1+ia2
               ibac =b +ia1+ic3
 204           CONTINUE
               id1=-nu
               id2=-nu2
               id3=-nu3
               DO 203 D=1,LD
                  BED=B.EQ.D
                  CED=C.EQ.D
                  id1=id1+nu
                  id2=id2+nu2
                  id3=id3+nu3
                  iabcd=id3+iabc
                  ibcda=id2+ibca
                  icdab=id1+icab
                  idabc=d  +iabc_
                  X=V(IABCD)
                  V(IABCD)=v(IBCDA)
                  V(IBCDA)=V(ICDAB)
                  V(ICDAB)=V(IDABC)
                  V(IDABC)=X
                  IF (BED.OR.AEB.AND.CED.OR.AEC)GO TO 203
                  iadcb=id1+iacb
                  idcba=d  +icba
                  icbad=id3+icba_
                  ibadc=id2+ibac
                  X=V(IADCB)
                  V(IADCB)=v(IDCBA)
                  V(IDCBA)=v(ICBAD)
                  V(ICBAD)=v(IBADC)
                  V(IBADC)=X
 203           CONTINUE
 202        CONTINUE
 201     CONTINUE
 200  CONTINUE
      GO TO 1000
 21   CONTINUE
      DO 210 A=1,NU
         ia=a-1
         ia1=ia*nu
         ia2=ia*nu2
         ia3=ia*nu3
         DO 211 B=1,A
            AEB=A.EQ.B
            ib=b-1
            ib1=ib*nu
            ib2=ib*nu2
            ib3=ib*nu3
            DO 212 C=1,A
               AEC=A.EQ.C
               ic=c-1
               ic1=ic*nu
               ic2=ic*nu2
               ic3=ic*nu3
               iabc =a  +ib1+ic2
               iabc_=ia1+ib2+ic3
               icab =c  +ia2+ib3
               ibca =b  +ic1+ia3
               IF (AEC)GO TO 214
               iacb =a  +ic2+ib3
               ibac =b  +ia1+ic3
               icba =c  +ib1+ia2
               icba_=ic1+ib2+ia3
 214           CONTINUE
               LD=B
               IF (AEB)LD=C
               id1=-nu
               id2=-nu2
               id3=-nu3
               DO 213 D=1,LD
                  BED=B.EQ.D
                  CED=C.EQ.D
                  id1=id1+nu
                  id2=id2+nu2
                  id3=id3+nu3
c                  id=d-1
c                  id1=id*nu
c                  id2=id*nu2
c                  id3=id*nu3
                  iabcd=id3+iabc
                  idabc=d  +iabc_
                  icdab=id1+icab
                  ibcda=id2+ibca
                  X=V(IABCD)
                  V(IABCD)=V(IDABC)
                  V(IDABC)=V(ICDAB)
                  V(ICDAB)=V(IBCDA)
                  V(IBCDA)=X
                  IF (BED.OR.AEB.AND.CED.OR.AEC)GO TO 213
                  iadcb=id1+iacb
                  ibadc=id2+ibac
                  icbad=id3+icba
                  idcba=d  +icba_
                  X=V(IADCB)
                  V(IADCB)=V(IBADC)
                  V(IBADC)=V(ICBAD)
                  V(ICBAD)=V(IDCBA)
                  V(IDCBA)=X
 213           CONTINUE
 212        CONTINUE
 211     CONTINUE
 210  CONTINUE
      GO TO 1000
 22   CONTINUE
      DO 220 A=1,NU
         ia=a-1
         ia1=ia*nu
         ia2=ia*nu2
         ia3=ia*nu3
         DO 221 B=1,A
            AEB=A.EQ.B
            ib=b-1
            ib1=ib*nu
            ib2=ib*nu2
            ib3=ib*nu3
            DO 222 D=1,A
               AED=A.EQ.D
               id=d-1
               id1=id*nu
               id2=id*nu2
               id3=id*nu3
               iabd=a  +ib1+id3
               ibda=b  +id1+ia2
               idba=d  +ib2+ia3
               iadb=ia1+id2+ib3
               IF (AED)GO TO 224
               iabd_=a  +ib2+id3
               idab =id1+ia2+ib3
               idba_=d  +ib1+ia3
               ibad =b  +ia1+id2
 224           CONTINUE
               LC=B
               IF (AEB)LC=D
               ic1=-nu
               ic2=-nu2
               ic3=-nu3
               DO 223 C=1,LC
                  BEC=B.EQ.C
                  CED=C.EQ.D
                  ic1=ic1+nu
                  ic2=ic2+nu2
                  ic3=ic3+nu3
c                  ic=c-1
c                  ic1=ic*nu
c                  ic2=ic*nu2
c                  ic3=ic*nu3
                  iabcd=ic2+iabd
                  ibdac=ic3+ibda
                  idcba=ic1+idba
                  icadb=c  +iadb
                  X=V(IABCD)
                  V(IABCD)=V(IBDAC)
                  V(IBDAC)=V(IDCBA)
                  V(IDCBA)=V(ICADB)
                  V(ICADB)=X
                  IF (BEC.OR.AEB.AND.CED.OR.AED)GO TO 223
                  iacbd=ic1+iabd_
                  icdab=c  +idab
                  idbca=ic2+idba_
                  ibadc=ic3+ibad
                  X=V(IACBD)
                  V(IACBD)=V(ICDAB)
                  V(ICDAB)=V(IDBCA)
                  V(IDBCA)=V(IBADC)
                  V(IBADC)=X
 223           CONTINUE
 222        CONTINUE
 221     CONTINUE
 220  CONTINUE
      GO TO 1000
 23   CONTINUE
      DO 230 A=1,NU
         ia=a-1
         ia1=ia*nu
         ia2=ia*nu2
         ia3=ia*nu3
         DO 231 B=1,A
            AEB=A.EQ.B
            ib=b-1
            ib1=ib*nu
            ib2=ib*nu2
            ib3=ib*nu3
            DO 232 D=1,A
               AED=A.EQ.D
               id=d-1
               id1=id*nu
               id2=id*nu2
               id3=id*nu3
               iabd=a  +ib1+id3
               iadb=ia1+id2+ib3
               idba=d  +ib2+ia3
               ibda=b  +id1+ia2
               IF (AED)GO TO 234
               iabd_=a  +ib2+id3
               ibad =b  +ia1+id2
               idba_=d  +ib1+ia3
               idab =id1+ia2+ib3
 234           CONTINUE
               LC=B
               IF (AEB)LC=D
               ic1=-nu
               ic2=-nu2
               ic3=-nu3
               DO 233 C=1,LC
                  BEC=B.EQ.C
                  CED=C.EQ.D
                  ic1=ic1+nu
                  ic2=ic2+nu2
                  ic3=ic3+nu3
                  iabcd=ic2+iabd
                  icadb=c  +iadb
                  idcba=ic1+idba
                  ibdac=ic3+ibda
                  X=V(IABCD)
                  V(IABCD)=V(ICADB)
                  V(ICADB)=V(IDCBA)
                  V(IDCBA)=V(IBDAC)
                  V(IBDAC)=X
                  IF (BEC.OR.AEB.AND.CED.OR.AED)GO TO 233
                  iacbd=ic1+iabd_
                  ibadc=ic3+ibad
                  idbca=ic2+idba_
                  icdab=c  +idab
                  X=V(IACBD)
                  V(IACBD)=V(IBADC)
                  V(IBADC)=V(IDBCA)
                  V(IDBCA)=V(ICDAB)
                  V(ICDAB)=X
 233           CONTINUE
 232        CONTINUE
 231     CONTINUE
 230  CONTINUE
      GO TO 1000
 1000 continue
      RETURN
      END
      SUBROUTINE TRANT3(V,NU,ID)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION V(NU,NU,NU)
      GO TO (1,2,3,4,5),ID
    1 CONTINUE
      DO 100 A=1,NU
      DO 100 B=1,NU
      DO 100 C=1,B
      X=V(A,B,C)
      V(A,B,C)=V(A,C,B)
      V(A,C,B)=X
  100 CONTINUE
      GO TO 1000
    2 CONTINUE
      DO 200 A=1,NU
      DO 200 B=1,A
      DO 200 C=1,NU
      X=V(A,B,C)
      V(A,B,C)=V(B,A,C)
      V(B,A,C)=X
  200 CONTINUE
      GO TO 1000
    3 CONTINUE
      DO 300 A=1,NU
      DO 300 B=1,NU
      DO 300 C=1,A
      X=V(A,B,C)
      V(A,B,C)=V(C,B,A)
      V(C,B,A)=X
  300 CONTINUE
      GO TO 1000
    4 CONTINUE
      DO 400 B=1,NU
      DO 400 C=1,B
      DO 400 A=1,C
      X=V(A,B,C)
      V(A,B,C)=V(B,C,A)
      V(B,C,A)=V(C,A,B)
      V(C,A,B)=X
      IF(B.EQ.C.OR.C.EQ.A)GO TO 400
      X=V(B,A,C)
      V(B,A,C)=V(A,C,B)
      V(A,C,B)=V(C,B,A)
      V(C,B,A)=X
  400 CONTINUE
      GO TO 1000
    5 CONTINUE
      DO 500 A=1,NU
      DO 500 C=1,A
      DO 500 D=1,C
      X=V(C,D,A)
      V(C,D,A)=V(A,C,D)
      V(A,C,D)=V(D,A,C)
      V(D,A,C)=X
      IF (A.EQ.C.OR.C.EQ.D)GO TO 500
      X=V(A,D,C)
      V(A,D,C)=V(C,A,D)
      V(C,A,D)=V(D,C,A)
      V(D,C,A)=X
  500 CONTINUE
      GO TO 1000
 1000 CONTINUE
      RETURN
      END
