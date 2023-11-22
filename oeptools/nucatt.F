      subroutine nucatt (
     &     natm, nprm, maxshl, maxang, 
     &     nshlatm, nprmshl, angfct, alpha, 
     &     crdnuc, crdatm, attint)
c
c The subroutines evaluates nuclear attraction integrals in 
c primitive functions representation (<mu|1/(r-R)|nu>). It is
c Ajith's TWODENS/repprdint () routine modified for the INTGRT
c module. The repprdint () uses subroutines borrowed from the
c VMOL(?) property integral package.
c
c Igor Schweigert, Oct 2002 
c $Id: nucatt.FPP,v 1.1.1.1 2008/05/29 19:35:40 taube Exp $
c 
      implicit none
c
c     .. input parameters ..
c
      integer
     &     natm, nprm, maxshl, maxang, i,j,
     &     nshlatm (natm), nprmshl (maxshl, natm), angfct (maxang,3)
c
      double precision
     &     crdnuc (3), crdatm (3, natm),
     &     alpha (nprm), attint (nprm, nprm)
c
c     .. local variables ..
c
      integer
     &     iprm, iatm, ishl, jprm, jatm, jshl, k, 
     &     ip (20), ioffst (15)
c
      double precision
     &     crdatmi (3), crdatmj (3), fac (9,9)
c      
c     .. initialize the factorial array ..
c
      call setrhf(fac, ioffst, ip)
c
c     .. clear the nuclear attraction integrals array ..
c
      call dzero(attint, nprm*nprm)
c
c     .. the first loop over atoms and angular momentum shells  ..
c
      iprm = 1
      do iatm = 1, natm
         do ishl = 1, nshlatm (iatm)
c
c           .. the second loop over atoms and angular momentum shells  ..
c
            jprm = 1
            do jatm = 1, natm
               do jshl = 1, nshlatm (jatm)
c                  
c                 .. coordinates of atoms ..
c
                  do k = 1, 3
                     crdatmi(k)  = crdatm(k, iatm)
                     crdatmj(k)  = crdatm(k, jatm)
                  enddo
c
c                 .. evaluate the nuc.att.int. for all the primitives in the shell ..
c                  
                  call evlnucatt (
     &                 iprm, nprmshl (ishl, iatm), 
     &                 jprm, nprmshl (jshl, jatm), nprm, 
     &                 angfct(ishl,1), angfct(ishl,2), angfct(ishl,3),
     &                 angfct(jshl,1), angfct(jshl,2), angfct(jshl,3),  
     &                 crdatmi, crdatmj, crdnuc,
     &                 alpha, fac, attint)
c
                  jprm = jprm + nprmshl (jshl, jatm)
               enddo
            enddo
c           
            iprm = iprm + nprmshl (ishl, iatm)
         enddo
      enddo
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c
      subroutine evlnucatt (
     &     iprm, niprm, jprm, njprm, nprm, li, mi, ni, lj, mj, nj,
     &     crd1, crd2, crd3, alpha, fac, attint)
c                                                                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c This routine evaluates the nuclear attraction intergrals for the 
c two sets of primitives. It's based on the Ajith's 
c TWODENS/evalrepuli () routine. IVS 10/18/01
c
      implicit none
c
c .. the input parameters ..
c
      integer
     &     iprm, niprm, jprm, njprm, nprm, 
     &     li, mi, ni, lj, mj, nj
c
      double precision
     &     crd1 (3), crd2 (3), crd3 (3), 
     &     alpha (nprm), attint (nprm, nprm), fac (9,9)
c
c .. local variables ..
c
      integer
     &     i, j, k, itm, jtm, iiq, jjq, il, im, in, jl, jm, jn,
     &     ibtand, ibtshr, 
     &     ilmn (27), jlmn (27)
c
      double precision
     &     sum, v, gama, efac, crd12 (3), 
     &     aaa (27), bbb (27)
c
c     .. bitwise operators ..
c
      ibtand(i,j) = iand(i,j)
      ibtshr(i,j) = ishft(i,-j)
c
c     .. loop over all the primitives in the ith and jth shells
c
      do i = iprm, iprm + niprm - 1
         do j = jprm, jprm + njprm - 1 
c
c           .. get the common center of the gaussians ..
c
            call civpt(
     &           crd1, alpha (i), crd2, alpha (j),
     &           crd12, gama, efac)
c
c           .. initialize the arrays needed for the integral evaluation ..
c
            call rhftce(aaa, crd1, crd12, fac, li, mi, ni, itm, ilmn)
            call rhftce(bbb, crd2, crd12, fac, lj, mj, nj, jtm, jlmn)
c
c           .. distance from the gaussians center and the nucleus ..
c
            do k=1,3
               crd12(k) =crd12 (k) - crd3 (k) 
            enddo
c
c           .. if you really want to understand what's going on in the 
c           .. two loops below, go to ~rjbprogs/Src/vprops and gods have
c           .. mercy for your soul ..
c
            sum = 0.0d0
            do iiq = 1, itm
               
               if (dabs (aaa (iiq)) .lt. 1.d-8 ) goto 101
c     
               il = ibtand(ibtshr(ilmn(iiq),20),2**10 - 1)
               im = ibtand(ibtshr(ilmn(iiq),10),2**10 - 1)
               in = ibtand(ilmn(iiq),2**10 - 1)
c     
               do jjq = 1, jtm
c                                           .. attention here! ..
                  if ( dabs (aaa (iiq) * bbb (jjq)) .lt. 1.d-8 )
     &                 goto 102
                  
                  jl = ibtand(ibtshr(jlmn(jjq),20),2**10-1) + il
                  jm = ibtand(ibtshr(jlmn(jjq),10),2**10-1) + im
                  jn = ibtand(jlmn(jjq),2**10 - 1) + in
                  
                  call evalreplp(jl, jm, jn, gama, v, crd12)
                  
                  sum = aaa(iiq)*bbb(jjq)*v + sum
 102              continue
               enddo
 101           continue 
            enddo
c
c           .. store the result in the global array ..
c     
            attint (i, j) = sum*efac
         enddo
      enddo
c
      return 
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c
      SUBROUTINE EVALREPLP(L,M,N,GA,V,D)
c                                                                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PRTBUG/ IPRTFS
C
      DIMENSION V(10),D(10)
C
      COMMON /VARS/ FNU(20),DP(20),FN(15),FD(15),PI,DFTR(11)
C
      DATA FQR/1.3333333333333D0/
      DATA SQRPI/1.772453851D0/
C
C Set the factorial arrays
C      
      PI = DATAN(1.0D0+00)*4.0D+00
      FN(1) = 1.0D0
      FD(1) = 1.0D0
      DO I = 2, 15
         FN(I) = (I - 1)*FN(I - 1)
         FD(I) = 1.0D0/FN(I)
      ENDDO
C
      NT=1
      ITOT=L+M+N
      ITOTH=ITOT/2
      PRE=(2.*PI/GA)*FN(L+1)*FN(M+1)*FN(N+1)
C
      IF (2*ITOTH - ITOT) 1, 2, 1
C
 1    PRE=-PRE
 2    DEL=.25/GA
      X=GA*(D(1)**2+D(2)**2+D(3)**2)
      XX = X + X
C
      CALL FMC(12,X,EXPMX,FMCH)
C
      FNU(13)=FMCH
      FNU(12)=(EXPMX+XX*FNU(13))/23.
      FNU(11)=(EXPMX+XX*FNU(12))/21.
      FNU(10)=(EXPMX+XX*FNU(11))/19.
      FNU(9)=(EXPMX+XX*FNU(10))/17.
      FNU(8)=(EXPMX+XX*FNU(9))/15.
      FNU(7)=(EXPMX+XX*FNU(8))/13.
      FNU(6)=(EXPMX+XX*FNU(7))/11.
      FNU(5)=(EXPMX+XX*FNU(6))/9.
      FNU(4)=(EXPMX+XX*FNU(5))/7.
      FNU(3)=(EXPMX+XX*FNU(4))/5.
      FNU(2)=(EXPMX+XX*FNU(3))/3.
      FNU(1)=EXPMX+XX*FNU(2)
C
      DP(1)=1.0D0
      DP(2)=DEL
      DP(3)=DEL**2
      DP(4)=DP(3)*DEL
      DP(5)=DP(4)*DEL
      DP(6)=DP(5)*DEL
      DP(7)=DP(6)*DEL
      V(1) =PRE*AAINER(0,0,0,L,M,N,D)
C$$      WRITE(6,*) V(1)
C
      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c
      SUBROUTINE SETRHF(FAC, IOFFST, IP)
c                                                                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DIMENSION IP(20), FAC(9, 9), IOFFST(15)
C
      IBTAND(I,J) = IAND(I,J)
      IBTOR(I,J)  = IOR(I,J)
      IBTXOR(I,J) = IEOR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTNOT(I)   = NOT(I)
c
      IOFFST(1) = 0
      IOFFST(3) = 1
      IOFFST(6) = 4 
      IOFFST(10) = 10
      IOFFST(15) = 20
C
c      CALL SETLMN
C
      FAC(1,1)=1.0
      FAC(2,1)=1.0
      FAC(2,2)=1.0
C
      DO  I=3,9
         FAC(I,1)=1.0
         FAC(I,I)=1.0
         JE=I-1
         DO J=2,JE
            FAC(I,J)=FAC(I-1,J-1)+FAC(I-1,J)
         ENDDO
      ENDDO
C
      IP(1)=0
      IP(2)=IBTSHL(1,20)
      IP(3)=IBTSHL(1,10)
      IP(4)=1
      IP(5)=IBTSHL(2,20)
      IP(6)=IBTSHL(2,10)
      IP(7)=2
      IP(8)=IBTOR(IBTSHL(1,20),IBTSHL(1,10))
      IP(9)=IBTOR(IBTSHL(1,20),1)
      IP(10)=IBTOR(IBTSHL(1,10),1)
      IP(11)=IBTSHL(3,20)
      IP(12)=IBTSHL(3,10)
      IP(13)=3
      IP(14)=IBTOR(IBTSHL(2,20),IBTSHL(1,10))
      IP(15)=IBTOR(IBTSHL(2,20),1)
      IP(16)=IBTOR(IBTSHL(1,20),IBTSHL(2,10))
      IP(17)=IBTOR(IBTSHL(2,10),1)
      IP(18)=IBTOR(IBTSHL(1,20),2)
      IP(19)=IBTOR(IBTSHL(1,10),2)
      IP(20)=IBTOR(IBTSHL(1,20),IBTOR(IBTSHL(1,10),1))
C
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c
      SUBROUTINE CIVPT(A,ALPHA,B,BETA,ECOORD,EEXP,EFACT)
c                                                                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
C     THIS COMPUTES THE CENTRE, EXPONENT, AND MULTIPLYING FACTOR OF
C     A SINGLE GAUSSIAN WHICH CAN REPLACE THE PRODUCT OF TWO GAUSSIA
C     CENTRES A AND B, AND EXPONENTS ALPHA AND BETA.
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(3), B(3), ECOORD(3)
      EEXP=ALPHA+BETA
      T = 1.0/EEXP
      DO 1 I=1,3
1     ECOORD(I) = ( ALPHA*A(I) + BETA*B(I) )*T
      ABSQD=0.0
      DO 2 I=1,3
2     ABSQD=ABSQD+(B(I)-A(I))*(B(I)-A(I))
      EFACT = EXP(-ALPHA*BETA*ABSQD*T)
      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c
      SUBROUTINE RHFTCE(C, A, E, FAC, LL, MM, NN, ITM, LMN)
c                                                                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
c The subroutine is slightly modified since in the INTGRT module
c the gaussian powers are given explicitly. All the lines containing
c LMNVAL are commented out. IVS 10/18/02
c
C$$$      SAVE IP

      DIMENSION LMN(27)
      DIMENSION FAC(9,9),XH(9),YH(9),ZH(9)
c
c      COMMON /HIGHL/ LMNVAL(3,84), ANORM(84)
c
C
      DIMENSION C(27)
      DIMENSION A(3),E(3)
C
      IBTAND(I,J) = IAND(I,J)
      IBTOR(I,J)  = IOR(I,J)
      IBTXOR(I,J) = IEOR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTNOT(I)   = NOT(I)
C
      AEX = E(1)-A(1)
      AEY = E(2)-A(2)
      AEZ = E(3)-A(3)
C
      Q = 1.00D0
      X = 1.00D0
C
c     L = LMNVAL(1,ITYP)+1
c     M = LMNVAL(2,ITYP)+1
c     N = LMNVAL(3,ITYP)+1
c
      l=ll+1
      m=mm+1
      n=nn+1
C
      DO L1=1,L
         XH(L1)=X
         X=X*AEX
      ENDDO
C      
      Y=1.0D0
      DO  M1=1,M
         YH(M1)=Y
         Y=Y*AEY
      ENDDO
C
      Z=1.0D0
      DO N1=1,N
         ZH(N1)=Z
         Z=Z*AEZ
      ENDDO
C
      ITM=0
C
      DO L1=1,L
C
         X=FAC(L,L1)*XH(L1)*Q
C
         DO  M1=1,M
C
           Y=FAC(M,M1)*YH(M1)*X
C
            DO N1=1,N
C
               Z=FAC(N,N1)*ZH(N1)*Y
               ITM=ITM+1
C
               LMN(ITM)=IBTOR(IBTSHL((L-L1),20),
     &                  IBTOR(IBTSHL((M-M1),10),(N-N1)))
                 C(ITM)=Z
C
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c
      SUBROUTINE FMC(MVAR,XVAR,EXPMX,FMCH)
c                                                                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      M=MVAR
      X=XVAR
      EXPMX= EXP(-X)
      IF(X - 20.0) 10,10,20
C     X LESS THAN OR EQUAL TO 20.0
10    A=M
      A = A+0.5
      TERM=1.0/A
      PTLSUM=TERM
      DO 11 I=2,200
      A=A+1.0
      TERM = TERM*X/A
      PTLSUM=PTLSUM+TERM
      IF(TERM/PTLSUM-1.E-10) 12, 11, 11
11    CONTINUE
      WRITE(6,999) M, X
12    FMCH=0.5*PTLSUM*EXPMX
      RETURN
C     X GREATER THAN 20.0.
20    A=M
      XD=1.0/X
      B=A+0.5
      A=A-0.5
      APPROX = 0.88622692*(SQRT(XD)*XD**M)
      IF (M) 21, 23, 21
21    DO 22 I=1,M
      B=B-1.0
22    APPROX=APPROX*B
23    FIMULT=0.5*EXPMX*XD
      FIPROP=FIMULT/APPROX
      TERM=1.0
      PTLSUM=TERM
      NOTRMS=X
      NOTRMS=NOTRMS+M
      DO 24 I=2,NOTRMS
      TERM=TERM*A*XD
      PTLSUM=PTLSUM+TERM
      IF ( ABS(TERM*FIPROP/PTLSUM) - 1.E-10) 25,25,24
24    A=A-1.0
      WRITE(6,999) M, X
25    FMCH=APPROX-FIMULT*PTLSUM
      RETURN
999   FORMAT (24H0NO CONVERGENCE FOR FMCH, I6, E16.9)
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c
      FUNCTION AAINER(R,S,T,L,M,N,D)
c                                                                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION D(*)
      COMMON /VARS/FNU(20),DP(20),FN(15),FD(15),PI,DFTR(11)
      INTEGER R,S,T,U,V,W,UVWT,RSTT,UVWTH
      RSTT=R+S+T
      LMNT=L+M+N
      LMNRST=RSTT+LMNT
      LH=L/2
      MH=M/2
      NH=N/2
      AAINER=0.0
      LAST1 = LH+1
      LAST2 = MH + 1
      LAST3 = NH + 1
      DO 6 II=1,LAST1
      I = II-1
      IL=L-2*I+1
      ILR=R+IL
      ILRH=(ILR-1)/2
      FI=FN(ILR)*FD(IL)*FD(I+1)
      LAST4 = ILRH + 1
      DO 5 JJ=1,LAST2
      J = JJ-1
      JM=M-2*J+1
      JMS=S+JM
      JMSH=(JMS-1)/2
      FIJ=FN(JMS)*FD(JM)*FD(J+1)*FI
      LAST5 = JMSH + 1
      DO 4 KK=1,LAST3
      K = KK-1
      KN=N-2*K+1
      KNT=T+KN
      KNTH=(KNT-1)/2
      IJKT=I+J+K
      FIJK=FN(KNT)*FD(KN)*FD(K+1)*DP(IJKT+1)*FIJ
      LMRSIJ=LMNRST-2*IJKT
      LAST6 = KNTH + 1
      DO 3 IU=1,LAST4
      U = IU-1
      ILRU=ILR-2*U
      FU=FD(U+1)*FD(ILRU)
      IF( ABS(D(1))-1.E-10) 10,10,11
10    IF(ILRU-1) 12,12, 3
11    FU=FU*D(1)**(ILRU-1)
12    CONTINUE
      DO 2 IV=1,LAST5
      V = IV-1
      JMSV=JMS-2*V
      FUV=FU*FD(V+1)*FD(JMSV)
      IF( ABS(D(2))-1.E-10) 20,20,21
20    IF(JMSV-1) 22,22,2
21    FUV=FUV*D(2)**(JMSV-1)
22    CONTINUE
      DO 1 IW=1,LAST6
      W = IW-1
      KNTW=KNT-2*W
      FUVW=FUV*FD(W+1) *FD(KNTW)
      IF( ABS(D(3))-1.E-10) 30,30,31
30    IF(KNTW-1) 32,32,1
31    FUVW=FUVW*D(3)**(KNTW-1)
32    UVWT=U+V+W
      UVWTH=UVWT/2
      IF(2*UVWTH-UVWT) 33,34,33
33    FUVW=-FUVW
34    NUINDX=LMRSIJ-UVWT
      FUVW=FUVW*FNU(NUINDX  +1)*DP(UVWT+1)
      AAINER=FIJK*FUVW+AAINER
1     CONTINUE
2     CONTINUE
3     CONTINUE
4     CONTINUE
5     CONTINUE
6     CONTINUE
      RETURN
      END
