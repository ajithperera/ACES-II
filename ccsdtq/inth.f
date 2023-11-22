      subroutine drinth(no,nu,o1)
      implicit double precision (a-h,o-z)
      logical print
      common/flags/iflags(100)
      dimension o1(1)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON/NEWOPT/NOPT(6)
      common/itrat/icycle,mx,icn
      common/totmem/mem 
      print=iflags(1).gt.10
      if(nopt(1).gt.3)then
      i1=1           !ti
      i2=i1+nu3      !t2pp
      i3=i2+no2u2    !t4
      i4=i3+nu4      !v
      it=i4+nu4
      if(icycle.eq.1.and.print)write(6,99)mem,it
      call inth2p4 (no,nu,o1,o1(i2),o1(i3),o1(i4))
      else
      i1=1           !ti
      i2=i1+nu3      !t2pp
      i3=i2+no2u2    !t4
      it=i3+nu4
      if(icycle.eq.1.and.print)write(6,99)mem,it
      call inth2p43(no,nu,o1(i2),o1(i3),o1)
      endif
      i1=1           !ti
      i2=i1+nou2     !t2pp
      i3=i2+no2u2    !t2vo
      i4=i3+no2u2    !voe
      i5=i4+no2u2    !vo
      i6=i5+no2u2    !pz
      i7=i6+no4      !vt
      it=i7+no2u2
      if(icycle.eq.1.and.print)write(6,98)mem,it
      call inth4p2(no,nu,o1,o1(i2),o1(i3),o1(i4),o1(i5),o1(i6),o1(i7))
 98   format('Space usage in inth4p2: available - ',i8,'   used - ',i8)
 99   format('Space usage in inth2p4: available - ',i8,'   used - ',i8)
      return
      end
      SUBROUTINE INTH2P4(NO,NU,TI,T2PP,T4,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SDTQ23,SDTQ4
      INTEGER A
      COMMON/ITRAT/IT,MAXIT,ICCNV
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON/NEWOPT/NOPT(6)
      DIMENSION TI(NU,NU,NU),T2PP(NU,NU,NO,NO),T4(NU,NU,NU,NU),
     *V(NU,NU,NU,NU)
      DATA HALF/0.5D+0/,TWO/2.0D+0/
      call ienter(9)
      SDTQ23=NOPT(1).EQ.4.OR.NOPT(1).EQ.5
      SDTQ4=NOPT(1).GE.6
      CALL RO2PPH(0,NO,NU,TI,T2PP)
      CALL TRANMD(T2PP,NU,NU,NO,NO,34)
      CALL VECMUL(T2PP,NO2U2,HALF)
      I10=0
      IF (SDTQ4)I10=10
      IF (SDTQ23)THEN
      CALL RDOV4A(0,NO,NU,TI,V)
      CALL MTRANS(V,NU,12)
      ELSE
      IF (NOPT(1).EQ.3)THEN
         do 1 a=1,nu
      CALL RDV43(a,0,NO,NU,V(1,1,1,a))
 1    continue
      ELSE
      CALL RDOV4(I10,NO,NU,TI,V)
      ENDIF
      ENDIF
      IF(I10.EQ.10)THEN
      CALL MTRANS(V,NU,12)
      ENDIF
      DO 210 I=1,NO
      DO 210 J=1,I
      kk=it2(i,j)
      IASV=NO3+(kk-1)*nu+1
      IF (SDTQ4.OR.SDTQ23) THEN
      CALL RDVA(IASV,NU,T4)
      CALL MTRANS(T4,NU,6)
      ELSE
      CALL ZEROMA(T4,1,NU4)
      ENDIF
      DO 190 A=1,NU
      CALL MATMUL(T2PP(1,1,I,J),V(1,1,1,A),T4(1,1,1,A),NU,NU2,NU,0,0)
      CALL TRANMD(T4(1,1,1,A),NU,NU,NU,1,12)
      CALL MATMUL(T2PP(1,1,J,I),V(1,1,1,A),T4(1,1,1,A),NU,NU2,NU,0,0)
 190  CONTINUE
      CALL MTRANS(T4,NU,4)
      CALL WRVA(IASV,NU,T4)
  210 CONTINUE
      CALL VECMUL(T2PP,NO2U2,TWO)
      call iexit(9)
      RETURN
      END
      SUBROUTINE INTH2P43(NO,NU,T2PP,T4,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SDTQ23,SDTQ4
      INTEGER A
      COMMON/ITRAT/IT,MAXIT,ICCNV
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON/NEWOPT/NOPT(6)
      DIMENSION T2PP(NU,NU,NO,NO),T4(NU,NU,NU,NU),V(1)
      DATA HALF/0.5D+0/,TWO/2.0D+0/
      call ienter(9)
      CALL RO2PPH(0,NO,NU,V,T2PP)
      CALL TRANMD(T2PP,NU,NU,NO,NO,34)
      CALL VECMUL(T2PP,NO2U2,HALF)
      DO 210 I=1,NO
      DO 210 J=1,i
      DO 190 A=1,NU
      CALL RDV43(A,0,NO,NU,V)
      CALL MATMUL(T2PP(1,1,I,J),V,T4(1,1,1,A),NU,NU2,NU,1,0)
      CALL TRANT3(T4(1,1,1,A),NU,2)
      CALL MATMUL(T2PP(1,1,J,I),V,T4(1,1,1,A),NU,NU2,NU,0,0)
 190  CONTINUE
      kk=it2(i,j)
      IASV=NO3+(kk-1)*nu+1
      CALL MTRANS(T4,NU,4)
      CALL WRVA(IASV,NU,T4)
  210 CONTINUE
      CALL VECMUL(T2PP,NO2U2,TWO)
      call iexit(9)
      RETURN
      END
      SUBROUTINE INTH2P4P(NO,NU,TI,T2PP,T4,VO,VOE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SDTQ23,SDTQ4
      COMMON/ITRAT/IT,MAXIT,ICCNV
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON/NEWOPT/NOPT(6)
      DIMENSION TI(NU,NU,NU),T2PP(NU,NU,NO,NO),T4(NU,NU,NU,NU),
     *VO(NO,NU,NU,NO),VOE(NO,NU,NU,NO)
      DATA HALF/0.5D+0/,TWO/2.0D+0/
c      return
      SDTQ23=NOPT(1).EQ.4.OR.NOPT(1).EQ.5
      SDTQ4=NOPT(1).GE.6
      CALL RO2PPH(0,NO,NU,TI,T2PP)
      CALL TRANMD(T2PP,NU,NU,NO,NO,34)
      CALL RO2HPP(2,NO,NU,TI,VO)
      CALL RO2HPP(1,NO,NU,TI,VOE)
      CALL TRANMD(VOE,NO,NU,NU,NO,23)
      DO 210 I=1,NO
      DO 210 J=1,NO
      IASV=NO3+(I-1)*NO*NU+(J-1)*NU+1
      CALL RDVA(IASV,NU,T4)
      call zeroma(t4,1,nu4)
      CALL MATMUL(T2PP(1,1,1,I),VO(1,1,1,J),T4,NU2,NU2,NO,0,1)
      CALL MTRANS(T4,NU,6)
      CALL MATMUL(T2PP(1,1,1,J),VO(1,1,1,I),T4,NU2,NU2,NO,0,1)
      CALL MTRANS(T4,NU,6)
      CALL MTRANS(T4,NU,7)
      CALL MATMUL(T2PP(1,1,1,I),VOE(1,1,1,J),T4,NU2,NU2,NO,0,1)
      CALL MTRANS(T4,NU,8)
      CALL MATMUL(T2PP(1,1,1,J),VOE(1,1,1,I),T4,NU2,NU2,NO,0,1)
      CALL MTRANS(T4,NU,8)
      CALL MTRANS(T4,NU,7)
      IASV=NO3+(I-1)*NO*NU+(J-1)*NU+1
      CALL WRVA(IASV,NU,T4)
  210 CONTINUE
      CALL VECMUL(T2PP,NO2U2,TWO)
      RETURN
      END
      SUBROUTINE INTH4P2(NO,NU,TI,T2PP,T2VO,VOE,VO,PZ,VT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SDTQ23,SDTQ4
      COMMON/ITRAT/IT,MAXIT,ICCNV
      COMMON/NVTAP/NTITER,NQ2,NVT
      COMMON/NEWOPT/NOPT(6)
      DIMENSION TI(NU,NU,NU),T2PP(NU,NU,NO,NO),VO(NU,NU,NO,NO),
     *PZ(NO,NO,NO,NO),VT(NU,NU,NO,NO),T2VO(NO,NU,NU,NO),VOE(1)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DATA HALF/0.5D+0/
      call ienter(10)
      SDTQ23=NOPT(1).EQ.4.OR.NOPT(1).EQ.5
      SDTQ4=NOPT(1).GE.6
      no4u2=no4*nu2
      I0=0
      CALL RO2PPH(0,NO,NU,TI,T2PP)
      CALL RO2HPP(0,NO,NU,TI,T2VO)
      CALL TRANMD(T2VO,NO,NU,NU,NO,23)
      I11=1
      IF(SDTQ4)I11=11
      IF (SDTQ23) THEN
      CALL RO2PPH(-1,NO,NU,TI,VO)
      ELSE 
      CALL RO2PPH(I11,NO,NU,TI,VO)
      ENDIF
      CALL TRANMD(VO,NU,NU,NO,NO,34)
      CALL TRANMD(VO,NU,NU,NO,NO,12)
      I12=2
      IF(SDTQ4)I12=12
      IF (SDTQ23)THEN
      CALL RO2PPH(-2,NO,NU,TI,VOE)
      ELSE
      CALL RO2PPH(I12,NO,NU,TI,VOE)
      ENDIF
      CALL TRANMD(VOE,NU,NU,NO,NO,12)
      CALL TRANMD(VOE,NU,NU,NO,NO,34)
      I11=1
      IF(SDTQ4)I11=11
      IF (SDTQ23)THEN
      CALL RDOV4A(1,NU,NO,TI,PZ)
      ELSE
      CALL RDOV4(I11,NU,NO,TI,PZ)
      ENDIF
      call mtrans(pz,no,1)
      call vecmul(pz,no4,half)
      kkk=0
      do 123 i=1,no
      do 123 j=1,no
      if (nopt(1).gt.3)then
      do 223 k=1,no
      kk=no2*(j-1)+no*(i-1)+k
      call rdgen(nvt,kk,nou2,vt(1,1,1,k))
 223  continue
      else
      call zeroma(vt,1,no2u2)
      endif
      call tranmd(vt,nu,nu,no,no,34)
      call insitu(nu,nu,no,no,ti,vt,13)
c
      CALL MATMUL(T2VO(1,1,1,i),VO(1,1,1,j),VT,NOU,NOU,NU,I0,0)
      call insitu(no,nu,nu,no,ti,vt,13)
      call tranmd(vt,nu,nu,no,no,34)
      CALL MATMUL(T2PP(1,1,J,I),VOE,VT,NU,NO2U,NU,0,0)
      CALL MATMUL(T2PP(1,1,1,I),PZ(1,1,1,J),VT,NU2,NO2,NO,0,1)
      call tranmd(vt,nu,nu,no,no,12)
      DO 123 K=1,NO
      kkk=no2*(j-1)+no*(i-1)+k
      CALL WRGEN(NVT,KKK,NOU2,VT(1,1,1,k))
 123  continue
      call symvt(no,nu,t2pp,t2vo)
      call iexit(10)
      RETURN
      END
