      SUBROUTINE CPHFSX(A,BUF,IBUF,IVO,IVOS,ISYMO,NINTMX,NIREP
     X ,NSIZ1,NSIZO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 LABEL
      DIMENSION LABEL(52)
      DIMENSION IVOS(NSIZ1,NSIZO,NIREP),IVO(NSIZ1,NSIZO)
      DIMENSION BUF(NINTMX),IBUF(NINTMX),ISYMO(1),A(1)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/INFOA/NBASIS,NUMSCF,NX,NMO2,NOC,NVT,NVO     
      COMMON/INFSYM/NSYMHF,NSO(8),NOCS(8),NVTS(8),IDPR(8,8),NVOS(8)
     X ,NIJS(8),NIJS2(8),NRDS(8),NIJSR(8)     
      DATA ZERO,TWO/0.D0,2.D0/ 
      WRITE(6,*) ' We are in CPHFSX '
C     write(6,*) ' IVO '
      do 100 i=1,nsiz1
C     write(6,*) i,'   ',(IVO(I,J),J=1,NSIZO)
  100 continue
C     write(6,*) ' (0) IVO(6,1)= ',IVO(6,1)
      NTAPE=10
      NVOS2=0
      DO 1 IREP=1,NIREP
      NVOS2=NVOS2+NVOS(IREP)*NVOS(IREP)
C     write(6,*) ' IVOS for IREP = ',IREP,NVOS(IREP)
C     do 101 i=1,nsiz1
C     write(6,*) (IVOS(I,J,IREP),J=1,NSIZO)
C 101 continue
    1 CONTINUE
      DO 2 IJ=1,NVOS2
    2 A(IJ)=ZERO
      OPEN(UNIT=NTAPE,FORM='UNFORMATTED',FILE='HF2',ACCESS='SEQUENTIAL')
      READ(NTAPE) LABEL(1)
      READ(NTAPE) LABEL(1)
      READ(NTAPE) LABEL(1)
      READ(NTAPE) LABEL(1)
C     WRITE(6,*) LABEL
C     WRITE(6,*) ' # , INTEGRAL I , J , K , L '
      NTOT=0
    5 READ(NTAPE,END=50) BUF,IBUF,NUT
      IF(NUT.LE.0) GO TO 30
      DO 10 IX=1,NUT
      NTOT = NTOT + 1
      L = AND(IBUF(IX),IALONE)
      K = AND(ISHFT(IBUF(IX),-IBITWD),IALONE)
      J = AND(ISHFT(IBUF(IX),-2*IBITWD),IALONE)
      I = AND(ISHFT(IBUF(IX),-3*IBITWD),IALONE)
      IJS = IDPR(ISYMO(I),ISYMO(J))
      KLS = IDPR(ISYMO(K),ISYMO(L)) 
      IKS = IDPR(ISYMO(I),ISYMO(K))
      JLS = IDPR(ISYMO(J),ISYMO(L))
      ILS = IDPR(ISYMO(I),ISYMO(L))
      JKS = IDPR(ISYMO(J),ISYMO(K))
      VAL = BUF(IX)
      VAL2 = VAL + VAL
      VAL2 = VAL2 + VAL2
      VALH = VAL
C     write(6,*) ' IVO(6,1) ',IVO(6,1),' when ',I,J,K,L
C ..... VAL2  is twiced by occupation # = 2 for close shell ....
C .... IVO is used only for feeding proper (IJ,KL) values for VAL in cases
C .... 1/2 when I=J and/or K=L ,  1/2 when (I,J)=(K,L) and I><J, K><L
C .... The operation could be done earlier using only I,J,K,L values
C ... here we took advantage to compare (I,J) and (K,L) without
C ... lengthy IF statements
      IF(I.GT.NOC.AND.J.LE.NOC.AND.K.GT.NOC.AND.L.LE.NOC) THEN
      N1=IVO(I,J)
C     write(6,*) ' I,J,IVO(I,J) ',I,J,IVO(I,J)
      N2=IVO(K,L)
C ..... those N1 and N2 are different from N1 and N2 determined later
C ..... by IVOS and used only for the purpose above.
      IF(N1.EQ.N2) THEN
      VAL2= VAL2/TWO
      VALH= VALH/TWO
      END IF
      N3=IVO(I,L)
      N4=IVO(K,J)
C     write(6,*) I,J,K,L,VAL,N1,N2,N3,N4
      IF(IJS.EQ.KLS) THEN
      N1=IVOS(I,J,IJS)
      N2=IVOS(K,L,KLS)
      N12=(N2-1)*NVOS(IJS)+N1+NIJS2(IJS)
      N21=(N1-1)*NVOS(IJS)+N2+NIJS2(IJS)
C     write(6,*) ' IJS = ',IJS,' N1,N2 ',N1,N2
      A(N12)=A(N12)+VAL2
      A(N21)=A(N21)+VAL2
      END IF
      IF(ILS.EQ.JKS) THEN
      N3=IVOS(I,L,ILS)
      N4=IVOS(K,J,JKS)
C     write(6,*) ' ILS = ',ILS,' N3,N4 ',N3,N4
      N34=(N4-1)*NVOS(ILS)+N3+NIJS2(ILS)
      N43=(N3-1)*NVOS(ILS)+N4+NIJS2(ILS)
      A(N34)=A(N34)-VALH
      A(N43)=A(N43)-VALH
      END IF
      END IF
      IF(I.GT.NOC.AND.J.LE.NOC.AND.K.LE.NOC.AND.L.GT.NOC) THEN
      N1=IVO(I,J)
      N2=IVO(L,K)
      IF(N1.EQ.N2) THEN
      VAL2= VAL2/TWO
      VALH= VALH/TWO
      END IF
      N3=IVO(I,K)
      N4=IVO(L,J)
C     write(6,*) I,J,K,L,VAL,N1,N2,N3,N4
      IF(IJS.EQ.KLS) THEN
      N1=IVOS(I,J,IJS)
      N2=IVOS(L,K,KLS)
      N12=(N2-1)*NVOS(IJS)+N1+NIJS2(IJS)
      N21=(N1-1)*NVOS(IJS)+N2+NIJS2(IJS)
C     write(6,*) ' IJS = ',IJS,' N1,N2 ',N1,N2
      A(N12)=A(N12)+VAL2
      A(N21)=A(N21)+VAL2
      END IF
      IF(IKS.EQ.JLS) THEN
      N3=IVOS(I,K,IKS)
      N4=IVOS(L,J,JLS)
      N34=(N4-1)*NVOS(IKS)+N3+NIJS2(IKS)
      N43=(N3-1)*NVOS(IKS)+N4+NIJS2(IKS)
C     write(6,*) ' IKS = ',IKS,' N3,N4 ',N3,N4
      A(N34)=A(N34)-VALH
      A(N43)=A(N43)-VALH
      END IF
      END IF
      IF(I.LE.NOC.AND.J.GT.NOC.AND.K.GT.NOC.AND.L.LE.NOC) THEN
      N1=IVO(J,I)
      N2=IVO(K,L)
      IF(N1.EQ.N2) THEN
      VAL2= VAL2/TWO
      VALH= VALH/TWO
      END IF
      N3=IVO(K,I)
      N4=IVO(J,L)
C     write(6,*) I,J,K,L,VAL,N1,N2,N3,N4
      IF(IJS.EQ.KLS) THEN
      N1=IVOS(J,I,IJS)
      N2=IVOS(K,L,KLS)
      N12=(N2-1)*NVOS(IJS)+N1+NIJS2(IJS)
      N21=(N1-1)*NVOS(IJS)+N2+NIJS2(IJS)
C     write(6,*) ' IJS = ',IJS,' N1,N2 ',N1,N2
      A(N12)=A(N12)+VAL2
      A(N21)=A(N21)+VAL2
      END IF
      IF(IKS.EQ.JLS) THEN
      N3=IVOS(K,I,IKS)
      N4=IVOS(J,L,JLS)
      N34=(N4-1)*NVOS(IKS)+N3+NIJS2(IKS)
      N43=(N3-1)*NVOS(IKS)+N4+NIJS2(IKS)
C     write(6,*) ' IKS = ',IKS,' N3,N4 ',N3,N4
      A(N34)=A(N34)-VALH
      A(N43)=A(N43)-VALH
      END IF
      END IF 
      IF(I.LE.NOC.AND.J.GT.NOC.AND.K.LE.NOC.AND.L.GT.NOC) THEN
      N1=IVO(J,I)
      N2=IVO(L,K)
      IF(N1.EQ.N2) THEN
      VAL2= VAL2/TWO
      VALH= VALH/TWO
      END IF
      N3=IVO(L,I)
      N4=IVO(J,K)
C     write(6,*) I,J,K,L,VAL,N1,N2,N3,N4
      IF(IJS.EQ.KLS) THEN
      N1=IVOS(J,I,IJS)
      N2=IVOS(L,K,KLS)
      N12=(N2-1)*NVOS(IJS)+N1+NIJS2(IJS)
      N21=(N1-1)*NVOS(IJS)+N2+NIJS2(IJS)
C     write(6,*) ' IJS = ',IJS,' N1,N2 ',N1,N2
      A(N12)=A(N12)+VAL2
      A(N21)=A(N21)+VAL2
      END IF
      IF(ILS.EQ.JKS) THEN
      N3=IVOS(L,I,ILS)
      N4=IVOS(J,K,JKS)
      N34=(N4-1)*NVOS(ILS)+N3+NIJS2(ILS)
      N43=(N3-1)*NVOS(ILS)+N4+NIJS2(ILS)
C     write(6,*) ' ILS = ',ILS,' N3,N4 ',N3,N4
      A(N34)=A(N34)-VALH
      A(N43)=A(N43)-VALH
      END IF
      END IF 
      IF(I.GT.NOC.AND.J.GT.NOC.AND.K.LE.NOC.AND.L.LE.NOC) THEN
      N1=IVO(I,K)
      N2=IVO(J,L)
      N3=IVO(I,L)
      N4=IVO(J,K)
      IF(I.EQ.J) VALH=VALH/TWO
      IF(K.EQ.L) VALH=VALH/TWO
C     write(6,*) I,J,K,L,VAL,N1,N2,N3,N4
      IF(IKS.EQ.JLS) THEN
      N1=IVOS(I,K,IKS)
      N2=IVOS(J,L,JLS) 
      N12=(N2-1)*NVOS(IKS)+N1+NIJS2(IKS)
      N21=(N1-1)*NVOS(IKS)+N2+NIJS2(IKS)
C     write(6,*) ' IKS = ',IKS,' N1,N2 ',N1,N2
      A(N12)=A(N12)-VALH 
      A(N21)=A(N21)-VALH 
      END IF
      IF(ILS.EQ.JKS) THEN
      N3=IVOS(I,L,ILS)
      N4=IVOS(J,K,JKS)
      N34=(N4-1)*NVOS(ILS)+N3+NIJS2(ILS)
      N43=(N3-1)*NVOS(ILS)+N4+NIJS2(ILS)
C     write(6,*) ' ILS = ',ILS,' N3,N4 ',N3,N4
      A(N34)=A(N34)-VALH
      A(N43)=A(N43)-VALH
      END IF
      END IF 
      IF(I.LE.NOC.AND.J.LE.NOC.AND.K.GT.NOC.AND.L.GT.NOC) THEN
      N1=IVO(K,I)
      N2=IVO(L,J)
      N3=IVO(L,I)
      N4=IVO(K,J)
      IF(I.EQ.J) VALH=VALH/TWO
      IF(K.EQ.L) VALH=VALH/TWO
C     write(6,*) I,J,K,L,VAL,N1,N2,N3,N4
      IF(IKS.EQ.JLS) THEN
      N1=IVOS(K,I,IKS)
      N2=IVOS(L,J,JLS)
      N12=(N2-1)*NVOS(IKS)+N1+NIJS2(IKS)
      N21=(N1-1)*NVOS(IKS)+N2+NIJS2(IKS)
C     write(6,*) ' IKS = ',IKS,' N1,N2 ',N1,N2
      A(N12)=A(N12)-VALH 
      A(N21)=A(N21)-VALH 
      END IF
      IF(ILS.EQ.JKS) THEN
      N3=IVOS(L,I,ILS)
      N4=IVOS(K,J,JKS)
      N34=(N4-1)*NVOS(ILS)+N3+NIJS2(ILS)
      N43=(N3-1)*NVOS(ILS)+N4+NIJS2(ILS)
C     write(6,*) ' ILS = ',ILS,' N3,N4 ',N3,N4
      A(N34)=A(N34)-VALH
      A(N43)=A(N43)-VALH
      END IF
      END IF 
   10 CONTINUE
      GO TO  5
   30 CONTINUE
   50 CONTINUE 
C     WRITE(6,*) ' A+B matrix without denominator shift '
C     CALL OUTMXD(A,NSIZVO,NVO,NVO)          
      CLOSE(NTAPE)     
      RETURN 
      END
