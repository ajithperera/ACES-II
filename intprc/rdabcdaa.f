      SUBROUTINE RDABCDAA(W,LENGTH,ITOPADR,ILOWADR,BUF,IBUF,
     &                   IRREPA,IPDSZ,IPDIS,ISYM,IPW,
     &                   NUMIRW,ILNBUF,ISPIN)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION W,BUF,PH,ONE,ONEM
      LOGICAL GRAD
      DIMENSION W(LENGTH),BUF(ILNBUF),IBUF(ILNBUF)
      DIMENSION IRREPA(*),IPDSZ(*),IPDIS(*),ISYM(*)
      DIMENSION IPW(*),NUMIRW(*)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON/INFO/ NOCCO(2),NVRTO(2)
C
      DATA ONE,ONEM /1.0D0,-1.0D0/
C
      IUPKI(INTX)=AND(INTX,IALONE)
      IUPKJ(INTX)=AND(ISHFT(INTX,-IBITWD),IALONE)
      IUPKK(INTX)=AND(ISHFT(INTX,-2*IBITWD),IALONE)
      IUPKL(INTX)=AND(ISHFT(INTX,-3*IBITWD),IALONE)
      INDXA(I,J)=I+((J-1)*(J-2))/2
C
C POSITION HF2 FILE AT FIRST INTEGRAL RECORD
C
      CALL SETHF2(ISPIN,BUF,IBUF,ILNBUF,1)
C
      NORBA  =NVRTO(1)
      NORBB  =NVRTO(2)
      NO   =NOCCO(ISPIN)
      NV   =NVRTO(ISPIN)
      IOFFSET=ILOWADR-1
      ITOP   =ITOPADR-IOFFSET
      ICOUNT =0
      CALL ZERO(W,ITOP)
C
C READ IN BATCHES OF INTEGRALS
C
1     READ(25,END=103)BUF,IBUF,NUT
      IF(NUT.EQ.-1)GOTO 103
      DO 11 IK=1,NUT
       I=IUPKJ(IBUF(IK))
       J=IUPKI(IBUF(IK))
       K=IUPKL(IBUF(IK))
       L=IUPKK(IBUF(IK))
C
C ASSURE THAT J>I, L>K.
C
       IF(I.GT.J)THEN
        IT=J
        J=I
        I=IT
       ENDIF
       IF(K.GT.L)THEN
        IT=L
        L=K
        K=IT
       ENDIF
       I=I-NO
       J=J-NO
       K=K-NO
       L=L-NO
       ITMP=K
       K=J
       J=ITMP
       ILOAD=MIN(I,J,K,L)
       IF(ILOAD.GT.0)THEN
        ICOUNT =ICOUNT+1
C
C WE HAVE AN ABCD INTEGRAL - CALCULATE ITS ABSOLUTE INDEX

        I1=INDXA(MIN(I,J),MAX(I,J))
        I2=INDXA(MIN(K,L),MAX(K,L))
        IF(I.NE.J.AND.K.NE.L)THEN
         PH=ONE
         IF(I.GT.J)PH=PH*ONEM
         IF(K.GT.L)PH=PH*ONEM
         IRREPIJ=DIRPRD(IRREPA(I+NO),IRREPA(J+NO))
         IOFF=IPW(IRREPIJ)
         IF(IOFF.LE.ITOPADR)THEN
          NN=NUMIRW(IRREPIJ+NIRREP)
          IADR1=IOFF+ISYM(I1)-IPDSZ(IRREPIJ)
     &              +NN*(ISYM(I2)-1-IPDIS(IRREPIJ))-IOFFSET
          IADR2=IOFF+ISYM(I2)-IPDSZ(IRREPIJ)
     &             +NN*(ISYM(I1)-1-IPDIS(IRREPIJ)) -IOFFSET
          IF(IADR1.LE.ITOP.AND.IADR1.GT.0)W(IADR1)=W(IADR1)+PH*BUF(IK)
          IF(IADR2.LE.ITOP.AND.IADR2.GT.0.AND.IADR1.NE.IADR2)THEN
           W(IADR2)=W(IADR2)+PH*BUF(IK)
          ENDIF
         ENDIF
        ENDIF
C
        if(j.ne.l.and.i.ne.k)then
        I1=INDXA(MIN(I,L),MAX(I,L))
        I2=INDXA(MIN(J,K),MAX(J,K))
        IF(I.NE.L.AND.J.NE.K)THEN 
         PH=ONE
         IF(I.GT.L)PH=PH*ONEM
         IF(K.GT.J)PH=PH*ONEM
         IRREPIL=DIRPRD(IRREPA(I+NO),IRREPA(L+NO))
         IOFF=IPW(IRREPIL)
         IF(IOFF.LE.ITOPADR)THEN
          NN=NUMIRW(IRREPIL+NIRREP)
          IADR1=IOFF+ISYM(I1)-IPDSZ(IRREPIL)
     &              +NN*(ISYM(I2)-1-IPDIS(IRREPIL)) -IOFFSET
          IADR2=IOFF+ISYM(I2)-IPDSZ(IRREPIL)
     &             +NN*(ISYM(I1)-1-IPDIS(IRREPIL))-IOFFSET
          IF(IADR1.LE.ITOP.AND.IADR1.GT.0)W(IADR1)=W(IADR1)+PH*BUF(IK)
          IF(IADR2.LE.ITOP.AND.IADR2.GT.0.AND.IADR1.NE.IADR2)THEN
           W(IADR2)=W(IADR2)+PH*BUF(IK)
          ENDIF
         ENDIF
         endif
        ENDIF
       ENDIF
11    CONTINUE
      IF(NUT.NE.-1)GOTO 1
C
103   CONTINUE
C
      WRITE(6,1000)ICOUNT
C
1000  FORMAT(T3,'@RDABCDAA-I, Processed ',I12,' ABCD integrals.')
C
      RETURN
      END
