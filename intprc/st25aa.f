      SUBROUTINE ST25AA(W1,W2,BUF,IBUF,NUMIRW1,ISYM,IPW1,
     &                 IPDIS,IPDSZ,NSIZE1,
     &                 NWDIS,NOCC,NVRT,ITYPE,ISPIN,IUHF,
     &                 NMO,IRREPA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C THIS ROUTINE PERFORMS AN IN-CORE SORT AND FORMS THE
C SYMMETRY PACKED LIST FOR THE <IJ//KA> AND <AB//CI>
C INTEGRALS (AAAA AND BBBB CASES ONLY)
C
CEND
C
C CODED JULY/90  JG
C
      INTEGER DIRPRD,A,B,C,POP,VRT,AND
      CHARACTER*4 SPCASE(3)
      CHARACTER*4 INTYPE(6)
      CHARACTER*8 INAME
      CHARACTER*80 FNAME
      DIMENSION W1(NSIZE1),W2(1),BUF(ILNBUF),IBUF(ILNBUF)
      DIMENSION NUMIRW1(1),IPW1(8),ISYM(1),IRREPA(NMO,2)
      DIMENSION IPDSZ(8),IPDIS(8),IOFFO(8),IOFFV(8)
      COMMON /FLAGS/IFLAGS(100)
      COMMON/MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/FILES/LUOUT,MOINTS
      COMMON /BUFLEN/ ILNBUF
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /SYM/POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /SHIFT/ ISHIFT,NDRGEO
      DATA INTYPE /'HHHH','PHHH','PPHH','PHPH','PPPH','PPPP'/
      DATA SPCASE /'AA  ','BB  ','AB  '/
C
      INDX(I,J)=I+(J*(J-1)/2)
      INDX2(I,J,N)=I+(J-1)*N
      IUPKI(INT)=AND(INT,IALONE)
      IUPKJ(INT)=AND(ISHFT(INT,-IBITWD),IALONE)
      IUPKK(INT)=AND(ISHFT(INT,-2*IBITWD),IALONE)
      IUPKL(INT)=AND(ISHFT(INT,-3*IBITWD),IALONE)
c
      IF(IFLAGS(1).GT.1)THEN
       WRITE(LUOUT,2009) INTYPE(ITYPE),SPCASE(ISPIN)
2009   FORMAT(T3,'@ST25AA-I, Processing integral type ',
     & A4,' spin case ',A2,'.')
      ENDIF
C
C  ZERO THE OUTPUT ARRAY
C
      CALL ZERO(W1,NSIZE1)
C
C  CALCULATE OFFSETS FOR ORBITALS
C
      IOFFO(1)=0
      IOFFV(1)=0
      DO 1 IRREP=1,NIRREP-1
      IOFFO(IRREP+1)=IOFFO(IRREP)+POP(IRREP,ISPIN)
      IOFFV(IRREP+1)=IOFFV(IRREP)+VRT(IRREP,ISPIN)
1     CONTINUE
C
C OPEN INTEGRAL FILE
C
      INAME=INTYPE(ITYPE)//SPCASE(ISPIN)
      CALL GFNAME(INAME,FNAME,ILENGTH)
      OPEN(UNIT=15,FILE=FNAME(1:ILENGTH),
     &     FORM='UNFORMATTED',ACCESS='SEQUENTIAL',STATUS='OLD')
C
C GET INTEGRALS ONE BY ONE AND PLACE THEM IN THE OUTPUT ARRAY
C NOTE THAT THE INTEGRALS <IJ/KA> AND <AB/CI> ARE SAVED WITH
C I<=K;J,A AND A<=B;CI RESPECTIVELY
C
      IF(ITYPE.EQ.2) THEN
C
C PROCESS <IJ/KA> INTEGRALS
C
1000   READ(15) BUF,IBUF,NUT 
C
CDIR$ IVDEP
C
*VOCL LOOP,NOVREC
       DO 1001 INT=1,NUT
        I=IUPKI(IBUF(INT))
        J=IUPKJ(IBUF(INT))
        K=IUPKK(IBUF(INT))
        A=IUPKL(IBUF(INT))
C
        I1=INDX(I,K)
        I2=INDX2(J,A,NOCC)
        IRREPIK=DIRPRD(IRREPA(I,ISPIN),IRREPA(K,ISPIN))
        IOFF1=IPW1(IRREPIK)
        NN=NUMIRW1(IRREPIK+NIRREP)
        IADR1=IOFF1+ISYM(I1+NWDIS)-IPDSZ(IRREPIK)
     &            +NN*(ISYM(I2)-1-IPDIS(IRREPIK)) 
        W1(IADR1)=BUF(INT)
1001   CONTINUE
       IF(NUT.EQ.ILNBUF) GOTO 1000
       DO 1002 IRREP=1,NIRREP
       IX=0
       NUMSYW=NUMIRW1(IRREP)
       DO 1003 IRREPAA=1,NIRREP 
       IRREPK=DIRPRD(IRREP,IRREPAA)
       IOFFA=IOFFV(IRREPAA)
       NVRTA=VRT(IRREPAA,ISPIN)
       IOFFK=IOFFO(IRREPK)
       NOCCK=POP(IRREPK,ISPIN)
       DO 1003 A=IOFFA+1,IOFFA+NVRTA
       DO 1003 K=IOFFK+1,IOFFK+NOCCK
       DO 1003 IRREPJ=1,NIRREP
       IOFFJ=IOFFO(IRREPJ)
       NOCCJ=POP(IRREPJ,ISPIN)
       IRREPI=DIRPRD(IRREP,IRREPJ)
       IOFFI=IOFFO(IRREPI)
       NOCCI=POP(IRREPI,ISPIN)
       IRREPIK=DIRPRD(IRREPI,IRREPK)
       IRREPJK=DIRPRD(IRREPJ,IRREPK)
       IF(IRREPI.EQ.IRREPJ) THEN
        DO 1004 J=IOFFJ+2,IOFFJ+NOCCJ
        DO 1004 I=IOFFI+1,J-1
        I1=INDX(MIN(I,K),MAX(I,K))
        I2=INDX2(J,A,NOCC)
        I3=INDX(MIN(J,K),MAX(J,K))
        I4=INDX2(I,A,NOCC)
        IOFF1=IPW1(IRREPIK)
        NN=NUMIRW1(IRREPIK+NIRREP)
        IADR1=IOFF1+ISYM(I1+NWDIS)-IPDSZ(IRREPIK)
     &             +NN*(ISYM(I2)-1-IPDIS(IRREPIK))
        IOFF2=IPW1(IRREPJK)
        NN=NUMIRW1(IRREPJK+NIRREP)
        IADR2=IOFF2+ISYM(I3+NWDIS)-IPDSZ(IRREPJK)
     &             +NN*(ISYM(I4)-1-IPDIS(IRREPJK))
        IX=IX+1
        W2(IX)=W1(IADR1)-W1(IADR2) 
1004    CONTINUE
       ELSE IF(IRREPI.LT.IRREPJ) THEN
        DO 1005 J=IOFFJ+1,IOFFJ+NOCCJ
        DO 1005 I=IOFFI+1,IOFFI+NOCCI
        I1=INDX(MIN(I,K),MAX(I,K))
        I2=INDX2(J,A,NOCC)
        I3=INDX(MIN(J,K),MAX(J,K))
        I4=INDX2(I,A,NOCC)
        IOFF1=IPW1(IRREPIK)
        NN=NUMIRW1(IRREPIK+NIRREP)
        IADR1=IOFF1+ISYM(I1+NWDIS)-IPDSZ(IRREPIK)
     &             +NN*(ISYM(I2)-1-IPDIS(IRREPIK))
        IOFF2=IPW1(IRREPJK)
        NN=NUMIRW1(IRREPJK+NIRREP)
        IADR2=IOFF2+ISYM(I3+NWDIS)-IPDSZ(IRREPJK)
     &             +NN*(ISYM(I4)-1-IPDIS(IRREPJK))
        IX=IX+1
        W2(IX)=W1(IADR1)-W1(IADR2) 
1005    CONTINUE 
       ENDIF
1003   CONTINUE
       CALL PUTLST(W2,1,NUMSYW,2,IRREP,6+ISPIN+ISHIFT)
1002   CONTINUE
C
       ELSE IF (ITYPE.EQ.5) THEN   
C
C PROCESS <AB//CI> INTEGRALS
C
2000   READ(15) BUF,IBUF,NUT 
C
CDIR$ IVDEP
C
*VOCL LOOP,NOVREC
       DO 2001 INT=1,NUT
        I=IUPKI(IBUF(INT))
        A=IUPKJ(IBUF(INT))
        B=IUPKK(IBUF(INT))
        C=IUPKL(IBUF(INT))
C
C CONTRIBUTIONS OF THIS INTEGRAL TO <AB/CI> AND <CB/AI>
C
        I1=INDX(A,C)
        I2=INDX2(B,I,NVRT)
        IRREPAC=DIRPRD(IRREPA(A+NOCC,ISPIN),IRREPA(C+NOCC,ISPIN))
        IOFF1=IPW1(IRREPAC)
        NN=NUMIRW1(IRREPAC+NIRREP)
        IADR1=IOFF1+ISYM(I1+NWDIS)-IPDSZ(IRREPAC)
     &            +NN*(ISYM(I2)-1-IPDIS(IRREPAC)) 
        W1(IADR1)=BUF(INT)
2001   CONTINUE
       IF(NUT.EQ.ILNBUF) GOTO 2000
       DO 2002 IRREP=1,NIRREP
       IX=0
       NUMSYW=NUMIRW1(IRREP)
       DO 2003 IRREPI=1,NIRREP 
       IRREPC=DIRPRD(IRREP,IRREPI)
       IOFFI=IOFFO(IRREPI)
       NOCCI=POP(IRREPI,ISPIN)
       IOFFC=IOFFV(IRREPC)
       NVRTC=VRT(IRREPC,ISPIN)
       DO 2003 I=IOFFI+1,IOFFI+NOCCI
       DO 2003 C=IOFFC+1,IOFFC+NVRTC
       DO 2003 IRREPBB=1,NIRREP
       IOFFB=IOFFV(IRREPBB)
       NVRTB=VRT(IRREPBB,ISPIN)
       IRREPAA=DIRPRD(IRREP,IRREPBB)
       IOFFA=IOFFV(IRREPAA)
       NVRTA=VRT(IRREPAA,ISPIN)
       IRREPAC=DIRPRD(IRREPAA,IRREPC)
       IRREPBC=DIRPRD(IRREPBB,IRREPC)
       IF(IRREPAA.EQ.IRREPBB) THEN
        DO 2004 B=IOFFB+2,IOFFB+NVRTB
        DO 2004 A=IOFFA+1,B-1
        I1=INDX(MIN(A,C),MAX(A,C))
        I2=INDX2(B,I,NVRT)
        I3=INDX(MIN(B,C),MAX(B,C))
        I4=INDX2(A,I,NVRT)
        IOFF1=IPW1(IRREPAC)
        NN=NUMIRW1(IRREPAC+NIRREP)
        IADR1=IOFF1+ISYM(I1+NWDIS)-IPDSZ(IRREPAC)
     &             +NN*(ISYM(I2)-1-IPDIS(IRREPAC))
        IOFF2=IPW1(IRREPBC)
        NN=NUMIRW1(IRREPBC+NIRREP)
        IADR2=IOFF2+ISYM(I3+NWDIS)-IPDSZ(IRREPBC)
     &             +NN*(ISYM(I4)-1-IPDIS(IRREPBC))
        IX=IX+1
        W2(IX)=W1(IADR1)-W1(IADR2) 
2004    CONTINUE
       ELSE IF(IRREPAA.LT.IRREPBB) THEN
        DO 2005 B=IOFFB+1,IOFFB+NVRTB
        DO 2005 A=IOFFA+1,IOFFA+NVRTA
        I1=INDX(MIN(A,C),MAX(A,C))
        I2=INDX2(B,I,NVRT)
        I3=INDX(MIN(B,C),MAX(B,C))
        I4=INDX2(A,I,NVRT)
        IOFF1=IPW1(IRREPAC)
        NN=NUMIRW1(IRREPAC+NIRREP)
        IADR1=IOFF1+ISYM(I1+NWDIS)-IPDSZ(IRREPAC)
     &             +NN*(ISYM(I2)-1-IPDIS(IRREPAC))
        IOFF2=IPW1(IRREPBC)
        NN=NUMIRW1(IRREPBC+NIRREP)
        IADR2=IOFF2+ISYM(I3+NWDIS)-IPDSZ(IRREPBC)
     &             +NN*(ISYM(I4)-1-IPDIS(IRREPBC))
        IX=IX+1
        W2(IX)=W1(IADR1)-W1(IADR2) 
2005    CONTINUE 
       ENDIF
2003   CONTINUE
       CALL PUTLST(W2,1,NUMSYW,2,IRREP,26+ISPIN+iSHIFT)
2002  CONTINUE
      ENDIF
      IF(IUHF.EQ.0) THEN
       CLOSE(UNIT=15,STATUS='KEEP')
      ELSE
       CLOSE(UNIT=15,STATUS='DELETE')
      ENDIF
      RETURN
      END  
