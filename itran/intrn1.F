      SUBROUTINE INTRN1(EVEC,AOINT,BUF1,RHF,NBASIS,
     &                  NBAS,NMO,NFIRST,NFIRST2,
     &                  NSTART,NEND,ISIZE,IOFFMO,
     &                  NSIZE,ISPIN,LAST,IREORD,IRREPX)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL RHF,DOALL,ABIJ,NOABCD,LAST
      CHARACTER*80 FNAME
      INTEGER DIRPRD,DISZAO,DISZMO
C
      DIMENSION EVEC(1),AOINT(1),BUF1(1)
      DIMENSION NBAS(8),NMO(8),NFIRST(8),NFIRST2(8),NSTART(8),NEND(8),
     &          ISIZE(8),IOFFMO(8,2),IREORD(1)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/AOOFST/INDOCC(8,2)   
      COMMON/SYMINF/NDUMMY,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/DOINTS/DOALL,ABIJ,NOABCD
C
      DATA AZERO,ONE/0.D0,1.D0/
C
      IF(.NOT.RHF.AND.(ISPIN.EQ.1.OR.(.NOT.DOALL))) THEN
       CALL GFNAME('HALFT   ',FNAME,ILENGTH)
      ENDIF
C
C LOOP OVER ALL IRREPS
C
      IOFF2=1
      IOFF3=1
C
      DO 1000 IRREP=1,NIRREP
C
       IF (LAST)  IRREPX = IRREP
C
       IF (IRREP .EQ. IRREPX)  THEN
       Write(6,*) irrep, irrepx
       IOFF2I=IOFF2
       IOFF3I=IOFF3
       NAOI=NBAS(IRREP)
       NMOI=NMO(IRREP)
       NFIRSTI=NFIRST(IRREP)
       NFIRST2I=NFIRST2(IRREP)
       NSTARTI=NSTART(IRREP)
       NENDI=NEND(IRREP)
       IF(NENDI.LT.NSTARTI) GO TO 1000
       DISZAO=NAOI*(NAOI+1)/2
       DISZMO=NMOI*(NMOI+1)/2
       NUMDIS=NMOI*(NENDI+1-NSTARTI)

#ifdef _DEBUG_LVLM1
       Write(6,*)
       write(*,"(t3,a,a,2(1x,i2))"), "@INTRN1: The 1-index",
     &      "transformed Ints for irrep at entry 2 intrn1", 
     &       Irrep, ispin
       Ioff = 0
       Do KL = 1, DISZAO
          Write(6,*), "Printing the distribution No. ", KL
          Write(*,'(4(1X,F12.7)))'), (Aoint(Ioff + I), I=1, NAOI*NMOI)
          Ioff = Ioff + NAOI*NMOI
       End Do
#endif

C
C LOOP OVER ALL OCCUPIED ORBITALS OR ALL ORBITALS IN CASE OF
C A FULL TRANSFORMATION
C
       IOFF1=INDOCC(IRREP,ISPIN)
C
       DO 100 IFIRST=NSTARTI,NENDI
C
C  TRANSFORM SECOND INDEX  (Mu I | SIGMA RHO) --> (P I | SIGMA RHO )
C
#ifdef _DEBUG_LVM1
         Write(6,*)
         Write(6,*), "Printing the distribution No. ", ifirst,
     &                Ioff2, ioff1
         Write(*,'(4(1X,F12.7)))'), (Aoint(Ioff2-1 + I), I=1, 
     &         diszao*NAOI)
#endif
  
        CALL XGEMM('T','N',NMOI,DISZAO,NAOI,ONE,EVEC(IOFF1),NAOI,
     &             AOINT(IOFF2),NAOI,AZERO,BUF1,NMOI)
C
#ifdef _DEBUG_LVL0
         Write(6,*) Nstarti, nendi
         Write(6,*), "Printing the distribution No. after T: ",ioff3
         Write(*,'(4(1X,F12.7)))'), (Buf1(I), I=1, diszao*nmoi)
         Write(6,*)
#endif
  
        CALL TRANSP(BUF1,AOINT(IOFF3),DISZAO,NMOI)
C
        IOFF2=IOFF2+NAOI*DISZAO
        IOFF3=IOFF3+NMOI*DISZAO
C
100    CONTINUE
C
C  NOW WE HAVE THE (PI,SIGMA RHO) INTEGRALS ORDERED AS   X X P I 
C
C     FOR UHF, WRITE INTEGRALS TO DISK SINCE WE WILL USE THEM AGAIN
C
#ifdef _DEBUG_LVM1
       Write(6,*)
       write(*,"(t3,a,a,2(1x,i2))"), "@INTRN1: The 1-index",
     &      "transformed Ints for irrep after stage 1 intrn1",
     &       Irrep, ispin
       Ioff = ioff3i 
       Do KL = 1, NUMDIS
          Write(6,*), "Printing the distribution No. ", KL
          Write(*,'(4(1X,F12.7)))'), (Aoint(Ioff + I), I=1, diszmo)
          Ioff = Ioff + Diszmo
       End Do
#endif

       IF(.NOT.RHF.AND.(ISPIN.EQ.1.OR.(.NOT.DOALL))) THEN
C
        OPEN(UNIT=79,FORM='UNFORMATTED',STATUS='UNKNOWN',
     &       ACCESS='SEQUENTIAL',FILE=FNAME(1:ILENGTH))
        REWIND(UNIT=79)
C
        IOFFW=IOFF3I
        DO 101 I=1,NUMDIS
         WRITE(79) (AOINT(IND),IND=IOFFW,IOFFW+DISZAO-1)
         IOFFW=IOFFW+DISZAO
101     CONTINUE
C
        REWIND(UNIT=79)
C
       ENDIF
C
C NOW TRANSFORM NEW LHS TO RS, SPIN CASE ISPIN,ISPIN
C
       IOFF4=IOFF3I
       IOFF5=IOFF3I
C
       IX=1+NAOI*NAOI
C
       DO 200 IDIS=1,NUMDIS
C
        CALL EXPND2(AOINT(IOFF4),BUF1,NAOI)
C     
        CALL GHTRAN(BUF1,BUF1,EVEC,IOFF1,IOFF1,1,1,BUF1(IX),
     &              0,1,NAOI,NAOI,NMOI,NMOI)
C
        CALL SQUEZ2(BUF1,AOINT(IOFF5),NMOI)
C
        IOFF4=IOFF4+DISZAO 
        IOFF5=IOFF5+DISZMO
C
200    CONTINUE
C
C DUMP INTEGRALS TO HF2
C
#ifdef _DEBUG_LVLM1
       Write(6,*)
       write(*,"(t3,a,a,2(1x,i2))"), "@INTRN1: The final 4-index",
     &      "transformed Ints for irrep", Irrep, ispin
       Ioff = ioff3i-1
       Do KL = 1, NUMDIS
          Write(6,*), "Printing the distribution No.", KL
          Write(*,'(4(1X,F12.7)))'), (Aoint(Ioff + I), I=1, diszmo)
          Write(6,*) 
          Write(6,'(8(1X,I5)))') (IREORD(I), I=1, Nbasis)
          Call checksum2("INTR1", AOint(ioff), diszmo)
          Ioff = Ioff + Diszmo
       End Do
#endif
      ILNBUF=600
      IOFF=1+ILNBUF
      IOFF11=IOFFMO(IRREP,ISPIN)
      IF(DOALL) THEN
       CALL DHF11(AOINT(IOFF3I),BUF1,BUF1(IOFF),ILNBUF,
     &            IRREP,NSTARTI,NENDI,NMOI,IOFF11,IOFF11,
     &            IREORD,RHF,ISPIN)
      ELSE
       CALL DHF21(AOINT(IOFF3I),BUF1,BUF1(IOFF),ILNBUF,
     &            IRREP,NSTARTI,NENDI,NFIRSTI,NMOI,
     &            IOFF11,IOFF11,IREORD,RHF,ISPIN)
       
      ENDIF
C
      IF(.NOT.RHF.AND.(ISPIN.EQ.1.OR.(.NOT.DOALL))) THEN
       IOFFW=IOFF3I
       DO 1101 I=1,NUMDIS
        READ(79) (AOINT(IND),IND=IOFFW,IOFFW+DISZAO-1)
        IOFFW=IOFFW+DISZAO
1101   CONTINUE
C
       CLOSE(UNIT=79,STATUS='DELETE')
C
C NOW TRANSFORM NEW LHS TO RS, SPIN CASE ISPIN,ISPIN
C
       IOFF4=IOFF3I
       IOFF5=IOFF3I
C
       IX=1+NAOI*NAOI
C
       IOFF1=INDOCC(IRREP,3-ISPIN)
C

       DO 1200 IDIS=1,NUMDIS
C
        CALL EXPND2(AOINT(IOFF4),BUF1,NAOI)
C     
        CALL GHTRAN(BUF1,BUF1,EVEC,IOFF1,IOFF1,1,1,BUF1(IX),
     &              0,1,NAOI,NAOI,NMOI,NMOI)
C
        CALL SQUEZ2(BUF1,AOINT(IOFF5),NMOI)
C
        IOFF4=IOFF4+DISZAO 
        IOFF5=IOFF5+DISZMO
C
1200    CONTINUE

C
#ifdef _DEBUG_LVM1
       Write(6,*)
       write(*,"(t3,a,1x,i2)"), "@INTRN1: The final 4-index 
     &      transformed Ints for irrep", Irrep, ispin 
       Ioff = Ioff3i 
       Do KL = 1, NUMDIS
          Write(6,*)
          Write(*,'(4(1X,F12.7)))'), (Aoint(Ioff + I), I=1, diszmo)
          Call checksum2("INTR1", AOint(ioff), diszmo)
          Ioff = Ioff + Diszmo
       End Do
#endif
C
C DUMP INTEGRALS TO HF2
C
       ILNBUF=600
       IOFF=1+ILNBUF
       IOFF11=IOFFMO(IRREP,ISPIN)
       IOFF22=IOFFMO(IRREP,3-ISPIN)
       IF(DOALL) THEN
        CALL DHF11(AOINT(IOFF3I),BUF1,BUF1(IOFF),ILNBUF,
     &             IRREP,NSTARTI,NENDI,NMOI,IOFF11,IOFF22,
     &             IREORD,RHF,2+ISPIN)
       ELSE
        CALL DHF21(AOINT(IOFF3I),BUF1,BUF1(IOFF),ILNBUF,
     &             IRREP,NSTARTI,NENDI,NFIRST2I,NMOI,
     &             IOFF11,IOFF22,IREORD,RHF,2+ISPIN)
       ENDIF
C
      ENDIF
C
      ENDIF
1000  CONTINUE
C
      RETURN
      END
