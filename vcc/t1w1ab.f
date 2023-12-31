
      SUBROUTINE T1W1AB(Z,ZT,T,TT,W,MAXSIZE,TA,TB,POP1,POP2,VRT1,VRT2,
     &                DISSYZ,DISSYWA,DISSYWB,NUMSYZ,NUMSYWA,NUMSYWB,
     &                NTASIZ,NTBSIZ,LISTT,LISTZ,LISTWA,LISTWB,IRREP,
     &                IUHF,TMP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DISSYZ,DISSYWA,DISSYWB,DIRPRD,POP1,POP2,VRT1,VRT2
      LOGICAL ROHF4,ITRFLG
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,
     &        SING1,QCISD,UCC,CC2
      DIMENSION Z(DISSYZ,NUMSYZ),W(DISSYWA,1),TA(NTASIZ)
      DIMENSION T(DISSYZ,NUMSYZ),TB(NTBSIZ)
      DIMENSION ZT(NUMSYZ,DISSYZ),TT(NUMSYZ,DISSYZ)
      DIMENSION TMP(1),POP1(8),POP2(8),VRT1(8),VRT2(8)
C
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON/SWITCH/MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,
     &              SING1,QCISD,UCC,CC2
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /ROHF/ ROHF4,ITRFLG
C
      DATA AZERO,ONE,ONEM /0.0D0,1.0D0,-1.0D0/
C
C     ZERO FIRST OUTPUT ARRAY (NECCESARY BECAUSE THIS IS NOT NECCESARILY
C     DONE IN THE MATRIX MULTIPLICATION
C
      CALL ZERO(Z,NUMSYZ*DISSYZ)
C
C     IF THERE ARE NO WA INTEGRALS THERE IS NOTHING TO DO
C
      IF(MIN(NUMSYWA,DISSYWA).NE.0) THEN
C
C GET T2 AMPLITUDES AND FORM TAU AMPLITUDES FIRST
C    
       IF (CC2) THEN
          CALL ZERO(T,NUMSYZ*DISSYZ)
       ELSE
          CALL GETLST(T,1,NUMSYZ,1,IRREP,LISTT)
       ENDIF

       IF(.NOT.ROHF4)THEN
        CALL FTAU(T,TA,TB,DISSYZ,NUMSYZ,POP1,POP2,VRT1,VRT2,
     &            IRREP,3,ONE)
       ENDIF
C
C     PROCESS AS MANY Am DISTRIBUTIONS AT ONCE AS POSSIBLE
C
C GET INTEGRALS <Ef|Am> FROM LISTWA
C
       NINCOR=MAXSIZE/(DISSYWA*IINTFP)
       NLEFT =NUMSYWA
       NFIRST=1
       NPASS=0
1      NREAD =MIN(NLEFT,NINCOR)
       CALL GETLST(W,NFIRST,NREAD,1,IRREP,LISTWA)
C
C   PERFORM MULTIPLICATION
C                   +
C           T(Ef,Ij) * W(Ef,Am) = Z(Ij,Am)
C
       CALL XGEMM ('T','N',NUMSYZ,NREAD,DISSYZ,ONE,T,
     &             DISSYZ,W,DISSYWA,AZERO,ZT(1,NFIRST),NUMSYZ)
       NFIRST=NFIRST+NREAD
       NLEFT =NLEFT-NREAD
       NPASS=NPASS+1
       IF(NLEFT.NE.0)GOTO 1
C
c       write(6,1001)npass
c1001   format(T3,'@T1W1AB-I, First phase required ',I5,' passes.')
C
C     DO THE SECOND PART OF THE MULTIPLICATION
C
      CALL ZERO(TT,NUMSYZ*DISSYZ)
      IOFF=1
      JOFFZ=1
      JOFFT=1
      DO 90 IRREPJ=1,NIRREP
C
       NOCCJ=POP2(IRREPJ)
       NVRTJ=VRT2(IRREPJ)
C
       IRREPI=DIRPRD(IRREPJ,IRREP)
C
       NVRTI=VRT1(IRREPI)
C
       IF(NVRTI.NE.0.AND.NOCCJ.NE.0.AND.NVRTJ.NE.0) THEN
C
        CALL XGEMM('N','T',NUMSYZ*NVRTI,NVRTJ,NOCCJ,ONE,ZT(1,JOFFZ),
     &             NUMSYZ*NVRTI,TB(IOFF),NVRTJ,AZERO,TT(1,JOFFT),      
     &             NUMSYZ*NVRTI)
C
       ENDIF
C
C      UPDATE POINTERS
C
        IOFF=IOFF+NVRTJ*NOCCJ
        JOFFZ=JOFFZ+NVRTI*NOCCJ
        JOFFT=JOFFT+NVRTJ*NVRTI
90     CONTINUE
C
       CALL TRANSP(TT,Z,DISSYZ,NUMSYZ) 
       CALL GETLST(T,1,NUMSYZ,1,IRREP,LISTZ)
       CALL VADD(T,T,Z,NUMSYZ*DISSYZ,ONEM)
       IF(IUHF.EQ.1) CALL PUTLST(T,1,NUMSYZ,1,IRREP,LISTZ)
      ENDIF
C
      IF(IUHF.EQ.0) THEN
C
C RHF CASE:
C
C IN RHF THIS IS SIMPLY A TRANSPOSITION
C
       IF(MIN(NUMSYWA,DISSYWA).EQ.0) RETURN
       CALL SYMTR1(IRREP,POP1,POP2,DISSYZ,Z,TMP,TMP(1+DISSYZ),
     &             TMP(1+2*DISSYZ))
       CALL SYMTR3(IRREP,VRT1,VRT2,DISSYZ,NUMSYZ,Z,TMP,
     &             TMP(1+NUMSYZ),TMP(1+2*NUMSYZ))
       CALL VADD(T,T,Z,NUMSYZ*DISSYZ,ONEM)
       CALL PUTLST(T,1,NUMSYZ,1,IRREP,LISTZ)
    
      ELSE
C
C      IF THERE ARE NO INTEGRALS SKIP MULTIPLICATION
C
       IF(MIN(NUMSYWB,DISSYWB).EQ.0) RETURN
C
C    DECIDE ABOUT ALGORITHM
C
       CALL ZERO(ZT,NUMSYZ*DISSYWB)
C
C   GET T2 AMPLITUDES AND FORM TAU AMPLITUDES
C     
        IF (CC2) THEN
           CALL ZERO(T,NUMSYZ*DISSYZ)
        ELSE
           CALL GETLST(T,1,NUMSYZ,1,IRREP,LISTT) 
        ENDIF 

        IF(.NOT.ROHF4)THEN
         CALL FTAU(T,TA,TB,DISSYZ,NUMSYZ,POP1,POP2,VRT1,VRT2,
     &             IRREP,3,ONE)
        ENDIF
C
C     PROCESS AS MANY Ej DISTRIBUTIONS AT ONCE AS POSSIBLE
C
C GET INTEGRALS <Ab|Ej> FROM LISTWA
C
       NINCOR=MAXSIZE/(DISSYWB*IINTFP)
       NLEFT =NUMSYWB
       NFIRST=1
       NPASS=0
2      NREAD =MIN(NLEFT,NINCOR)
       CALL GETLST(W,NFIRST,NREAD,1,IRREP,LISTWB)
C
C   PERFORM MULTIPLICATION
C                   +
C           T(Ab,Mn) * W(Ab,Ej) = Z(Mn,Ej)
C
       CALL XGEMM ('T','N',NUMSYZ,NREAD,DISSYZ,ONE,T,
     &             DISSYZ,W,DISSYWB,AZERO,ZT(1,NFIRST),NUMSYZ)
       NFIRST=NFIRST+NREAD
       NLEFT =NLEFT-NREAD
       NPASS=NPASS+1
       IF(NLEFT.NE.0)GOTO 2
C
c       write(6,1002)npass
c1002   format(T3,'@T1W1AB-I, Second phase required ',I5,' passes.')
C
C   DO THE SECOND PART OF THE MULTIPLICATION 
C
        CALL ZERO(TT,NUMSYZ*DISSYZ)
        CALL SYMTR1(IRREP,POP1,VRT2,NUMSYZ,ZT,TMP,TMP(1+NUMSYZ),
     &               TMP(1+2*NUMSYZ))
C
        IOFF=1
        JOFFT=1
        JOFFZ=1
        DO 190 IRREPJ=1,NIRREP
C
        NOCCJ=POP1(IRREPJ)
        NVRTJ=VRT1(IRREPJ)
        IRREPI=DIRPRD(IRREPJ,IRREP)
        NVRTI=VRT2(IRREPI)
        IF(NVRTI.NE.0.AND.NOCCJ.NE.0.AND.NVRTJ.NE.0) THEN  
C
         CALL XGEMM('N','T',NUMSYZ*NVRTI,NVRTJ,NOCCJ,ONE,ZT(1,JOFFZ),
     &              NUMSYZ*NVRTI,TA(IOFF),NVRTJ,ONE,TT(1,JOFFT),
     &              NUMSYZ*NVRTI)
        ENDIF
       JOFFT=JOFFT+NVRTJ*NVRTI
       JOFFZ=JOFFZ+NVRTI*NOCCJ
       IOFF=IOFF+NOCCJ*NVRTJ
190    CONTINUE
       CALL SYMTRA(IRREP,VRT2,VRT1,NUMSYZ,TT,ZT)
       CALL TRANSP(ZT,T,DISSYZ,NUMSYZ)
       CALL GETLST(Z,1,NUMSYZ,2,IRREP,LISTZ)
       CALL VADD(Z,Z,T,NUMSYZ*DISSYZ,ONEM)
       CALL PUTLST(Z,1,NUMSYZ,1,IRREP,LISTZ)
      ENDIF
      RETURN
      END
