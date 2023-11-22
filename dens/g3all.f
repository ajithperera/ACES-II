      SUBROUTINE G3ALL(T1,T2,G,CC,MBPT4,FACT,ISPIN,ISAME,
     &                 DISSYT,NUMSYT,LISTT1,LISTT2,LISTG,IRREP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ISAME,MBPT4,CC,CIS,EOM,LTRP
      INTEGER DISSYT,DIRPRD,POP,VRT
      DIMENSION T1(DISSYT,NUMSYT),T2(DISSYT,NUMSYT),G(NUMSYT,NUMSYT)
      COMMON/SYM/POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /EXCITE/ CIS,EOM
      COMMON /LTRIP/ LTRP
C
      DATA AZERO,ONE,ONEM,TWO /0.0D0,1.D0,-1.D0,2.D0/
C
C      PICK UP FIRST THE APPROBIATE T AMPLITUDES
C
      IF(.NOT.CC) THEN
C
      CALL GETLST(T1,1,NUMSYT,2,IRREP,LISTT1)
C
C GET SECOND SET OF AMPLITUDES BUT NOT FOR MBPT(3) (ISAME = TRUE )
C
      IF(.NOT.ISAME)CALL GETLST(T2,1,NUMSYT,1,IRREP,LISTT2)
C
C FOR MBPT4 BUILD   [2*T2 - T1]
C
      IF(MBPT4) THEN
       CALL SSCAL(NUMSYT*DISSYT,TWO,T2,1)
       CALL SAXPY(NUMSYT*DISSYT,ONEM,T1,1,T2,1)
      ENDIF
C
C     PERFORM THE MULTIPLICATION IF WE WANT THE INTERMEDIATE.
C
      CALL XGEMM('T','N',NUMSYT,NUMSYT,DISSYT,FACT,T1,DISSYT,
     &           T2,DISSYT,AZERO,G,NUMSYT)
C
      ELSE
       CALL GETLST(G,1,NUMSYT,2,IRREP,LISTT1)
      ENDIF
      IF(EOM .OR. LTRP)THEN
       IX=NUMSYT+1
       CALL GETLST(G(1,IX),1,NUMSYT,1,IRREP,LISTG)
       CALL SAXPY (NUMSYT*NUMSYT,ONE,G(1,IX),1,G,1)
      ENDIF
C
C  SYMMETRIZE THE CALCULATED GAMMA INTERMEDIATE. THIS IS ONLY NECESSARY
C  WHEN THE TWO AMPLITUDES LISTS ARE DIFFERENT (ISAME = FALSE )
C
      IF(.NOT.ISAME) THEN
       CALL SYMMET2(G,NUMSYT) 
      ENDIF
C
C     SAVE THE RESULT ON FILE
C
      CALL PUTLST(G,1,NUMSYT,1,IRREP,LISTG)
C
      RETURN
      END