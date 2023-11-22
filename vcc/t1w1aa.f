      
      SUBROUTINE T1W1AA(Z,ZT,T,TT,W,MAXSIZE,TA,POP,VRT,
     &                DISSYZ,DISSYW,NUMSYZ,NUMSYW,NVRTSQ,
     &                NTASIZ,LISTT,LISTZ,LISTW,IRREP,
     &                IUHF,ISPIN,TMP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DISSYZ,DISSYW,DIRPRD,POP,VRT
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,
     &        SING1,QCISD,UCC,CC2
      LOGICAL ITRFLG,ROHF4
      DIMENSION Z(DISSYZ,NUMSYZ),W(DISSYW,1),TA(NTASIZ)
      DIMENSION T(DISSYZ,NUMSYZ)
      DIMENSION ZT(NUMSYZ,DISSYZ),TT(NUMSYZ,DISSYZ)
      DIMENSION TMP(1),POP(8),VRT(8)
C
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON/SWITCH/MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,
     &              SING1,QCISD,UCC,CC2
      COMMON /ROHF/ ROHF4,ITRFLG
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA AZERO,ONE,ONEM /0.0D0,1.0D0,-1.0D0/
C
C     ZERO FIRST OUTPUT ARRAY (NECCESARY BECAUSE THIS IS NOT NECCESARILY
C     DONE IN THE MATRIX MULTIPLICATION
C
      CALL ZERO(Z,NUMSYZ*DISSYZ)
C
C     IF THERE ARE NO W INTEGRALS THERE IS NOTHING TO DO
C
      IF(MIN(NUMSYW,DISSYW).NE.0) THEN
C
C   GET T2 AMPLITUDES AND FORM TAU AMPLITUDES FIRST, For CC2, TAu
C   is only T1 * T1.
C
       IF (CC2) THEN
           CALL ZERO(T,NUMSYZ*DISSYZ)
       ELSE
           CALL GETLST(T,1,NUMSYZ,1,IRREP,LISTT)
       ENDIF 

       IF(.NOT.ROHF4)THEN
        CALL FTAU(T,TA,TA,DISSYZ,NUMSYZ,POP,POP,VRT,VRT,
     &            IRREP,ISPIN,ONE)
       ENDIF
C
C     PROCESS AS MANY AM DISTRIBUTIONS AT ONCE AS POSSIBLE
C
C   GET INTEGRALS W(E<F,AM) FROM LISTW
C
       NINCOR=MAXSIZE/(DISSYW*IINTFP)
       NLEFT =NUMSYW
       NFIRST=1
       NPASS=0
1      NREAD =MIN(NLEFT,NINCOR)
       CALL GETLST(W,NFIRST,NREAD,1,IRREP,LISTW)
C
C   PERFORM MULTIPLICATION
C                     +
C           T(E<F,I<J) * W(E<F,AM) = Z(I<J,AM)
C
       CALL XGEMM ('T','N',NUMSYZ,NREAD,DISSYZ,ONE,T,
     &             DISSYZ,W,DISSYW,AZERO,ZT(1,NFIRST),NUMSYZ)
       NFIRST=NFIRST+NREAD
       NLEFT =NLEFT-NREAD
       NPASS=NPASS+1
       IF(NLEFT.NE.0)GOTO 1
C
c       write(6,1001)npass
c1001   format(T3,'@T1W1AA-I, First phase required ',I5,' passes.')
C
C     DO THE SECOND PART OF THE MULTIPLICATION
C
      CALL ZERO(TT,NUMSYZ*NVRTSQ)
      IOFF=1
      JOFFZ=1
      JOFFT=1
      DO 90 IRREPJ=1,NIRREP
C
       NOCCJ=POP(IRREPJ)
       NVRTJ=VRT(IRREPJ)
C
       IRREPI=DIRPRD(IRREPJ,IRREP)
C
       NVRTI=VRT(IRREPI)
C
       IF(NVRTI.NE.0.AND.NOCCJ.NE.0.AND.NVRTJ.NE.0) THEN
C
        CALL XGEMM('N','T',NUMSYZ*NVRTI,NVRTJ,NOCCJ,ONE,ZT(1,JOFFZ),
     &             NUMSYZ*NVRTI,TA(IOFF),NVRTJ,AZERO,TT(1,JOFFT),      
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
       CALL ASSYM(IRREP,VRT,NUMSYZ,NUMSYZ,ZT,TT)
       CALL TRANSP(ZT,T,DISSYZ,NUMSYZ) 
       CALL GETLST(Z,1,NUMSYZ,1,IRREP,LISTZ)
       CALL VADD(Z,Z,T,NUMSYZ*DISSYZ,ONEM)
       CALL PUTLST(Z,1,NUMSYZ,1,IRREP,LISTZ)
      ENDIF
C
      RETURN
      END