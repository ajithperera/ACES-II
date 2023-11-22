      SUBROUTINE GTAO(IFLG,NATOMS,NBAS,NBASP,NTANGM,ICNTR,
     & COEFF,COEFFN,NANGMOM,NMOMAO,NAOATM,XOCC,DORDR,IORDR,
     & EVALS,LEN,ISPIN,IPRINT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION COEFF(NBAS*NBASP),COEFFN(NBAS*NBASP),ICNTR(NBAS),
     & NANGMOM(NATOMS),NMOMAO(NATOMS*NTANGM),
     & NAOATM(NATOMS),DORDR(NBASP),XOCC(NBASP,2),EVALS(NBASP),
     & IORDR(NBASP)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /IPAR/ LUOUT
C
      NBAS2=NBAS*NBASP
C
      IF(IFLG.EQ.0) THEN
C Get list of centers to which atomic orbitals belong
         CALL GETREC(20,'JOBARC','CNTERBF0',NBAS,ICNTR)
C Get the matrix of alpha MO's in terms of AO's.
         CALL GETREC(20,'JOBARC','EVECAO_A',IINTFP*NBAS2,COEFF)
         IF (LEN .GT. 0) THEN 
            WRITE(6,*) 
            WRITE(6,"(2A)") "The alpha correlated natural orbitals",
     &                      " are read"
            CALL GETREC(20,'JOBARC','CORR_EVA',IINTFP*NBAS2,COEFF)
         ELSE
            CALL GETREC(20,'JOBARC','ORB_OCCA',NBASP*IINTFP,DORDR)
            CALL GETREC(20,"JOBARC","SCFEVLA0",NBASP*IINTFP,EVALS)
            DO IMO = 1, NBASP
               IORDR(IMO) = IMO
            ENDDO 
            CALL PIKSR2(NBASP,EVALS,IORDR)
            DO IMO = 1, NBASP
               XOCC(IMO,1) = DORDR(IORDR(IMO))
            ENDDO
         ENDIF
      ELSE 
C  Get the matrix of alpha MO's in terms of AO's.
         CALL GETREC(20,'JOBARC','EVECAO_B',IINTFP*NBAS2,COEFF)
         IF (LEN .GT. 0) THEN
            WRITE(6,"(2A)") "The beta  correlated natural orbitals",
     &                      " are read"
            CALL GETREC(20,'JOBARC','CORR_EVB',IINTFP*NBAS2,COEFF)
         ELSE
            CALL GETREC(20,'JOBARC','ORB_OCCB',NBASP*IINTFP,DORDR)
            CALL GETREC(20,"JOBARC","SCFEVLB0",NBASP*IINTFP,EVALS)
            DO IMO = 1, NBASP
               IORDR(IMO) = IMO
            ENDDO 
            CALL PIKSR2(NBASP,EVALS,IORDR)
            DO IMO = 1, NBASP
               XOCC(IMO,2) = DORDR(IORDR(IMO))
            ENDDO
         ENDIF
      ENDIF

C     Put higher angular momentum vectors in correct order
      ICNT=1
      IATMOFF=0
      DO 110 IATM=1,NATOMS
         DO 120 IANGMOM=1,NANGMOM(IATM)
            NCNT=1
            DO 125 IMOMAO=1,NMOMAO((IATM-1)*NTANGM+IANGMOM)
               IF(IANGMOM.EQ.1)THEN
                  DO 130 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=COEFF((IMO-1)*NBAS+ICNT)
  130             CONTINUE
               ENDIF
               IOFF=IATMOFF+NMOMAO((IATM-1)*NTANGM+1)
               IF(IANGMOM.EQ.2)THEN
                  DO 141 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     &                           COEFF((IMO-1)*NBAS+(NCNT-1)*3+IOFF+1)
  141             CONTINUE
               ENDIF
               IF(IANGMOM.EQ.3)THEN
                  DO 142 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     &                          COEFF((IMO-1)*NBAS+(NCNT-1)*3+IOFF+2)
  142             CONTINUE
               ENDIF
               IF(IANGMOM.EQ.4)THEN
                  DO 143 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     &                          COEFF((IMO-1)*NBAS+(NCNT-1)*3+IOFF+3)
  143             CONTINUE
               ENDIF
               IOFF=IATMOFF+NMOMAO((IATM-1)*NTANGM+1)+
     &                                    3*NMOMAO((IATM-1)*NTANGM+2)
               IF(IANGMOM.EQ.5)THEN
                  DO 151 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     &                           COEFF((IMO-1)*NBAS+(NCNT-1)*6+IOFF+1)
  151             CONTINUE
               ENDIF
               IF(IANGMOM.EQ.6)THEN
                  DO 152 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     &                           COEFF((IMO-1)*NBAS+(NCNT-1)*6+IOFF+2)
  152             CONTINUE
               ENDIF
               IF(IANGMOM.EQ.7)THEN
                  DO 153 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     &                           COEFF((IMO-1)*NBAS+(NCNT-1)*6+IOFF+3)
  153             CONTINUE
               ENDIF
               IF(IANGMOM.EQ.8)THEN
                  DO 154 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     &                           COEFF((IMO-1)*NBAS+(NCNT-1)*6+IOFF+4)
  154             CONTINUE
               ENDIF
               IF(IANGMOM.EQ.9)THEN
                  DO 155 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     &                           COEFF((IMO-1)*NBAS+(NCNT-1)*6+IOFF+5)
  155             CONTINUE
               ENDIF
               IF(IANGMOM.EQ.10)THEN
                  DO 156 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     &                           COEFF((IMO-1)*NBAS+(NCNT-1)*6+IOFF+6)
  156             CONTINUE
               ENDIF
               IOFF=IATMOFF+NMOMAO((IATM-1)*NTANGM+1)+
     &                    3*NMOMAO((IATM-1)*NTANGM+2)+
     &                    6*NMOMAO((IATM-1)*NTANGM+5)
               IF(IANGMOM.EQ.11)THEN
                  DO 160 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     &                           COEFF((IMO-1)*NBAS+(NCNT-1)*10+IOFF+1)
  160             CONTINUE
               ENDIF
               IF(IANGMOM.EQ.12)THEN
                  DO 161 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     &                           COEFF((IMO-1)*NBAS+(NCNT-1)*10+IOFF+2)
  161             CONTINUE
               ENDIF
               IF(IANGMOM.EQ.13)THEN
                  DO 162 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     &                           COEFF((IMO-1)*NBAS+(NCNT-1)*10+IOFF+3)
  162             CONTINUE
               ENDIF
               IF(IANGMOM.EQ.14)THEN
                  DO 163 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     &                           COEFF((IMO-1)*NBAS+(NCNT-1)*10+IOFF+4)
  163             CONTINUE
               ENDIF
               IF(IANGMOM.EQ.15)THEN
                  DO 164 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     &                           COEFF((IMO-1)*NBAS+(NCNT-1)*10+IOFF+5)
  164             CONTINUE
               ENDIF
               IF(IANGMOM.EQ.16)THEN
                  DO 165 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     &                           COEFF((IMO-1)*NBAS+(NCNT-1)*10+IOFF+6)
  165             CONTINUE
               ENDIF
               IF(IANGMOM.EQ.17)THEN
                  DO 166 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     &                           COEFF((IMO-1)*NBAS+(NCNT-1)*10+IOFF+7)
  166             CONTINUE
               ENDIF
               IF(IANGMOM.EQ.18)THEN
                  DO 167 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     &                           COEFF((IMO-1)*NBAS+(NCNT-1)*10+IOFF+8)
  167             CONTINUE
               ENDIF
               IF(IANGMOM.EQ.19)THEN
                  DO 168 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     &                           COEFF((IMO-1)*NBAS+(NCNT-1)*10+IOFF+9)
  168             CONTINUE
               ENDIF
               IF(IANGMOM.EQ.20)THEN
                  DO 169 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     &                           COEFF((IMO-1)*NBAS+(NCNT-1)*10+IOFF+10)
  169             CONTINUE
               ENDIF
C Offset for g functions
c               IOFF=IATMOFF+NMOMAO((IATM-1)*NTANGM+1)+
c     &                    3*NMOMAO((IATM-1)*NTANGM+2)+
c                          6*NMOMAO((IATM-1)*NTANGM+5)+
c                         10*NMOMAO((IATM-1)*NTANGM+11)
C
            ICNT=ICNT+1
            NCNT=NCNT+1
  125       CONTINUE
  120    CONTINUE
      IATMOFF=IATMOFF+NAOATM(IATM)
  110 CONTINUE
C
C     Write the MO's to output in groups of four.
      IPRINT = 0
      IF(IPRINT.EQ.1)THEN
C
      WRITE(LUOUT,1010)
 1010 FORMAT(/'               OCCUPIED MOLECULAR ORBITALS')
C
      IF(IFLG.EQ.0)THEN
      WRITE(LUOUT,1020)
 1020 FORMAT(/'               ******* ALPHA BLOCK *******')
      ELSE
      WRITE(LUOUT,1025)
 1025 FORMAT(/'               ******* BETA BLOCK *******')
      ENDIF
C
      INBASP=NBASP/4
      DO 30 ICNT=1,INBASP
         IB=4*(ICNT-1)
         IE=4*(ICNT-1)+3
         WRITE(LUOUT,1030) (K,K=IB+1,IE+1)
         DO 30 J=1,NBAS
            WRITE(LUOUT,1040) ICNTR(J),(COEFFN(I*NBAS+J),I=IB,IE)
   30 CONTINUE
      IF(INBASP*4.NE.NBASP)THEN
         KNBASP=INBASP*4
         JNBASP=NBASP-KNBASP
         IF(JNBASP.EQ.1) WRITE(LUOUT,1031) (K,K=KNBASP+1,NBASP)
         IF(JNBASP.EQ.2) WRITE(LUOUT,1032) (K,K=KNBASP+1,NBASP)
         IF(JNBASP.EQ.3) WRITE(LUOUT,1033) (K,K=KNBASP+1,NBASP)
         DO 40 J=1,NBAS
          WRITE(LUOUT,1040) ICNTR(J),(COEFFN(I*NBAS+J),I=KNBASP,NBASP-1)
   40    CONTINUE
      ENDIF
C
      ENDIF
C
 1030 FORMAT(/'CENTER',6X,'MO',I4,8X,'MO',I4,8X,'MO',I4,8X,'MO',I4)
 1031 FORMAT(/'CENTER',6X,'MO',I4)
 1032 FORMAT(/'CENTER',6X,'MO',I4,8X,'MO',I4)
 1033 FORMAT(/'CENTER',6X,'MO',I4,8X,'MO',I4,8X,'MO',I4)
 1040 FORMAT(I4,4F14.5)
      RETURN
      END
