      SUBROUTINE GENINT(ICORE,MAXCOR,IUHF,INTPCK)
C
C DRIVER FOR INTERMEDIATE FORMATION.  
C
C   INTPCK = 0  , DUMP PACKED INTEGRALS ONLY FOR W INTERMEDIATES (NO F).
C   INTPCK = 1  , DUMP ONLY QUADRATIC PART OF W INTERMEDIATES AND F.
C   INTPCK = 2  , DUMP FULL W AND F INTERMEDIATES (For standard CCSD or beyond)
C   INTPCK = 3  , DUMP FIRST ORDER W AND F ONLY.  (For Linearized CC methods)
C
CEND
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION TCPU,TSYS
      DOUBLE PRECISION C1,C2,C3,C4,FACT
C
      DIMENSION ICORE(MAXCOR)
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,LAMBDA,NONHF,ROHF4,ITRFLG,UCC
      EQUIVALENCE (METHOD,IFLAGS(2))
C
      COMMON/SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC
      COMMON/NHFREF/NONHF
      COMMON/ROHF/ROHF4,ITRFLG
      COMMON/FLAGS/IFLAGS(100)
      COMMON/SYM/POP(8,2),VRT(8,2),NT1AA,NT1BB,NF1AA,NF1BB,
     &           NF2AA,NF2BB
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      EQUIVALENCE(METHOD,IFLAGS(2))
C
C W INTERMEDIATES. SET COEFFICIENTS BASED ON CALCULATION TYPE.
C
      IF(IFLAGS(1).GE.10)THEN
         WRITE(6,100)
 100     FORMAT(T3,'@GENINT-I, Calculating intermediates.')
         CALL TIMER(1)
      ENDIF
C
      IF(INTPCK.EQ.2)THEN
         C1=1.0
         C2=0.0
         C3=0.0
         C4=-0.5
      ELSEIF(INTPCK.EQ.1)THEN
         C1=0.0
         C2=0.0
         C3=0.0
         C4=-0.5 
      ELSEIF(INTPCK.EQ.0.OR.INTPCK.EQ.3)THEN
         C1=1.0
         C2=0.0
         C3=0.0
         C4=0.0
      ENDIF
C
      FACT=1.D0
      IF(UCC) THEN
         FACT=0.5
         C4=C4*FACT
      ENDIF
C
      IF(METHOD.GT.9.AND.SING1.AND..NOT.QCISD.AND.INTPCK.NE.3)C3=1.0
C
C ---------------The first part of the W(mn,ij) begins Here -------------
C
C The <mn||ij> + 1/4 Sum Tau(ij,ef)<mn|ef> contribution to the W(mn,ij)
C intermediate. The other piece is calculated in T1INW2 (see below).
C
      CALL QUAD1(ICORE,MAXCOR,INTPCK,IUHF,FACT)
C
C ---------------The fist part of the W(mn,ij) ends Here -----------------
C
C--------The W(mbej) and parts of F(ae) intermediate begins Here ---------
C
      DO 333 ILIST = 54,59,2-IUHF
         CALL ZERSYM(ICORE,ILIST) 
 333  CONTINUE
C
C The T1RING is responsible for calculating Sum T(jf)*<mb||ef> - 
C Sum T(nb)<mn||ej> contributions to W(mbej) intermediates. Also, calculates 
C the Sum T(mf)<ma||fe>  contributions to the F(ae) intermediate.
C
      INCRF=0
      IF ((METHOD.GT.9.AND.SING1.AND.(.NOT.QCISD)).OR.
     &   (ROHF4.AND.INTPCK.NE.0)) THEN
C
C Non-zero T1 contributions to W intermediates for RHF/UHF ref., only if 
C CALCLEVEL is beyond CCSD, but not QCISD and SING1 which is true if there is
C non zero T1's available from a prior iteration or ROHF-MBPT(4) and 
C 
         CALL ZERLST(ICORE,NF2AA,1,1,1,92)
         IF (IUHF.NE.0) CALL ZERLST(ICORE,NF2BB,1,1,2,92)
         LAMBDA=.FALSE.
         CALL T1RING(ICORE,MAXCOR,IUHF,LAMBDA)
         INCRF = 1
      ENDIF
C
C The following calls calculate <mb||ej>-(1/2)*SUM [T(jn,fb)<mn||ef>] 
C contribution to W(mbej) intermediate. 

      CALL DWMBEJ(ICORE,MAXCOR,'ABAB',IUHF,C1,C2,C3,C4)
      CALL DWMBEJ(ICORE,MAXCOR,'ABBA',IUHF,C1,C2,C3,C4)
C
      IF(IUHF.EQ.1)THEN
         CALL DWMBEJ(ICORE,MAXCOR,'AAAA',IUHF,C1,C2,C3,C4)
         CALL DWMBEJ(ICORE,MAXCOR,'BBBB',IUHF,C1,C2,C3,C4)
         CALL DWMBEJ(ICORE,MAXCOR,'BABA',IUHF,C1,C2,C3,C4)
         CALL DWMBEJ(ICORE,MAXCOR,'BAAB',IUHF,C1,C2,C3,C4)
      ENDIF
C
      IF(IFLAGS(1).GE.10)THEN
       CALL TIMER(1)
       WRITE(6,101)TIMENEW
101    FORMAT(T3,'@GENINT-I, W intermediates required ',F9.3,
     &           ' seconds.')
      ENDIF
C
C ----------------------W(mb,ej) and parts of F(ae) ends here-----------
C
      IF (ROHF4.AND.INTPCK.EQ.0) RETURN
C
C  Linearized CC or ROHF-MBPT(4)
C
      IF(INTPCK.EQ.3.OR.ROHF4) THEN
C
         IF(NONHF)THEN
C
            CALL ZERLST(ICORE,NF1AA,1,1,1,91)
            CALL ZERLST(ICORE,NT1AA,1,1,1,93)
            IF(IUHF.NE.0)THEN
               CALL ZERLST(ICORE,NF1BB,1,1,2,91)
               CALL ZERLST(ICORE,NT1BB,1,1,2,93)
            ENDIF
C 
            IF(ROHF4) THEN
               CALL NHFFINT(ICORE,MAXCOR,IUHF,.TRUE.)
            ELSE
               CALL NHFFINT(ICORE,MAXCOR,IUHF,.FALSE.)
            ENDIF
C
            IF(.NOT.SING1)THEN
               CALL ZERLST(ICORE,NT1AA,1,1,1,90)
               CALL ZERLST(ICORE,NT1BB,1,1,2,90)
            ENDIF
            IF(.NOT.ROHF4)RETURN
C

         ELSEIF(SING1)THEN
C
            CALL ZERLST(ICORE,NF1AA,1,1,1,91)
            CALL ZERLST(ICORE,NF2AA,1,1,1,92)
            CALL ZERLST(ICORE,NT1AA,1,1,1,93)
C
            IF(IUHF.NE.0)THEN
               CALL ZERLST(ICORE,NF1BB,1,1,2,91)
               CALL ZERLST(ICORE,NF2BB,1,1,2,92)
               CALL ZERLST(ICORE,NT1BB,1,1,2,93)
            ENDIF
         ENDIF
C
         IF(.NOT.ROHF4)RETURN
C
      ENDIF
C
C ---- Pieces of F(AE), F(MI) AND F(ME) intermediates begin Here------------
C ---- Note that part of the F(AE) is already calculated in T1RING----------
C
      IF(INTPCK.EQ.0) RETURN  ! No need for any of the F intermediates.

      IF(IFLAGS(1).GE.10)CALL TIMER(1)
C
      IF(ROHF4) THEN
        NINCRF=INCRF
      ELSE
        NINCRF=0
      ENDIF
C
C Calculate 1/2 Sum Tau(in,ef)<mn|ef> contribution to the F(mi) intermediate
C      
      CALL QUAD2(ICORE,MAXCOR,IUHF,NINCRF,FACT)
C
C Calculate 1/2 Sum Tau(mn,af)<mn|ef> contribution to the F(ae) intermediate
C The other piece (Sum T(mf) <ma||fe>) is calculated in T1RING (see above).

      CALL QUAD3(ICORE,MAXCOR,IUHF,INCRF,FACT)
C
      IF((METHOD.GE.9.AND.SING1).OR.ROHF4)THEN
C 
C Calculate F(me) + Sum T(nf)<mn||ef> contribution to the F(me). That 
C is all for the total F(me) intermediate. 
C
         CALL MAKFME(ICORE,MAXCOR,IUHF,1)
         IF(IUHF.NE.0) CALL MAKFME(ICORE,MAXCOR,IUHF,2)
C 
      ENDIF
C
      IF(NONHF.AND..NOT.ROHF4) THEN
C
C Calculate (1-delta) f(ae) - 1/2 Sum f(me) t(am) and (1-delta)f(mi) - 1/2 Sum
C T(ie) f(me). These terms have non zero contributions only for NON-HF methods.
C ROHF-MBPT(4) does not need these contributions. 
C
        CALL NHFFINT(ICORE,MAXCOR,IUHF,SING1)
C
      ENDIF
C
      IF(IFLAGS(1).GE.10)THEN
       CALL TIMER(1)
       WRITE(6,102)TIMENEW
102    FORMAT(T3,'@GENINT-I, F intermediates required ',F9.3,
     &           ' seconds.')
      ENDIF
C
C------The F(AE) AND F(ME) and parts of F(MI) intermediates end Here -------
C
C----- Rest of the W(mnij) and F(MI) and part of W(abef) contribution to
C-----------------------begins Here -----------------------------------
C
      IF((METHOD.GT.9.AND.SING1.AND.(.NOT.QCISD)).OR.ROHF4) THEN
C
C Calculate Sum T(ej)<mn||ie> contribution to W(mn,ij) intermediate. Also,
C add the Sum T(ne)<mn||ie> contriution to the F(MI) intermediate. At this 
C point F(MI) intermediate is complete. 
C
       CALL T1INW2(ICORE,MAXCOR,IUHF)
C
C Calculate the Sum T(bm) T(ef,ij) <am||ef> contribution to T2. Note that 
C this is not a contributio to the W(ab,ef) intermediate, instead it is a 
C contribution to the T2(ij,ab). This can avoid storing W(ab,ef) on the disk
C at this point. 
C
       CALL T1INW1(ICORE,MAXCOR,IUHF)
C
      ENDIF
C-------------------------W(mn,ij) and F(MI) is complete----------------
C
      IF(IFLAGS(1).GE.10)THEN
       CALL TIMER(1)
       WRITE(6,103)TIMENEW
103    FORMAT(T3,'@GENINT-I, additional terms for CCSD required ',F9.3,
     &           ' seconds.')
      ENDIF
C
      RETURN
      END
