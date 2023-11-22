
C THIS PROGRAM TRANSFORMS THE T2 AMPLITUDES FROM THE MOLECULAR TO
C  ATOMIC ORBITAL BASIS.  CURRENTLY, IT ASSUMES THAT THE NUMBER OF
C  ACTIVE MOLECULAR ORBITALS IS THE SAME AS THE NUMBER OF BASIS FUNCTIONS.
C
C J.F. STANTON AND J. GAUSS, GAINESVILLE, 1991.

      PROGRAM BCKTRN
      IMPLICIT INTEGER (A-Z)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,INCORE
      LOGICAL SPHHRM,MOEQAO,ALASKA
      COMMON / / ICORE(1)
      COMMON /ISTART/ I0,ICRSIZ
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /INFO /  NOCCO(2),NVRTO(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /AOPOPS/ AOPOP(8),MOPOP(8),NAO,NAOSPH,NMO
      COMMON /AOOFST/ INDOCC(8,2),INDVRT(8,2)
      COMMON /BASTYP/ SPHHRM,MOEQAO
      COMMON /METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
C
      NNP1O2(I)=(I*(I+1))/2
C
      CALL CRAPSI(ICORE,IUHF,0)
C
C Check to see whether we are doing a gradient calculation using ALASKA. For
C ALASKA gradient calculations, MO gammas need only to be backtransformed to AO
C basis (not to Cartesian AO basis), 02/2002 Ajith Perera.
C
      ALASKA = (IFLAGS(56).EQ.4)
      CALL SETMET(ALASKA)
C
      IF (.NOT.MBPT2.AND.IFLAGS(35).NE.0) CALL INCOR(I0,ICRSIZ,IUHF)
      MAXCOR=ICRSIZ
C
C  FILL /AOPOPS/ WITH THE NUMBER OF BASIS FUNCTIONS PER IRREP
C
      CALL FAOPOP(ALASKA)
C
C  CALCULATE SEVERAL MEMORY REQUIREMENTS
C
      CALL MEMSIZ(N1A,N1B,N1C,N1D,N2A,N2B,N2C,N3A,N3B,N3C,N4A,
     &            N4B,N4C,NBUF,IUHF,N1AA,N2AA,N3AA,N4AA)
C
C  REORDER AND SYMMETRY PACK EIGEN VECTORS
C  ALPHA EIGENVECTORS FIRST ( RHF AND UHF)
C  
C
C SECOND ORDER
C
      IF(MBPT2) THEN
      I000=I0
      I010=I000+IINTFP*NAO*NMO
      I020=I010+IINTFP*NAO*NMO
      CALL GETREC(20,'JOBARC','SCFEVECA',NAOSPH*NMO*IINTFP,
     &            ICORE(I000))
C
C CALL ROUTINE WHICH FORMS THE SYMMETRY-ORDERED EIGENVECTOR MATRIX
C  AND THE SYMMETRY OFFSETS (STORED IN COMMON BLOCK AOOFST).
C
      IOFF0=0
      IF(SPHHRM .AND. .NOT. ALASKA)CALL SPH2CRT(ICORE(I000),
     &                             ICORE(I010),MAXCOR,IUHF,1)
      CALL SYPKEV(ICORE(I000),ICORE(I010),NAO,TOTLENA,IOFF0,1)
      I010=I000+IINTFP*TOTLENA
C
C  NOW REORDER BETA EIGENVECTORS ( UHF ONLY)
C
      IF(IUHF.NE.0)THEN
       IOFF0=IOFF0+TOTLENA
       CALL GETREC(20,'JOBARC','SCFEVECB',NAOSPH*NMO*IINTFP,
     &             ICORE(I010))
       I030=I020+IINTFP*NAO*NMO
       IF(SPHHRM .AND. .NOT. ALASKA)CALL SPH2CRT(ICORE(I010),
     &                              ICORE(I020),MAXCOR,IUHF,2)
       CALL SYPKEV(ICORE(I010),ICORE(I020),NAO,TOTLENB,IOFF0,2)
       I020=I010+IINTFP*TOTLENB
      ELSE
       CALL ICOPY(NIRREP,INDOCC(1,1),1,INDOCC(1,2),1)
       CALL ICOPY(NIRREP,INDVRT(1,1),1,INDVRT(1,2),1)
       I020=I010
      ENDIF
      ELSE
       I000=I0
       I010=I000+IINTFP*NAO*NMO
       I020=I010+IINTFP*NAO*NMO
       CALL GETREC(20,'JOBARC','SCFEVECA',NAOSPH*NMO*IINTFP,
     &             ICORE(I000))
       IF(SPHHRM .AND. .NOT. ALASKA)CALL SPH2CRT(ICORE(I000),
     &                                   ICORE(I010),MAXCOR,IUHF,1)
       CALL REORDC(ICORE(I000),ICORE(I010),NAO,POP(1,1),VRT(1,1),
     &             AOPOP,NMO)
       IF(IUHF.EQ.1) THEN
        I020=I010+IINTFP*NAO*NMO
        CALL GETREC(20,'JOBARC','SCFEVECB',NAOSPH*NMO*IINTFP,
     &              ICORE(I010))
        IF(SPHHRM .AND. .NOT. ALASKA)CALL SPH2CRT(ICORE(I010),
     &                                    ICORE(I020),MAXCOR,IUHF,2)
        CALL REORDC(ICORE(I010),ICORE(I020),NAO,POP(1,2),VRT(1,2),
     &              AOPOP,NMO)
       ELSE
        I020=I010
       ENDIF
        CALL GETOFS(ICORE(I000),ICORE(I010),IUHF,NAO)
      ENDIF
C
      I020REL=I020-I0+1
C
C SECOND-ORDER CODE
C
      IF(MBPT2) THEN
C
C RESORT THE T2-AMPLITUDES FROM DIRAC TO MULLIKEN  ORDER, SPIN ADAPT
C THE AMPLITUDES IN THE CASES OF RHF
C
      CALL PREHF(ICORE(I020),MAXCOR-I020REL+1,IUHF)
C
C INITIALIZE THE NEW MOINTS FILE WHICH HOLDS THE TRANSFORMED GAMMAS
C
      CALL SETLST
      MEMLEFT=I0+MAXCOR-I020
      IF (IFLAGS(35).NE.0) CALL INCOR(I020,MEMLEFT,IUHF)
CJDW 12/28/93. Redefinition of MAXSIZE and address of SCR. SCR now starts
C              at I050. Previously it started at I020, which sometimes
C              caused BUF2 and BUF3 to be trashed following the call to
C              TRANX.
Cwas  MAXSIZE=(MAXCOR-I020REL+1)/IINTFP
      I030=I020+N1A*IINTFP
      I040=I030+N1B*IINTFP
      I050=I040+N1C*IINTFP
      IEND=I050+NBUF*IINTFP-I0+1
      INCORE=.TRUE.
      IF(IEND.GE.MAXCOR)THEN
       INCORE=.FALSE.
       N1A=N1D
       I030=I020+N1D*IINTFP
       I040=I030+N1B*IINTFP
       I050=I040+N1C*IINTFP
       IEND=I050+NBUF*IINTFP-I0+1
       IF(IEND.GE.MAXCOR)THEN
        CALL INSMEM('BCKTRN',IEND,MAXCOR)
       ENDIF
      ENDIF
      I050REL = I050 - I0 + 1
      MAXSIZE = (MAXCOR-I050REL+1)/IINTFP
      CALL BXAAAA(ICORE(I000),ICORE(I020),ICORE(I030),ICORE(I040),
     &            ICORE(I050),N1A,N1B,N1C,NBUF,IUHF,ICORE(I050),
     &            MAXSIZE,INCORE)
C
      IF(NIRREP.GE.2)THEN
C
       I030=I020+N2A*IINTFP
       I040=I030+N2B*IINTFP
       I050=I040+N2C*IINTFP
       IEND=I050+NBUF*IINTFP-I0+1
       IF(IEND.GE.MAXCOR) CALL INSMEM('BXAABB',IEND,MAXCOR)
       CALL BXAABB(ICORE(I000),ICORE(I020),ICORE(I030),
     &             ICORE(I040),ICORE(I050),N2A,N2B,N2C,NBUF,IUHF)
C
C TRY FIRST ALLOCATION FOR FULL INCORE
C
       I030=I020+N3A*IINTFP
       I040=I030+N3B*IINTFP
       I050=I040+N3C*IINTFP
       IEND=I050+NBUF*IINTFP-I0+1
       IF(IEND.GE.MAXCOR) THEN
C
C THERE IS NOT SUFFICIENT MEMORY FOR FULL INCORE, ALLOCATE
C MEMORY FOR OUT OF CORE
C 
        I030=I020
        I040=I030+N3B*IINTFP
        I050=I040+N3C*IINTFP
        ISTART1=I050+NBUF*IINTFP
        MEMLEFT=(MAXCOR-(ISTART1-I0))
        IF(MEMLEFT.LT.N3AA) THEN
C
C THERE IS NOT SUFFICIENT MEMORY FOR OUT OF CORE, ERROR EXIT
C
         CALL INSMEM('BXABAB', N3AA*IINTFP, MEMLEFT)
C
        ENDIF
C
        CALL BXABAB1(ICORE(I000),ICORE(ISTART1),ICORE(I030),
     &               ICORE(I040),ICORE(I050),MEMLEFT,N3B,N3C,
     &               NBUF,IUHF)
C 
       ELSE
C
C FULL INCORE ALGORITHM
C
        CALL BXABAB(ICORE(I000),ICORE(I020),ICORE(I030),
     &              ICORE(I040),ICORE(I050),N3A,N3B,N3C,NBUF,IUHF)
C
       ENDIF
C
       ENDIF
C
       IF(NIRREP.GE.4)THEN
        I030=I020+N4A*IINTFP
        I040=I030+N4B*IINTFP
        I050=I040+N4C*IINTFP
        IEND=I050+NBUF*IINTFP-I0+1
        IF(IEND.GE.MAXCOR) THEN
C
C THERE IS NOT SUFFICIENT MEMORY FOR FULL INCORE, ALLOCATE
C MEMORY FOR OUT OF CORE
C
         I030=I020
         I040=I030+N4B*IINTFP
         I050=I040+N4C*IINTFP
         ISTART1=I050+NBUF*IINTFP
         MEMLEFT=(MAXCOR-(ISTART1-I0))
         IF(MEMLEFT.LT.N4AA) THEN 
C
C THERE IS NOT SUFFICIENT MEMORY FOR OUT OF CORE, ERROR EXIT  
C
          CALL INSMEM('BXABCD', N4AA*IINTFP, MEMLEFT)
C
         ENDIF

         CALL BXABCD1(ICORE(I000),ICORE(ISTART1),ICORE(I030), 
     &                ICORE(I040),ICORE(I050),MEMLEFT,N4B,N4C,
     &                NBUF,IUHF)
C
        ELSE 
C
         CALL BXABCD(ICORE(I000),ICORE(I020),ICORE(I030),
     &               ICORE(I040),ICORE(I050),N4A,N4B,N4C,NBUF,IUHF)
        ENDIF
C
       ENDIF
C
      ELSE
C
C  ALL METHODS WHICH ARE NOT MBPT(2)
C
       MAXSIZE=(MAXCOR-I020REL+1)/IINTFP
       I030=I020+N1A*IINTFP
       I040=I030+N1B*IINTFP
       IEND=I040+N1C*IINTFP-I0+1
       IF(IEND.LE.MAXCOR)THEN
          INCORE=.TRUE.
       ELSE
          INCORE=.FALSE. 
          N1A=N1C
          I030=I020+N1A*IINTFP
          I040=I030+N1B*IINTFP
          IEND=I040+N1C*IINTFP-I0+1
C       
C No longer needed. See comments in BXAAAA2 and 
C TRNLST, Ajith Perera, 11/2001. 
C 
cSSS    IF (.NOT. MOEQAO) THEN
cSSS    Write(6, *) "@BCKTRN Not enough memory for incore algorithm"
cSSS    CALL INSMEM('BXAAAA2', IEND, MAXCOR)
cSSS    ENDIF
C 
       ENDIF 
C
       CALL BXAAAA2(ICORE(I000),ICORE(I020),ICORE(I030),ICORE(I040),
     &               N1A,N1B,N1C,IUHF,ICORE(I020),MAXSIZE,INCORE,ALASKA)
C
       IF(NIRREP.GE.2)THEN
        I030=I020+N2A*IINTFP
        I040=I030+N2B*IINTFP
        IEND=I040+N2C*IINTFP-I0+1
        IF(IEND.GE.MAXCOR) CALL INSMEM('BXAABB2',IEND,MAXCOR)
        CALL BXAABB2(ICORE(I000),ICORE(I020),ICORE(I030),ICORE(I040),
     &               N2A,N2B,N2C,IUHF,ALASKA)
C
        I030=I020+N3A*IINTFP
        I040=I030+N3B*IINTFP
        IEND=I040+N3C*IINTFP-I0+1
        IF(IEND.GE.MAXCOR) CALL INSMEM('BXABAB2',IEND,MAXCOR)
        CALL BXABAB2(ICORE(I000),ICORE(I020),ICORE(I030),ICORE(I040),
     &               N3A,N3B,N3C,IUHF,ALASKA)
C
       ENDIF
C
       IF(NIRREP.GE.4)THEN
        I030=I020+N4A*IINTFP*2
        I040=I030+N4B*IINTFP*2
        IEND=I040+N4C*IINTFP*2-I0+1
        IF(IEND.GE.MAXCOR) CALL INSMEM('BXABCD2',IEND,MAXCOR)
        CALL BXABCD2(ICORE(I000),ICORE(I020),ICORE(I030),ICORE(I040),
     &               N4A,N4B,N4C,IUHF,ALASKA)
       ENDIF
C
      ENDIF
      CALL SUMGAM
      IF (IFLAGS(35).NE.0) THEN
         CALL ACES_AUXCACHE_FLUSH
         CALL ACES_AUXCACHE_RESET
      END IF
      CALL ACES_IO_REMOVE(51,'GAMLAM') 
      call aces_fin  
      STOP
      END
