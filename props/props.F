      SUBROUTINE PROPS(ICORE,MAXCOR,IUHF)
      IMPLICIT INTEGER (A-Z)
C
      DOUBLE PRECISION ICORE
      DIMENSION ICORE(MAXCOR)
      CHARACTER*3  TYPE(3),FRONT
      CHARACTER*8  LABEL
      CHARACTER*80 LABEL1
      LOGICAL NONHF, RELAXED
      DOUBLE PRECISION ONE,TWO,SCFOCC
      INTEGER DIRPRD,JUNK(8),NOCCO(2),NVRTO(2)
      INTEGER NDRPOP(8),NDRVRT(8),POPFUL(8,2),VRTFUL(8,2),NBFIRR(8)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /INFO/   NOCA,NOCB,NVRTA,NVRTB
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/LUOUT,MOINTS 
C 
      DATA TYPE /'SCF','REL','EXT'/
      DATA ONE /1.0D0/,TWO/2.D0/
C
      IONE=1
      ITWO=2
      LUOUT=6
      CALL GETREC(20,'JOBARC','NBASTOT ',IONE,NORBS)
C
      IF(IFLAGS(2).GT.0)THEN
        CALL GETREC(20,'JOBARC','NBASCOMP',IONE,NCOMP)
      ELSE
        NCOMP = NORBS
      ENDIF
C-----------------------------------------------------------------------
C
CJDW/AP 3/19/96. Try to get props working properly for dropped core.
C
C     The next block code sets the values for NCOMP and NORBS. These
C     two variables are both the total number of molecular orbitals,
C     even in the dropped MO case. Note that 'NBASTOT ' contains the
C     number of active MOs in correlated calculations. It is written
C     in vtran and its meaning is therefore different from its original
C     meaning in vmol2ja.
C
C     Note that NOCCO, NVRTO are local arrays and are used to set NORBS
C     and NCOMP. Entries in /INFO/ are NOT changed.
C
C     The relaxed density in the MO basis has correlated order. In the
C     dropped MO case, the full matrix is stored and is ordered as
C     though no MOs were dropped. The density is transformed to the
C     AO basis using MOs from SCFEVCA0/SCFEVCB0. Although these records
C     are for SCF ordered MOs, full correlated order MOs are written
C     in vtran. Note how these records are read in MO2AO2. Formerly,
C     MO2AO3 read SCFEVECA/SCFEVECB, which was sensible, given that
C     these records are for correlated order MOs. However, for dropped
C     MOs, these records do not contain the dropped MOs. Hence, we mod-
C     ify MO2AO3 so SCFEVCA0/SCFEVCB0 are read. Of course, this depends
C     on vtran always writing these records, and one day we should
C     perhaps create records which contain all MOs in correlated order.
C
C     Note that vtran also changes 'NUMBASIR' !
C
C-----------------------------------------------------------------------
      CALL GETREC(20,'JOBARC','NUMBASIR',NIRREP,NBFIRR)
      CALL IZERO(NDRPOP,8)
      CALL IZERO(NDRVRT,8)
      CALL GETREC(20,'JOBARC','NUMDROPA',1,NDROP)
      IF(NDROP.NE.0) then
       CALL GETREC(20,'JOBARC','NDROPPOP',NIRREP,NDRPOP)
       CALL GETREC(20,'JOBARC','NDROPVRT',NIRREP,NDRVRT)
      ENDIF
C
      CALL IZERO(NOCCO,2)
      CALL IZERO(NVRTO,2)
      CALL GETREC(20,'JOBARC','SYMPOPOA',NIRREP,JUNK)
      DO 151 IRREP=1,NIRREP
       POPFUL(IRREP,1) = JUNK(IRREP)   + NDRPOP(IRREP)
       VRTFUL(IRREP,1) = NBFIRR(IRREP) + NDRVRT(IRREP) + NDRPOP(IRREP)
     &                 - POPFUL(IRREP,1)
       NOCCO(1)     = NOCCO(1)      + POPFUL(IRREP,1)
       NVRTO(1)     = NVRTO(1)      + VRTFUL(IRREP,1)
  151 CONTINUE
C
      IF(IUHF.NE.0)THEN
       CALL GETREC(20,'JOBARC','SYMPOPOB',NIRREP,JUNK)
       DO 251 IRREP=1,NIRREP
        POPFUL(IRREP,2) = JUNK(IRREP)   + NDRPOP(IRREP)
        VRTFUL(IRREP,2) = NBFIRR(IRREP) + NDRVRT(IRREP) + NDRPOP(IRREP)
     &                  - POPFUL(IRREP,2)
        NOCCO(2)     = NOCCO(2)      + POPFUL(IRREP,2)
        NVRTO(2)     = NVRTO(2)      + VRTFUL(IRREP,2)
  251  CONTINUE
      ELSE
       NOCCO(2)=NOCCO(1)
       NVRTO(2)=NVRTO(1)
       CALL ICOPY(NIRREP,POPFUL(1,1),1,POPFUL(1,2),1)
       CALL ICOPY(NIRREP,VRTFUL(1,1),1,VRTFUL(1,2),1)
      ENDIF
C
      NCOMP = NOCCO(1) + NVRTO(1)
      NORBS = NOCCO(1) + NVRTO(1)
C
      IF(IFLAGS(1).GT.0)THEN
       write(6,*) ' data ' , ncomp,norbs,noca,nvrta
       write(6,*) ' popful ',popful
       write(6,*) ' vrtful ',vrtful
       write(6,*) ' nocco, nvrto ', nocco,nvrto
       write(6,*) ' ndrpop, ndrvrt', ndrpop,ndrvrt
       write(6,*) ' nbfirr ',nbfirr
      ENDIF
C
C-----------------------------------------------------------------------
C
      OPEN(UNIT=10,FILE='VPOUT',FORM='UNFORMATTED',STATUS='OLD')
      NORBS2=NORBS*NORBS
      IONE=1
      ICBOT=1
      ITOP=2
      IGRN=0
      ISKIP=1
      IF(IFLAGS(18).EQ.0)ICBOT=2
      IF(IFLAGS(18).EQ.3)ITOP=2
c      IF(IFLAGS(91).NE.0)THEN
c       CALL GETREC(20,'JOBARC','EXTNTOT ',IONE,NROOT)
c       ITOP=3
c       ISKIP=2
c       ICBOT=1
c      ENDIF
      IF(IFLAGS(2).EQ.0)ITOP=1
C
C  CHECK THE FLAG FOR NON HARTREE-FOCK-CCSD RUNS, IN THESE
C  CASES THE DIRECT AND NON-DIRECT PART OF THE RELAXED 
C  DENSITY HAVE TO BE READ IN SEPARATELY
C    To accomodate property calculations with "RESPONSE" density, we need
C  to know the "DENSITY" flag. The above comment is somewhat misleading
C  since except for QRHF, the RELAXED density matrix already has both
C  direct and indirect parts.  (see below). Ajith Perera, 09/2004
C
      ITHRU=1
      NONHF=(IFLAGS(77).NE.0)
      RELAXED=(IFLAGS(19).EQ.0)
C
C CHECK IF THIS A POST-SCF GRNFNC RUN.  IF SO, THE SCF DENSITY IS NOT
C  ON DISK, AND WE'LL HAVE TO FIGURE IT OUT FROM THE RELAXED DENSITY.
C  KLUDGY, BY GRNFNC IS GOING TO BE DEEP-SIXED SOON ENOUGH.
C
      IF(IFLAGS(11).GE.3)IGRN=1
      IF(ITOP.EQ.2.AND.IFLAGS(11).GE.3)IGRN=1
      DO 1000 ICTYP=ICBOT,MAX(ICBOT,ITOP),ISKIP
       I000=1
       I010=I000+NORBS2
       I020=I010+NORBS2
       I030=I020+NORBS2
       I040=I030+NORBS2
       I050=I040+NORBS2
       I060=I050+NORBS2
15     FRONT=TYPE(ICTYP)
C
C LOAD DENSITY MATRIX INTO BOTTOM OF CORE, AND TRANSFORM TO THE
C  AO BASIS IF THIS IS FOR A CORRELATED DENSITY.
C
       IF(FRONT.EQ.'SCF')THEN
        CALL GETREC(20,'JOBARC',FRONT//'DENSA',
     &              IINTFP*NCOMP*NCOMP,ICORE)
       ELSEIF(FRONT.EQ.'REL')THEN
        CALL GETREC(20,'JOBARC',FRONT//'DENSA',
     &              IINTFP*NCOMP*NCOMP,ICORE)
        IF(NONHF.AND.RELAXED) THEN
         CALL GETREC(20,'JOBARC','NHF'//'DENSA',
     &               IINTFP*NCOMP*NCOMP,ICORE(I010))
         CALL SAXPY(NCOMP*NCOMP,ONE,ICORE(I010),1,ICORE,1)
        ENDIF
        IF(FRONT.EQ.'REL'.AND.NCOMP.NE.NORBS)THEN
           CALL EXPDEN(ICORE,ICORE(I010),NORBS,NCOMP,'IDROPA',
     &                 ICORE(I020),NOCA,IUHF)
        ENDIF
       ELSEIF(FRONT.EQ.'EXT')THEN
        LABEL='EXTDNA  '
        WRITE(LABEL(7:7),'(I1)')ITHRU
        CALL GETREC(20,'JOBARC',LABEL,NCOMP*NCOMP,ICORE)
       ENDIF
       IF(ICTYP.EQ.2)THEN
        CALL MO2AO2(ICORE,ICORE(I040),ICORE(I020),ICORE(I030),
     &            NORBS,1)
       ELSE
c YAU : old
c         CALL ICOPY(IINTFP*NORBS2,ICORE,1,ICORE(I040),1)
c YAU : new
          CALL DCOPY(NORBS2,ICORE,1,ICORE(I040),1)
c YAU : end
       ENDIF
C
C IF THIS IS UHF, PULL THE BETA DENSITY AND ADD IT.
C
       IF(IUHF.EQ.1)THEN
        IF(FRONT.EQ.'SCF')THEN
         CALL GETREC(20,'JOBARC',FRONT//'DENSB',
     &               IINTFP*NCOMP*NCOMP,ICORE)
        ELSEIF(FRONT.EQ.'REL')THEN
         CALL GETREC(20,'JOBARC',FRONT//'DENSB',
     &               IINTFP*NCOMP*NCOMP,ICORE)
         IF(NONHF.AND.RELAXED) THEN
            CALL GETREC(20,'JOBARC','NHF'//'DENSB',
     &                  IINTFP*NCOMP*NCOMP,ICORE(I010))
            CALL SAXPY(NCOMP*NCOMP,ONE,ICORE(I010),1,ICORE,1)
         ENDIF
         IF(FRONT.EQ.'REL'.AND.NCOMP.NE.NORBS)THEN
          CALL EXPDEN(ICORE,ICORE(I010),NORBS,NCOMP,'IDROPB',
     &                ICORE(I020),NOCB,IUHF)
         ENDIF
        ELSEIF(FRONT.EQ.'EXT')THEN
         LABEL='EXTDNB  '
         WRITE(LABEL(7:7),'(I1)')ITHRU
         CALL GETREC(20,'JOBARC',LABEL,NCOMP*NCOMP,ICORE)
        ENDIF
        IF(ICTYP.EQ.2)THEN
         CALL MO2AO2(ICORE,ICORE(I050),ICORE(I020),ICORE(I030),
     &               NORBS,2)
        ELSE
c YAU : old
c          CALL ICOPY(IINTFP*NORBS2,ICORE,1,ICORE(I050),1)
c YAU : new
           CALL DCOPY(NORBS2,ICORE,1,ICORE(I050),1)
c YAU : end
        ENDIF
C
C PUT ALPHA-BETA DENSITY INTO ICORE(I030) FOR SPIN DENSITY CALCULATIONS,
C  AND PUT FULL DENSITY (ALPHA+BETA) INTO ICORE(I040).  THEN COPY THIS
C  TO BOTTOM OF CORE.
C
COLD    CALL IZERO(ICORE(I030),NORBS2*IINTFP) 
        CALL ZERO(ICORE(I030),NORBS2) 
        CALL VADD(ICORE(I030),ICORE(I040),ICORE(I050),NORBS2,-1.D0)
        CALL VADD(ICORE(I040),ICORE(I040),ICORE(I050),NORBS2,1.D0)
       ENDIF
c YAU : old
c      CALL ICOPY(IINTFP*NORBS2,ICORE(I040),1,ICORE,1)
c YAU : new
       CALL DCOPY(NORBS2,ICORE(I040),1,ICORE,1)
c YAU : end
C
C NOW WE HAVE THE FULL DENSITY MATRIX AT THE BOTTOM OF CORE.
C   DO THE DIPOLE MOMENT.
C
C Get the number of real (Ghost atoms are not included) atoms
C
       CALL GETREC(20,'JOBARC','NATOMS  ',IONE,NATOM)
c
       IF(ICTYP.EQ.1)WRITE(LUOUT,2000)
       IF(ICTYP.EQ.2)WRITE(LUOUT,2001)
       IF(ICTYP.EQ.3)WRITE(LUOUT,2002)ITHRU
2000   FORMAT(T3,'@DRVPRP-I, Properties computed from ',
     &       'the SCF density matrix follow.')
2001   FORMAT(T3,'@DRVPRP-I, Properties computed from ',
     &       'the correlated density matrix follow.')
2002   FORMAT(T3,'@DRVPRP-I, Properties for excited state ',
     &        '#',I2,' follow.')
       CALL DRVDIP(ICORE,MAXCOR,NORBS)
       CALL DRVSEC(ICORE,MAXCOR,NORBS)
       CALL DRVQUD(ICORE,MAXCOR,NORBS)
       CALL DRVOCT(ICORE,MAXCOR,NORBS)
       CALL DRVREL(ICORE,MAXCOR,NORBS)
c
c calculate multipole moments in alternative way for checking
c
       if (ictyp .eq. 2)
     $      call calc_tmom(icore(i000), icore(i050), icore(i060),
     $      norbs, (norbs*(norbs+1))/2)
C
C DO SPIN DENSITY IF THIS IS A UHF CALCULATION.  DO THIS FOR
C   ALL ATOMS IN VMOL FILE FOR NOW.
C 
       I000 = 1 
       I040 = I000 + 4*NORBS2
       I050 = I040 +  NATOM
       I060 = I050 +  NATOM 
       IF (I060 .GT. MAXCOR) CALL INSMEM("@-PROPS",I060,MAXCOR)
       IF(IUHF.EQ.1)CALL DRVSPD(ICORE(I000),ICORE(I040),ICORE(I050),
     $                          NORBS, NATOM)
       IF(IUHF.EQ.1)CALL DRVAHF(ICORE(I000),ICORE(I040),NORBS, NATOM)
       CALL DRVEFG(ICORE(I000),ICORE(I040),NORBS, NATOM)
C       CALL DRVPLD(ICORE(I000), ICORE(I040), NORBS, NATOM)
        CALL DRVPOT(ICORE(I000), ICORE(I040), NORBS, NATOM)
c       IF(FRONT.EQ.'EXT'.AND.ITHRU.LT.NROOT)THEN
c        ITHRU=ITHRU+1
c        GOTO 15
c       ENDIF
1000  CONTINUE

#ifdef _DEBUG_LVLM
C This is added as a demonstration of what Zack was trying to do
      Call built_scf_orbdens(Icore(I000),Maxcor,Ncomp,Norbs,
     &                       Nirrep,Popful,Vrtful,Nocco,
     &                       Nvrto,Iuhf)
#endif 
C The orbital kintic energies at CC for Y. Park (01/2022). This employ
C CC natural orbitals. 

      IF (Itop .EQ. 2) Call Built_orb_kint_enegs(Icore(I000),
     &                                           Maxcor,Ncomp,
     &                                           Norbs,Nirrep,
     &                                           Popful,Vrtful,
     &                                           Nocco,Nvrto,Iuhf)
C
C PERFORM A MULLIKEN POPULATION ANALYSIS
C
      DO 1100 ICTYP=ICBOT,MAX(ICBOT,ITOP),ISKIP
C
C ALLOCATE MEMORY
C
       IDMOA=1
       IDMOB=1+IUHF*NCOMP*NCOMP
       IDAO=IDMOB+NCOMP*NCOMP
       ISCR=IDAO+NORBS*NORBS
C
       MXCOR=MAXCOR-((1+IUHF)*NCOMP*NCOMP+NORBS*NORBS)
C
       FRONT=TYPE(ICTYP)
       IF(ICTYP.EQ.1) THEN
C
C CONSTRUCT HF-MO DENSITY FROM OCCUPATION
C
C     The LENGTH logic is as follows. If 'SCFEVECA' exists, it
C     means that vtran has run, ie we are doing a correlated
C     calculation. If so, we prepare SCF density matrix in
C     correlated order so we can MOs in correlated order in MO2AO3.
C     If not, we form SCF density in SCF order.
C
        CALL GETREC(0,'JOBARC','SCFEVECA',LENGTH,ICORE)
        IF(LENGTH.GE.0) THEN
         SCFOCC=TWO
         IF(IUHF.EQ.1) SCFOCC=ONE
         CALL ZERO(ICORE(IDMOA),NCOMP*NCOMP)
         CALL FILLDHF(SCFOCC,ICORE(IDMOA),NCOMP,NOCCO(1))
         IF(IUHF.EQ.1) THEN
          CALL ZERO(ICORE(IDMOB),NCOMP*NCOMP)
          CALL FILLDHF(SCFOCC,ICORE(IDMOB),NCOMP,NOCCO(2))
         ENDIF
         ITYPE=1
        ELSE
         INBASIR=ISCR
         IOCCA=INBASIR+NIRREP
         IOCCB=IOCCA
         SCFOCC=TWO
         IF(IUHF.EQ.1) SCFOCC=ONE
         CALL ZERO(ICORE(IDMOA),NCOMP*NCOMP)
         CALL GETREC(20,'JOBARC','NUMBASIR',NIRREP,ICORE(INBASIR))
         CALL GETREC(20,'JOBARC','OCCUPYA ',NIRREP,ICORE(IOCCA))
         CALL FILLDHF2(SCFOCC,ICORE(IDMOA),NCOMP,ICORE(IOCCA),
     &                 ICORE(INBASIR))
         IF(IUHF.EQ.1) THEN
          CALL ZERO(ICORE(IDMOB),NCOMP*NCOMP)
          CALL GETREC(20,'JOBARC','OCCUPYB ',NIRREP,ICORE(IOCCB))
          CALL FILLDHF2(SCFOCC,ICORE(IDMOB),NCOMP,ICORE(IOCCB),
     &                  ICORE(INBASIR))
         ENDIF
         ITYPE=2
        ENDIF
        LABEL1='SCF'
        ILENGTH=3

C If HFSTABILITY=FOLLOW and a HF calculation then the above
C logic is creating a problem (because hfstab=follow requires running 
C xvtran).

        CALL GETREC(-20,'JOABRC',"STAB_TAB",1,ITEST)

        IF (IFLAGS(2) .EQ. 0 .AND. ITEST .EQ. 2) THEN
           CALL GETREC(20,'JOBARC',"SCF"//'DENSA',IINTFP*NCOMP*NCOMP,
     &                 ICORE(IDMOA))
           IF(IUHF.EQ.1) THEN
              CALL GETREC(20,'JOBARC',"SCF"//'DENSB',IINTFP*NCOMP*NCOMP,
     &                    ICORE(IDMOB))
           ENDIF
           ITYPE=2
        ENDIF 

       ELSE
C
C READ IN CORRELATED DENSITIES
C
        CALL GETREC(20,'JOBARC',FRONT//'DENSA',IINTFP*NCOMP*NCOMP,
     &              ICORE(IDMOA))
        IF(IUHF.EQ.1) THEN
         CALL GETREC(20,'JOBARC',FRONT//'DENSB',IINTFP*NCOMP*NCOMP,
     &               ICORE(IDMOB))
        ENDIF
        ITYPE=1
        IF(IFLAGS(2).EQ.1) THEN
         LABEL1='MBPT(2)'
         ILENGTH=7
        ELSE IF(IFLAGS(2).EQ.2) THEN
         LABEL1='MBPT(3)'
         ILENGTH=7
        ELSE IF(IFLAGS(2).EQ.3) THEN
         LABEL1='SDQ-MBPT(4)'
         ILENGTH=11
        ELSE IF(IFLAGS(2).EQ.4) THEN
         LABEL1='MBPT(4)'
         ILENGTH=7
        ELSE IF(IFLAGS(2).EQ.7) THEN
         LABEL1='UCCSD(4)'
         ILENGTH=8
        ELSE IF(IFLAGS(2).EQ.8) THEN
         LABEL1='CCD'
         ILENGTH=3
        ELSE IF(IFLAGS(2).EQ.9) THEN
         LABEL1='UCC(4)'
         ILENGTH=6
        ELSE IF(IFLAGS(2).EQ.10) THEN
         LABEL1='CCSD'
         ILENGTH=4
        ELSE IF(IFLAGS(2).EQ.11) THEN
         LABEL1='CCSD+T(CCSD)'
         ILENGTH=12
        ELSE IF(IFLAGS(2).EQ.21) THEN
         LABEL1='QCISD(T)'
         ILENGTH=8
        ELSE IF(IFLAGS(2).EQ.22) THEN
         LABEL1='CCSD(T)'
         ILENGTH=7
        ELSE IF(IFLAGS(2).EQ.23) THEN
         LABEL1='QCISD'
         ILENGTH=5
        ELSE IF(IFLAGS(2) .EQ. 42) THEN
         LABEL1='aCCSD(T)'
         ILENGTH=8
        ENDIF

       ENDIF
C
       CALL POPANAL(NORBS,NCOMP,ICORE(IDMOA),ICORE(IDAO),ICORE(ISCR),
     &              MXCOR,IUHF,1,LABEL1,ILENGTH,ITYPE)
       IF(IUHF.EQ.1) THEN
        CALL POPANAL(NORBS,NCOMP,ICORE(IDMOB),ICORE(IDAO),ICORE(ISCR),
     &               MXCOR,IUHF,2,LABEL1,ILENGTH,ITYPE)
       ENDIF
1100  CONTINUE
C
      call aces_fin
      STOP
      END
