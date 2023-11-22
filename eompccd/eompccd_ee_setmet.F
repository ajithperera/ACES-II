
#include "flags.h"

      SUBROUTINE EOMPCCD_EE_SETMET(ICORE,IUHF,pCCD)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*8 STRING
      LOGICAL CIS,EOMCC,FULDIAG,INCORE,READGUES,DOUBLE,ESPROP,NONSTD
      LOGICAL EMINFOL,EVECFOL,VPROP,LAMTHERE,GABCD,CISD
      LOGICAL MBPT2,CC, SS, SD, DS, DD, NODAVID, TRIPLET,CCD,rCCD,
     +        drCCD,RPA
      LOGICAL LCCD, LCCSD,CC2,CCSD
      LOGICAL DAVID,BLOCK_DAVID,LANCZOS
      LOGICAL pCCD

      CHARACTER*11 STRING1(6)
      CHARACTER*6  STRING2(4)
                                
      DIMENSION ICORE(*)
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /METH/ CIS,RPA,EOMCC,CISD,FULDIAG,INCORE,READGUES
      COMMON /EXTRAP/ MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /REFTYPE/ MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
      COMMON /PROPGRAD/ ESPROP,IOPTROOT,IOPTSYM
      COMMON /CALCINFO/ NROOT(8)
      COMMON /GUESS/ DOUBLE,NONSTD
      COMMON /GUESS2/ IMAP(100,8)
      COMMON /CNVRGE/ EMINFOL,EVECFOL
      COMMON /LAMSTATE/ LAMTHERE
      COMMON /SPINSTATE/TRIPLET
      COMMON /VDINTPRT/NTPERT,NPERT(8),KPERT(8),IDIPSM(3),
     &                 IYZPERT,IXZPERT,IXYPERT,ITRANSX,
     &                 ITRANSY,ITRANSZ,NUCIND
      COMMON /VDINTLEN/ LENS(8),LENT(8)
      COMMON /INTPROG/ VPROP
      COMMON /ABCD/ GABCD
      COMMON /LISTPROJ/ LISTH0, ICOLPR1, ICOLPR2
      COMMON /PROJECT/ IPROJECT, IPATTERN, NCALC, ICALC, IWINDOW(8)
      COMMON /DRVHBAR/ SS, SD, DS, DD
      COMMON /TDALIST/ LISTETDA, LISTVTDA
      COMMON /PARTEOM/ NODAVID
      COMMON /ROOTSEARCH/DAVID,BLOCK_DAVID,LANCZOS 

      DATA STRING1 /"NO", "SEARCH_ONLY","PRJCT_ALL","PRJCT_NOISE",
     &              "CORE","NTO"/
      DATA STRING2 /"LOWEST", "CORE", "LUMO", "HOMO"/

      CALL IZERO(NROOT,8)
      LISTETDA = 95
      LISTVTDA = 94
C
       WRITE(6,*)
       WRITE(6,*)'      **********************************************'
       WRITE(6,*)'      *                                            *'
       WRITE(6,*)'      * PCCD EXCITATION ENERGY EOMCC CALCULATION   *'
       WRITE(6,*)'      *                                            *'
       WRITE(6,*)'      *            CODED BY A. PERERA              *'
       WRITE(6,*)'      *                     Z. WINDOM              *'
       WRITE(6,*)'      **********************************************'
       WRITE(6,*)
C
C
C PROCESS INPUT FOR EXCITATION ENERGY CALCULATION. IFLAGS(87) refers to EXCITE keyword 
C ILFAGS(2) is CALCLEVEL.

      CIS=.FALSE.
      IF(IFLAGS(87).EQ.1.OR.IFLAGS(87).EQ.5.OR.IFLAGS(87).EQ.6) THEN
       CIS=.TRUE.
      ENDIF

      CISD =IFLAGS(87).EQ.6
      RPA  = (IFLAGS(2).EQ. 48 .AND. IFLAGS(87).EQ.17) .OR. 
     &       (IFLAGS(2).EQ. 49 .AND. IFLAGS(87).EQ.18) .OR.
     &       (IFLAGS(2).EQ. 48 .AND. IFLAGS(87).EQ.15) .OR.
     &       (IFLAGS(2).EQ. 49 .AND. IFLAGS(87).EQ.13) .OR.
     &       (IFLAGS(2).EQ. 49 .AND. IFLAGS(87).EQ.19) .OR.
     &       (IFLAGS(2).EQ. 49 .AND. IFLAGS(87).EQ.20) 

      EOMCC=IFLAGS(87).EQ.3  .OR. IFLAGS(87) .EQ. 7  .OR.
     &      IFLAGS(87).EQ.8  .OR. IFLAGS(87) .EQ. 9  .OR. 
     &      IFLAGS(87).EQ.4  .OR. IFLAGS(87) .EQ. 10 .OR. 
     &      IFLAGS(87).EQ.11 .OR. IFLAGS(87) .EQ. 12 .OR. 
     &      IFLAGS(87).EQ.14 .OR. IFLAGS(87) .EQ. 16 
C
C  determine reference state (IFLAGS2(117) referes to EOMREF)
C
      IF (EOMCC) THEN

        CCSD  = IFLAGS2(117) .EQ. 1
        MBPT2 = IFLAGS2(117) .EQ. 2
        CCD   = IFLAGS2(117) .EQ. 3
        RCCD  = (IFLAGS2(117).EQ. 10 .OR. IFLAGS(2) .EQ. 48)
        DRCCD = (IFLAGS2(117).EQ. 11 .OR. IFLAGS2(117) .EQ. 15
     &          .OR. IFLAGS(2) .EQ. 49)
        LCCD  = IFLAGS2(117) .EQ. 6
        LCCSD = IFLAGS2(117) .EQ. 7 
    
        CC2   = IFLAGS2(117) .EQ. 9 
C
        CC = .NOT. (MBPT2 .OR. LCCD .OR. LCCSD .OR. RCCD .OR. DRCCD)

        MBPT2 =  MBPT2 .OR. RCCD .OR. DRCCD 
C 10/2021. 
        CALL GETREC(0,"JOBARC","PCCD_RUN",ILENGTH,JUNK)
        pCCD = .FALSE.
        IF (ILENGTH .GT. 0) pCCD = .TRUE.
        IF (pCCD .AND. IFLAGS2(117) .NE. 1) THEN
           CCD = .TRUE.
        ELSE 
           CCSD = .TRUE.
        ENDIF 

        IF (CC .AND. .NOT. 
     &     (IFLAGS2(117) .EQ. 1 .OR. IFLAGS2(117) .EQ.9 .OR.PCCD)) THEN
          WRITE(6,*) ' SOMETHING STRANGE IN SETMET '
        ENDIF

      ELSE IF (RPA) THEN

        RCCD  = (IFLAGS2(117) .EQ. 10 .OR. IFLAGS(2) .EQ. 48)
        DRCCD = (IFLAGS2(117) .EQ. 11 .OR. IFLAGS2(117) .EQ. 15 
     &          .OR. IFLAGS(2) .EQ. 49)
      ENDIF 
C
      ESPROP=IFLAGS(91).NE.0
C
      IF(CIS)THEN
        IF(IFLAGS(87).EQ.1) THEN
          WRITE(6,100)
  100     FORMAT(T3,'TDA excitation energies will be evaluated.')
        ELSE IF(IFLAGS(87).EQ.5) THEN
          WRITE(6,101)
  101     FORMAT(T3,'CIS excitation energies will be evaluated.')
        ELSEIF(IFLAGS(87).EQ.6) THEN
          WRITE(6,102)
  102     FORMAT(T3,'CIS(D) excitation energies will be evaluated.')
        ENDIF

      ELSEIF(EOMCC)THEN

        IF (CCSD) THEN
        WRITE(6,103)
  103   FORMAT(T3,'EOM-CCSD excitation energies will be evaluated.')
        ELSE IF (CC2) THEN
        WRITE(6,104)
  104   FORMAT(T3,'EOM-CC2 excitation energies will be evaluated.')
        ELSE IF (LCCD) THEN
        WRITE(6,105)
  105   FORMAT(T3,'EOM-LCCD excitation energies will be evaluated.')
        ELSE IF (LCCSD) THEN
        WRITE(6,106)
  106   FORMAT(T3,'EOM-LCCSD excitation energies will be evaluated.')
        ELSE IF (DRCCD .AND. (IFLAGS(87).EQ.13)) THEN
        WRITE(6,107)
  107   FORMAT(T3,'EOM(S)-DRCCD excitation energies will be evaluated.')
        ELSE IF (DRCCD .AND.(IFLAGS(87).EQ.14)) THEN
        WRITE(6,108)
  108   FORMAT(T3,'EOM(D)-DRCCD excitation energies will be evaluated.')
        ELSE IF (RCCD .AND.(IFLAGS(87).EQ.15)) THEN
        WRITE(6,109)
  109   FORMAT(T3,'EOM(S)-RCCD excitation energies will be evaluated.')
        ELSE IF (RCCD .AND.(IFLAGS(87).EQ.16)) THEN
        WRITE(6,110)
  110   FORMAT(T3,'EOM(D)-RCCD excitation energies will be evaluated.')
        ELSE IF (pCCD .AND.  IFLAGS2(117) .NE. 1) THEN
        WRITE(6,115)
  115   FORMAT(T3,'EOM-pCCD/CCD excitation energies will be evaluated.')
        ELSE IF (pCCD .AND. IFLAGS2(117) .EQ. 3) THEN
        WRITE(6,116)
  116   FORMAT(T3,'EOM-pCCD/CCSD excitation energies will be evaluated.'
     +         )
        ENDIF 
      ELSEIF (RPA) THEN

          IF (IFLAGS(2) .EQ. 48 .AND. IFLAGS(87) .EQ. 17) THEN
          WRITE(6,111)
  111     FORMAT(T3,'EOM(Sf)-rCCD (RPA) excitation energies will be'
     +           ' evaluated.')
          ENDIF 
          IF (IFLAGS(2) .EQ. 49 .AND. IFLAGS(87) .EQ. 18) THEN
          WRITE(6,112)
  112     FORMAT(T3,'EOM(Sf)-drCCD (dRPA) excitation energies will be'
     +           ' evaluated.')
          ENDIF 
         IF (IFLAGS(2) .EQ. 48 .AND. IFLAGS(87) .EQ. 15) THEN
          WRITE(6,113)
  113     FORMAT(T3,'EOM(S)-rCCD excitation energies will be'
     +           ' evaluated.')
          ENDIF
          IF (IFLAGS(2) .EQ. 49 .AND. IFLAGS(87) .EQ. 13) THEN
          WRITE(6,114)
  114     FORMAT(T3,'EOM(S)-drCCD excitation energies will be'
     +           ' evaluated.')
          ENDIF

      ENDIF 

c check first if this is some kind of continuation of a single root
c calculation
c
      ione=1
      nbas=nocco(1)+nvrto(1)
      call getrec(-1,'JOBARC','PRINSPIN',IONE,ICHECK)
      if(icheck.ne.0)then
       call fixsym(isym,nbas,icore,icore(1+nbas))
       nroot(isym)=1
       ntot=1
      else
       CALL GETREC(-1,'JOBARC','EESYMINF',NIRREP,NROOT)

C DO not set IPATTERN if IFLAGS2(119)=4 (i.e. ee_search=double), Ajith Perera,
C 10/2021

       IF (IFLAGS2(119) .NE. 4) THEN
          IPATTERN = IFLAGS2(119)
       ENDIF 
       IF (IPATTERN .EQ. 2) CALL INITNUMR(ICORE,IUHF,NROOT)
       NTOT=0
       DO 10 I=1,NIRREP
        NTOT=NTOT+NROOT(I)
   10  CONTINUE
      endif
C
C SEE IF SPECIAL INPUT EXISTS
C
      NONSTD=.FALSE.
      OPEN(UNIT=30,FILE='ZMAT',FORM='FORMATTED')
1     READ(30,'(A)',END=99)STRING

      IF(INDEX(STRING,'%EXCITE').EQ.0.AND.
     &   INDEX(STRING,'%excite').EQ.0)GOTO 1

      NONSTD=INDEX(STRING,'*').NE.0
  99  CONTINUE
c
c hard coded stuff for now until keyword input available
c
      MAXEXP=IFLAGS(97)
      NREDUCE=MIN(12,IFLAGS(97)-3)
      NTOL=IFLAGS(98)
      READGUES=.TRUE.
      INCORE=.FALSE.
      FULDIAG=.FALSE.
      DOUBLE=IFLAGS(11).EQ.0.AND..NOT.NONSTD.AND.ICHECK.EQ.0
C
C   AT PRESENT DOUBLE IS NOT USED, AND PUT TO .FALSE. IN NEWGES.
C   LET US PUT IT TO .FALSE. NOW.
C 
      DOUBLE = .FALSE.
C
C SET UP IMAP VECTOR - THERE ARE TWO POSSIBLE TYPES OF NONSTANDARD
C INPUT AND THEY MUST BE DEALT WITH HERE!
C
      IF(.NOT.NONSTD)THEN
       IOMAX=0
       DO 30 IRREP=1,NIRREP
        DO 31 I=1,NROOT(IRREP)+3
         IMAP(I,IRREP)=I
31      CONTINUE
30     CONTINUE
C
C SETUP FOR MORE COMPLICATED INPUT
C
      ELSE
       I000=1
       I010=I000+IINTFP*2*NBAS
       I011=I010+3*NBAS*NBAS
       I012=I011+3*NBAS*NBAS
       I020=I012+3*NBAS*NBAS
       I030=I020+2*NBAS
       I040=I030+2*NBAS
       CALL PARSEINP(NBAS,ICORE(I000),ICORE(I010),ICORE(I011),
     &               ICORE(I012),ICORE(I020),ICORE(I030),ICORE(I040))
      ENDIF
C
      NTOTAL = 0
      DO IRREP = 1, NIRREP
        NTOTAL = NTOTAL + NROOT(IRREP)
        IF(NROOT(IRREP).NE.0)THEN
         IOMAX=MAX(IOMAX,NROOT(IRREP)+3)
        ELSE
         IOMAX=0
        ENDIF
      ENDDO
C
C
      IF(DOUBLE)CALL PUTREC(20,'JOBARC','STRTGUES',IOMAX,ICORE)
C
      CLOSE(UNIT=30,STATUS='KEEP')
C
C WRITE BLANK JOBARC RECORDS FOR DIAGONAL DENSITIES
C
      IONE=1
      CALL PUTREC(20,'JOBARC','EXTNTOT ',IONE,NTOT)
      CALL PUTREC(20,'JOBARC','ESTATPOP',NIRREP,NROOT)
C
C SG 1/24/97 Set IRREPX to 1, since it does not matter what it is.
C
      IRREPX = 1
      CALL GETAOINF(IUHF,IRREPX)
      CALL GETREC(-1,'JOBARC','NSYMPERT',27,NTPERT) 
      CALL GETREC(-1,'JOBARC','SYMMLENG',16,LENS)
      VPROP=NTPERT.EQ.0.OR.IFLAGS(91).EQ.2
C
      IF(VPROP)THEN
       WRITE(6,6000)
      ELSE
       WRITE(6,6001)
      ENDIF
C
C CHECK TO SEE IF THIS IS A GRADIENT OR FREQUENCY CALCULATION.  IF
C SO, THEN THE LAMBDA STATE IS NOT AVAILABLE AND WE SHOULD SKIP THE
C CALCULATION OF TRANSITION PROPERTIES AND GROUND STATE PROPERTIES.
C
      IF(IFLAGS(87).GT.1)THEN
       IF(IFLAGS2(h_IFLAGS2_geom_opt).NE.0.OR.IFLAGS(54).NE.0)THEN
        LAMTHERE=.FALSE.
       ELSE
        LAMTHERE=.TRUE.
       ENDIF
      ENDIF
C
      GABCD=IFLAGS(100).EQ.0
      IF(CIS) GABCD=.FALSE.
C
C If more than one root is desired, a gradient is only calculated for
C the last one.
C
      WRITE(6,6004)NTOTAL
      IF(NTOTAL.NE.1)THEN
       IF(NTOTAL.EQ.0)THEN
        WRITE(6,6002)
        CALL ERREX
       ELSE
        IF(IFLAGS(91).GT.1)THEN
         WRITE(6,6003)  
        ENDIF
       ENDIF
      ENDIF
C
C  LOGIC TO DETERMINE COMMON BLOCK /DRVHBAR/
C
      SS = .TRUE.
      SD = .TRUE.
      DS = .TRUE.
C
C EOM triples are introduced. DD is false only for partitioned
C and EOM-BWPT2 methods, 03/2014. Ajith Perera

      IF (IFLAGS(87) .EQ. 7 .OR. IFLAGS(87) .EQ. 8) THEN
        DD = .FALSE.
      ELSE
        DD = .TRUE.
      ENDIF
C
C CHECK HBARLISTS
C
      IF (DD) THEN
        IF (IFLAGS2(122) .EQ. 2) THEN
          WRITE(6,*) ' HBARABCD SHOULD CORRESPOND TO BARE INTEGRALS'
          CALL ERREX
        ENDIF
        IF (CC .AND. IFLAGS2(123) .EQ. 2) THEN
          WRITE(6,*) ' T1*ABCD SHOULD NOT BE INCLUDED IN HBARABCI'
          CALL ERREX
        ENDIF
      ENDIF
C
C IF IUHF = 1, BUT IT IS REALLY A CLOSED SHELL CALCULATION
C (POP(I,1) = POP(I,2))
C THEN WE WANT TO CALCULATE ONLY TRIPLET COMPONENTS. WE CAN RESTRICT OURSELVES
C TO TRIPLETS BY A SUITABLE CHOICE OF STARTING VECTORS.
C
       TRIPLET =.FALSE.
       IF (IUHF .EQ. 1) THEN
         IDIF = 0
         DO IRREP = 1, NIRREP
           IDIF = IDIF + ABS(POP(IRREP,1) - POP(IRREP,2))
         ENDDO
         TRIPLET = (IDIF .EQ. 0)
         TRIPLET = .FALSE. 
         IF (TRIPLET) THEN
           WRITE(6,*) ' TDA search is restricted to triplet states'
         ENDIF
       ENDIF
C
C  determine nodavid (nodavid for p-EOMEE,EOM-BWPT2) 
C
      NODAVID = (IFLAGS(87) .EQ. 7 .OR. IFLAGS(87) .EQ. 8)
      IF (NONSTD) NODAVID = .FALSE.
      IF (NODAVID) THEN
        NREDUCE = 0
        MAXEXP  = 0
      ENDIF

C ee_prjct options { "NO","SEARCH_ONLY","PRJCT_ALL","PRJCT_NOISE",
C                    "CORE",NTO}
C ee_search_options{"LOWEST","CORE","LUMO","HOMO",DOUBLE}

C
C   LOGIC TO DETERMINE COMMON BLOCK /PROJECT/
C
      IPROJECT = Iflags2(h_IFLAGS2_eom_prjct)
      LISTH0 = 472
      ICOLPR1 = 5
      ICOLPR2 = 6

C Do not set IPATTERN for double excited state searches since this 
C does not involve any projections (not at the moment), Ajith Perera,10/21

      IF (Iflags2(h_IFLAGS2_ee_search) .NE. 4) THEN
         IPATTERN =  Iflags2(h_IFLAGS2_ee_search)
      ENDIF 
C
C
C Projection must be ensured when IPROJECT=5 (NTO)
C
      IF (IPATTERN .EQ. 0 .AND. .NOT. (IPROJECT .EQ. 5)) IPROJECT = 0
C
C If the core window is set, it is assuemed that core-valence separation
C is used. The IPROJECT is set to 4 (to indicate) core. When the core
C window is not set but if ee_search is set to core, the code follows
C Marcel's implementation and IPROJECT is set to 2. To best of my 
C knowledge Marcel scheme is very much like core-valence separation
C business. However, I do not believe Marcel made lot of noise about
C this around 1996 (contrast this to 2013 noise from Koch et. al.).
C Ajith Perera.

      CALL GETREC(0,"JOBARC","EEWINDOW",LENGTH,JUNK)
      IF (LENGTH .GT. 0) THEN
         CALL GETREC(20,"JOBARC","EEWINDOW",NIRREP,IWINDOW)
         IPROJECT = 4
      ELSE
         IWINDOW(1)=0
         IF (IPROJECT .NE. 1) IPROJECT = 2
      ENDIF

C If the user specified the core site he wants, do not try to automatically
C do it. 

      IF (.NOT. NONSTD) THEN 
         CALL DETNCALC(ICORE,IUHF)
      ELSE
         NCALC=1
      ENDIF 

      IF (LENGTH .GT. 0) THEN 
      Write(6,"(a,8I2)") "  The core window size      : ", 
     +                   (iwindow(i),i=1,Nirrep)
      Write(6,"(a,I2)") "  The number of calculations: ", ncalc 
       
      ENDIF 
C
C  write out some stuff, related to projection
C
      IF (IPROJECT .NE. 0) THEN
      write(6,*) ' Projection used in eom-cc calculation: ', 
     &             String1(iproject+1)
      write(6,*) ' Excitation pattern: ', 
     &             String2(ipattern+1)

      ELSE

      write(6,"(a,I3)") '  Number of different eom-cc calculations:',
     &                     ncalc
      WRITE(6,*)
      ENDIF
C
C DETERMINE ROOT SEARCHING PROCEDURE  EVECFOL / EMINFOL
C
      IF (NONSTD.OR.ICHECK.NE.0) THEN
        EVECFOL = .TRUE.
      ELSE
        IF (IPATTERN .EQ. 0) EVECFOL = .FALSE.
        IF (IPATTERN .EQ. 1) EVECFOL =  (IPROJECT .EQ. 0)
        IF (IPATTERN .GT. 1) EVECFOL = .TRUE.
      ENDIF
C
      EMINFOL = .NOT. EVECFOL
      EVECFOL = .TRUE.
C
      IF (EVECFOL) THEN
        WRITE(6,*) ' Root following procedure: overlap guess vector'
      ELSE
        WRITE(6,*) ' Root following procedure: Minimum energy'
      ENDIF

C Set the root search algorithm. The defalut is standard Davidson 
C algorithm 

      DAVID       = .TRUE.
      BLOCK_DAVID = .FALSE.
      LANCZOS     = .FALSE.

      IF (IFLAGS2(h_IFLAGS2_ee_root_find) .EQ. 1) THEN
         DAVID       = .FALSE.
         BLOCK_DAVID = .TRUE.
      ELSEIF  (IFLAGS2(h_IFLAGS2_ee_root_find) .EQ. 2) THEN
         DAVID       = .FALSE.
         LANCZOS     = .TRUE.
      ENDIF 
C
4000  FORMAT(T3,'@SETMET-F, %EXCITE record missing from ZMAT file.')
5000  FORMAT(T3,'@SETMET-I, Maximum overlap method used for root ',
     &          'searches.')
5001  FORMAT(T3,'@SETMET-I, Energy minimization method used for root ',
     &          'searches.')
6000  FORMAT(T3,'@SETMET-I, VPROP property integrals are used.')
6001  FORMAT(T3,'@SETMET-I, VDINT property integrals are used.')
6002  FORMAT(T3,'@SETMET-F, No excited states specified.')
6003  FORMAT(T3,'@SETMET-I, Gradients and properties are ',
     &          'calculated for only the last state.')
6004  FORMAT(T3,'@SETMET-I, A total of ',I3,' final states will be ',
     &          'studied.')
C
      RETURN
C
      END
