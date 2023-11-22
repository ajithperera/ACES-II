C
C THIS PROGRAM SOLVES FOR EXCITATION ENERGIES USING THE PARTITIONED
C COUPLED-CLUSTER EQUATION OF MOTION APPROACH.  OTHER OPTIONS INCLUDE
C EOM-SINGLES AND CI-SINGLES METHODS.
C
C PROGRAMMED BY J.F. STANTON, GAINESVILLE, 1992
C

#include "flags.h"

      PROGRAM VEE
      IMPLICIT INTEGER (A-Z)
      LOGICAL CIS,EOMCC,CISD,FULDIAG,INCORE,READGUES,DOUBLE,NONSTD
      LOGICAL ESPROP,RPA,VPROP,LAMTHERE, NODAVID, TRIPLET
      LOGICAL SS,SD,DS,DD,CC,MBPT2,TRANABCI,CCD,RCCD,DRCCD
      LOGICAL LCCD,LCCSD,CC2,ADC2
      LOGICAL EOM_exite_exist
      LOGICAL DAVID,BLOCK_DAVID,LANCZOS
      DOUBLE PRECISION TDAVID, TMULT, POLTOT, R
      COMMON / / ICORE(1)
      COMMON /ISTART/ I0,ICRSIZ
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /INFO /  NOCCO(2),NVRTO(2)
      COMMON /PROPGRAD/ ESPROP,IOPTROOT,IOPTSYM
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /LAMSTATE/ LAMTHERE
      COMMON /EXTRAP/ MAXEXP,NREDUCE,NTOL ,NSIZEC
      COMMON /RESTART_COM/ IRES
      COMMON /POLAR/ POLTOT(3,3)
      COMMON /RMAT/ R(10000)
      COMMON /GUESS/ DOUBLE,NONSTD
      COMMON /METH/CIS,RPA,EOMCC,CISD,FULDIAG,INCORE,READGUES
      COMMON /VDINTPRT/NTPERT,NPERT(8),KPERT(8),IDIPSM(3),
     &                 IYZPERT,IXZPERT,IXYPERT,ITRANSX,
     &                 ITRANSY,ITRANSZ,NUCIND
      COMMON/EIGPROB/ISIDE
      COMMON /INTPROG/ VPROP
      COMMON/PROJECT/IPROJECT, IPATTERN, NCALC, ICALC, IWINDOW(8)
      COMMON/LISTPROJ/LISTH0, ICOLPR1, ICOLPR2
      COMMON/TIMSUB/TDAVID, TMULT
      COMMON/PARTEOM/NODAVID
      COMMON /SPINSTATE/TRIPLET
      COMMON/DRVHBAR/SS, SD, DS, DD
      COMMON /REFTYPE/ MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
      COMMON/ROOTSEARCH/DAVID,BLOCK_DAVID,LANCZOS
C 
#include "parcc.com"
C
      IONE=1
      TDAVID = 0.0D0
      TMULT = 0.0D0

      CALL CRAPSI(ICORE,IUHF,0)
      CALL ACES_IO_REMOVE(54,'DERGAM')
      CALL SETMET(ICORE(I0),IUHF)

      Coulomb = .False.
      if (CC .OR. CCD) Then
         call parread(iuhf)
         if (ispar) then
           If ((iflags(87) .EQ. 11)) Then
              Write(6,"(a,a)") " Excitation energies for triple",
     &        " excitations pCC methods have not been verified"
              call aces_exit(1)
           Endif
           write(6,*) ' Parameterized EOM-CC calculation is performed'
           write(6,2010) paralpha
           write(6,2011) parbeta
           write(6,2012) pargamma
           write(6,2013) pardelta
           write(6,2014) parepsilon
 2010      format(' PCCSD   alpha parameter : ', F14.6)
 2011      format(' PCCSD    beta parameter : ', F14.6)
 2012      format(' PCCSD   gamma parameter : ', F14.6)
 2013      format(' PCCSD   delta parameter : ', F14.6)
 2014      format(' PCCSD epsilon parameter : ', F14.6)
           if (coulomb) Write(6,"(a,a)") " The Coulomb integrals are ",
     $                    "used in W(mbej) intermediate."
           write(6,*)
           Fae_scale    = (Paralpha - 1.0D0)
           Fmi_scale    = (Parbeta  - 1.0D0)
           Wmnij_scale  = Pargamma
           Wmbej_scale  = Pardelta
           Gae_scale    = Paralpha 
           Gmi_scale    = Parbeta 
         else
           write(6,*) ' Regular EOM-CC calculation is performed'
           write(6,*)
           Fae_scale    = 0.0D0
           Fmi_scale    = 0.0D0
           Wmnij_scale  = 1.0D0
           Wmbej_scale  = 1.0D0
           Gae_scale    = 1.0D0
           Gmi_scale    = 1.0D0
         endif
      endif

      IF (IFLAGS(35).NE.0) CALL INCOR(I0,ICRSIZ,IUHF)
      MAXCOR=ICRSIZ
      CALL ZERO(POLTOT,9)
C
C ADC2 check; CALC=MBPT(2),EXCITE=EOMEE,EOMREF=ADC2
C
      ADC2 = .FALSE.
      ADC2 =  (IFLAGS(2) .EQ. 1 .AND. IFLAGS(87) .EQ. 3 .AND.
     &         IFLAGS2(117) .EQ. 10)

      IF(ESPROP)CALL NUCDIP(ICORE(I0))
      IF ((IFLAGS2(h_IFLAGS2_geom_opt).EQ.0 .AND.
     &     IFLAGS(54).EQ.0) .OR. IFLAGS(1).GE.1)
     &   CALL ORBANAL(ICORE(I0),MAXCOR,IUHF)
C
      CALL GETREC(-1,'JOBARC','RESTART ',IONE,IRES)
      IF(IFLAGS(87).NE.3)IRES=0
CMN      IF(IRES.NE.1)THEN
        CALL MAKLST(ICORE(I0),MAXCOR,IUHF)
CMN      ENDIF
      CALL INITTDA(ICORE(I0),IUHF)
CMN
C  START CALCULATIONS. POSSIBLY NCALC CYCLES 
C
      IFIRST = 1
      ILAST = NCALC
      DO 500 JCALC = IFIRST, ILAST
        WRITE(6,*)
        WRITE(6,"(A,I3,A,I3,A)") "  Currently doing the site ", 
     +                        JCALC, " out of ", NCALC, " sites!"
        WRITE(6,"(a)") "        ------------------------"
        WRITE(6,"(a)") 
        ICALC = JCALC
        IF (IFLAGS(35).NE.0) CALL INCOR(I0,ICRSIZ,IUHF)
        MAXCOR=ICRSIZ
CMN END
      IF(CIS)THEN
       If (Ispar) CALL RESTORE_CC_WMBEJ(ICORE(I0),MAXCOR,IUHF)
       CALL MAKESS(ICORE(I0),MAXCOR,IUHF)
       CALL DRVTDA(ICORE(I0),MAXCOR/IINTFP,IUHF)
       CALL RESETSS(ICORE(I0),MAXCOR,IUHF)
      ELSEIF (RPA) THEN
       NMULT = 1
       IF ((RCCD .OR. DRCCD) .AND. IUHF .EQ. 0) NMULT = 2
       DO IMULT = 1, NMULT
          CALL EOM_RCC_DRIVER(ICORE(I0),MAXCOR,IUHF,IMULT)
          CALL EOM_RCC_DIAGS(ICORE(I0),MAXCOR/IINTFP,IUHF,IMULT)
       ENDDO 
      ELSE 
       IF(READGUES)THEN
        CIS=.TRUE.
        If (Ispar) CALL RESTORE_CC_WMBEJ(ICORE(I0),MAXCOR,IUHF)
        CALL MAKESS(ICORE(I0),MAXCOR,IUHF)
C
C READGUES is set to tru  and NONSTD is TRUE only if %EXCITE*, So
C if GUESS vector is to be read from ZMAT then use the %EXCITE*.
C
        IF (.NOT. NONSTD) THEN
         CALL DRVTDA(ICORE(I0),MAXCOR/IINTFP,IUHF)
        ENDIF 

        CALL RESETSS(ICORE(I0),MAXCOR,IUHF)
        CIS=.FALSE.
       ENDIF
CSSS#ifdef _DCC_FLAG
       If (Ispar) Then
          CALL GETREC(0,"JOBARC","LAMBDA  ",LENGTH1,JUNK)
          CALL GETREC(0,"JOBARC","HBAR    ",LENGTH2,JUNK)

C If lambda is done then the FAE ad FMI need to be modified. If hbar,lhbar
C route is taken then this step can be skipped. The lambda record indicates
C that the code has taken the lambda route. 

          IF (LENGTH1 .GT. 0 .AND. LENGTH2 .LT. 0) 
     &    CALL PDCC_FIXF(ICORE(I0),MAXCOR,IUHF)
          CALL MODF(ICORE(I0),MAXCOR,IUHF,1)
          CALL PDCC_MODF(ICORE(I0),MAXCOR,IUHF,1)
       Else
CSSS#else  
          CALL MODF(ICORE(I0),MAXCOR,IUHF,1)
       Endif 
CSSS#endif 
C Open a file call EOM_excite to write the dominant contrbutions to excited
C states in %excite format.
C      
      Iunit= 556
      Inquire(File="EOM_excite",Exist=EOM_exite_exist)
      If (.NOT. EOM_exite_exist) Then
          Open(Unit=Iunit, File="EOM_excite",Status="New",
     &         Form="Formatted") 
      Else
          Open(Unit=Iunit, File="EOM_excite",Status="Old",
     &         Form="Formatted") 
               Close(Iunit,Status="delete")
          Open(Unit=Iunit, File="EOM_excite",Status="New",
     &         Form="Formatted") 
      Endif 
CMN
C  THE CODE HAS CHANGED!! IN PREVIOUS VERSIONS MAKESS
C  PUT A MINUS SIGN FOR UHF IN THE RELEVANT LISTS 54 AND 55
C  THIS IS NO LONGER REQUIRED
CMN END
C
       if (iuhf .eq. 0) CALL MAKESS(ICORE(I0),MAXCOR,IUHF)
       IF (CC2 .OR. ADC2) CALL MODF2(ICORE(I0),MAXCOR,IUHF,1)
C
C  CHECK IF HBARABCI AND HBARABCD INTEGRALS ARE CORRECT
C
       IF (DD .AND. CC .AND. IFLAGS2(122) .EQ. 2) THEN
         Write(6,*) ' VEE program is fed transformed ABCD integrals'
         Write(6,*) ' This is not supported in the current version'
         call errex
       ENDIF
       IF (DD .AND. CC .AND. IFLAGS2(123) .EQ. 2) THEN
       Write(6,*) ' VEE program is fed fully transformed ABCI integrals'
         Write(6,*) ' This is not supported in the current version'
         call errex
       ENDIF
C
C For MBPT(2), Linear CC, CCSD this does not do anything (as the
C logic indicates;IFLAGS2(123)=1=HBARABCI=OFF is the default setting)
C
       TRANABCI = (CC .AND. IFLAGS2(123) .EQ. 1 .AND. .NOT. DD)
       CALL MODHBAR(ICORE(I0), MAXCOR,IUHF, .FALSE., TRANABCI)

#ifdef _DEBUG_LVL0
      Write(6,*) "Enter checkhbar",TRANABCI
      Call checkhbar(ICORE(I0), MAXCOR/IINTFP, IUHF)
#endif
C
       IF (NODAVID) THEN
          IF (IFLAGS(87) .EQ. 8)
     &      CALL DOBWPT2(IUHF, ICORE(I0), MAXCOR/IINTFP)
          IF (IFLAGS(87) .EQ. 7)
     &      CALL DOPARTEOM(IUHF, ICORE(I0), MAXCOR/IINTFP)
       ELSE
            IF (DAVID) THEN
               CALL DOEOMEE_DAVID(ICORE(I0),MAXCOR,IUHF)
            ELSE IF (BLOCK_DAVID) THEN 
               CALL DOEOMEE_BLOCKDAVID(ICORE(I0),MAXCOR,IUHF)
            ELSE IF (LANCZOS) THEN 
               CALL DOEOMEE_LANCZOS(ICORE(I0),MAXCOR,IUHF)
            ENDIF
       ENDIF

CSSS       CALL Diagall(ICORE(I0),MAXCOR/IINTFP,IUHF)
C
C   clean up lists
C
       if (iuhf. eq. 0) CALL RESETSS(ICORE(I0),MAXCOR,IUHF)
CSSS#ifdef _DCC_FLAG
       If (Ispar) Then
          CALL MODF(ICORE(I0),MAXCOR,IUHF,-1)
          CALL PDCC_MODF(ICORE(I0),MAXCOR,IUHF,-1)
       Else
CSSS#else
          CALL MODF(ICORE(I0),MAXCOR,IUHF,-1)
       Endif 
cSSS#endif
      ENDIF
C
  500 CONTINUE
C
CMN      CALL PUTREC(20,'JOBARC','RESTART ',IONE,IONE)
      IF (IFLAGS(35).NE.0) THEN
         CALL ACES_AUXCACHE_FLUSH
         CALL ACES_AUXCACHE_RESET
      END IF
      CALL ACES_IO_REMOVE(54,'DERGAM')
      IF (.NOT.CIS .AND. IFLAGS(1).GE.2) THEN
        WRITE(6,*)
        IF (.NOT. NODAVID) THEN
          WRITE(6,1000) TDAVID
 1000     FORMAT(' TOTAL TIME LAPSED IN DAVIDSON: ',F12.4)
        ENDIF
        WRITE(6,1001) TMULT
 1001   FORMAT(' TOTAL TIME LAPSED IN HBARXC:   ',F12.4)
        WRITE(6,*)
      ENDIF
      call aces_fin
      STOP
      END
