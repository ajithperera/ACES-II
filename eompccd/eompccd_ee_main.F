#include "flags.h"
      SUBROUTINE EOMPCCD_EE_MAIN(ICORE,MAXCOR,IUHF)

      IMPLICIT INTEGER (A-Z)

      DIMENSION ICORE(MAXCOR)

      LOGICAL CIS,EOMCC,CISD,FULDIAG,INCORE,READGUES,DOUBLE,NONSTD
      LOGICAL ESPROP,RPA,VPROP,LAMTHERE, NODAVID, TRIPLET
      LOGICAL SS,SD,DS,DD,CC,MBPT2,TRANABCI,CCD,RCCD,DRCCD
      LOGICAL LCCD,LCCSD,CC2,ADC2
      LOGICAL EOM_exite_exist
      LOGICAL DAVID,BLOCK_DAVID,LANCZOS,PCCD

      DOUBLE PRECISION TDAVID,TMULT,POLTOT,R,DNULL

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

      DATA IONE,DNULL/1,0.0D0/

      I0     = IONE
      TDAVID = DNULL
      TMULT  = DNULL

      CALL EOMPCCD_EE_SETMET(ICORE(I0),IUHF,PCCD)
      CALL ZERO(POLTOT,9)
C
      IF(ESPROP)CALL NUCDIP(ICORE(I0))

      IF ((IFLAGS2(h_IFLAGS2_geom_opt).EQ.0 .AND. IFLAGS(54).EQ.0) 
     &    .OR. IFLAGS(1).GE.1) CALL ORBANAL(ICORE(I0),MAXCOR,IUHF)
C
      CALL GETREC(-1,'JOBARC','RESTART ',IONE,IRES)
      IF(IFLAGS(87).NE.3)IRES=0
      CALL EOMPCCD_EE_MAKLST(ICORE(I0),MAXCOR,IUHF)
      CALL INITTDA(ICORE(I0),IUHF)

C  START CALCULATIONS. POSSIBLY NCALC CYCLES 

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
        IF(CIS)THEN
           CALL MAKESS(ICORE(I0),MAXCOR,IUHF)
           CALL DRVTDA(ICORE(I0),MAXCOR/IINTFP,IUHF)
           CALL RESETSS(ICORE(I0),MAXCOR,IUHF)
        ELSE 
           IF(READGUES)THEN
             CIS=.TRUE.
             CALL MAKESS(ICORE(I0),MAXCOR,IUHF)

C READGUES is set to tru  and NONSTD is TRUE only if %EXCITE*, So
C if GUESS vector is to be read from ZMAT then use the %EXCITE*.
C
             IF (.NOT. NONSTD) THEN
                CALL DRVTDA(ICORE(I0),MAXCOR/IINTFP,IUHF)
             ENDIF 
             CALL RESETSS(ICORE(I0),MAXCOR,IUHF)
             CIS=.FALSE.
           ENDIF 

           CALL MODF(ICORE(I0),MAXCOR,IUHF,1)

C Open a file call EOM_excite to write the dominant contrbutions to excited
C states in %excite format.
C      
           Iunit= 556
           Inquire(File="EOM_excite",Exist=EOM_exite_exist)
           If (.NOT. EOM_exite_exist) Then
               Open(Unit=Iunit, File="EOM_excite",Status="New",
     &              Form="Formatted") 
           Else
               Open(Unit=Iunit, File="EOM_excite",Status="Old",
     &              Form="Formatted") 
                    Close(Iunit,Status="delete")
               Open(Unit=Iunit, File="EOM_excite",Status="New",
     &              Form="Formatted") 
           Endif 
C
C  CHECK IF HBARABCI AND HBARABCD INTEGRALS ARE CORRECT
C
           IF (DD .AND. CC .AND. IFLAGS2(122) .EQ. 2) THEN
           Write(6,*) ' VEE program is fed transformed ABCD integrals'
           Write(6,*) ' This is not supported in the current version'
           call errex
           ENDIF

           IF (DD .AND. CC .AND. IFLAGS2(123) .EQ. 2) THEN
           Write(6,*) ' VEE program is fed transformed ABCI integrals'
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
      Write(6,"(a)") " Hbar before entry to Davidson step"
      Call checkhbar(ICORE(I0), MAXCOR/IINTFP, IUHF)
#endif
           IF (DAVID) THEN
              CALL DOEOMEE_DAVID(ICORE(I0),MAXCOR,IUHF,PCCD)
           ELSE IF (BLOCK_DAVID) THEN 
              CALL DOEOMEE_BLOCKDAVID(ICORE(I0),MAXCOR,IUHF)
           ELSE IF (LANCZOS) THEN 
               CALL DOEOMEE_LANCZOS(ICORE(I0),MAXCOR,IUHF)
           ENDIF

           CALL MODF(ICORE(I0),MAXCOR,IUHF,-1)
    
        ENDIF 

  500 CONTINUE
C
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

      RETURN
      END
