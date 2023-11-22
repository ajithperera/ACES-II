      SUBROUTINE PCCD_GENINT(ICORE,MAXCOR,IUHF)
C
#include "flags.h" /* for h_IFLAGS_* and h_IFLAGS2_* */
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION C1,C2,C3,C4,FACT,C5
      LOGICAL PCCD,CCD,LCCD
C
      DIMENSION ICORE(MAXCOR)

      COMMON/CALC/PCCD,CCD,LCCD
      COMMON/NHFREF/NONHF
      COMMON/FLAGS/IFLAGS(100)
      COMMON/FLAGS2/IFLAGS2(500)
      COMMON/SYM/POP(8,2),VRT(8,2),NT1AA,NT1BB,NF1AA,NF1BB,
     &           NF2AA,NF2BB
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      C1   = 1.0D0
      C2   = 0.0D0
      C3   = 0.0D0
      IF (LCCD) THEN
         C4   =  0.0D0
      ELSE
         C4   = -0.50D0
      ENDIF 
      FACT = 1.0D0

      CALL PCCD_QUAD1(ICORE,MAXCOR,IUHF,FACT,LCCD)

#ifdef _DEBUG_LVLM
      call pccd_checkintms(icore,Maxcor,iuhf,0)
#endif 
      CALL ZERSYM(ICORE,56) 
      CALL ZERSYM(ICORE,58) 
      CALL PCCD_DWMBEJ(ICORE,MAXCOR,'ABAB',IUHF,C1,C2,C3,C4)
      CALL PCCD_DWMBEJ(ICORE,MAXCOR,'ABBA',IUHF,C1,C2,C3,C4)

#ifdef _DEBUG_LVLM
      call pccd_checkintms(icore,Maxcor,iuhf,0)
#endif 
      IF (.NOT. LCCD) THEN
         CALL PCCD_QUAD2(ICORE,MAXCOR,IUHF,FACT)
         CALL PCCD_QUAD3(ICORE,MAXCOR,IUHF,FACT)
      ELSE
        CALL ZERLST(ICORE,NF1AA,1,1,1,91)
        CALL ZERLST(ICORE,NF2AA,1,1,1,92)
        CALL ZERLST(ICORE,NT1AA,1,1,1,93)
        IF (IUHF .NE. 0) THEN
            CALL ZERLST(ICORE,NF1BB,1,1,2,91)
            CALL ZERLST(ICORE,NF2BB,1,1,2,92)
            CALL ZERLST(ICORE,NT1BB,1,1,2,93)
        ENDIF
      ENDIF 

#ifdef _DEBUG_LVLM
      call pccd_checkintms(icore,Maxcor,iuhf,0)
#endif 
C
C THe last argument says that there are no singles in pCCD or OCCD 
C variants. 

      IF(NONHF) THEN
        CALL PCCD_NHFFINT(ICORE,MAXCOR,IUHF,.FALSE.)
      ENDIF

#ifdef _DEBUG_LVL0
      call pccd_checkintms(icore,Maxcor,iuhf,0)
#endif 
      RETURN
      END
