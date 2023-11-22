      SUBROUTINE PCCD_DUMPIT(IOO,IVV,IVO,XOV,IUHF,FLAG)
C
C DUMP OUT I AND X INTERMEDIATES
C
      IMPLICIT INTEGER (A-Z)
      CHARACTER*5 SPCASE
      CHARACTER*2 FLAG
      DIMENSION IOO(*),IVV(*),IVO(*),XOV(*),SPCASE(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON/SYM2/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
C
      DATA SPCASE/'alpha','beta '/

#ifdef _DEBUG_LVL0
C
      WRITE(6,1000)
1000  FORMAT(T3,"@DUMPIT-I, D matrices, I or X intermediates") 
C
      IF (FLAG .EQ. "OO" .OR. FLAG .EQ. "AL") THEN
      DO 10 ISPIN=1,1+IUHF
       WRITE(6,1001)SPCASE(ISPIN)
1001   FORMAT(T3,'Occupied-occupied D/I/X intermediate for spin ',A,':')
       CALL PRVECR(IOO(1+IINTFP*(ISPIN-1)*NFMI(1)),NFMI(ISPIN))
       call checksum("G(oo)   :", IOO(1+IINTFP*(ISPIN-1)*NFMI(1)),
     +                            NFMI(ISPIN))
10    CONTINUE
      ENDIF 
C
      IF (FLAG .EQ. "VV" .OR. FLAG .EQ. "AL") THEN
      DO 20 ISPIN=1,1+IUHF
       WRITE(6,2001)SPCASE(ISPIN)
2001   FORMAT(T3,'Virtual-virtual  D/I/X intermediate for spin ',A,':')
       CALL PRVECR(IVV(1+IINTFP*(ISPIN-1)*NFEA(1)),NFEA(ISPIN))
       call checksum("G(vv)   :", IVV(1+IINTFP*(ISPIN-1)*NFEA(1)),
     +                            NFEA(ISPIN))
20    CONTINUE
      ENDIF 
C
      IF (FLAG .EQ. "OV" .OR. FLAG .EQ. "AL") THEN
      DO 30 ISPIN=1,1+IUHF
       WRITE(6,3001)SPCASE(ISPIN)
3001   FORMAT(T3,'Occupied-Virtual D/I/X intermediate for spin ',A,':')
       CALL PRVECR(XOV(1+IINTFP*(ISPIN-1)*NT(1)),NT(ISPIN))
       call checksum("G(vo)   :", XOV(1+IINTFP*(ISPIN-1)*NT(1)),
     +                            NT(ISPIN))
30    CONTINUE
      ENDIF 
C
      IF (FLAG .EQ. "VO" .OR. FLAG .EQ. "AL") THEN
      DO 40 ISPIN=1,1+IUHF
       WRITE(6,4001)SPCASE(ISPIN)
4001   FORMAT(T3,'Virtual-Occupied D/I/X intermediate for spin ',A,':')
       CALL PRVECR(IVO(1+IINTFP*(ISPIN-1)*NT(1)),NT(ISPIN))
       call checksum("G(ov)   :", IVO(1+IINTFP*(ISPIN-1)*NT(1)),
     +                            NT(ISPIN))
40    CONTINUE
      ENDIF 
#endif
C
      RETURN
      END
