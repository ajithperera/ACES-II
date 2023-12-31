      SUBROUTINE ZEROV1(WMAT, MAXCOR, IRREPXR, IUHF, LISTL1, LISTOFF)
C
C In order to zero-out a block of a Matrix
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)

      DIMENSION WMAT(*)
C
      I0L1AA  = 1
      I0L1BB  = I0L1AA + IRPDPD(IRREPXR,9)
      IEND    = I00L1B + IRPDPD(IRREPXR,10)

      CALL GETLST(WMAT(I0L1AA), 1, 1, 1,  LISTOFF, LIST1)
      LEN1AA = IRPDPD(IRREPXR,9)
      CALL ZERO(WMAT(I0L1AA), LEN1AA)
      CALL  PUTLST(WMAT(I0L1AA), 1, 1, 1,  LISTOFF, LIST1)
C
      IF (IUHF .NE. 0) THEN
         CALL GETLST(WMAT(I0L1BB), 1, 1, 1,  1+LISTOFF, LIST1)
         LEN1AA = IRPDPD(IRREPXR,9)
         CALL ZERO(WMAT(I0L1BB), LEN1AA)
         CALL  PUTLST(WMAT((I0L1BB)), 1, 1, 1,  1+LISTOFF, LIST1)
      ENDIF 
   
C 
      RETURN
      END

