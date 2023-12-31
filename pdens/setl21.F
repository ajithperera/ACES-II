      SUBROUTINE SETL21(WMAT, MAXCOR, IRREPXR, IUHF, LISTL1, LISTL2, 
     &                  LISTOFF, SKIP_ONEP)
C
C In order to zero-out a block of a Matrix
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL SKIP_ONEP
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)

      DIMENSION WMAT(*)
C
      I0L1AA  = I000
      I0L1BB  = I0L1AA + IRPDPD(IRREPXR,9)
      IEND    = I00L1B + IRPDPD(IRREPXR,10)

      IF (.NOT. SKIP_ONEP) THEN
      CALL GETLST(WMAT(I0L1AA), 1, 1, 1,  LISTOFF, LISTL1)
      LEN1AA = IRPDPD(IRREPXR,9)
#ifdef _DEBUG_LVL0
      Write(6,*)
      Print*, "The lambda1 (AA) vector readin"
      call output(WMAT(I0L1AA), 1, LEN1AA, 1, 1, LEN1AA, 1, 1)
#endif
      CALL ONE(WMAT(I0L1AA), LEN1AA)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Print*, "The lambda1 (AA) vector putin"
      call output(WMAT(I0L1AA), 1, LEN1AA, 1, 1, LEN1AA, 1, 1)
      Write(6,*)
#endif
      CALL  PUTLST(WMAT(I0L1AA), 1, 1, 1,  LISTOFF, LISTL1)
C
      IF (IUHF .NE. 0) THEN
         CALL GETLST(WMAT(I0L1BB), 1, 1, 1,  1+LISTOFF, LISTL1)
#ifdef _DEBUG_LVL0
      Write(6,*)
      Print*, "The lambda1 (BB) vector readin"
      call output(WMAT(I0L1BB), 1, LEN1AA, 1, 1, LEN1AA, 1, 1)
      write(6,*)
#endif
         LEN1BB = IRPDPD(IRREPXR,9)
         CALL ONE(WMAT(I0L1BB), LEN1BB)
         CALL  PUTLST(WMAT((I0L1BB)), 1, 1, 1,  1+LISTOFF, LISTL1)
#ifdef _DEBUG_LVL0
      Print*, "The lambda1 (BB) vector putin"
      call output(WMAT(I0L1BB), 1, LEN1AA, 1, 1, LEN1AA, 1, 1)
      Write(6,*)
#endif
      ENDIF 
C endif for SKIP_ONEP
      ENDIF
C 
      LENL2AA   = IDSYMSZ(IRREPXR,ISYTYP(1,LISTL2),ISYTYP(2,LISTL2))
      LENL2AA10 = IDSYMSZ(IRREPXR,ISYTYP(1,LISTL2+10),
     &                            ISYTYP(2,LISTL2+10))
      IF (IUHF .NE. 0) THEN
          LENL2BB   = IDSYMSZ(IRREPXR,ISYTYP(1,LISTL2+1),
     &                                ISYTYP(2,LISTL2+1))
          LENL2BB10 = IDSYMSZ(IRREPXR,ISYTYP(1,LISTL2+10+1),
     &                                ISYTYP(2,LISTL2+10+1))
          LENL2AB0  = IDSYMSZ(IRREPXR,ISYTYP(1,LISTL2+2),
     &                                ISYTYP(2,LISTL2+2))
          LENL2AB1  = IDSYMSZ(IRREPXR,ISYTYP(1,LISTL2+3),
     &                                ISYTYP(2,LISTL2+3))
          LENL2AB2  = IDSYMSZ(IRREPXR,ISYTYP(1,LISTL2+4),
     &                                ISYTYP(2,LISTL2+4))
          LENL2AB3  = IDSYMSZ(IRREPXR,ISYTYP(1,LISTL2+5),
     &                                ISYTYP(2,LISTL2+5))
          LENL2AB10 = IDSYMSZ(IRREPXR,ISYTYP(1,LISTL2+2+10),
     &                                ISYTYP(2,LISTL2+2+10))
      ENDIF
      
      I0L2AA   = IEND
      I0L2BB   = I0L2AA   + LENL2AA
      I0L2AB0  = I0L2BB   + LENL2BB
      I0L2AB1  = I0L2AB0  + LENL2AB0
      I0L2AB2  = I0L2AB1  + LENL2AB1
      I0L2AB3  = I0L2AB2  + LENL2AB2
      I0L2AA10 = I0L2AB3  + LENL2AB3
      I0L2BB10 = I0L2AA10 + LENL2AB3
      I0L2AB10 = I0L2BB10 + LENL2AB3
      IEND     = I0L2AB10 + LENL2AB3
 

C
      CALL GETALL(WMAT(I0L2AA),   LENL2AA,   IRREPXR, LISTL2)
      CALL GETALL(WMAT(I0L2AA10), LENL2AA10, IRREPXR, LISTL2+10)
C
      CALL MONE(WMAT(I0L2AA)  , LENL2AA)
      CALL MONE(WMAT(I0L2AA10), LENL2AA10)
C
      CALL PUTALL(WMAT(I0L2AA),   LENL2AA,   IRREPXR, LISTL2)
      CALL PUTALL(WMAT(I0L2AA10), LENL2AA10, IRREPXR, LISTL2+10)
C
      IF (IUHF .NE. 0) THEN
C
         CALL GETALL(WMAT(I0L2BB),    LENL2BB ,   IRREPXR, LISTL2+1)
         CALL GETALL(WMAT(I0L2BB10),  LENL2BB10,  IRREPXR, 
     &               LISTL2+1+10)
C
         CALL GETALL(WMAT(I0L2AB0),  LENL2AB0,  IRREPXR, LISTL2+2)
         CALL GETALL(WMAT(I0L2AB10), LENL2AB10, IRREPXR, 
     &               LISTL2+2+10)
C
         CALL GETALL(WMAT(I0L2AB1), LENL2AB1, IRREPXR, LISTL2+3)
         CALL GETALL(WMAT(I0L2AB2), LENL2AB2, IRREPXR, LISTL2+4)
         CALL GETALL(WMAT(I0L2AB3), LENL2AB3, IRREPXR, LISTL2+5)
C
         CALL MONE(WMAT(I0L2BB) ,   LENL2BB )
         CALL MONE(WMAT(I0L2BB10) , LENL2BB10)
C
         CALL ONE(WMAT(I0L2AB0),  LENL2AB0)
         CALL ONE(WMAT(I0L2AB10), LENL2AB10)
C
         CALL ONE(WMAT(I0L2AB1), LENL2AB1)
         CALL ONE(WMAT(I0L2AB2), LENL2AB2)
         CALL ONE(WMAT(I0L2AB3), LENL2AB3)
C
         CALL PUTALL(WMAT(I0L2BB) , LENL2BB , IRREPXR, LISTL2+1)
         CALL PUTALL(WMAT(I0L2BB10) , LENL2BB , IRREPXR, LISTL2+1+10)
C
         CALL PUTALL(WMAT(I0L2AB0), LENL2AB0, IRREPXR, LISTL2+2)
         CALL PUTALL(WMAT(I0L2AB10), LENL2AB10, IRREPXR, LISTL2+2+10)
C
         CALL PUTALL(WMAT(I0L2AB1), LENL2AB1, IRREPXR, LISTL2+3)
         CALL PUTALL(WMAT(I0L2AB2), LENL2AB2, IRREPXR, LISTL2+4)
         CALL PUTALL(WMAT(I0L2AB3), LENL2AB3, IRREPXR, LISTL2+5)
      ENDIF
C
#ifdef _DEBUG_LVL0
      Print*, "The lambda2 AA,BB,AB vectors putin"
      call output(WMAT(I0L2AA),  1, LENL2AA , 1, 1, LENL2AA , 1, 1)
      call output(WMAT(I0L2BB),  1, LENL2BB , 1, 1, LENL2BB , 1, 1)
      call output(WMAT(I0L2AB0), 1, LENL2AB0, 1, 1, LENL2AB0, 1, 1)
      call output(WMAT(I0L2AB1), 1, LENL2AB1, 1, 1, LENL2AB1, 1, 1)
      call output(WMAT(I0L2AB2), 1, LENL2AB2, 1, 1, LENL2AB2, 1, 1)
      call output(WMAT(I0L2AB3), 1, LENL2AB3, 1, 1, LENL2AB3, 1, 1)
      call output(WMAT(I0L2AA10), 1, LENL2AA10, 1, 1, LENL2AA10, 1, 1)
      call output(WMAT(I0L2BB10), 1, LENL2BB10, 1, 1, LENL2BBA0, 1, 1)
      call output(WMAT(I0L2AB10), 1, LENL2AB10, 1, 1, LENL2AB10, 1, 1)
      Write(6,*)
#endif

 
      RETURN
      END

