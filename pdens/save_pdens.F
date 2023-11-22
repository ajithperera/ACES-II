      SUBROUTINE SAVE_PDENS(DOOA,DOOB,DVVA,DVVB,DVOA,DVOB,DOVA,
     &                      DOVB,IUHF)
C
C THIS SUBROUTINE SAVES THE 1ST ORDER PERTURB DENSITY MATRIX ON THE
C GAMLAM FILE. THIS IS ONLY REQUIRED FOR SECOND DERIVATIVE CALCULATIONS.
C
C THE LISTS WRITTEN ARE:
C
C  DOOA    1,160
C  DOOB    2,160    UHF AND ROHF ONLY
C  DVVA    3,160
C  DVVB    4,160    UHF AND ROHF ONLY
C  DVOA    5,160  
C  DVOB    6,160    UHF AND ROHF ONLY
C  DOVA    7,160
C  DOVB    8,160    UHF AND ROHF ONLY
C

      INTEGER POP,VRT
      LOGICAL ANTI_SYMM
      DIMENSION DOOA(1),DOOB(1),DVVA(1),DVVB(1),DVOA(1),
     &          DVOB(1),DOVA(1),DOVB(1)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
C
      IENTER = 0
      IOFF = 0
C
      DO 100 ISPIN=1,IUHF+1
C
       CALL UPDMOI(3,NF1(ISPIN),ISPIN,160,IENTER,IOFF)
       CALL UPDMOI(3,NF2(ISPIN),2+ISPIN,160,IENTER,IOFF)
       CALL UPDMOI(3,NT(ISPIN),4+ISPIN,160,IENTER,IOFF)
       CALL UPDMOI(3,NT(ISPIN),6+ISPIN,160,IENTER,IOFF)
C
100   CONTINUE
C
#if _DEBUG_LVLM
      Print*, "-------The A and B Density Mats in SAVE_PDENS ------"
      Print*, "occ-occ A and B"
      call output(dooa, 1, nf1(1), 1, 1, nf1(1), 1, 1)
      Write(6,*)
      call output(doob, 1, nf1(2), 1, 1, nf1(2), 1, 1)
      Print*, "vrt-vrt A and B"
      call output(dvva, 1, nf2(1), 1, 1, nf2(1), 1, 1)
      Write(6,*)
      call output(dvvb, 1, nf2(2), 1, 1, nf2(2), 1, 1)
      Write(6,*)
      Print*, "occ-vrt A and B"
      call output(dova, 1, nt(1), 1, 1, nt(1), 1, 1)
      Write(6,*)
      call output(dovb, 1, nt(2), 1, 1, nt(2), 1, 1)
      Write(6,*)
      Print*, "vrt-occ A and B"
      call output(dvoa, 1, nt(1), 1, 1, nt(1), 1, 1)
      Write(6,*) 
      call output(dvob, 1, nt(2), 1, 1, nt(2), 1, 1)

#endif
      CALL PUTLST(DOOA,1,1,1,1,160)
      CALL PUTLST(DVVA,1,1,1,3,160)
      CALL PUTLST(DVOA,1,1,1,5,160)
      CALL PUTLST(DOVA,1,1,1,7,160)
C
      IF (IUHF.NE.0) THEN
C
       CALL PUTLST(DOOB,1,1,1,2,160)
       CALL PUTLST(DVVB,1,1,1,4,160)
       CALL PUTLST(DVOB,1,1,1,6,160)
       CALL PUTLST(DOVB,1,1,1,8,160)
C
      ENDIF
C
      RETURN
      END
