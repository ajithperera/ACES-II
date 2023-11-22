C
C  THIS ROUTINE CREATES THE LISTS FOR THE DERIVATIVES OF THE
C  MO INTEGRALS
C
      SUBROUTINE INIDI(IRREPX,IUHF)
      IMPLICIT INTEGER (A-Z)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      COMMON/PERTINF/INIT,MPERT(8)   
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
C
C FOR ALL METHODS WE NEED <IJ||AB>^(chi)
C
C NOW CREATE <IJ||AB>^(chi) LISTS
C
C  THE STRUCTURE IS THE FOLLOWING
C
c      IENTER=1
c      IOFF=-1
      call aces_io_remove(53,'DERINT')
      IENTER=0
      IOFF=0
C
      IBOT=1
      IF(IUHF.EQ.0) IBOT=3
      DO 40 ISPIN=IBOT,3
       RPOINT=2+ISPIN
       LPOINT=ISPIN
       IF(ISPIN.EQ.3) THEN
        RPOINT=14
        LPOINT=13
       ENDIF
       LIST=313+ISPIN
       CALL INIPCK(IRREPX,LPOINT,RPOINT,LIST,IENTER,IOFF,1)
       do iGrp = 1, nirrep
          call aces_list_memset(iGrp,LIST,0)
       end do
       IENTER=0
       IOFF=0
40    CONTINUE
C
C FOR SECOND ORDER WE NEED IN ADDITION <IJ||KA>^chi
C AND <AB||CI>^chi
C
      IBOT=1
      IF(IUHF.EQ.0) IBOT=4
      DO 140 ISPIN=IBOT,4
       RPOINT=ISYTYP(2,6+ISPIN)
       LPOINT=ISYTYP(1,6+ISPIN)
       LIST=306+ISPIN
       CALL INIPCK(IRREPX,LPOINT,RPOINT,LIST,IENTER,IOFF,1)
       do iGrp = 1, nirrep
          call aces_list_memset(iGrp,LIST,0)
       end do
       RPOINT=ISYTYP(2,26+ISPIN)
       LPOINT=ISYTYP(1,26+ISPIN)
       LIST=326+ISPIN
       CALL INIPCK(IRREPX,LPOINT,RPOINT,LIST,IENTER,IOFF,1)
       do iGrp = 1, nirrep
          call aces_list_memset(iGrp,LIST,0)
       end do
140   CONTINUE
C
C CREATE THE <IA||JB>^chi LISTS
C
      IBOT=1
      IF(IUHF.EQ.0) IBOT=4
      DO 240 ISPIN=IBOT,4
       LIST=322+ISPIN-1+IUHF
       RPOINT=ISYTYP(2,LIST-300)
       LPOINT=ISYTYP(1,LIST-300)
       CALL INIPCK(IRREPX,LPOINT,RPOINT,LIST,IENTER,IOFF,1)
       do iGrp = 1, nirrep
          call aces_list_memset(iGrp,LIST,0)
       end do
240   CONTINUE
C
C CREATE THE <Aj|Ib> AND <aJ|iB> LISTS
C
       ITOP=2
       IF(IUHF.EQ.0) ITOP=1
       DO 241 ISPIN=1,ITOP
        LIST=320+ISPIN
        RPOINT=ISYTYP(2,LIST-300)
        LPOINT=ISYTYP(1,LIST-300)
        CALL INIPCK(IRREPX,LPOINT,RPOINT,LIST,IENTER,IOFF,1)
        do iGrp = 1, nirrep
           call aces_list_memset(iGrp,LIST,0)
        end do
241    CONTINUE
C
C <IJ||KL>^(chi)
C
      IBOT=1
      IF(IUHF.EQ.0) IBOT=3
      DO 340 ISPIN=IBOT,3
       LIST=310+ISPIN
       RPOINT=ISYTYP(2,LIST-300)
       LPOINT=ISYTYP(1,LIST-300)
       CALL INIPCK(IRREPX,LPOINT,RPOINT,LIST,IENTER,IOFF,1)
       do iGrp = 1, nirrep
          call aces_list_memset(iGrp,LIST,0)
       end do
340   CONTINUE
C
C FOR HIGHER ORDERS OF PERTURBATION THEORY AS WELL AS 
C COUPLED-CLUSTER METHODS, WE NEED ADDITIONAL FILES
C
      IF(.NOT.MBPT2) THEN
C
C <AB||CD>^(chi)
C
      ENDIF
C
C ALL DONE, RETURN
C
      RETURN
      END