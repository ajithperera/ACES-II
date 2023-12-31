      SUBROUTINE NEWT1(ICORE,MAXCOR,IUHF,INCREM)
C
C THIS ROUTINE PICKS UP THE FINAL T1 OR T1 INCREMENTS, DENOMINATOR
C  WEIGHTS THEM AND THEN OVERWRITES THE T1 INCREMENT LIST WITH THE NEW VALUES.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR)
      DOUBLE PRECISION ONE
      LOGICAL NONHF
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/ POP1(8),POP2(8),VRT1(8),VRT2(8),NT(2),NF1AA,
     &             NF1BB,NF2AA,NF2B
      COMMON /NHFREF/ NONHF
      DATA ONE /1.0/
      I000=1
      DO 5 ISPIN=1,1+IUHF
       LSTDEN=63+ISPIN
       LSTINC=2+ISPIN
       NSIZE=NT(ISPIN)
       I010=I000+NSIZE*IINTFP
       I020=I010+NSIZE*IINTFP
       I030=I020+NSIZE*IINTFP
       IF(I030.GT.MAXCOR)CALL INSMEM('NEWT1',I020,MAXCOR)
       CALL GETLST(ICORE,1,1,1,LSTINC,90)
       CALL GETLST(ICORE(I010),1,1,1,9,LSTDEN)
       CALL VECDIV(ICORE(I000),ICORE(I010),ICORE(I000),NSIZE)
C
C IF THIS IS SOME SORT OF ITERATIVE SCHEME, WE WANT TO STICK THE
C  FIRST ORDER T1 INTO THIS STUFF
C
        IF(INCREM.NE.0 .AND. NONHF)THEN
         CALL GETLST(ICORE(I020),1,1,1,ISPIN+2,93)
         CALL VECDIV(ICORE(I020),ICORE(I010),ICORE(I020),NSIZE)
         CALL SAXPY (NSIZE,ONE,ICORE(I020),1,ICORE(I000),1)
        ENDIF
        CALL PUTLST(ICORE(I000),1,1,1,LSTINC,90)

#ifdef _DEBUG_LVLM
        write(*,"(a)") "NEWT1 statistics"
        call checksum ("@-NEWT1", ICORE(I000), NSIZE)
#endif
5     CONTINUE
      RETURN
      END
