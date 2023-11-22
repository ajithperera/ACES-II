
      SUBROUTINE RCL_NEWT2(ICORE,MAXCOR,IUHF)
C
C THIS ROUTINE PICKS UP THE FINAL T2 OR T2 INCREMENTS, DENOMINATOR
C  WEIGHTS THEM AND THEN OVERWRITES THE T2 INCREMENT LIST WITH THE NEW VALUES.
C  GOES ONE IRREP AT A TIME OVER ALL SPIN CASES.  FOR RHF CASES,
C  IT ALSO CALLS A ROUTINE WHICH FORMS THE AA AMPLITUDES FROM THE 
C  AB VALUES.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ELAST,FACT
      DIMENSION ICORE(MAXCOR)
      LOGICAL LINCC,CICALC
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/ POP1(8),POP2(8),VRT1(8),VRT2(8),NTAA,NTBB,NF1AA,
     &             NF1BB,NF2AA,NF2BB
      COMMON /CORENG/ ELAST

      IBOT=1
      IF(IUHF.EQ.0)IBOT=3
      I000=1
      DO 5 ISPIN=IBOT,3
       LSTDEN=47+ISPIN
       LSTINC=60+ISPIN
       DO 10 IRREP=1,NIRREP
        NDSSYM_T2=IRPDPD(IRREP,ISYTYP(2,LSTINC))
        DISSYM_T2=IRPDPD(IRREP,ISYTYP(1,LSTINC))
        NDSSYM_D2=IRPDPD(IRREP,ISYTYP(2,LSTDEN))
        DISSYM_D2=IRPDPD(IRREP,ISYTYP(1,LSTDEN))

        NSIZE_T2=NDSSYM_T2*DISSYM_T2
        NSIZE_D2=NDSSYM_D2*DISSYM_D2

        I010=I000+NSIZE_T2*IINTFP
        I020=I010+NSIZE_D2*IINTFP
        I030=I020+NSIZE_T2*IINTFP
        IEND=I030
        IF(IEND.GT.MAXCOR)CALL INSMEM('NEWT2',IEND,MAXCOR)
        LEFT=MAXCOR - IEND 

        CALL GETLST(ICORE,1,NDSSYM_T2,1,IRREP,LSTINC)
        CALL GETLST(ICORE(I010),1,NDSSYM_D2,1,IRREP,LSTDEN)

        IF (ISPIN .LT. 3) THEN
           CALL RCC_EXPAND(ICORE(I010),ICORE(I020),NDSSYM_T2,
     &                     DISSYM_T2,NDSSYM_D2,DISSYM_D2,
     &                     IRREP,ISPIN,ICORE(IEND),LEFT)

#ifdef _DEBUG_LVLM
           if (ispin .eq. 1) then
           Write(6,*) "L2(AB,IJ) in"
           call output(icore(i000),1,DISSYM_T2,1,NDSSYM_T2,DISSYM_T2,
     &                 NDSSYM_T2,1)
           endif 
#endif 
           CALL VECPRD(ICORE(I000),ICORE(I020),ICORE(I000),NSIZE_T2)

#ifdef _DEBUG_LVLM
CSSS           Write(6,*) "AAAA/BBBB denominator"
CSSS           call output(icore(i010),1,DISSYM_T2,1,NDSSYM_T2,DISSYM_T2,
CSSS     &                 NDSSYM_T2,1)
           if (ispin .eq. 1) then
           Write(6,*) "L2(AB,IJ) out"
           call output(icore(i000),1,DISSYM_T2,1,NDSSYM_T2,DISSYM_T2,
     &                 NDSSYM_T2,1)
           endif 
#endif 
        ELSE 

#ifdef _DEBUG_LVLM
           Write(6,*) "L2(Ab,Ij) in"
           call output(icore(i000),1,DISSYM_D2,1,NDSSYM_D2,DISSYM_D2,
     &                 NDSSYM_D2,1)
#endif 

           CALL VECPRD(ICORE(I000),ICORE(I010),ICORE(I000),NSIZE_T2)

#ifdef _DEBUG_LVLM
CSSS           Write(6,*) "ABAB denominator"
CSSS           call output(icore(i010),1,DISSYM_D2,1,NDSSYM_D2,DISSYM_D2,
CSSS     &                 NDSSYM_D2,1)
           Write(6,*) "L2(Ab,Ij) out"
           call output(icore(i000),1,DISSYM_T2,1,NDSSYM_T2,DISSYM_T2,
     &                 NDSSYM_T2,1)
#endif 
        ENDIF 

        CALL PUTLST(ICORE(I000),1,NDSSYM_T2,1,IRREP,LSTINC)
#ifdef _DEBUG_LVL0
        If (Ispin .Eq.1) call checksum("new-T2(AAAA): ",ICORE(I000),
     &                                  NSIZE_T2)
        If (Ispin .Eq.2) call checksum("new-T2(BBBB): ",ICORE(I000),
     &                                  NSIZE_T2)
        If (Ispin .Eq.3) call checksum("new-T2(ABAB): ",ICORE(I000),
     &                                  NSIZE_T2)
#endif 
10     CONTINUE
5     CONTINUE

      IF(IUHF.EQ.0)THEN

       IF (DRCCD) THEN

          DO IRREP=1,NIRREP
             NUMAB=IRPDPD(IRREP,ISYTYP(2,63))
             DSZAB=IRPDPD(IRREP,ISYTYP(1,63))
             NUMAA=IRPDPD(IRREP,ISYTYP(2,61))
             DSZAA=IRPDPD(IRREP,ISYTYP(1,61))
             NSIZ1=DSZAB*NUMAA
             NSIZ2=DSZAB*NUMAB

             I000=1
             I010=I000+IINTFP*NSIZ1
             I020=I010+IINTFP*NSIZ2
             CALL GETLST(ICORE(I010),1,NUMAB,1,IRREP,63)

CSSS        CALL ASSYM(IRREP,POP1,DSZAB,DSZAB,ICORE(I000),ICORE(I010))
CSSS        CALL SQSYM(IRREP,VRT1,DSZAA,DSZAB,NUMAA,ICORE(I010),ICORE(I000))

            CALL PUTLST(ICORE(I010),1,NUMAA,1,IRREP,61)

#ifdef _DEBUG_LVL0
       call checksum("new-T2(ABAB):",ICORE(I010),NSIZ1)
#endif 
          ENDDO

       ELSEIF (RCCD) THEN

          LSTDEN=50 
          LSTINC=61
          DO IRREP=1,NIRREP
             NDSSYM=IRPDPD(IRREP,ISYTYP(2,LSTINC))
             DISSYM=IRPDPD(IRREP,ISYTYP(1,LSTINC))
             NSIZE = NDSSYM * DISSYM

             I000=1
             I010=I000+NSIZE*IINTFP
             I020=I010+NSIZE*IINTFP
             IEND=I020
             IF(IEND.GT.MAXCOR)CALL INSMEM('NEWT2',IEND,MAXCOR)

             CALL GETLST(ICORE(I000),1,NDSSYM,1,IRREP,LSTINC)
             CALL GETLST(ICORE(I010),1,NDSSYM,1,IRREP,LSTDEN)

             CALL VECPRD(ICORE(I000),ICORE(I010),ICORE(I000),NSIZE)
             CALL PUTLST(ICORE(I000),1,NDSSYM,1,IRREP,LSTINC) 
#ifdef _DEBUG_LVL0
      call checksum("new-T2(ABAB): ",ICORE(I000),NSIZE)
#endif 
          ENDDO 
       ENDIF 
      ENDIF

      RETURN
      END
