      SUBROUTINE S2HHS2_2L(ICORE,MAXCOR,IUHF, ISPIN)
C
C   THIS SUBROUTINE CALCULATES THE CONTRIBUTION OF S2 TO
C   S2 INVOLVING THE HH BLOCK OF THE EFFECTIVE ONE-PARTICLE 
C   INTEGRALS.
C
C CONTRACTION :
C
C     Z(iJ,ap) = -SUM [Hbar(m,j) (2S(mi,ap) - S(im,ap)) [ RHF CASE ]
C                 m 
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      LOGICAL SKIP
C
      DIMENSION ICORE(MAXCOR), NUMSZS(8)
      DIMENSION IOFF_HBAR(8),IOFF_S2IN(8),IOFF_S2EX(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SINFO/ NS(8), SIRREP
      COMMON/SLISTS/LS1IN, LS1OUT, LS2IN(2,2), LS2OUT(2,2)
C
      DATA ONE,IONE  /1.0,1/
      DATA ZILCH/0.0/
      DATA ONEM/-1.0/

      Write(6,"(a)") "-----------------S2HHS2_2R-------------------"
C 
C S2(M=ISPIN,AJ=MIXSPIN)
C
      DO 5 XIRREP = 1, NIRREP

C Here consider ISPIN=MIXSPIN (see below).

         IRREPA = DIRPRD(XIRREP,SIRREP)
         NUMSZS(XIRREP) = VRT(IRREPA,ISPIN) * NS(SIRREP)
 5    CONTINUE      
C
      DO MIXSPIN = 1, 1+IUHF
         LISTS2IN = LS2IN(ISPIN,MIXSPIN + 1 - IUHF)
         CALL GETLEN_HHP(LENS,POP(1,MIXSPIN),POP(1,ISPIN),
     $                   VRT(1,MIXSPIN),NS)

       I000=1
       I010=I000+LENS*IINTFP
C
C AAA and BBB blocks 
C
       IF ((IUHF.NE.0) .AND. (MIXSPIN.EQ.ISPIN)) THEN
          CALL GETEXP2_HHP(ICORE(I000),LENS,NUMSZS,
     +                     IRPDPD(1,20+ISPIN),
     +                     LISTS2IN,POP(1,MIXSPIN), 
     +                     IRPDPD(1,ISYTYP(1,LISTS2IN)))
       ELSE
          CALL GETALLS2_HHP(ICORE(I000),LENS,VRT(1,MIXSPIN),NS(1),1,
     $                      LISTS2IN)
       ENDIF

#ifdef _DEBUG_LVL0
      call checksum("->s2    :",icore(i000),lens)
#endif 
C
C TRANSPOSE THE S BRA INDICES S2IN(MI,AP) -> S2IN(IM,AP). 
C
      IOFF = I000
      DO XIRREP  = 1, NIRREP
         IRREPA  = DIRPRD(XIRREP,SIRREP)
         
         NUMDIS = VRT(IRREPA,MIXSPIN) * NS(SIRREP) 
         IF ((IUHF.NE.0) .AND. (MIXSPIN.EQ.ISPIN)) THEN
            DISSIZ = IRPDPD(XIRREP,20+ISPIN)
         ELSE
            DISSIZ = IRPDPD(XIRREP,ISYTYP(1,LISTS2IN))
         ENDIF 
         MAXSIZ = MAX(NUMDIS,DISSIZ)

         I020   = I010  + MAXSIZ * IINTFP
         I030   = I020  + MAXSIZ * IINTFP
         IEND   = I030  + MAXSIZ 

         IF (IEND .GT. MAXCOR) CALL INSMEM("S2HHPS2",IEND,MAXCOR)
   
         CALL SYMTR3(XIRREP,POP(1,MIXSPIN),POP(1,ISPIN),DISSIZ,
     +               NUMDIS,ICORE(IOFF),ICORE(I010),ICORE(I020),
     +               ICORE(I030))
         IOFF = IOFF + NUMDIS * DISSIZ * IINTFP 
       ENDDO 

C TRANSPOSE, S2IN(IM,AP) -> S2IN(AP,IM). 

      IEND = I010 + LENS 
      IF (IEND .GT. MAXCOR) CALL INSMEM("S2HHPS2",IEND,MAXCOR)

      CALL DCOPY(LENS,ICORE(I000),1,ICORE(I010),1)
      IOFF = I010
      JOFF = I000
      DO XIRREP  = 1, NIRREP
         IRREPA  = DIRPRD(XIRREP,SIRREP)

         NUMDIS = VRT(IRREPA,MIXSPIN) * NS(SIRREP)
         IF ((IUHF.NE.0) .AND. (MIXSPIN.EQ.ISPIN)) THEN
            DISSIZ = IRPDPD(XIRREP,20+ISPIN)
         ELSE
            DISSIZ = IRPDPD(XIRREP,ISYTYP(1,LISTS2IN)) 
         ENDIF 
         CALL TRANSP(ICORE(IOFF),ICORE(JOFF),NUMDIS,DISSIZ)
         IOFF = IOFF + NUMDIS * DISSIZ * IINTFP
         JOFF = JOFF + NUMDIS * DISSIZ * IINTFP
       ENDDO
C
C ICORE(I000) CONTAINS ALL [SPIN ADAPTED] S2 COEFFICIENTS [2S(MI,AP) -
C S(IM,AP)], ORDERED AS (AP,IM)

       NDIM = MAX(LENS,NFMI(ISPIN))
       I020 = I010 + NDIM * IINTFP
       IEND = I020 + NDIM * IINTFP
       IF (IEND .GT. MAXCOR) CALL INSMEM("S2HHPS2",IEND,MAXCOR)
C
C Read the entire Hbar(M,J) list and do the multiplication to obtain
C S2(IJ,AP)
C
       CALL GETLST(ICORE(I020),1,1,1,ISPIN,91)
#ifdef _DEBUG_LVL0
       call checksum("fmi     :",icore(i020),nfmi(Ispin))
#endif 
C TRANSPOS HBAR(M,J) -> HBAR(J,M)

        IOFF = I020
        JOFF = I010
        KOFF = 0
        DO IRREP_J = 1, NIRREP
           IRREP_M = DIRPRD(IRREP_J,1)
           NM = POP(IRREP_M,MIXSPIN)
           NJ = POP(IRREP_J,MIXSPIN)
           CALL TRANSP(ICORE(IOFF),ICORE(JOFF),NJ,NM)
           IOFF_HBAR(IRREP_J) = I010 + KOFF 
           IOFF = IOFF + NM*NJ
           JOFF = JOFF + NM*NJ
           KOFF = KOFF + NM*NJ
        ENDDO 

C S2IN(API,M)*H(J,M) -> S2EX(API,J)

       CALL ZERO(ICORE(I020),LENS)
       IS2IN = I000
       IS2EX = I020 
       DO XIRREP = 1, NIRREP 
          IRREPA = DIRPRD(XIRREP,SIRREP)
          NDIM_AP = VRT(IRREPA,MIXSPIN)*NS(SIRREP)

          DO IRREPM = 1, NIRREP
             IRREPI = DIRPRD(IRREPM,XIRREP)
             IRREPJ = IRREPM
             IOFF_S2IN(IRREPM) = IS2IN
             IOFF_S2EX(IRREPJ) = IS2EX
             IS2IN = IS2IN + NDIM_AP*POP(IRREPI,ISPIN)*
     &                               POP(IRREPM,MIXSPIN)*IINTFP
             IS2EX = IS2EX + NDIM_AP*POP(IRREPJ,MIXSPIN)*
     &                               POP(IRREPI,ISPIN)*IINTFP
          ENDDO

          DO IRREPM = 1, NIRREP
             IRREPJ = DIRPRD(IRREPM,1)
             IRREPI = DIRPRD(IRREPM,XIRREP)

             IHBAR = IOFF_HBAR(IRREPM)
             IS2IN = IOFF_S2IN(IRREPM)
             IS2EX = IOFF_S2EX(IRREPJ)

             NROW = POP(IRREPI,ISPIN)*NDIM_AP
             NCOL = POP(IRREPJ,MIXSPIN)
             NSUM = POP(IRREPM,MIXSPIN)

             IF (NROW*NCOL*NSUM .NE. 0) THEN
             CALL XGEMM('N','N',NROW,NCOL,NSUM,ONEM,
     &                   ICORE(IS2IN),NROW,
     &                   ICORE(IHBAR),NCOL,ONE,
     &                   ICORE(IS2EX),NROW)
             ENDIF
          ENDDO
       ENDDO

C AT THE MOMENT S2EX IS S2EX(AP,IJ). TRANSPOSE IT TO S2EX(IJ,AP)

      CALL DCOPY(LENS,ICORE(I020),1,ICORE(I000),1)
      IOFF = I000
      JOFF = I020
      DO XIRREP  = 1, NIRREP
         IRREPA  = DIRPRD(XIRREP,SIRREP)
         NROW = VRT(IRREPA,MIXSPIN) * NS(SIRREP)
         IF ((IUHF.NE.0) .AND. (MIXSPIN.EQ.ISPIN)) THEN
            NCOL = IRPDPD(XIRREP,20+ISPIN)
         ELSE
            NCOL = IRPDPD(XIRREP,ISYTYP(1,LISTS2IN))
         ENDIF 
         CALL TRANSP(ICORE(IOFF),ICORE(JOFF),NCOL,NROW)
         IOFF = IOFF + NROW * NCOL * IINTFP
         JOFF = JOFF + NROW * NCOL * IINTFP
      ENDDO
C
C S2EX(IJ,AP) -> S2EX(JI,AP)
C
      IOFF = I020
      I030 = IEND
      DO XIRREP  = 1, NIRREP 
         IRREPA = DIRPRD(XIRREP,SIRREP)

         NUMDIS = VRT(IRREPA,MIXSPIN) * NS(SIRREP)
         IF ((IUHF.NE.0) .AND. (MIXSPIN.EQ.ISPIN)) THEN
            DISSIZ = IRPDPD(XIRREP,20+ISPIN)
         ELSE
            DISSIZ = IRPDPD(XIRREP,ISYTYP(1,LISTS2IN))
         ENDIF 
         MAXSIZ = MAX(NUMDIS,DISSIZ)

         I040   = I030  + MAXSIZ * IINTFP
         I050   = I040  + MAXSIZ * IINTFP
         IEND   = I050  + MAXSIZ

         IF (IEND .GT. MAXCOR) CALL INSMEM("@-S2HHPS2",I060,MAXCOR)

         CALL SYMTR3(XIRREP,POP(1,ISPIN),POP(1,MIXSPIN),DISSIZ,
     +               NUMDIS,ICORE(IOFF),ICORE(I030),ICORE(I040),
     +               ICORE(I050))
         IOFF = IOFF + NUMDIS * DISSIZ * IINTFP
       ENDDO

C Add this increment to S2 increment lists. 
C
      LISTS2EX = LS2OUT(ISPIN, MIXSPIN + 1 - IUHF)

      IF ((IUHF.NE.0) .AND. (MIXSPIN.EQ.ISPIN)) THEN
          CALL GETEXP2_HHP(ICORE(I010),LENS,NUMSZS,IRPDPD(1,20+ISPIN),
     &                 LISTS2EX,POP(1,MIXSPIN),
     &                 IRPDPD(1,ISYTYP(1,LISTS2EX)))
       ELSE
          CALL GETALLS2_HHP(ICORE(I010),LENS,VRT(1,MIXSPIN),NS(1),
     &                      1,LISTS2EX)
       ENDIF

#ifdef _DEBUG_LVL0
      call checksum("->hs2   :",icore(i010),lens)
      call checksum("s2->s2  :",icore(i020),lens)
#endif
       CALL SAXPY (LENS,ONE,ICORE(I020),1,ICORE(I010),1)

#ifdef _DEBUG_LVL0
      call checksum("s2->fs2 :",icore(i010),lens)
#endif
C
       IF ((IUHF.NE.0).AND. (MIXSPIN.EQ.ISPIN)) THEN
          CALL PUTSQZ_HHP(ICORE(I010),LENS,NUMSZS,IRPDPD(1,20+ISPIN),
     &                    LISTS2EX,POP(1,MIXSPIN), 
     &                    IRPDPD(1,ISYTYP(1,LISTS2EX)),ICORE(I020),
     &                    MAXCOR-I020+1)
       ELSE
          CALL PUTALLS2_HHP(ICORE(I010),LENS,VRT(1,MIXSPIN),NS(1),1,
     $                      LISTS2EX)
       ENDIF
C
C ENDDO for MIXSPIN LOOP
C
      ENDDO
C  
      RETURN
      END
