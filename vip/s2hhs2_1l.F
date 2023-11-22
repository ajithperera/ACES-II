      SUBROUTINE S2HHS2_1L(ICORE,MAXCOR,IUHF, ISPIN)
C
C   THIS SUBROUTINE CALCULATES THE CONTRIBUTION OF S2 TO
C   S2 INVOLVING THE HH BLOCK OF THE EFFECTIVE ONE-PARTICLE 
C   INTEGRALS.
C
C CONTRACTION :
C
C     Z(iJ,ap) = -SUM [Hbar(m,i) (2S(Jm,ap) - S(mJ,ap)) [ RHF CASE ]
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

      Write(6,"(a)") "-----------------S2HHS2_1R------------------"
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
C TRANSPOSE S2IN(JM,AP) -> S2IN(AP,JM). 

      IEND=I010+LENS*IINTFP
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
C ICORE(I000) CONTAINS ALL [SPIN ADAPTED] S2 COEFFICIENTS [2S(JM,AP) -
C S(MJ,AP)], ORDERED AS (AP,JM)

       NDIM = MAX(LENS,NFMI(ISPIN))
       I020 = I010 + NDIM * IINTFP
       IEND = I020 + NDIM * IINTFP
       IF (IEND .GT. MAXCOR) CALL INSMEM("S2HHPS2",IEND,MAXCOR)
C
C Read the entire Hbar(M,I) list and do the multiplication to obtain
C S2(IJ,AP)
C
       CALL GETLST(ICORE(I020),1,1,1,ISPIN,91)
#ifdef _DEBUG_LVL0
       call checksum("fmi     :",icore(i020),nfmi(Ispin))
#endif 
C TRANSPOSE HBAR(M,I) -> HBAR(I,M)

       IOFF = I020
       JOFF = I010
       KOFF = 0
       DO IRREP_I = 1, NIRREP
          IRREP_M = DIRPRD(IRREP_I,1)
          NM = POP(IRREP_M,ISPIN)
          NI = POP(IRREP_I,ISPIN)
          IOFF_HBAR(IRREP_I) = I010 + KOFF
          CALL TRANSP(ICORE(IOFF),ICORE(JOFF),NI,NM)
          IOFF = IOFF + NM*NI
          JOFF = JOFF + NM*NI
          KOFF = KOFF + NM*NI
       ENDDO 

C S2IN(AP,JM)*HBAR(I,M) -> S2EX(AP,JI)

       CALL ZERO(ICORE(I020),LENS)
       IS2IN = I000
       IS2EX = I020
       DO XIRREP = 1, NIRREP 
          IRREPA = DIRPRD(XIRREP,SIRREP)
          NDIM_AP = VRT(IRREPA,MIXSPIN)*NS(SIRREP)

          DO IRREPM = 1, NIRREP
             IRREPJ = DIRPRD(IRREPM,XIRREP)
             IRREPI = IRREPM
             IOFF_S2IN(IRREPM) = IS2IN
             IOFF_S2EX(IRREPI) = IS2EX
             IS2IN = IS2IN + NDIM_AP*POP(IRREPJ,MIXSPIN)*
     &                               POP(IRREPM,ISPIN)*IINTFP
             IS2EX = IS2EX + NDIM_AP*POP(IRREPJ,MIXSPIN)*
     &                               POP(IRREPI,ISPIN)*IINTFP
          ENDDO 

          DO IRREPM = 1, NIRREP
             IRREPI = DIRPRD(IRREPM,1)
             IRREPJ = DIRPRD(IRREPM,XIRREP)

             IHBAR = IOFF_HBAR(IRREPM)
             IS2IN = IOFF_S2IN(IRREPM)
             IS2EX = IOFF_S2EX(IRREPI)

             NROW = POP(IRREPJ,MIXSPIN)*NDIM_AP
             NCOL = POP(IRREPI,ISPIN)
             NSUM = POP(IRREPM,ISPIN)

             IF (NROW*NCOL*NSUM .NE. 0) THEN
                CALL XGEMM('N','N',NROW,NCOL,NSUM,ONEM,
     &                      ICORE(IS2IN),NROW,
     &                      ICORE(IHBAR),NCOL,ONE,
     &                      ICORE(IS2EX),NROW)
             ENDIF
          ENDDO
       ENDDO

C AT THE MOMENT S2EX IS S2EX(AP,JI). TRANSPOSE IT TO S2EX(JI,AP)

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
#ifdef _DEBUG_LVLM
        if (iuhf.eq.0) call check_ab2aa(icore(i010),lens,icore(iend),
     &                                  maxcor,1)
#endif
C
C ENDDO for MIXSPIN LOOP
C
      ENDDO
C  
      RETURN
      END
