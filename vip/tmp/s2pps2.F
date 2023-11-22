      SUBROUTINE S2PPS2(ICORE,MAXCOR,IUHF, ISIDE, ISPIN)
C
C   THIS SUBROUTINE CALCULATES THE CONTRIBUTION OF S2 TO
C   S2 INVOLVING THE PP BLOCK OF THE EFFECTIVE ONE-PARTICLE 
C   INTEGRALS.
C
C CONTRACTION :
C
C  Z(ji,ap) = SUM [Hbar(a,e) * S2(ji,ep)]  S2 is spin adapted RHF
C             e 
C  Ispin = spin-label of e and Mixspin = spin-lable of j,e
C
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
      DATA ONE  /1.0/
      DATA ZILCH/0.0/
      DATA ONEM/-1.0/

       Write(6,"(a)") "-----------------S2PPS2---------------------"

      DO 5 XIRREP = 1, NIRREP
         IRREPE = DIRPRD(XIRREP,SIRREP)
         NUMSZS(XIRREP) = VRT(IRREPE,ISPIN) * NS(SIRREP)
 5    CONTINUE      
C
      DO MIXSPIN = 1, 1+IUHF
         LISTS2IN = LS2IN(ISPIN, MIXSPIN + 1 - IUHF)
         CALL GETLEN_HHP(LENS,POP(1,MIXSPIN),POP(1,ISPIN),
     $                   VRT(1,MIXSPIN),NS)

       I000=1
       I010=I000+LENS*IINTFP
C
C AAA and BBB blocks 
C
       IF ((IUHF.NE.0) .AND. (MIXSPIN.EQ.ISPIN)) THEN
          CALL GETEXP2_HHP(ICORE(I000), LENS,NUMSZS,
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
       NDIM = MAX(NFEA(MIXSPIN),LENS)
       I020 = I010 + NDIM * IINTFP
       IEND = I020 + NDIM * IINTFP
       IF (IEND .GE. MAXCOR) CALL INSMEM("S2PPS2",IEND,MAXCOR)
C
C Read the entire Hbar(A,E) list and do the multiplication to obtain
C S2(JI,EP) * Hbar(A,E) = S2(JI,AP)
C
       CALL GETLST(ICORE(I010), 1,1,1,MIXSPIN,92)
#ifdef _DEBUG_LVL0
       call checksum("Hbae    :",icore(i010),nfea(Ispin))
#endif
C Offser of HBAR(E,A) (NOTE THAT HBAR(A,E) IS STORED AS HBAR(E,A))

       IOFF_HBAR(1) = I010
       DO IRREP_A = 1, NIRREP-1
           IRREP_E = DIRPRD(IRREP_A,1)
           NA = VRT(IRREP_A,MIXSPIN)
           NE = VRT(IRREP_E,MIXSPIN)
           IOFF_HBAR(IRREP_A+1) = IOFF_HBAR(IRREP_A)+NA*NE
        ENDDO

C S2IN(IJP,E)*HBAR(E,A) --> S2EX(IJP,A) = S2EX(JI,AP)

      CALL ZERO(ICORE(I020),LENS)
      IOFF = I000
      JOFF = I020
      DO XIRREP = 1, NIRREP
         NROW   = IRPDPD(XIRREP,ISYTYP(1,LISTS2IN))*NS(SIRREP)
         IRREPE = DIRPRD(XIRREP,SIRREP)
         NSUM   = VRT(IRREPE,MIXSPIN)

         DO IRREPA = 1, NIRREP
            NCOL = VRT(IRREPA,MIXSPIN)
            IF (IRREPA .EQ. IRREPE) THEN

               IOFF_EA  = IOFF_HBAR(IRREPE) 
               IOFF_JIP = IOFF 
               IOFF_JIA = JOFF

#ifdef _DEBUG_LVLM
         Write(6,'(a,6i2)'), "NCOL,NSUM,NROW,IRREPE,IRREPA,XIRREP=", 
     +            NCOL,NSUM,NROW,IRREPE,IRREPA,XIRREP
         Write(6,'(a,6i4)'), "IOFF_EA,IOFF_JIP,IOFF_JIA=",
     +                        IOFF_EA,IOFF_JIP,IOFF_JIA
#endif 
               CALL VV_CONTRCT(ICORE(IOFF_JIP),ICORE(IOFF_EA),
     +                         ICORE(IOFF_JIA),NROW,NCOL,
     +                         NSUM)
            ENDIF 

            IOFF_EA  = IOFF_EA  + NSUM * NCOL * IINTFP 
            IOFF_JIP = IOFF_JIP + NROW * NSUM * IINTFP 
            IOFF_JIA = IOFF_JIA + NROW * NCOL * IINTFP 

         ENDDO 
         IOFF = IOFF + NROW * NSUM * IINTFP
         JOFF = JOFF + NROW * NSUM * IINTFP
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
     &                   LISTS2EX,POP(1,MIXSPIN), 
     &                   IRPDPD(1,ISYTYP(1,LISTS2EX)),ICORE(I020),
     &                   MAXCOR-I020+1)
      ELSE
         CALL PUTALLS2_HHP(ICORE(I010),LENS,VRT(1,MIXSPIN),NS(1),1,
     $                     LISTS2EX)
      ENDIF
C
C ENDDO for MIXSPIN LOOP
C
      ENDDO
C  
      RETURN
      END

