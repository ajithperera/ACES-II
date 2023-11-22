      SUBROUTINE S1HPS2(ICORE,MAXCOR,IUHF, ISIDE, ISPIN)
C
C   THIS SUBROUTINE CALCULATES THE CONTRIBUTION OF S1 TO
C   S2 INVOLVING THE HP BLOCK OF THE EFFECTIVE ONE-PARTICLE 
C   INTEGRALS.
C
C CONTRACTION :
C
C     Z(JI,AP) = Hbar(J,A) S1(I,P)   [ only for iside = 2 ]
C                
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      logical check
C
      DIMENSION ICORE(MAXCOR), NUMSZS(8)
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
C
      IF (ISIDE .EQ. 1) THEN
         RETURN
      ENDIF
      Write(6,"(a)") "-----------------S1HPS2-------------------"
C
      DO 5 XIRREP = 1, NIRREP
         IRREPA = DIRPRD(XIRREP,SIRREP)
         NUMSZS(XIRREP) = VRT(IRREPA,ISPIN) * NS(SIRREP)
 5    CONTINUE      
C
      DO 25 MIXSPIN = 1, 1+IUHF
         LISTS2EX = LS2OUT(ISPIN, MIXSPIN + 1 - IUHF)
         CALL GETLEN_HHP(LENS, POP(1,MIXSPIN), POP(1,ISPIN),
     $                   VRT(1,MIXSPIN), NS)
       I000 = 1
       I010 = I000 + LENS * IINTFP
       I020 = I010 + LENS * IINTFP
       I030 = I020 + NT(MIXSPIN) * IINTFP
       I040 = I030 + POP(SIRREP,ISPIN) * NS(SIRREP) * IINTFP
       IF (I040 .GT. MAXCOR) THEN
          CALL INSMEM("S1HPS2",I040,MAXCOR)
          CALL ERREX
       ENDIF
C
C DO IN CORE ALGORITHM
C
        CALL ZERO(ICORE(I000), LENS)
        CALL GETLST(ICORE(I020), 1, 1, 1, MIXSPIN, 93)
        CALL GETLST(ICORE(I030), 1, 1, 1, ISPIN, LS1IN)
#ifdef _DEBUG_LVL0
      ndim=pop(sirrep,ispin)*ns(sirrep)
      call checksum("->s1    :",icore(i030),ndim)
      call checksum("->faj   :",icore(i020),nt(mixspin))
#endif
C TRA
C
C  Hbar(A,J) stored as AJ
C  S(AJ,IP) = SUM Hbar(AJ,1) S(1,IP) 
C             
        NROW = NT(MIXSPIN)
        NCOL = POP(SIRREP,ISPIN) * NS(SIRREP)
        NSUM = 1
        CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(I020),NROW,
     &              ICORE(I030),NSUM,ONE,ICORE(I000),NROW)
C
C REORDER THE S COEFFICIENT ACCORDING TO S(AJ,IP) -> S(IJ,AP)
C
       CALL SSTGEN(ICORE(I000),ICORE(I010),LENS,VRT(1,MIXSPIN),
     $             POP(1,MIXSPIN),POP(1,ISPIN),NS,ICORE(I020),1,
     $             '3214')
C#ifdef _NOSKIP
C
C Transpose the first two indices. 
C
       IF ((IUHF.EQ.0) .OR. (MIXSPIN.NE.ISPIN)) THEN
          IOFF = I010
          DO XIRREP  = 1, NIRREP
             IRREPA  = DIRPRD(XIRREP,SIRREP) 
             IRREPJI = XIRREP
             NUMDSS  = VRT(IRREPA,MIXSPIN) * NS(SIRREP)
             DISSYS  = IRPDPD(IRREPJI,14)

             I030 = I020 + MAX(NUMDSS,DISSYS)
             I040 = I030 + MAX(NUMDSS,DISSYS)
             IEND = I040 + MAX(NUMDSS,DISSYS)
             write(6,"(6(1x,i4))") xirrep,irrepa,irrepji,numdss,
     +                             dissys

             IF (IEND .GT. MAXCOR) THEN
                CALL INSMEM("@-S1HPS2",IEND,MAXCOR)
             ENDIF
             CALL SYMTR3(IRREPJI,POP(1,ISPIN),POP(1,MIXSPIN),  
     +                   DISSYS,NUMDSS,ICORE(IOFF),ICORE(I020),
     +                   ICORE(I030),ICORE(I040))
             IOFF = IOFF + NUMDSS*DISSYS*IINTFP
          ENDDO 
       ELSE
          CALL DSCAL(LENS, -1.0D0, ICORE(I010), 1)
       ENDIF 
C#endif 
#ifdef _DEBUG_LVL0
      call checksum("s1->s2  :",icore(i010),lens)
      call output(icore(i010),1,1,1,lens,1,lens,1)
#endif

C Antisymmetrize The S2(JI,AP)

       IF ((IUHF.NE.0) .AND. (MIXSPIN.EQ.ISPIN)) THEN
          CALL ASSYMALL_HHP(ICORE(I010),LENS,NUMSZS,IRPDPD(1,20+ISPIN),
     $        POP(1,ISPIN),ICORE(I020),MAXCOR-I020+1)
          CALL GETEXP2_HHP(ICORE(I000),LENS,NUMSZS, IRPDPD(1,20+ISPIN),
     $       LISTS2EX,POP(1,ISPIN),IRPDPD(1,ISYTYP(1, LISTS2EX)))
       ELSE
          CALL GETALLS2_HHP(ICORE(I000),LENS,VRT(1,MIXSPIN),NS(1),1,
     $                      LISTS2EX)
       ENDIF
#ifdef _DEBUG_LVL0
       if ((check) .and. (mspin .eq. ISPIN)) then
          write(6,*) ' correction in s1uhps2 '
          call printall(icore(i010), lens, irpdpd(1, 18+ispin),
     $       numszs)
          write(6,*) ' input s2 '
          call printall(icore(i000), lens, irpdpd(1, 18+ispin),
     $       numszs)
       endif
#endif 
C
#ifdef _DEBUG_LVL0
      call checksum("->s2    :",icore(i000),lens)
      call output(icore(i000),1,1,1,lens,1,lens,1)
#endif 
       CALL SAXPY(LENS,ONE,ICORE(I010),1,ICORE(I000),1)
C
#ifdef _DEBUG_LVL0
       if ((check) .and. (mspin .eq. ISPIN)) then
          write(6,*) ' corrected s-vector'
          call printall(icore(i000), lens, irpdpd(1, 18+ispin),
     $       numszs)
       endif
#endif 
C
       IF ((IUHF.NE.0).AND. (MIXSPIN.EQ.ISPIN)) THEN
          CALL PUTSQZ_HHP(ICORE(I000),LENS,NUMSZS,IRPDPD(1,20+ISPIN),
     $       LISTS2EX,POP(1,ISPIN),IRPDPD(1,ISYTYP(1,LISTS2EX)),
     $       ICORE(I020),MAXCOR-I020+1)

       ELSE
          CALL PUTALLS2_HHP(ICORE(I000),LENS,VRT(1,MIXSPIN),NS(1),1,
     $         LISTS2EX)
       ENDIF       
#ifdef _DEBUG_LVL0
      call checksum("s1->s2f :",icore(i000),lens)
      call output(icore(i000),1,1,1,lens,1,lens,1)
#endif 
C
 25   CONTINUE
C
      RETURN
      END
