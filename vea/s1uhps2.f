      SUBROUTINE S1UHPS2(ICORE,MAXCOR,IUHF, ISIDE, ISPIN)
C
C   THIS SUBROUTINE CALCULATES THE CONTRIBUTION OF S1 TO
C   S2 INVOLVING THE HP BLOCK OF THE EFFECTIVE ONE-PARTICLE 
C   INTEGRALS.
C
C CONTRACTION :
C
C     Z(ab,ip) = U(ia) S1(b,p)   [ only for iside = 2 ]
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
C
      DO 5 XIRREP = 1, NIRREP
         MIRREP = DIRPRD(XIRREP, SIRREP)
         NUMSZS(XIRREP) = POP(MIRREP, ISPIN) * NS(SIRREP)
 5    CONTINUE      
C
      DO 25 SINSPIN = 1, 1+IUHF
         MSPIN = SINSPIN
         LISTS2EX = LS2OUT(ISPIN, SINSPIN + 1 - IUHF)
       CALL GETLEN(LENS, VRT(1,SINSPIN), VRT(1,ISPIN),
     $      POP(1,SINSPIN), NS)
       I000 = 1
       I010 = I000 + LENS * IINTFP
       I020 = I010 + LENS * IINTFP
       I030 = I020 + NT(MSPIN) * IINTFP
       I040 = I030 + VRT(SIRREP, ISPIN) * NS(SIRREP) * IINTFP
       IF (I040 .GT. MAXCOR) THEN
          WRITE(6,*) ' INSUFFICIENT MEMORY, S1UHPS2'
          CALL ERREX
       ENDIF
C
C   DO IN CORE ALGORITHM
C
          CALL ZERO(ICORE(I000), LENS)
          CALL GETLST(ICORE(I020), 1, 1, 1, MSPIN, 93)
          CALL GETLST(ICORE(I030), 1, 1, 1, ISPIN, LS1IN)
C
C DO MULTIPLICATION
C
C  S(AI, BP) = SUM U(AI, 1) S(1, BP)
C               ..
C
          NROW = NT(MSPIN)
          NCOL = VRT(SIRREP, ISPIN) * NS(SIRREP)
          NSUM = 1
          CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(I020),NROW,
     &       ICORE(I030),NSUM,ONE,ICORE(I000),NROW)
C
C REORDER THE S COEFFICIENT ACCORDING TO S(AM, BP) -> S(AB, MP)
C
       CALL SSTGEN(ICORE(I000), ICORE(I010), LENS, VRT(1,MSPIN),
     $    POP(1,MSPIN), VRT(1,ISPIN), NS, ICORE(I020), 1, '1324')
C
       IF ((IUHF.NE.0) .AND. (SINSPIN.EQ.ISPIN)) THEN
          CALL ASSYMALL(ICORE(I010), LENS, NUMSZS, IRPDPD(1, 18+ISPIN),
     $        VRT(1, ISPIN), ICORE(I020), MAXCOR-I020+1)
          CALL GETEXP2(ICORE(I000), LENS, NUMSZS, IRPDPD(1, 18+ISPIN),
     $       LISTS2EX, VRT(1, ISPIN), IRPDPD(1, ISYTYP(1, LISTS2EX)))
       ELSE
          CALL GETALLS2(ICORE(I000), LENS,POP(1,SINSPIN),NS(1),1,
     $       LISTS2EX)
       ENDIF
       check = .false.
       if ((check) .and. (mspin .eq. ISPIN)) then
          write(6,*) ' correction in s1uhps2 '
          call printall(icore(i010), lens, irpdpd(1, 18+ispin),
     $       numszs)
          write(6,*) ' input s2 '
          call printall(icore(i000), lens, irpdpd(1, 18+ispin),
     $       numszs)
       endif
C
       CALL SAXPY(LENS, ONE, ICORE(I010), 1, ICORE(I000), 1)
C
       if ((check) .and. (mspin .eq. ISPIN)) then
          write(6,*) ' corrected s-vector'
          call printall(icore(i000), lens, irpdpd(1, 18+ispin),
     $       numszs)
       endif
C
       IF ((IUHF.NE.0).AND. (MSPIN.EQ.ISPIN)) THEN
          CALL PUTSQZ(ICORE(I000), LENS, NUMSZS, IRPDPD(1, 18+ISPIN),
     $       LISTS2EX, VRT(1, ISPIN), IRPDPD(1, ISYTYP(1, LISTS2EX)),
     $       ICORE(I020), MAXCOR-I020+1)
          if (check) then
          CALL GETEXP2(ICORE(I000), LENS, NUMSZS, IRPDPD(1, 18+ISPIN),
     $       LISTS2EX, VRT(1, ISPIN), IRPDPD(1, ISYTYP(1, LISTS2EX)))
          write(6,*) ' after putsqz and inverse operation: '
             call printall(icore(i000), lens, irpdpd(1, 18+ispin),
     $       numszs)
          endif
       ELSE
          CALL PUTALLS2(ICORE(I000),LENS,POP(1,MSPIN), NS(1),1,
     $       LISTS2EX)
       ENDIF       
C
 25   CONTINUE
C
      RETURN
      END
