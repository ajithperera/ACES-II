C
C *****************************************************************
C  GENERAL PROCEDURES RELATED TO GETTING LENGTHS OF S-VECTORS
C *****************************************************************
C         
      SUBROUTINE GETLEN(LEN, POP1, POP2, POP3, POP4)
C
C THIS ROUTINE GETS THE TOTAL LENGTH OF A GENERAL FOUR-INDEX ARRAY.
C
C INPUT: THE POPULATION OVER THE FOUR INDICES.
C OUTPUT: LEN, THE LENGTH OF THE ARRAY.
C
      IMPLICIT INTEGER (A-Z)
C
      DIMENSION POP1(8), POP2(8), POP3(8), POP4(8)
C

      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREP0(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/INFO/NOCCO(2),NVRTO(2)
C
      LEN=0
      DO 10 XIRREP = 1, NIRREP
         DO 20 IRREP3 = 1, NIRREP
            IRREP4 = DIRPRD(IRREP3, XIRREP)
            N34 = POP3(IRREP3)*POP4(IRREP4)
            DO 30 IRREP2 = 1, NIRREP
               IRREP1=DIRPRD(IRREP2,XIRREP)
               N12 =POP1(IRREP1)*POP2(IRREP2)
               LEN = LEN + N12*N34
 30         CONTINUE
 20      CONTINUE
 10      CONTINUE
         RETURN
         END