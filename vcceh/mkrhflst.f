C
      SUBROUTINE MKRHFLST(SCR, MAXCOR, IRREPX, IUHF, POP1, SYTYPL, 
     &                    SYTYPR, INEWFIL, LISTIN, LISTTAR)
C
C This subroutine creates antisymmetrized lists for RHF calculations
C for ABAB lists. This routine is used in T1QIAE create antisymmetrized
C lists used in RHF calculations.
C
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION SCR
      DIMENSION SCR(MAXCOR), POP1(8)
C      
      COMMON /INFO/ NOCA,NOCB,NVRTA,NVRTB
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /FLAGS/ IFLAGS(100)
C
C Create the pointers for the target list. Call the INIPCK
C with INEWFIL = 0 since these files are already exist.
C
      CALL INIPCK(IRREPX, SYTYPL, SYTYPR, LISTTAR, INEWFIL, 0, 1)
C
      DO 10 IRREPR = 1, NIRREP
C
         IRREPL = DIRPRD(IRREPR, IRREPX)
C     
         DISSYW = IRPDPD(IRREPL, ISYTYP(1, LISTIN))
         NUMSYW = IRPDPD(IRREPR, ISYTYP(2, LISTIN))
C
C Allocate memory for the input lists
C
         I000 = 1
         I010 = I000 + DISSYW*NUMSYW
C
         CALL GETLST(SCR(I000), 1, NUMSYW, 2, IRREPR, LISTIN)
C
C Antisymmetrizes the integral just read in. Notice the in place
C antisymmetrization.
C
         CALL ASSYM2A(IRREPL, POP1, DISSYW, NUMSYW, SCR(I000), 
     &                SCR(I010 + NUMSYW), SCR(I010 + 2*NUMSYW))
C
C Put the antisymmetrized integrals in the newly created lists.
C
         CALL PUTLST(SCR(I000), 1, NUMSYW, 2, IRREPR, LISTTAR)
C
 10   CONTINUE
C
      RETURN
      END
