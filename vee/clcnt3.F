      SUBROUTINE CLCNT3(CORE,MAXCOR,IUHF,ISIDE,IRREPX,EL1R3,EL2R3,EL3R2)
C-----------------------------------------------------------------------
C     This subroutine calculates noniterative triple excitation correct-
C     ions to the EOM-CCSD excitation energy from R increments and L_1
C     and L_2 amplitudes. 
C
C     List expectations :
C
C     401 R(A<B;I<J); 402 R(a<b;i<j); 403 R(Ab;Ij).
C     404 L(A<B;I<J); 405 L(a<b;i<j); 406 L(Ab;Ij).
C     406-406 R_2 increments stored as above.
C     410 Parts 1 and 2 : R_1
C     410 Parts 3 and 4 : L_1
C     410 Parts 5 and 6 : R_1 increments.
C
C     Nowadays this routine is only called with ISIDE=1 (i.e. R increm-
C     ents are contracted with L amplitudes.
C-----------------------------------------------------------------------
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     Arguments.
C-----------------------------------------------------------------------
      DOUBLE PRECISION CORE,EL1R3,EL2R3,EL3R2
      INTEGER MAXCOR,IUHF,ISIDE,IRREPX
C-----------------------------------------------------------------------
C     Common blocks.
C-----------------------------------------------------------------------
      INTEGER IINTLN,IFLTLN,IINTFP,IALONE,IBITWD,
     &        NSTART,NIRREP,IRREPA,IRREPB,DIRPRD,
     &        IRPDPD,ISYTYP,ID,
     &        POP,VRT,NTAA,NTBB,NF1AA,NF2AA,NF1BB,NF2BB
C-----------------------------------------------------------------------
C     Local variables.
C-----------------------------------------------------------------------
      DOUBLE PRECISION ESPAD,E1A,E1B,E2AA,E2BB,E2AB,E1TOT,E2TOT
      INTEGER DISSIZ,NDIS,NSIZ,LNAA,LNBB,LNAB,NEED,
     &        I000,I010,I020,I030,I040,I050,I060,I070,I080,I090,I100,
     &        IRREPR,IRREPL
C-----------------------------------------------------------------------
C     Functions.
C-----------------------------------------------------------------------
      DOUBLE PRECISION SDOT
      INTEGER IDSYMSZ
C-----------------------------------------------------------------------
      DIMENSION CORE(1)
C-----------------------------------------------------------------------
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     &                NF1BB,NF2BB
C-----------------------------------------------------------------------
C
      IF(ISIDE.NE.1)THEN
        WRITE(6,*) ' CLCNT3-F, Invalid value of ISIDE ',ISIDE
        CALL ERREX
      ENDIF
C-----------------------------------------------------------------------
C
      IF(ISIDE.EQ.1)THEN
       IF(IUHF.NE.0)THEN
C-----------------------------------------------------------------------
C      ISIDE=1, UHF.
C-----------------------------------------------------------------------
        I000 = 1
        I010 = I000 + IRPDPD(IRREPX, 9)
        I020 = I010 + IRPDPD(IRREPX, 9)
        I030 = I020 + IRPDPD(IRREPX,10)
        I040 = I030 + IRPDPD(IRREPX,10)
        I050 = I040 + IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
        I060 = I050 + IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
        I070 = I060 + IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
        I080 = I070 + IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
        I090 = I080 + IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
        I100 = I090 + IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
        NEED = I100 * IINTFP
C
        IF(NEED.GT.MAXCOR)THEN
         CALL INSMEM('CLCNT3',NEED,MAXCOR)
        ENDIF
C
        CALL GETLST(CORE(I000),1,1,1,3,410)
        CALL GETLST(CORE(I010),1,1,1,5,410)
        CALL GETLST(CORE(I020),1,1,1,4,410)
        CALL GETLST(CORE(I030),1,1,1,6,410)
C
        CALL GETALL(CORE(I040),1,IRREPX,404)
        CALL GETALL(CORE(I050),1,IRREPX,407)
C
        CALL GETALL(CORE(I060),1,IRREPX,405)
        CALL GETALL(CORE(I070),1,IRREPX,408)
C
        CALL GETALL(CORE(I080),1,IRREPX,406)
        CALL GETALL(CORE(I090),1,IRREPX,409)
C
        LNAA = IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
        LNBB = IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
        LNAB = IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
C
        E1A = 0.0D+00
        E1B = 0.0D+00
        E2AA = 0.0D+00
        E2BB = 0.0D+00
        E2AB = 0.0D+00

CSSS        call checksum("E1A-1", CORE(I000), IRPDPD(IRREPX, 9))
CSSS        call checksum("E1A-1", CORE(I010), IRPDPD(IRREPX, 9))
C
        E1A = E1A + SDOT(IRPDPD(IRREPX, 9),CORE(I000),1,CORE(I010),1)
        E1B = E1B + SDOT(IRPDPD(IRREPX,10),CORE(I020),1,CORE(I030),1)
C
        E2AA = E2AA + SDOT(LNAA,CORE(I040),1,CORE(I050),1)
        E2BB = E2BB + SDOT(LNBB,CORE(I060),1,CORE(I070),1)
        E2AB = E2AB + SDOT(LNAB,CORE(I080),1,CORE(I090),1)
C
#ifdef _DEBUG_LVL0
        WRITE(6,1000) e1a,e1b,e2aa,e2bb,e2ab
#endif
C
C     L1 * Q1 Hbar * R3, L2 * Q2 Hbar * R3
C
        E1TOT = E1A + E1B
        E2TOT = E2AA + E2BB + E2AB
#ifdef _DEBUG_LVL0
        WRITE(6,1010) E1A,E1B,E1TOT,E2AA,E2BB,E2AB,E2TOT
#endif
C
       ELSE
C-----------------------------------------------------------------------
C      ISIDE=1, RHF.
C-----------------------------------------------------------------------
        I000 = 1
        I010 = I000 + IRPDPD(IRREPX, 9)
        I020 = I010 + IRPDPD(IRREPX, 9)
        I030 = I020 + IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
        I040 = I030 + IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
        NEED = I040 * IINTFP
C
        IF(NEED.GT.MAXCOR)THEN
         CALL INSMEM('CLCNT3',NEED,MAXCOR)
        ENDIF
C
        CALL GETLST(CORE(I000),1,1,1,3,410)
CSSS        call sumblk(core(i000),irpdpd(irrepx,9))
        CALL GETLST(CORE(I010),1,1,1,5,410)
CSSS        call sumblk(core(i010),irpdpd(irrepx,9))
C
        LNAB = IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
C
c       CALL GETALL(CORE(I020),1,IRREPX,406)
c       CALL GETALL(CORE(I030),1,IRREPX,409)
        CALL GETALL(CORE(I020),LNAB,IRREPX,406)
CSSS        call sumblk(core(i020),LNAB)
        CALL GETALL(CORE(I030),LNAB,IRREPX,409)
CSCC        call sumblk(core(i030),LNAB)
C
        LNAB = IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
C
        E1A =  0.0D+00
        E1B =  0.0D+00
        E2AA = 0.0D+00
        E2BB = 0.0D+00
        E2AB = 0.0D+00
C
        E1A  = E1A  + SDOT(IRPDPD(IRREPX, 9),CORE(I000),1,CORE(I010),1)
        E2AB = E2AB + SDOT(LNAB,CORE(I020),1,CORE(I030),1)
C
#ifdef _DEBUG_LVL0
        WRITE(6,1000) e1a,e1b,e2aa,e2bb,e2ab
#endif
C
        ESPAD = 0.0D+00
        DO 30 IRREPR=1,NIRREP
C
        IRREPL = DIRPRD(IRREPX,IRREPR)
C
        DISSIZ = IRPDPD(IRREPL,13)
        NDIS   = IRPDPD(IRREPR,14)
        NSIZ   = DISSIZ * NDIS
C
        IF(NSIZ.EQ.0) GOTO 30
        I000 = 1
        I010 = I000 + NSIZ * IINTFP
        I020 = I010 + NSIZ * IINTFP
        I030 = I020 + MAX(DISSIZ,NDIS)
        I040 = I030 + MAX(DISSIZ,NDIS)
        NEED = I040 * IINTFP
C
        IF(NEED.LE.MAXCOR)THEN
         CALL GETLST(CORE(I000),1,NDIS,1,IRREPR,406)
         CALL GETLST(CORE(I010),1,NDIS,1,IRREPR,409)
         CALL SPINAD3(IRREPL,VRT(1,1),DISSIZ,NDIS,CORE(I000),
     &                CORE(I020),CORE(I030))
         ESPAD = ESPAD + SDOT(NSIZ,CORE(I000),1,CORE(I010),1)
        ELSE
         CALL INSMEM('CLCNT3',NEED,MAXCOR)
        ENDIF
   30   CONTINUE
C
#ifdef _DEBUG_LVL0
        WRITE(6,2010) ESPAD
#endif
C
        E1B   = E1A
        E1TOT = E1A + E1B
        E2AA  = (ESPAD - E2AB)/2.0D+00
        E2BB  = E2AA
        E2TOT = ESPAD
#ifdef _DEBUG_LVL0
        WRITE(6,1010) E1A,E1B,E1TOT,E2AA,E2BB,E2AB,E2TOT
#endif
C
       ENDIF
C
      ELSE
C
       IF(IUHF.NE.0)THEN
C-----------------------------------------------------------------------
C      ISIDE=2, UHF.
C-----------------------------------------------------------------------
        I000 = 1
        I010 = I000
        I020 = I010
        I030 = I020
        I040 = I030
        I050 = I040 + IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
        I060 = I050 + IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
        I070 = I060 + IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
        I080 = I070 + IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
        I090 = I080 + IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
        I100 = I090 + IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
        NEED = I100 * IINTFP
C
        IF(NEED.GT.MAXCOR)THEN
         CALL INSMEM('CLCNT3',NEED,MAXCOR)
        ENDIF
C
        CALL GETALL(CORE(I040),1,IRREPX,401)
        CALL GETALL(CORE(I050),1,IRREPX,407)
C
        CALL GETALL(CORE(I060),1,IRREPX,402)
        CALL GETALL(CORE(I070),1,IRREPX,408)
C
        CALL GETALL(CORE(I080),1,IRREPX,403)
        CALL GETALL(CORE(I090),1,IRREPX,409)
C
        LNAA = IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
        LNBB = IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
        LNAB = IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
C
        E1A  = 0.0D+00
        E1B  = 0.0D+00
        E2AA = 0.0D+00
        E2BB = 0.0D+00
        E2AB = 0.0D+00
C
        E2AA = E2AA + SDOT(LNAA,CORE(I040),1,CORE(I050),1)
        E2BB = E2BB + SDOT(LNBB,CORE(I060),1,CORE(I070),1)
        E2AB = E2AB + SDOT(LNAB,CORE(I080),1,CORE(I090),1)
C
#ifdef _DEBUG_LVL0
        WRITE(6,1000) e1a,e1b,e2aa,e2bb,e2ab
#endif
C 
        E2TOT = E2AA + E2BB + E2AB
#ifdef _DEBUG_LVL0
        WRITE(6,3010) E2AA,E2BB,E2AB,E2TOT
#endif`
C
       ELSE
C-----------------------------------------------------------------------
C      ISIDE=2, RHF.
C-----------------------------------------------------------------------
        I000 = 1
        I010 = I000
        I020 = I010
        I030 = I020
        I040 = I030
        I050 = I040
        I060 = I050
        I070 = I060
        I080 = I070
        I090 = I080 + IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
        I100 = I090 + IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
C
        CALL GETALL(CORE(I080),1,IRREPX,403)
        CALL GETALL(CORE(I090),1,IRREPX,409)
C
        LNAB = IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
C
        E1A = 0.0D+00
        E1B = 0.0D+00
        E2AA = 0.0D+00
        E2BB = 0.0D+00
        E2AB = 0.0D+00
C
        E2AB = E2AB + SDOT(LNAB,CORE(I080),1,CORE(I090),1)
C
#ifdef _DEBUG_LVL0
        write(6,1000) e1a,e1b,e2aa,e2bb,e2ab
#endif
C
        ESPAD = 0.0D+00
        DO 130 IRREPR=1,NIRREP
C
        IRREPL = DIRPRD(IRREPX,IRREPR)
C
        DISSIZ = IRPDPD(IRREPL,13)
        NDIS   = IRPDPD(IRREPR,14)
        NSIZ   = DISSIZ * NDIS
C
        IF(NSIZ.EQ.0) GOTO 130
        I000 = 1
        I010 = I000 + NSIZ * IINTFP
        I020 = I010 + NSIZ * IINTFP
        I030 = I020 + MAX(DISSIZ,NDIS)
        I040 = I030 + MAX(DISSIZ,NDIS)
        NEED = I040 * IINTFP
C
        IF(NEED.LT.MAXCOR)THEN
C
         CALL GETLST(CORE(I000),1,NDIS,1,IRREPR,403)
         CALL GETLST(CORE(I010),1,NDIS,1,IRREPR,409)
         CALL SPINAD3(IRREPL,VRT(1,1),DISSIZ,NDIS,CORE(I000),
     &                CORE(I020),CORE(I030))
         ESPAD = ESPAD + SDOT(NSIZ,CORE(I000),1,CORE(I010),1)
C
        ELSE
         CALL INSMEM('CLCNT3',NEED,MAXCOR)
        ENDIF
  130   CONTINUE
C
#ifdef _DEBUG_LVL0
        WRITE(6,2010) ESPAD
        WRITE(6,*)  
#endif
        E2AA  = (ESPAD - E2AB)/2.0D+00
        E2BB  = E2AA
        E2TOT = ESPAD
#ifdef _DEBUG_LVL0
        WRITE(6,3010) E2AA,E2BB,E2AB,E2TOT
#endif 
C
       ENDIF
      ENDIF
C
      IF(ISIDE.EQ.1)THEN
       EL1R3 = E1TOT
       EL2R3 = E2TOT
      ELSE
       EL3R2 = E2TOT
      ENDIF
C
      RETURN
 1000 FORMAT(' e1a  ',F20.12,/,' e1b  ',F20.12,/,' e2aa ',F20.12,/,
     &       ' e2bb ',F20.12,/,' e2ab ',F20.12)
 1010 FORMAT('    Spin case         <L1 Q1 Hbar*R3>   ',/,
     &       '        AA       ',F20.12,/,
     &       '        BB       ',F20.12,/,
     &       '       Total     ',F20.12,/,
     &       '                        ',/,
     &       '    Spin case         <L2 Q2 Hbar*R3>   ',/,
     &       '       AAAA      ',F20.12,/,
     &       '       BBBB      ',F20.12,/,
     &       '       ABAB      ',F20.12,/,
     &       '       Total     ',F20.12)
 1020 FORMAT(' @CLCNT3-F, Insufficient memory. Need ',I12,' Got ',I12)
 2010 FORMAT(' @CLCNT3-I, Spin adapted energy is ',F20.12)
 3010 FORMAT('    Spin case         <L3 Q3 Hbar*R2>   ',/,
     1       '       AAAA      ',F20.12,/,
     1       '       BBBB      ',F20.12,/,
     1       '       ABAB      ',F20.12,/,
     1       '       Total     ',F20.12)
      END
