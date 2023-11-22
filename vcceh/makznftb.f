C     
      SUBROUTINE MAKZNFTB(ZFN, ICORE, MAXCOR, NFVOZ, LISTT1, LISTT2, 
     &                    LISTW1, LISTW2, IUHF, IRREPX, ISPIN)
C
C This routine drives the formation of the Z(F,N)) intermediate
C in the calculation of three-body contribution to the quadratic
C term. 
C
C   F(F,N) = t(E;M)*Hbar(MN,EF) + t(e,m)*Hbar(Nm,Fe)  [ISPIN = 1]
C   F(e,m) = t(e;m)*Hbar(mn,ef) + t(E;M)*Hbar(Mn,Ef)  [ISPIN = 2]
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,ZFN
      DIMENSION ICORE(MAXCOR), ZFN(NFVOZ)
C      
      COMMON /INFO/ NOCA,NOCB,NVRTA,NVRTB
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                    DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
C   
C Common blocks used in the quadratic term
C 
      COMMON /QADINTI/ INTIMI, INTIAE, INTIME
      COMMON /QADINTG/ INGMNAA, INGMNBB, INGMNAB, INGABAA,
     &                 INGABBB, INGABAB, INGMCAA, INGMCBB,
     &                 INGMCABAB, INGMCBABA, INGMCBAAB,
     &                 INGMCABBA
      COMMON /APERTT1/ IAPRT1AA             
      COMMON /BPERTT1/ IBPRT1AA
      COMMON /APERTT2/ IAPRT2AA1, IAPRT2BB1, IAPRT2AB1, IAPRT2AB2, 
     &                 IAPRT2AB3, IAPRT2AB4, IAPRT2AA2, IAPRT2BB2,
     &                 IAPRT2AB5 

      DATA ONE /1.0/ 
      DATA ONEM /-1.0/
      DATA ZILCH /0.0/
C
C We only need to consider one DPD irrep; the irrepx.
C 
      NSIZT1  = IRPDPD(IRREPX, 8 + ISPIN)
      NSIZT2  = IRPDPD(IRREPX, 11 - ISPIN)
C
C I000 Holds first T1 in equation above.
C I001 Holds second T1 in equation above.
C I002 Holds the Hbar integral lists.
C
      I000 = 1
C     
      IF (IUHF .EQ. 0) THEN
         I001 = I000
      ELSE
         I001 = I000 + NSIZT1*IINTFP
      ENDIF
C
      I002 = I001 + NSIZT2*IINTFP
C     
      CALL GETLST(ICORE(I000), 1, 1, 1, ISPIN, LISTT1)
C
      IF(IUHF .NE. 0) CALL GETLST(ICORE(I001), 1, 1, 1, 3-ISPIN,
     &                            LISTT2)
C
C Irrep of multiplication is IRREPFN
C      
      IRREPFN = IRREPX
      IRREPEM = IRREPX
C
      NWDSZ1 = IRPDPD(IRREPEM, ISYTYP(1, LISTW1))
      NWDIS1 = IRPDPD(IRREPFN, ISYTYP(2, LISTW1))
      NWDSZ2 = IRPDPD(IRREPEM, ISYTYP(1, LISTW2))
      NWDIS2 = IRPDPD(IRREPFN, ISYTYP(2, LISTW2))
C     
      I003 = I002 + NWDSZ1*NWDIS1*IINTFP
C
      IF(I003 .GT. MAXCOR) CALL INSMEM('MAKZNFTB', I003, MAXCOR)
C
      CALL GETLST(ICORE(I002), 1, NWDIS1, 2, IRREPFN, LISTW1)
C
      CALL XGEMM('N', 'N', 1, NFVOZ, NSIZT1, ONE, ICORE(I000),
     &            1, ICORE(I002), NWDSZ1, ZILCH, ZFN, 1)
C
      I003 = I002 + NWDSZ2*NWDIS2*IINTFP
C
      IF(I003 .GT. MAXCOR) CALL INSMEM('MAKZNFTB', I003, MAXCOR)
C     
      CALL GETLST(ICORE(I002), 1, NWDIS2, 2, IRREPFN, LISTW2)
C
      CALL XGEMM('N','N', 1, NFVOZ, NSIZT2, ONE, ICORE(I001),
     &            1, ICORE(I002), NWDSZ2, ONE, ZFN, 1)
 100  CONTINUE
C
      RETURN
      END
