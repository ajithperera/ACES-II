C
      SUBROUTINE T1QIME(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF, ISPIN)
C
C This routine drives the calculation of the I(ME) intermediate.
C
C  I(m,e) =  Sum t(n,f)*Hbar(mn,ef)
C               n,f
C Spin integrated formulas are as given below. 
C
C UHF
C
C  F(E,M) = t(F,N)*Hbar(MN,EF) + t(f,n)*Hbar(Mn,Ef)  [AA]
C  F(e,m) = t(f,n)*Hbar(mn,ef) + t(F,N)*Hbar(mN,eF)  [BB]
C
C RHF
C
C  F(E,M) = t(F,N)*Hbar(MN,EF) + t(f,n)*Hbar(Mn,Ef) 
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,SDOT
      CHARACTER*6 SPCASE(2)
      DIMENSION ICORE(MAXCOR)
C      
      COMMON /INFO/ NOCA,NOCB,NVRTA,NVRTB
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
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
C
      DATA ONE /1.0/ 
      DATA ONEM /-1.0/
      DATA ZILCH /0.0/
      DATA SPCASE /'AA =  ', 'BB =  '/
C
      LISTT1 = IAPRT1AA
      LISTT2 = IAPRT1AA
      LISTW1 = 18 + ISPIN
      LISTW2 = 16 + ISPIN
      LISTF  = INTIME
      IF(IUHF .EQ. 0) LISTW2 = 18
C
C We only need to consider one DPD irrep; the irrepx.
C
      NTOTTAR = IRPDPD(IRREPX, 8 + ISPIN)
      NSIZT1  = IRPDPD(IRREPX, 8 + ISPIN)
      NSIZT2  = IRPDPD(IRREPX, 11 - ISPIN)
C
C I000 Holds F target in equation above.
C I010 Holds first T in equation above.
C I020 Holds second T in equation above.
C I030 Holds individual IRREPX of blocked integral list.
C
      I000 = 1
      I010 = I000 + NTOTTAR*IINTFP
C
      IF (IUHF .EQ. 0) THEN
         I020 = I010
      ELSE 
         I020 = I010 + NSIZT1*IINTFP
      ENDIF
C
      I030 = I020 + NSIZT2*IINTFP
C      
      CALL GETLST(ICORE(I010), 1, 1, 1, ISPIN, LISTT1)
      IF(IUHF .NE. 0) CALL GETLST(ICORE(I020), 1, 1, 1, 3 - ISPIN,
     &                            LISTT2)
C
C Irrep of multiplication is IRREPFN
C      
      IRREPFN = IRREPX
      IRREPEM = IRREPFN
C
      NWDSZ1 = IRPDPD(IRREPFN, ISYTYP(1, LISTW1))
      NWDIS1 = IRPDPD(IRREPEM, ISYTYP(2, LISTW1))
      NWDSZ2 = IRPDPD(IRREPFN, ISYTYP(1, LISTW2))
      NWDIS2 = IRPDPD(IRREPEM, ISYTYP(2, LISTW2))
      NSIZF  = NSIZT1
C     
      I040 = I030 + NWDSZ1*NWDIS1*IINTFP
C
      IF(I040 .GT. MAXCOR) CALL INSMEM('T1QIME', I040, MAXCOR)
C
      CALL GETLST(ICORE(I030), 1, NWDIS1, 2, IRREPEM, LISTW1)
C
      CALL XGEMM('N', 'N', 1, NTOTTAR, NSIZT1, ONE, ICORE(I010),
     &            1, ICORE(I030), NWDSZ1, ZILCH, ICORE(I000), 1)
C
      I040 = I030 + NWDSZ2*NWDIS2*IINTFP
C
      IF(I040 .GT. MAXCOR) CALL INSMEM('MAKFME', I040, MAXCOR)
C     
      CALL GETLST(ICORE(I030), 1, NWDIS2, 2, IRREPEM, LISTW2)
C
      CALL XGEMM('N','N', 1, NTOTTAR, NSIZT2, ONE, ICORE(I020),
     &            1, ICORE(I030), NWDSZ2, ONE, ICORE(I000), 1)
 100  CONTINUE
C
      CALL PUTLST(ICORE(I000), 1, 1, 1, ISPIN, LISTF)
C
      IF (IFLAGS(1) .GE. 20) THEN
C
         IF (ISPIN .EQ. 1) THEN
            NSIZE = IRPDPD(IRREPX, 9)
         ELSE
            NSIZE = IRPDPD(IRREPX, 10)
         ENDIF
C
         CALL HEADER('Checksum @-T1QIME', 0, LUOUT)
            
         WRITE(LUOUT, *) SPCASE(ISPIN), SDOT(NSIZE, ICORE(I000),
     &                    1, ICORE(I000), 1)
      ENDIF
C      
      RETURN
      END
