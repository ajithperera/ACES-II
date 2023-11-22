C 
      SUBROUTINE DBLRESORT(ICORE, MAXCOR, IRREPX, IUHF)
C
C This subroutine drives the formation of the resorted T(Ab,Ij)
C and T(AB,IJ) amplitudes and  write them to disk.
C
C For RHF, the (AI,bj), (Aj,bI) and (AJ,BI) packed lists are written.
C For UHF, the (bj,AI), (bI,Aj) and (ai,bj) packed lists are also
C written.
C
      IMPLICIT INTEGER (A-Z)
      CHARACTER*4 SPCASE(2),SPCAS2,RERDTP
      CHARACTER*11 STRNG(9)
      DOUBLE PRECISION SDOT
      DIMENSION ICORE(MAXCOR),LIST(6)
C
      COMMON /MACHSP/ NSTART,NIRREP,IINTFP,IALONE,IBITWD
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /FLAGS/IFLAGS(100)
      COMMON /FILES/LUOUT, MOINTS
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

      DATA SPCASE /'AABB','BBAA'/
C
      NNM1O2(I) = (I*(I-1))/2
      PRTLEVL = IFLAGS(1)
C
      NOCA  = NOCCO(1)
      NOCB  = NOCCO(2)
      NVRTA = NVRTO(1)
      NVRTB = NVRTO(2)
C
      LIST(1) = IAPRT2AB2
      LIST(2) = IAPRT2AB1
      LIST(3) = IAPRT2AB4
      LIST(4) = IAPRT2AB3
      LIST(5) = IAPRT2AA1
      LIST(6) = IAPRT2BB1
C
      STRNG(1) = ' (AI,bj) = '
      STRNG(2) = ' (bj,AI) = '
      STRNG(3) = ' (Aj,bI) = '
      STRNG(4) = ' (bI,Aj) = '
      STRNG(5) = ' (AJ,BI) = '
      STRNG(6) = ' (aj,bi) = '
      STRNG(7) = ' (Ab,Ij) = '
      STRNG(8) = ' (AB,IJ) = '
      STRNG(9) = ' (ab,ij) = '
C
      REFLST = IAPRT2AB5
      RERDTP = 'AJBI'
C
      ISZTOT  = IDSYMSZ(IRREPX, 13, 14)
      NSCRSZ  = NOCA*NOCB + NVRTA*NVRTB + NVRTA*NOCA + NVRTB*NOCB
      ISZTAR  = IDSYMSZ(IRREPX, 9, 10)
      ISZTAR2 = IDSYMSZ(IRREPX, 11, 12)
C
      I000 = 1
      I010 = I000 + IINTFP*MAX(ISZTOT, ISZTAR, ISZTAR2)
      I020 = I010 + IINTFP*MAX(ISZTOT, ISZTAR, ISZTAR2)
      I030 = I020 + NSCRSZ
C
      IF(I030 .GT. MAXCOR) CALL INSMEM('DBLRESORT', I030, MAXCOR)
C
C Write (AI,bj) and (Aj,bI) orderings and (bj,AI) and (bI,Aj) if UHF.
C
      DO 20 I = 1, 1 + IUHF
C
         CALL GETALL(ICORE(I000), ISZTOT, IRREPX, REFLST)
C
         IF (PRTLEVL .GE. 20) THEN
            CALL HEADER('Checksum of original and transformed lists', 
     &                   0, LUOUT)
            WRITE(LUOUT, *) STRNG(7), SDOT(ISZTOT, ICORE(I000), 1,
     &                      ICORE(I000), 1)
         ENDIF
C
         CALL ALTSYMPCK1(ICORE(I000), ICORE(I010), ISZTOT, ISZTAR,
     &                   ICORE(I020), IRREPX, SPCASE(I))
C     
         CALL PUTALL(ICORE(I010), ISZTAR, IRREPX, LIST(I))
C
         IF (PRTLEVL .GE. 20) THEN
            WRITE(LUOUT, *) STRNG(I), SDOT(ISZTAR, ICORE(I010), 1,
     &                      ICORE(I010), 1)
         ENDIF
C
         CALL ALTSYMPCK0(ICORE(I010), ICORE(I000), ISZTAR, ISZTAR2,
     &                   ICORE(I020), IRREPX, SPCASE(I))
C
         CALL PUTALL(ICORE(I000), ISZTAR2, IRREPX, LIST(2+I))
C
         IF (PRTLEVL .GE. 20) THEN
            WRITE(LUOUT, *) STRNG(2+I), SDOT(ISZTAR2, ICORE(I000), 1,
     &                      ICORE(I000), 1)
         ENDIF
C
 20   CONTINUE
C
C Write (AI,BJ) integrals (ai,bj) for UHF if ISPIN = 2
C
      SPCAS2 = 'AAAA'
C
      DO 30 ISPIN = 1, 1 + IUHF
C
         ISZTOT = IDSYMSZ(IRREPX, ISPIN, ISPIN + 2)
         ISZTAR = IDSYMSZ(IRREPX, 8 + ISPIN, 8 + ISPIN)
         NSCRSZ = NNM1O2(NVRTO(ISPIN)) + NNM1O2(NOCCO(ISPIN)) +
     &            NVRTO(ISPIN)*NOCCO(ISPIN)
C
         I000 = 1
         I010 = I000 + IINTFP*MAX(ISZTOT, ISZTAR)
         I020 = I010 + IINTFP*ISZTAR
         I030 = I020 + NSCRSZ
C
         IF(I030 .GT. MAXCOR) CALL INSMEM('DBLRESORT', I030, MAXCOR)
C
         CALL GETALL(ICORE(I000), ISZTOT, IRREPX, (REFLST - 3 + ISPIN))
C
         IF (PRTLEVL .GE. 20) THEN
            CALL HEADER('Checksum of original and transformed lists', 
     &                   0, LUOUT)
            
            WRITE(LUOUT, *) STRNG(7+ISPIN), SDOT(ISZTOT, ICORE(I000),
     &                      1, ICORE(I000), 1)
         ENDIF
C
         CALL ALTSYMPCK2(ICORE(I000), ICORE(I010), ISZTOT, ISZTAR, 
     &                   ICORE(I020), IRREPX, SPCAS2, RERDTP)
C
         CALL PUTALL(ICORE(I010), ISZTAR, IRREPX, LIST(4 + ISPIN))
C
         IF(PRTLEVL .GE. 20) THEN
            WRITE(LUOUT, *) STRNG(4+ISPIN), SDOT(ISZTAR, ICORE(I010),
     &                      1, ICORE(I010), 1)
         ENDIF
C
         SPCAS2 = 'BBBB'
C
30    CONTINUE
C
      RETURN
      END
