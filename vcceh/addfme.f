C
      SUBROUTINE ADDFME(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
C
C This subroutine computes the part of the I intermediates which 
C is due Hbar(M,E) elements and then augments the I intermediates
C with these values.
C
C  - (1/2) SUM Hbar(E,M)*T1(A,M) (For case F(A,E))
C                         M
C  + (1/2) SUM Hbar(E,M)*T1(E,I) (For case F(M,I))
C                         E
C Contribution from the bare Hbar elements is absent for 
C the quadratic contribution. Spin interated formulas 
C are straightforward and not explicitly given.
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION HALF,MHALF,ZILCH,ONE,SDOT
      CHARACTER*6 SPCASE(2)
      DIMENSION ICORE(MAXCOR),IOFFT(8,2), ITAROFF1(8),
     &          ITAROFF2(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FILES/ LUOUT, MOINTS
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
      DATA HALF  /0.5/
      DATA MHALF /-0.5/
      DATA ZILCH /0.0/
      DATA ONE  /1.0/
      DATA SPCASE /'AA =  ', 'BB =  '/
C
C Get T1(Alpha) amplitudes
C
      LIST1 = IAPRT1AA
C
      CALL GETPERT1(ICORE, MAXCOR, MXCOR, IUHF, IRREPX, LIST1,
     &              IOFFSET, IOFFT)
C
C Loop over spin cases
C
      DO 5 ISPIN = 1, 1 + IUHF
C     
C Determine the size of the symmetry-packed Hbar(M,E) elements.
C
         SIZFOK = IRPDPD(1, 8 + ISPIN)
         SIZTAR = IRPDPD(IRREPX, 18 + ISPIN)
C
C Do I(E,A) intermediates.
C
         I000 = 1
         I010 = I000 + SIZFOK*IINTFP
         I020 = I010 + SIZTAR*IINTFP
C
         IOFFFME = I000
         ITAROFF1(1) = I010
C
         CALL GETLST(ICORE(I000), 1, 1, 1, ISPIN, 93)
         CALL ZERO(ICORE(I010), SIZTAR)
C     
C Set up an offset array for the target
C
         DO 5000 IRREPAO = 2, NIRREP
C
            IRREP = IRREPAO - 1
            IRREPE = DIRPRD(IRREP, IRREPX)

            ITAROFF1(IRREPAO) = ITAROFF1(IRREPAO - 1) + VRT(IRREP,
     &                         ISPIN)*VRT(IRREPE, ISPIN)*IINTFP
 5000    CONTINUE
C     
         DO 10 IRREPM = 1, NIRREP
C
            IRREPHBM = IRREPM
            IRREPHBE = IRREPHBM
            IRREPT1M = IRREPHBM
            IRREPT1A = DIRPRD(IRREPT1M, IRREPX)
C
            NVRTE = VRT(IRREPHBE, ISPIN)
            NOCCM = POP(IRREPT1M, ISPIN)
            NVRTA = VRT(IRREPT1A, ISPIN)
C
            IOFFT1 = IOFFT(IRREPT1A, ISPIN)
            IOFFTAR = ITAROFF1(IRREPT1A)
C         
            CALL XGEMM('N', 'T', NVRTE, NVRTA, NOCCM, MHALF, 
     &                  ICORE(IOFFFME), NVRTE, ICORE(IOFFT1), NVRTA,
     &                  ZILCH, ICORE(IOFFTAR), NVRTE)
C     
            IOFFFME = IOFFFME + NOCCM*NVRTE*IINTFP
C
 10      CONTINUE
C
         IF (IFLAGS(1) .GE. 20) THEN
            NSIZE = SIZTAR
            CALL HEADER('Checksum @-ADDFME-IAE', 0, LUOUT)
            WRITE(LUOUT, *) SPCASE(ISPIN), SDOT(NSIZE, 
     &                      ICORE(I010), 1, ICORE(I010), 1)
         ENDIF
C
         I030 = I020 + SIZTAR*IINTFP
C
C Up date the I(A,E) contribution on the disk. Notice that
C we dont have to add bare Hbar(A,E) to the intermediate.
C
         CALL GETLST(ICORE(I020), 1, 1, 1, ISPIN, INTIAE)
         CALL SAXPY(SIZTAR, ONE, ICORE(I010), 1, ICORE(I020), 1)
         CALL PUTLST(ICORE(I020), 1, 1, 1, ISPIN, INTIAE)
C
C         IF (IFLAGS(1) .GE. 20) THEN
C            NSIZE = SIZTAR
C            CALL HEADER('Checksum @-ADDFME-IAE', 0, LUOUT)
C            WRITE(LUOUT, *) SPCASE(ISPIN), SDOT(NSIZE, 
C     &                      ICORE(I020), 1, ICORE(I020), 1)
C         ENDIF
C     
C Do F(M,I) intermediates.
C     
         SIZTAR = IRPDPD(IRREPX, 20 + ISPIN)
C
         I000 = 1
         I010 = I000 + SIZFOK*IINTFP
         I020 = I010 + SIZTAR*IINTFP
C
         IOFFFME = I000
         ITAROFF2(1) = I010
C
C Set an offset array for the target
C     
         DO 5100 IRREPI = 2, NIRREP
C
            IRREP = IRREPI - 1
            IRREPM = DIRPRD(IRREP, IRREPX)

            ITAROFF2(IRREPI) = ITAROFF2(IRREPI - 1) + POP(IRREP,
     &                         ISPIN)*POP(IRREPM, ISPIN)*IINTFP
 5100    CONTINUE
C
         CALL GETLST(ICORE(I000), 1, 1, 1, ISPIN, 93)
         CALL ZERO(ICORE(I010), SIZTAR)
C
         DO 20 IRREPE = 1, NIRREP
C
            IRREPHBE = IRREPE
            IRREPHBM = IRREPHBE
            IRREPT1E = IRREPHBE
            IRREPT1I = DIRPRD(IRREPT1E, IRREPX)
C
            NVRTE = VRT(IRREPHBE, ISPIN)
            NOCCM = POP(IRREPHBM, ISPIN)
            NOCCI = POP(IRREPT1I, ISPIN)
C
            IOFFT1  = IOFFT(IRREPT1E, ISPIN)
            IOFFTAR = ITAROFF2(IRREPT1I)
C
            IF (NVRTE .NE. 0) THEN 
            CALL XGEMM('T', 'N', NOCCM, NOCCI, NVRTE, HALF,   
     &                  ICORE(IOFFFME), NVRTE, ICORE(IOFFT1), NVRTE,
     &                  ZILCH, ICORE(IOFFTAR), NOCCM)
C
            IOFFFME = IOFFFME + NOCCM*NVRTE*IINTFP
            ELSE
            CALL ZERO(ICORE(IOFFTAR),NOCCM*NOCCI)
            ENDIF 

 20      CONTINUE
C     
         I030 = I020 + SIZTAR*IINTFP
C
C Up date the I(M,I) contribution on the disk. Notice that
C we dont have to add bare Hbar(M,I) to the intermediate.
C
         CALL GETLST(ICORE(I020), 1, 1, 1, ISPIN, INTIMI)
         CALL SAXPY(SIZTAR, ONE, ICORE(I010), 1, ICORE(I020), 1)
         CALL PUTLST(ICORE(I020), 1, 1, 1, ISPIN, INTIMI)
C     
C Update the F(MI) intermediate lists
C
         IF (IFLAGS(1) .GE. 20) THEN
C     
            NSIZE = SIZTAR
C
            CALL HEADER('Checksum @-ADDFME-IMI', 0, LUOUT)
            
            WRITE(LUOUT, *) SPCASE(ISPIN), SDOT(NSIZE, 
     &                      ICORE(I020), 1, ICORE(I020), 1)
         ENDIF
C
 5    CONTINUE
C
      RETURN
      END
