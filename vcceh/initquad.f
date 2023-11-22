C
      SUBROUTINE INITQUAD(ICORE, MAXCOR, IUHF)
C
C This routine initiate the calculation of quadratic contribution
C to EOM-CCSD second-order properties. This assign list numbers
C required in quadratic term.
C
      IMPLICIT INTEGER(A-Z) 
C     
      DIMENSION ICORE(MAXCOR)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
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
      COMMON /AOPARTT1/ IAOTT1AA1, IAOTT1BB1, IAOTT1AB1, IAOTT1AA2, 
     &                  IAOTT1BB2, IAOTT1AB2, IAOTT1AA3, IAOTT1BB3, 
     &                  IAOTT1AB3
C
C I(MI), I(AE) and I(ME) intermediates
C
      INTIMI = 430
      INTIAE = 431
      INTIME = 432
C
C G(MN,IJ) intermediate. Left index is MN and right index is IJ.
C
      INGMNAA = 435
      INGMNBB = 436
      INGMNAB = 437
C
C G(AB,EF) intermediate. Left index is AB and right index is EF.
C No attempt is made to store these integrals on a disk. 
C
      INGABAA = 438
      INGABBB = 439
      INGABAB = 440
C
C G(MB,EJ) intermediate. Left index is EM and right index is BJ.
C 
      INGMCAA   = 441 
      INGMCBB   = 442
      INGMCABAB = 443
      INGMCBABA = 444
      INGMCBAAB = 445
      INGMCABBA = 446
C
C Perturb T1 (T1(ALPHA) and T1(BETA)) amplitudes. 
C need separate lists for both perturb amplitudes.
C  
      IAPRT1AA = 447
      IBPRT1AA = 448
C
C Perturb T2 (T2(ALPHA) and T2(BETA)) amplitudes. 
C In the actual calculation we can overwrite alpha perturb
C T2 amplitudes by beta perturb T2 amplitudes.
C
      IAPRT2AA1 = 449
      IAPRT2BB1 = 450
      IAPRT2AB1 = 451
      IAPRT2AB2 = 452
      IAPRT2AB3 = 453
      IAPRT2AB4 = 454
C
      IAPRT2AA2 = 455
      IAPRT2BB2 = 456
      IAPRT2AB5 = 457
C
C Perturbed T1 and T1 product (T1 tau) for AO ladder, 05/2020
C
      IF (IFLAGS(93) .EQ. 2) THEN

         IAOTT1AA1 = 470
         IAOTT1BB1 = 471
         IAOTT1AB1 = 472

         IAOTT1AA2 = 473
         IAOTT1BB2 = 474
         IAOTT1AB2 = 475

         IAOTT1AA3 = 476
         IAOTT1BB3 = 477
         IAOTT1AB3 = 478

      ENDIF 

      RETURN
      END
