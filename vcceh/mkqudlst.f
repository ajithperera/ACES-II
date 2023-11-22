C
      SUBROUTINE MKQUDLST(ICORE, MAXCOR, IUHF, IRREPX)
C
C Create lists used in EOM-CCSD property calculations (quadratic term)
C
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22), ISYTYP(2,500), ID(18)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FILES/LUOUT,MOINTS
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
C Create pointers for occupied-occupied, virtual-virtual and
C occupied-virtual singles lists. First determine the maximum size of
C each the lists.
C
c      INEWFIL = 1
      call aces_io_remove(54,'DERGAM')
      INEWFIL = 0
C
      DO 10 ISPIN = 1, IUHF + 1
C
         NFVV = IRPDPD(IRREPX, 18 + ISPIN)
         NFOO = IRPDPD(IRREPX, 20 + ISPIN)
         NFVO = IRPDPD(IRREPX, 8  + ISPIN)
C
C The very first time list 400 is addressed, INEWFIL = 1, 0 otherwise.
C
         CALL UPDMOI(1, NFVV, ISPIN, INTIAE, INEWFIL, 0)
C
         INEWFIL = 0
C
         CALL UPDMOI(1, NFOO, ISPIN, INTIMI, INEWFIL, 0)
         CALL UPDMOI(1, NFVO, ISPIN, INTIME, INEWFIL, 0)
C
C Create pointers for perturb T1(Alpha) and T1(Beta) amplitudes.
C
         CALL UPDMOI(1, NFVO, ISPIN, IAPRT1AA, INEWFIL, 0)
         CALL UPDMOI(1, NFVO, ISPIN, IBPRT1AA, INEWFIL, 0)
 10   CONTINUE
C
C Create list pointers for G(MN,IJ), G(AB, EF), G(MB,EJ) intermediates.
C First determine the dimension of these intermediates.
C First do the AAAA and BBBB lists.
C
      DO 20 ISPIN = 1, 1 + IUHF
C
         LISTGMN = (INGMNAA - 1) + ISPIN
         LISTGAB = (INGABAA - 1) + ISPIN
         LISTGMB = (INGMCAA - 1) + ISPIN
C
         SYTYPGMN = 2 + ISPIN
         SYTYPGIJ = 2 + ISPIN
         SYTYPGAB = ISPIN
         SYTYPGEF = ISPIN
         SYTYPGEM = 8 + ISPIN
         SYTYPGBJ = 8 + ISPIN
C
         CALL INIPCK(IRREPX, SYTYPGMN, SYTYPGIJ, LISTGMN, INEWFIL,
     &               0, 1)
C
         CALL INIPCK(IRREPX, SYTYPGAB, SYTYPGEF, LISTGAB, INEWFIL,
     &               0, 1)
C
         CALL INIPCK(IRREPX, SYTYPGEM, SYTYPGBJ, LISTGMB, INEWFIL,
     &               0, 1)
C
 20   CONTINUE
C
C Now do MOIO pointers for ABAB lists for G(MN,IJ), G(AB,EF) and
C G(MB,EJ) intermediate. In addition to that this also create
C pointers for BAAB list for G(MB,EJ) intermediate.
C
      SYTYPGMN  = 14
      SYTYPGIJ  = 14
      SYTYPGAB  = 13
      SYTYPGEF  = 13
      SYTYPGEM1 = 9
      SYTYPGBJ1 = 10
      SYTYPGEM2 = 11
      SYTYPGBJ2 = 11
C
      CALL INIPCK(IRREPX, SYTYPGMN, SYTYPGIJ, INGMNAB, INEWFIL, 0, 1)
C
      CALL INIPCK(IRREPX, SYTYPGAB, SYTYPGEF, INGABAB, INEWFIL, 0, 1)
C
      CALL INIPCK(IRREPX, SYTYPGEM1, SYTYPGBJ1, INGMCABAB, INEWFIL,
     &            0, 1)
C
      CALL INIPCK(IRREPX, SYTYPGEM2, SYTYPGBJ2, INGMCBAAB, INEWFIL,
     &            0, 1)
C
C There are two additional lists left for G(MB,EF) intermediate
C for UHF cases, namely BABA and ABBA spin cases.
C
      IF (IUHF .NE. 0) THEN
C
         SYTYPGEM1 = 10
         SYTYPGBJ1 = 9
         SYTYPGEM2 = 12
         SYTYPGBJ2 = 12
C
         CALL INIPCK(IRREPX, SYTYPGEM1, SYTYPGBJ1, INGMCBABA, INEWFIL,
     &               0, 1)
C
         CALL INIPCK(IRREPX, SYTYPGEM2, SYTYPGBJ2, INGMCABBA, INEWFIL,
     &               0, 1)
C
      ENDIF
C
C Create list pointers for perturb T2 lists. Fist do the AAAA and
C BBBB spin cases. There are no  BBBB lists for RHF calculations.
C
      DO 60 ISPIN = 1, IUHF + 1
C
         LISTT2A = (IAPRT2AA1 - 1) + ISPIN
         LISTT2B = (IAPRT2AA2 - 1) + ISPIN
C
         SYTYPTAJ = 8 + ISPIN
         SYTYPTBI = 8 + ISPIN
         SYTYPTAB = ISPIN
         SYTYPTIJ = 2 + ISPIN
C
         CALL INIPCK(IRREPX, SYTYPTAJ, SYTYPTBI, LISTT2A, INEWFIL,
     &               0, 1)
C
         CALL INIPCK(IRREPX, SYTYPTAB, SYTYPTIJ, LISTT2B, INEWFIL,
     &               0, 1)
C
 60   CONTINUE
C
C Now do the ABAB T2 lists which are required for both RHF and
C UHF calculations.
C
      SYTYPTAI = 9
      SYTYPTBJ = 10
      SYTYPTAJ = 11
      SYTYPTBI = 12
      SYTYPTAB = 13
      SYTYPTIJ = 14
C
      CALL INIPCK(IRREPX, SYTYPTAI, SYTYPTBJ, IAPRT2AB2, INEWFIL, 0, 1)
C
      CALL INIPCK(IRREPX, SYTYPTAJ, SYTYPTBI, IAPRT2AB4, INEWFIL, 0, 1)
C
      CALL INIPCK(IRREPX, SYTYPTAB, SYTYPTIJ, IAPRT2AB5, INEWFIL, 0, 1)
C
C Set the MOIO pointers ABAB T2 list which requires only for UHF
C calculations.
C
      IF (IUHF .NE. 0) THEN
C
         SYTYPTBJ = 10
         SYTYPTAI = 9
         SYTYPTBI = 12
         SYTYPTAJ = 11
C
         CALL INIPCK(IRREPX, SYTYPTBJ, SYTYPTAI, IAPRT2AB1, INEWFIL,
     &               0, 1)
C
         CALL INIPCK(IRREPX, SYTYPTBI, SYTYPTAJ, IAPRT2AB3, INEWFIL,
     &               0, 1)
C
      ENDIF
C
C Create two scratch lists as a temporary storge area for
C ABAB and BABA intermediate contributions which is ordered as
C (Ej,bM) [ISPIN = 1] and (eJ,Bm) [ISPIN = 2]. In the original CC
C code these contributions are stored in the lists which are structured as
C (EM,bj) and (em,BJ) by changing the list structure with a NEWTYP call.
C The subroutine NEWTYP works only for the totally symmetric perturbations
C and also a dangerous routine since it changes the structure of
C existing lists. The following two scratch lists are used to store
C the intermediate. For further information see MKT1GMBEJC1 and
C MKT1GMBEJC2 subroutines.
C
      CALL INIPCK(IRREPX, 11, 12, 460, INEWFIL, 0, 1)
      CALL INIPCK(IRREPX, 12, 11, 461, INEWFIL, 0, 1)

C The AO ladders added on 05/2019.
  
      IF (IFLAGS(93) .EQ. 2) THEN

         SYTYPTAB = 13
         SYTYPTIJ = 14
         IRREPT2  = 1

         CALL INIPCK(IRREPT2, SYTYPTAB, SYTYPTIJ, IAOTT1AB1, INEWFIL, 
     &               0, 1)

         SYTYPTAB = 15
         CALL INIPCK(IRREPT2, SYTYPTAB, SYTYPTIJ, IAOTT1AB2, INEWFIL, 
     &               0, 1)
         CALL INIPCK(IRREPT2, SYTYPTAB, SYTYPTIJ, IAOTT1AB3, INEWFIL, 
     &               0, 1)

         IF (IUHF .NE. 0) THEN

            SYTYPTAB = 1
            SYTYPTIJ = 3
         
            CALL INIPCK(IRREPT2, SYTYPTAB, SYTYPTIJ, IAOTT1AA1, INEWFIL,
     &                  0, 1)

            SYTYPTAB = 15

            CALL INIPCK(IRREPT2, SYTYPTAB, SYTYPTIJ, IAOTT1AA2, INEWFIL,
     &                  0, 1)
            CALL INIPCK(IRREPT2, SYTYPTAB, SYTYPTIJ, IAOTT1AA3, INEWFIL,
     &                  0, 1)
C
            SYTYPTAB = 2
            SYTYPTIJ = 4

            CALL INIPCK(IRREPT2, SYTYPTAB, SYTYPTIJ, IAOTT1BB1, INEWFIL,
     &                  0, 1)

            SYTYPTAB = 15

            CALL INIPCK(IRREPT2, SYTYPTAB, SYTYPTIJ, IAOTT1BB2, INEWFIL,
     &                  0, 1)
            CALL INIPCK(IRREPT2, SYTYPTAB, SYTYPTIJ, IAOTT1BB3, INEWFIL,
     &                  0, 1)
         ENDIF

      ENDIF 

C This will finish the setting up the pointers for lists used
C in the quadratic term.
C
      RETURN
      END
