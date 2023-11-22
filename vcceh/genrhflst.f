C
      SUBROUTINE GENRHFLST(ICORE, MAXCOR, IUHF)
C
C This subroutine creates several lists which are required in
C RHF calculations and calls only for RHF calculations. 
C The description of the lists  being created is given 
C prior to each subroutine calls which creat the particular lists.
C
      IMPLICIT INTEGER(A-Z)
C
      DIMENSION ICORE(MAXCOR)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /FILES/ LUOUT, MOINTS
C
C For RHF runs the list number 56 contains symmetry adapted
C Hbar integrals (2*Hbar(Mb,Ej) - Hbar(Mb,Je)). The following
C call create the Hbar(Mb,Ej) non symmetry adapted integrals. It also
C creates the Hbar(MB,EJ) lists (list No. 54).
C
      CALL QRESET(ICORE, MAXCOR, IUHF)
C
C Create the Habr(CI,AB) (list No. 27) which is required in
C in the calculation of T1 to singles contributions (T1QIAE)
C
C
      CALL MKRHFLST(ICORE, MAXCOR, 1, IUHF, VRT(1, 1), 19, 9,
     &              0, 30, 27)
C
      RETURN
      END
