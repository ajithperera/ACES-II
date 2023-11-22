      SUBROUTINE PCCD_OV(DOO,DVV,DVO,XIA,WORK,MAXCOR,IUHF,
     &                   NONHF_TERMS_EXIST)
C
C DRIVER FOR THE FOLLOWING CONTRIBUTION TO THE Htau_OV
C
C       X(ia) =   SUM E,M (<IE//AM> + <IM//AE>) D(E,M)
C       X(ia) = - SUM M,N <IM//NA> D(MN)
C       X(ia) =   SUM E,F <AE//FI/> D(EF)
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION DOO,DVO,DVV
      LOGICAL NONHF_TERMS_EXIST
      DIMENSION WORK(MAXCOR),DOO(1),DVO(1),DVV(1),XIA(1)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPP(255,2),DIRPRD(8,8)
      COMMON /SYM/  POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)

      Call Dzero(Xia,Nt(1)+Iuhf*Nt(2))

      Call Pccd_xint1(Xia,Doo,Work,Maxcor,Iuhf)
      CALL Pccd_xint2(Xia,Dvv,Work,Maxcor,Iuhf)
      If (Nonhf_terms_exist) Call Pccd_xint3(Xia,Dvo,Work,
     +                                       Maxcor,Iuhf)
#ifdef _DEBUG_LVLM
      call checksum("Xov    :",Xia,Nt(1))
#endif
      Return
      End
