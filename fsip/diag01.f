      SUBROUTINE DIAG01(HEFF,SCR,MAXCOR,IUHF,NBAS)
C
C DIAGONALIZE EFFECTIVE HAMILTONIAN FOR 0,1 SECTOR (ELECTRON
C REMOVAL) FOCK SPACE CALCULATIONS.  DIAGONAL PART OF ONE-BODY
C HBAR IS ADDED HERE. 
C
CEND
      IMPLICIT INTEGER  (A-H,O-Z)
      DOUBLE PRECISION HEFF(1),SCR(1),EREF,EFINAL
      CHARACTER*5 SPCASE(2)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /FSSYMPOP/ FSDPDAN(8,22),FSDPDNA(8,22),FSDPDAA(8,22),
     &                  FSDPDIN(8,22),FSDPDNI(8,22),FSDPDII(8,22),
     &                  FSDPDAI(8,22),FSDPDIA(8,22)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /FSSYM/ POPA(8,2),POPI(8,2),VRTA(8,2),VRTI(8,2),
     &               NTAN(2),NTNA(2),NTAA(2),
     &               NTIN(2),NTNI(2),NTII(2),
     &               NTAI(2),NTIA(2),
     &               NF1AN(2),NF1AA(2),NF1IN(2),NF1II(2),NF1AI(2),
     &               NF2AN(2),NF2AA(2),NF2IN(2),NF2II(2),NF2AI(2)
      COMMON /INFO/ NOCCO(2),NVRTO(2) 
      COMMON /REFSTATE/ EREF,efinal
C
      DATA SPCASE /'alpha',' beta'/
C
      IOFF0=1
      DO 10 ISPIN=1,1+IUHF
       CALL FSGETT1(HEFF(IOFF0),ISPIN+2,91,'AA',20+ISPIN)
       WRITE(6,1000)SPCASE(ISPIN)
1000   FORMAT(T3,'--- Spin case ',A,' ---')
       WRITE(6,501)
501    FORMAT(T5,'Symmetry',T20,'Root',T30,'E (IP)',T44,'E(final)',
     &        T58,'Dominant MO')
       DO 20 IRREP=1,NIRREP
        NDIM =POPA(IRREP,ISPIN)
        NSIZE=NDIM*NDIM
        IOFF=IOFF0
        IF(NDIM.NE.0)THEN
         DO 30 IROW=1,NDIM
c          WRITE(6,1002)(HEFF(I),I=IOFF,IOFF+NDIM-1)
1002      FORMAT((8F10.5))
          IOFF=IOFF+NDIM
30       CONTINUE
         I000=1
         I010=I000+NSIZE
         I020=I010+NDIM
         I030=I020+NDIM
         I035=I030+NSIZE
         I040=I035+NDIM*NDIM 
         CALL SCOPY(NSIZE,HEFF(IOFF0),1,SCR(I000),1)
CSSS
CSSS         CALL RG(NDIM,NDIM,SCR(I000),SCR(I010),SCR(I020),1,
CSSS     &           SCR(I030),SCR(I035),SCR(I040),IERR)
CSSS
        LWORK = (MAXCOR - I040 + 1)
        CALL MN_GEEV(NDIM,NDIM,SCR(I000),SCR(I010),SCR(I020),
     &                SCR(I030),SCR(I035),LWORK,IERR)

         ILOCEVC=I030
         ILOCEVL=I010
         DO 500 I=1,NDIM
          iposrel=isamax(ndim,scr(ilocevc),1)
          iposabs=iposrel+ilocevc-1
          write(6,1010)irrep,i,scr(ilocevl),eref+scr(ilocevl),
     &                 abs(scr(iposabs)),iposrel
1010      format(t8,i2,t20,i3,t24,f14.8,t39,f14.8,t57,f7.5,t65,'[',i3,
     &           ']')
          efinal=eref+scr(ilocevl)
          ILOCEVL=ILOCEVL+1
          ILOCEVC=ILOCEVC+NDIM
500      CONTINUE
1003     FORMAT(T3,2F20.10)
C
        ENDIF
        IOFF0=IOFF0+NSIZE
20     CONTINUE
C
10    CONTINUE
C
      RETURN
      END
