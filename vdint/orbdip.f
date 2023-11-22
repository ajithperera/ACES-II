      SUBROUTINE ORBDIP(DFX,DFY,DFZ,DIPGX,DIPGY,DIPGZ)
C
C 6-Jun-1985 hjaaj & tuh
C
C Purpose:
C   To add the orbital dipole gradients in DIPGX, DIPGY, DIPGZ
C
C Input:
C   The total dipole Fock matrices DFX, DFY, DFZ
C
C Output:
C   Dipole orbital gradients in DIPGX, DIPGY, DIPGZ
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      PARAMETER (D2 = 2.0D0)
      DIMENSION DFX(NOCCT,*), DFY(NOCCT,*), DFZ(NOCCT,*),
     *          DIPGX(*), DIPGY(*), DIPGZ(*)
C
C Used from common blocks:
C   INFVAR : NWOPT,JWOP(2,*)
C   INFORB : NOCCT
C
      COMMON /INFVAR/ NCONF,NWOPT,NVAR,JWOPSY,NWOP(8),NWOPH,NVARH,
     *                JWOP(2,3 000),KLWOP(60,150)
      COMMON /INFORB/ MULD2H(8,8), NRHF(8),NFRO(8),
     *       NISH(8),NASH(8),NSSH(8),NOCC(8),NORB(8),NBAS(8),
     *       NNORB(8),NNBAS(8), N2ORB(8),N2BAS(8),
     *       IISH(8),IASH(8),ISSH(8),IOCC(8),IORB(8),IBAS(8),
     *       IIISH(8),IIASH(8),IIORB(8),IIBAS(8),I2ORB(8),I2BAS(8),
     *       ICMO(8), NSYM,
     *       NISHT,NASHT,NSSHT,NOCCT,NORBT,NBAST,NCMOT,NRHFT,
     *       N2ISHX,NNASHX,N2ASHX,NNASHY,NNOCCX,N2OCCX,
     *       NNORBT,NNORBX,N2ORBT,N2ORBX,NNBAST,N2BAST,NNBASX,N2BASX,
     *       NAS1(8), NAS2(8), NAS3(8),
     *       NTINT,NT1AM,INTBUF,NO2,NT2,NV2,NT2IND,NNVIRX,NVIRT,NT2AM
C
      DO 100 IG = 1, NWOPT
         K = JWOP(1,IG)
         L = JWOP(2,IG)
         DIPGX(IG) = DIPGX(IG) + D2*DFX(K,L)
         DIPGY(IG) = DIPGY(IG) + D2*DFY(K,L)
         DIPGZ(IG) = DIPGZ(IG) + D2*DFZ(K,L)
         IF (L .LE. NOCCT) THEN
            DIPGX(IG) = DIPGX(IG) - D2*DFX(L,K)
            DIPGY(IG) = DIPGY(IG) - D2*DFY(L,K)
            DIPGZ(IG) = DIPGZ(IG) - D2*DFZ(L,K)
         END IF
  100 CONTINUE
      RETURN
      END
