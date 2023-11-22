      SUBROUTINE PCCD_ROTG(Grd,Grd_oo,Grd_vv,Grd_vo,Grd_ov,
     +                     Grad_stata,Grad_statb,B_like,Icore,
     +                     Maxcor,Lenoo,Lenvv,Lenvo,Nocc,Nvrt,
     +                     Nbas,Ispin,Iuhf,Ncycle,Tol)
C
C THIS ROUTINE ROTATES THE MOLECULAR ORBITALS BY THE T1-LIKE
C  AMPLITUDES ON LISTS 3,90 (AND 4,90 FOR UHF) CALCULATIONS.  THIS
C  ROUTINE GENERATES A NEW SET OF ORTHONORMALIZED MOLECULAR ORBITALS.
C
C Extended to work with dropmo. Ajith Perera, 07/2005.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      DOUBLE PRECISION GRD(NBAS,NBAS)
      DOUBLE PRECISION GRD_OO(LENOO)
      DOUBLE PRECISION GRD_VV(LENVV)
      DOUBLE PRECISION GRD_VO(LENVO)
      DOUBLE PRECISION GRD_OV(LENVO)
      DOUBLE PRECISION GRAD_STATA(6)
      DOUBLE PRECISION GRAD_STATB(6)
      DOUBLE PRECISION DSUM,DDOT,TOL
      CHARACTER*1 ISP(2)
      LOGICAL B_LIKE 
      LOGICAL DIIS
      DIMENSION ICORE(MAXCOR)
      DIMENSION NOCC_EXPND(8,2),NVRT_EXPND(8,2),
     &          NOCCO_EXPND(2),NVRTO_EXPND(2),
     &          NBF4IRREP(8)
      COMMON /SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF(4)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON/EXTRAPO/DIIS
C
      DATA ONE  /1.0/
      DATA ZILCH/0.0/
      DATA ONEM /-1.0/
      DATA ISP /'A','B'/
C
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "--------------------pCCD_rotg--------------------"
#endif 
      CALL GETREC(20, "JOBARC", 'NUMDROPA', 1, NDROPMO)
      NBAS2=NBAS*NBAS
      IF (NDROPMO.GT.0) THEN
         NBAS_DROP = NOCCO(1)+NVRTO(1)
         CALL FILL_DROP_MOS(NIRREP, IUHF, NOCC_EXPND, NVRT_EXPND,
     &                      NBAS_FULL, NOCCO_EXPND, NVRTO_EXPND)
         NBAS  = NBAS_FULL
         NBAS2 = NBAS_FULL*NBAS_FULL
      END IF
C
C NOCC_EXPND, POP = The # of occ per irrep in the full and reduced spaces.
C NVRT_EXPND, VRT = The # of vrt per irrep in the full and reduced spaces.
C NOCCO_EXPND, NOCCO = The # of occ in the full and reduced spaces.
C NVRTO_EXPND, NVRTO = The # of vrt in the full and reduced spaces.
C NBAS = The number of basis functions in the full space.
C
       IF (NDROPMO.GT.0) THEN
          NOCC=NOCCO_EXPND(ISPIN)
          NVRT=NVRTO_EXPND(ISPIN)
       END IF

       I000=1
       I010=I000+NBAS2*IINTFP
       I020=I010+NBAS2*IINTFP
       I030=I020+NOCC*NVRT*IINTFP
       I040=I030+NOCC*NVRT*IINTFP
       IOFF=IINTFP*(NBAS*NOCC)

C (Always read the MO vectors in the full space.)
C
       IF (NDROPMO.GT.0) THEN
          CALL GETREC(20,'JOBARC','SCFEVC'//ISP(ISPIN)//'0',
     &                NBAS2*IINTFP,ICORE(I000))
       ELSE
          CALL GETREC(20,'JOBARC','SCFEVEC'//ISP(ISPIN),
     &                NBAS2*IINTFP,ICORE(I000))
       END IF

#ifdef _DEBUG_LVL0
      If (Ispin .Eq. 1) then
         Write(6,"(a)") " The original alpha SCF vectors"
      Else 
         Write(6,"(a)") " The original beta SCF vectors"
      Endif 

CSSS      call checksum("pccd_rotg      :", ICORE(I000), NBAS2)
      call output(icore(i000),1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif
       CALL SCOPY (NBAS2,ICORE(I000),1,ICORE(I010),1)
CSSS       Call GETREC(-20,'JOBARC','ORBOPITR',Ione,Ncycle)

       MEMLEFT = MAXCOR-I040 
       CALL PCCD_NR_SEARCH(GRD,GRD_OO,GRD_VV,GRD_VO,GRD_OV,GRAD_STATA,
     &                     GRAD_STATB,B_LIKE,LENOO,LENVV,LENVO,NBAS,
     &                     NOCC,NVRT,ICORE(I040),MEMLEFT,Iuhf,ISPIN,
     &                     TOL,NCYCLE)

C This B_LIKE option do a Brueckner like fake orbital optimization 

       IF (B_LIKE) RETURN

#ifdef _DEBUG_LVL0
       write(6,"(a)") "The U (in pccd_rotg)"
       call output(Grd,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif 
       IF (DIIS) THEN

          Do I = 1, Nbas
             Grd(I,I) = Grd(I,I) + One
          Enddo

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "U=(1+G)"
      Call output(Grd,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif
          CALL PCCD_ROTATE(GRD,ICORE,MAXCOR,ISPIN,NBAS,NCYCLE)
C This call has been moved to pccd_rotgrd.F!
CSSS       CALL PCCD_ROTATE_AUX(ICORE,MAXCOR,NBAS,IUHF,NCYCLE)
       ELSE
          CALL PCCD_NODIIS(GRD,ICORE,MAXCOR,ISPIN,NBAS,NCYCLE)
       ENDIF 

      RETURN
      END

