      SUBROUTINE GETDENS3(DENS,CZAO,CSO,ZAOSO,SCR,SCR2,
     &                    NBASX,NBAS,NATOMS,
     &                    NIRREP,LDIM1,ITRIOF,IREPS,NBFIRR,IUHF,
     &                    IONEORTWO)
      IMPLICIT NONE
C-----------------------------------------------------------------------
      DOUBLE PRECISION DENS,CZAO,CSO,ZAOSO,SCR,SCR2
      INTEGER NBASX,NBAS,NATOMS,NIRREP,LDIM1,ITRIOF,IREPS,NBFIRR,IUHF,
     &        IONEORTWO
C-----------------------------------------------------------------------
      CHARACTER*6 BASTYP
      CHARACTER*80 BASNAM(100)
      INTEGER IATOM,NCOLS,NFULL,NLEFT,IBAS,J,IFULL
      INTEGER IRREP
      INTEGER IMAP,IOFF,IOFFMIN,IZATOM,NSMIN,NPMIN,NDMIN,NSFULL,NPFULL,
     &        NDFULL,NMIN,NTOT,NSTO3G,ISPIN,JFUN,IFUN,JATOM
      INTEGER IOCC,ICOL,P,Q,IJUNK,NBFIRRMIN,IOCCIND,IOFFSYM
      INTEGER IJUNK1,IJUNK2
      integer irowoff,icoloff
      INTEGER IPRINT
      DOUBLE PRECISION EJUNK
      DOUBLE PRECISION DENMIN,FACTOR
C-----------------------------------------------------------------------
      INTEGER INDX2,I,N
C-----------------------------------------------------------------------
      INTEGER IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      INTEGER NOCC
      INTEGER IFLAGS
C-----------------------------------------------------------------------
      DIMENSION CZAO(NBASX,NBAS*(IUHF+1)),CSO(NBAS,NBAS),
     &          ZAOSO(NBAS*NBASX),DENS(1),
     &          SCR(NBAS*NBAS),SCR2(NBAS*NBAS),
     &          ITRIOF(8),IREPS(9),NBFIRR(8)
      dimension denmin(10000)
c     DIMENSION IMAP(24,100),FACTOR(24,100),IOFF(100),IOFFMIN(100),
      DIMENSION IMAP(29,100),FACTOR(29,100),IOFF(100),IOFFMIN(100),
     &          IZATOM(100),NSMIN(100),NPMIN(100),NDMIN(100),
     &          NSFULL(100),NPFULL(100),NDFULL(100),NMIN(100),NTOT(100)
      DIMENSION BASTYP(100)
      DIMENSION NBFIRRMIN(8),IOCC(100),IOFFSYM(8)
C-----------------------------------------------------------------------
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /POPUL/  NOCC(16)
      COMMON /FLAGS/  IFLAGS(100)
C-----------------------------------------------------------------------
      INDX2(I,J,N)=I+(J-1)*N
C
C     This routine reads a density matrix in one basis set and from this
C     density creates an initial density matrix in the current basis.
C     The density matrix is expected to be in the ZMAT-ordered AO basis
C     set and alpha and beta density matrices are expected if this is an
C     open-shell calculation.
C
C     NSTO3G  --- The number of functions in the original basis set.
C     NATOMS  --- The number of atoms.
C     IZATOM  --- Array containing the atomic numbers of each atom.
C     NSMIN   --- Number of s functions in original basis set for each
C                 atom.
C     NPMIN   --- Number of p functions in original basis set for each
C                 atom.
C     NDMIN   --- Number of d functions in original basis set for each
C                 atom.
C     NSFULL  --- Number of s functions in current basis set for each
C                 atom.
C     NPFULL  --- Number of p functions in current basis set for each
C                 atom.
C     NDFULL  --- Number of d functions in current basis set for each
C                 atom.
C
      IPRINT = IFLAGS(1)
C-----------------------------------------------------------------------
C     Get atomic numbers.
C-----------------------------------------------------------------------
      CALL GETREC(20,'JOBARC','NATOMS'  ,1     ,NATOMS)
      CALL GETREC(20,'JOBARC','ATOMCHRG',NATOMS,IZATOM)
      IF(IPRINT.GE.10) write(6,*) ' @GETDENS3-I, IZATOM ',izatom
C
C-----------------------------------------------------------------------
C     Get basis name of each atom.
C-----------------------------------------------------------------------
      DO 1 IATOM=1,NATOMS
      CALL GETBASNAM(IATOM,IZATOM(IATOM),BASNAM(IATOM))
      WRITE(6,*) BASNAM(IATOM)
    1 CONTINUE
C
C-----------------------------------------------------------------------
C     Get numbers of contracted functions of current basis and splitting
C     type.
C-----------------------------------------------------------------------
      DO 2 IATOM=1,NATOMS
      CALL GETBASINF(IZATOM(IATOM),BASNAM(IATOM),BASTYP(IATOM),
     &               NSFULL(IATOM),NPFULL(IATOM),NDFULL(IATOM))
    2 CONTINUE
cccc      call errex
C
      IF(IONEORTWO.NE.1 .AND. IONEORTWO.NE.2)THEN
       write(6,*) ' @GETDENS3-F, Invalid value for IONEORTWO ',IONEORTWO
       call errex
      ENDIF
C
      IF(IONEORTWO.EQ.1)THEN
       OPEN(UNIT=71,FILE='EHTMOFILE',STATUS='OLD',FORM='FORMATTED',
     &      ACCESS='SEQUENTIAL')
      ENDIF
C
      IF(IONEORTWO.EQ.2)THEN
       OPEN(UNIT=71,FILE='MOFILE',STATUS='OLD',FORM='FORMATTED',
     &      ACCESS='SEQUENTIAL')
      ENDIF
C
C     Read size of STO-3G basis set and symmetry populations for the
C     STO-3G basis set.
C
      READ(71,*) NSTO3G,IJUNK1,IJUNK2,(NBFIRRMIN(I),I=1,NIRREP)
C
C-----------------------------------------------------------------------
C     Compute mapping and factor arrays for each atom.
C-----------------------------------------------------------------------
C
      DO 30 IATOM=1,NATOMS
C
      IF(IPRINT.GE.10)THEN
       write(6,*) ' @GETDENS3-I, Calling mkbasmap '
       write(6,*) izatom(iatom),nsfull(iatom),npfull(iatom),
     &                                        ndfull(iatom)
       write(6,'(A)') bastyp(iatom)
      ENDIF

      CALL MKBASMAP(IMAP(1,IATOM),FACTOR(1,IATOM),IZATOM(IATOM),
     &               NSMIN(IATOM), NPMIN(IATOM), NDMIN(IATOM),
     &              NSFULL(IATOM),NPFULL(IATOM),NDFULL(IATOM),
     &              BASTYP(IATOM))
   30 CONTINUE
C
C-----------------------------------------------------------------------
C     Calculate number of functions for each atom in both basis sets.
C-----------------------------------------------------------------------
      DO 26 IATOM=1,NATOMS
      NMIN(IATOM) = NSMIN(IATOM) +3*NPMIN(IATOM) +6*NDMIN(IATOM)
      NTOT(IATOM) = NSFULL(IATOM)+3*NPFULL(IATOM)+6*NDFULL(IATOM)
   26 CONTINUE
C
      IF(IPRINT.GE.10)THEN
       write(6,*) ' @GETDENS3-I, Mapping array '
       do 31 i=1,29
       write(6,*) (imap(i,j),j=1,natoms)
   31  continue
       write(6,*) ' @GETDENS3-I, Factor  array '
       do 32 i=1,29
       write(6,*) (factor(i,j),j=1,natoms)
   32  continue
      ENDIF
C
C-----------------------------------------------------------------------
C     Compute offsets in the two basis sets for each atom.
C-----------------------------------------------------------------------
C
      DO 35 IATOM=1,NATOMS
      IF(IATOM.EQ.1)THEN
       IOFF(IATOM) = 0
       IOFFMIN(IATOM) = 0
      ELSE
       IOFF(IATOM) = IOFF(IATOM-1) + NTOT(IATOM-1)
       IOFFMIN(IATOM) = IOFFMIN(IATOM-1) + NMIN(IATOM-1)
      ENDIF
   35 CONTINUE
C
C-----------------------------------------------------------------------
C     Read the STO-3G MOs.
C-----------------------------------------------------------------------
      NCOLS = 4
      NFULL = NSTO3G/NCOLS
      NLEFT = NSTO3G - NFULL*NCOLS

      IF(IPRINT.GE.10) write(6,*) nsto3g,ncols,nfull,nleft

      DO 300 ISPIN=1,IUHF+1
C
      CALL ZERO(DENMIN,NSTO3G*NSTO3G)
C
      DO 120 IFULL=1,NFULL
      DO 110 IBAS =1,NSTO3G
      READ(71,'(4F20.10)')
     &(DENMIN((IFULL-1)*NCOLS*NSTO3G + (J-1)*NSTO3G + IBAS),J=1,4)
  110 CONTINUE
  120 CONTINUE
C
      IF(NLEFT.GT.0)THEN
       DO 130 IBAS=1,NSTO3G
       IF(NLEFT.EQ.1) READ(71,'(1F20.10)')
     &  (DENMIN(NFULL*NCOLS*NSTO3G + (J-1)*NSTO3G + IBAS),J=1,1)
       IF(NLEFT.EQ.2) READ(71,'(2F20.10)')
     &  (DENMIN(NFULL*NCOLS*NSTO3G + (J-1)*NSTO3G + IBAS),J=1,2)
       IF(NLEFT.EQ.3) READ(71,'(3F20.10)')
     &  (DENMIN(NFULL*NCOLS*NSTO3G + (J-1)*NSTO3G + IBAS),J=1,3)
  130  CONTINUE
      ENDIF
C
      DO 140 IBAS=1,NSTO3G
      READ(71,*) IOCC(IBAS),EJUNK
  140 CONTINUE
C
      DO 150 IBAS=1,NSTO3G
      READ(71,*) IJUNK
  150 CONTINUE
C
C-----------------------------------------------------------------------
C     Compute offset vector so we know where each symmetry begins in
C     the STO-3G list.
C-----------------------------------------------------------------------
      IOFFSYM(1)=1
      IF(NIRREP.GT.1)THEN
       DO 160 IRREP=2,NIRREP
       IOFFSYM(IRREP) = IOFFSYM(IRREP-1) + NBFIRRMIN(IRREP-1)
  160  CONTINUE
      ENDIF
C
C-----------------------------------------------------------------------
C     Form new set of occupied orbitals.
C-----------------------------------------------------------------------
      CALL ZERO(CZAO,NBASX*NBAS)
      ICOL=0
      DO 190 IRREP  =1,NIRREP
      DO 180 IOCCIND=1,NOCC(IRREP+(ISPIN-1)*8)
      ICOL=ICOL+1
      DO 170 IATOM=1,NATOMS
      DO 165 IFUN =1,NTOT(IATOM)
C
      IF(IMAP(IFUN,IATOM) .NE. 999)THEN
       CZAO(IOFF(IATOM) + IFUN,ICOL) =
C
     & DENMIN( (IOFFSYM(IRREP)-1)*NSTO3G + (IOCCIND-1)*NSTO3G
     &                                   + IOFFMIN(IATOM)
     &                                   + IMAP(IFUN,IATOM) ) *
C
     & FACTOR(IFUN,IATOM)
      ENDIF

  165 CONTINUE
  170 CONTINUE
  180 CONTINUE
  190 CONTINUE
C
      IF(IPRINT.GE.1)THEN
       write(6,*) ' @GETDENS3-I, STO-3G orbitals '
       do 191 i=1,nsto3g
       write(6,'(7f9.5)') (denmin( (j-1)*nsto3g + i),j=1,7)
  191  continue
C
       write(6,*) ' @GETDENS3-I, Occupied orbitals in current basis '
       call output(czao,1,nbasx,1,icol,nbasx,icol,1)
       write(6,*) ' @GETDENS3-I, Spin and number of occ ',ispin,icol
      ENDIF
C
C-----------------------------------------------------------------------
C     Let us now try to orthonormalize these orbitals.
C     First transform to the SO basis. MOs in the SO basis are put in
C     CSO. 
C     *** CSO MUST be preserved ***
C-----------------------------------------------------------------------
      CALL GETREC(20,'JOBARC','ZMAT2CMP',NBAS*NBASX*IINTFP,ZAOSO)
      CALL XGEMM('N','N',NBAS,ICOL,NBASX,1.0D+00,
     &           ZAOSO,NBAS,CZAO,NBASX,0.0D+00,
     &           CSO,NBAS)
C
C     Occupied MOs in the SO basis are in CSO.
C
C     CZAO and ZAOSO may now be used for other purposes.
C
C     Orthonormalize one block at a time.
C
      IROWOFF=1
      ICOLOFF=1

      DO 192 IRREP=1,NIRREP
      IF(NOCC(IRREP+8*(ISPIN-1)) .EQ. 0) GOTO 192
C
C     Move the appropriate block of CSO to ZAOSO, bearing in mind the
C     symmetry blocking of CSO.
C
      CALL GETBLK2(CSO,ZAOSO,NBAS,ICOL,
     &             NBFIRR(IRREP),NOCC(IRREP+8*(ISPIN-1)),
     &             IROWOFF,ICOLOFF)
C
      CALL SORTHOI(ZAOSO,CZAO,SCR,SCR2,NBAS,.FALSE.,
     &             NOCC(1+8*(ISPIN-1)),IRREP)
C
C     SORTHOI returns orthornormal block of occupied orbitals in ZAOSO.
C     Form the contribution to the density matrix from ZAOSO.
C
      CALL MKDEN(ZAOSO,SCR,NOCC(IRREP+8*(ISPIN-1)),NBFIRR(IRREP),IUHF)
      IF(ISPIN.EQ.1)THEN
       CALL SQUEZ2(SCR,DENS(ITRIOF(IRREP)),NBFIRR(IRREP))
      ELSE
       CALL SQUEZ2(SCR,DENS(LDIM1+ITRIOF(IRREP)),NBFIRR(IRREP))
      ENDIF
C
      IROWOFF=IROWOFF+NBFIRR(IRREP)
      ICOLOFF=ICOLOFF+NOCC(IRREP+8*(ISPIN-1))
  192 CONTINUE
C
  300 CONTINUE
C
      CLOSE(UNIT=71,STATUS='KEEP')
      RETURN
      END
