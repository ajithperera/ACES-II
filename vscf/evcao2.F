
C This subroutine reads the MOs in the SO basis and writes them to
C the output file in the ZMAT-ordered AO basis. It also writes
C the AOBASMOS file. This file can be read to form initial guesses
C and determine occupations, particularly in finite difference
C vibrational frequency calculations.

c#define _DEBUG_EVCAO2

      SUBROUTINE EVCAO2(CSO,    CZAO, SOZAO,  EVAL,
     &                  ATMCHG, LBAS, ICENT,  NS,
     &                  NP,     ND,   NF,     NG,
     &                  ANGTYP, NOCC, IOCC,   COORD,
     &                  NBASX,  NBAS, NATOMS, IUHF)

C     CSO   --- used for SO basis MOs, NBAS rows, NBAS columns.
C     CZAO  --- used for ZMAT ordered AO basis, NBASX rows, NBAS columns. 
C     SOZAO --- matrix to transform from SO to ZMAT AO basis
C               (CZAO=SOZAO*CSO).
C     EVAL  --- eigenvalues.
C     ATMCHG --- atomic numbers of atoms in ZMAT order (value is 0 for
C                a dummy, 110 for a ghost atom) (INTEGER).
C     LBAS   --- angular momentum for the NBASX basis functions.
C     ICENT  --- center to which a basis function belongs (ZMAT order).
C     NS     --- number of s basis functions on a center.
C     [NP,ND,NF,NG analogous to NS for p,d,f,g]
C     ATMLBL --- atomic symbol of center to which a basis function belongs.
C     BASLAB --- angular symbol (eg S,X,XX,XXY) of a basis function.
C     ANGTYP --- integer code according to value of BASLAB for each
C                basis function.
C     NOCC   --- Occupation in terms of computational point group irreps.
C     IOCC   --- Occupation numbers of the SOs for each spin.
C     COORD  --- Cartesian coordinates (in Bohrs).
C     NBASX  --- number of AOs in Cartesian basis (supplied by caller).
C     NBAS   --- number of basis functions or MOs (supplied by caller).
C     NATOMS --- number of centers (includes dummy and ghost atoms)
C                (supplied by caller).
C     IUHF   --- 0 if RHF, 1 otherwise.
C
C     Known limitations : no logic for higher than g functions (numbers
C                         printed may be OK, but labels will not).
CEND
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DOUBLE PRECISION CSO,CZAO,SOZAO,EVAL,COORD
      INTEGER ATMCHG,LBAS,ICENT,NS,NP,ND,NG,ANGTYP,NOCC,IOCC,NBASX,NBAS,
     &        NATOMS,IUHF
#include "mxcbf.par"
#include "flags.h"
      logical geom_opt
      CHARACTER*4 BASLAB(mxcbf),ATMLBL(mxcbf)
C-----------------------------------------------------------------------
      CHARACTER*4 ATSYMB(110)
      CHARACTER*5 SPIN(2)
      integer ijunk(100)
      LOGICAL YESNO
C-----------------------------------------------------------------------
      PARAMETER(LUOUT=6)
C-----------------------------------------------------------------------
      DIMENSION CSO(NBAS*NBAS),CZAO(NBASX*NBAS),SOZAO(NBASX*NBAS),
     &          EVAL(NBAS*2),
     &          ATMCHG(NATOMS),LBAS(NBASX),ICENT(NBASX),
     &          NS(NATOMS),NP(NATOMS),ND(NATOMS),NF(NATOMS),NG(NATOMS),
     &          ANGTYP(NBASX),NOCC(8,2),
     &          IOCC(NBAS,2),COORD(3,NATOMS)
C-----------------------------------------------------------------------
#include "machsp.com"
#include "symm2.com"
#include "flags.com"
#include "flags2.com"
C-----------------------------------------------------------------------
      DATA SPIN /'ALPHA',' BETA'/
      DATA (ATSYMB(I),I= 1, 54)
     &    /'H   ','He  ',
     &     'Li  ','Be  ','B   ','C   ','N   ','O   ','F   ','Ne  ',
     &     'Na  ','Mg  ','Al  ','Si  ','P   ','S   ','Cl  ','Ar  ',
     &     'K   ','Ca  ',
     &     'Sc  ','Ti  ','V   ','Cr  ','Mn  ','Fe  ','Co  ','Ni  ',
     &                                               'Cu  ','Zn  ',
     &                   'Ga  ','Ge  ','As  ','Se  ','Br  ','Kr  ',
     &     'Rb  ','Sr  ',
     &     'Y   ','Zr  ','Nb  ','Mo  ','Tc  ','Ru  ','Rh  ','Pd  ',
     &                                               'Ag  ','Cd  ',
     &                   'In  ','Sn  ','Sb  ','Te  ','I   ','Xe  '/
      DATA (ATSYMB(I),I=55,110)
     &    /'Cs  ','Ba  ',
     &     'La  ','Ce  ','Pr  ','Nd  ','Pm  ','Sm  ','Eu  ','Gd  ',
     &     'Tb  ','Dy  ','Ho  ','Er  ','Tm  ','Yb  ','Lu  ',
     &            'Hf  ','Ta  ','W   ','Re  ','Os  ','Ir  ','Pt  ',
     &                                               'Au  ','Hg  ',
     &                   'Tl  ','Pb  ','Bi  ','Po  ','At  ','Rn  ',
     &     'Fr  ','Ra  ',
     &     'Ac  ','Th  ','Pa  ','U   ','Np  ','Pu  ','Am  ','Cm  ',
     &     'Bk  ','Cf  ','Es  ','Fm  ','Md  ','No  ','Lr  ',
     &     '    ','    ','    ','    ','    ','    ','GH  '       /
C-----------------------------------------------------------------------

#ifdef _TRACE_VSCF
      write(*,*) '@VSCF: Entering ', __FILE__
#endif /* _TRACE_VSCF */

c     Calculate vector of occupation numbers.
      CALL GETREC(20,'JOBARC','OCCUPYA0',NIRREP,NOCC(1,1))
      IF (IUHF.EQ.0) THEN
         CALL GETREC(20,'JOBARC','OCCUPYA0',NIRREP,NOCC(1,2))
      ELSE
         CALL GETREC(20,'JOBARC','OCCUPYB0',NIRREP,NOCC(1,2))
      END IF

      call izero(iocc(1,1),nbas)
      call izero(iocc(1,2),nbas)
      do ispin=1,iuhf+1
         do irrep=1,nirrep
            if (nocc(irrep,ispin).ge.1) then
               do j=1,nocc(irrep,ispin)
                  iocc(irpoff(irrep) + j,ispin) = 2-iuhf
               end do
            end if
         end do
      end do

C     Get atomic numbers of centers (0 for dummy atoms, 110 for ghosts),
C     angular momenta of basis functions, and centers to which functions
C     belong.
      CALL GETREC(20,'JOBARC','ATOMCHRG',NATOMS,ATMCHG)
      CALL GETREC(20,'JOBARC','ANMOMBF0',NBASX,LBAS)
      CALL GETREC(20,'JOBARC','CNTERBF0',NBASX,ICENT)
#ifdef _DEBUG_EVCAO2
      WRITE(*,*) '@EVCAO2: ATMCHG ',(ATMCHG(I),I=1,NATOMS)
      WRITE(*,*) '@EVCAO2: LBAS   ',(LBAS(I),  I=1,NBASX)
      WRITE(*,*) '@EVCAO2: ICENT  ',(ICENT(I), I=1,NBASX)
#endif /* _DEBUG_EVCAO2 */

C     Compute number of s,p,d,f,g for each center.
      CALL IZERO(NS,NATOMS)
      CALL IZERO(NP,NATOMS)
      CALL IZERO(ND,NATOMS)
      CALL IZERO(NF,NATOMS)
      CALL IZERO(NG,NATOMS)
      DO IATOM=1,NATOMS
         NUMS = 0
         NUMP = 0
         NUMD = 0
         NUMF = 0
         NUMG = 0
         DO IBAS=1,NBASX
            IF (ICENT(IBAS).EQ.IATOM) THEN
               IF(LBAS(IBAS).EQ.0) NUMS = NUMS + 1
               IF(LBAS(IBAS).EQ.1) NUMP = NUMP + 1
               IF(LBAS(IBAS).EQ.2) NUMD = NUMD + 1
               IF(LBAS(IBAS).EQ.3) NUMF = NUMF + 1
               IF(LBAS(IBAS).EQ.4) NUMG = NUMG + 1
            END IF
         END DO
         NS(IATOM) = NUMS
         NP(IATOM) = NUMP / 3
         ND(IATOM) = NUMD / 6
         NF(IATOM) = NUMF /10
         NG(IATOM) = NUMG /15
#ifdef _DEBUG_EVCAO2
         WRITE(*,*) '@EVCAO2: NUM values for atom ',IATOM
         WRITE(*,*) NUMS,NUMP,NUMD,NUMF,NUMG
         WRITE(*,*) '@EVCAO2: N values for atom ',IATOM
         WRITE(*,*) NS(IATOM),NP(IATOM),ND(IATOM),NF(IATOM),NG(IATOM)
#endif /* _DEBUG_EVCAO2 */
         IF (NP(IATOM)* 3 .NE. NUMP .OR. ND(IATOM)* 6 .NE. NUMD .OR.
     &       NF(IATOM)*10 .NE. NUMF .OR. NG(IATOM)*15 .NE. NUMG) THEN
            WRITE(*,*) '@EVCAO2: Unexpected numbers of functions. '
#ifdef _DEBUG_EVCAO2
            WRITE(*,*) '   The following pairs should all be equal:'
            WRITE(*,*) 'NUMP: ', NP(IATOM)* 3, NUMP
            WRITE(*,*) 'NUMD: ', ND(IATOM)* 6, NUMD
            WRITE(*,*) 'NUMF: ', NF(IATOM)*10, NUMF
            WRITE(*,*) 'NUMG: ', NG(IATOM)*15, NUMG
#endif /* _DEBUG_EVCAO2 */
            CALL ERREX
         END IF
      END DO

C     Determine angular label for each basis function.       
      IBAS=1
      DO 60 IATOM=1,NATOMS
       IF(NS(IATOM).GE.1)THEN
        DO 10 IS=1,NS(IATOM)
        BASLAB(IBAS) = 'S   '
        ANGTYP(IBAS) = 1
        IBAS = IBAS + 1
   10   CONTINUE
       ENDIF
C
       IF(NP(IATOM).GE.1)THEN
        DO 20 IP=1,NP(IATOM)
        BASLAB(IBAS  ) = 'X   '
        BASLAB(IBAS+1) = 'Y   '
        BASLAB(IBAS+2) = 'Z   '
        ANGTYP(IBAS  ) = 2
        ANGTYP(IBAS+1) = 3
        ANGTYP(IBAS+2) = 4
        IBAS = IBAS + 3
   20   CONTINUE
       ENDIF
C
       IF(ND(IATOM).GE.1)THEN
        DO 30 IS=1,ND(IATOM)
        BASLAB(IBAS  ) = 'XX  '
        BASLAB(IBAS+1) = 'XY  '
        BASLAB(IBAS+2) = 'XZ  '
        BASLAB(IBAS+3) = 'YY  '
        BASLAB(IBAS+4) = 'YZ  '
        BASLAB(IBAS+5) = 'ZZ  '
        ANGTYP(IBAS  ) = 5
        ANGTYP(IBAS+1) = 6
        ANGTYP(IBAS+2) = 7
        ANGTYP(IBAS+3) = 8
        ANGTYP(IBAS+4) = 9
        ANGTYP(IBAS+5) = 10
        IBAS = IBAS + 6
   30   CONTINUE
       ENDIF
C
       IF(NF(IATOM).GE.1)THEN
        DO 40 IS=1,NF(IATOM)
        BASLAB(IBAS   ) = 'XXX '
        BASLAB(IBAS+ 1) = 'XXY '
        BASLAB(IBAS+ 2) = 'XXZ '
        BASLAB(IBAS+ 3) = 'XYY '
        BASLAB(IBAS+ 4) = 'XYZ '
        BASLAB(IBAS+ 5) = 'XZZ '
        BASLAB(IBAS+ 6) = 'YYY '
        BASLAB(IBAS+ 7) = 'YYZ '
        BASLAB(IBAS+ 8) = 'YZZ '
        BASLAB(IBAS+ 9) = 'ZZZ '
        ANGTYP(IBAS  ) = 11
        ANGTYP(IBAS+1) = 12
        ANGTYP(IBAS+2) = 13
        ANGTYP(IBAS+3) = 14
        ANGTYP(IBAS+4) = 15
        ANGTYP(IBAS+5) = 16
        ANGTYP(IBAS+6) = 17
        ANGTYP(IBAS+7) = 18
        ANGTYP(IBAS+8) = 19
        ANGTYP(IBAS+9) = 20
        IBAS = IBAS + 10
   40   CONTINUE
       ENDIF
C
       IF(NG(IATOM).GE.1)THEN
        DO 50 IG=1,NG(IATOM)
        BASLAB(IBAS   ) = 'XXXX'
        BASLAB(IBAS+ 1) = 'XXXY'
        BASLAB(IBAS+ 2) = 'XXXZ'
        BASLAB(IBAS+ 3) = 'XXYY'
        BASLAB(IBAS+ 4) = 'XXYZ'
        BASLAB(IBAS+ 5) = 'XXZZ'
        BASLAB(IBAS+ 6) = 'XYYY'
        BASLAB(IBAS+ 7) = 'XYYZ'
        BASLAB(IBAS+ 8) = 'XYZZ'
        BASLAB(IBAS+ 9) = 'XZZZ'
        BASLAB(IBAS+10) = 'YYYY'
        BASLAB(IBAS+11) = 'YYYZ'
        BASLAB(IBAS+12) = 'YYZZ'
        BASLAB(IBAS+13) = 'YZZZ'
        BASLAB(IBAS+14) = 'ZZZZ'
        ANGTYP(IBAS   ) = 21
        ANGTYP(IBAS+ 1) = 22
        ANGTYP(IBAS+ 2) = 23
        ANGTYP(IBAS+ 3) = 24
        ANGTYP(IBAS+ 4) = 25
        ANGTYP(IBAS+ 5) = 26
        ANGTYP(IBAS+ 6) = 27
        ANGTYP(IBAS+ 7) = 28
        ANGTYP(IBAS+ 8) = 29
        ANGTYP(IBAS+ 9) = 30
        ANGTYP(IBAS+10) = 31
        ANGTYP(IBAS+11) = 32
        ANGTYP(IBAS+12) = 33
        ANGTYP(IBAS+13) = 34
        ANGTYP(IBAS+14) = 35
        IBAS = IBAS + 15
   50   CONTINUE
       ENDIF
C
   60 CONTINUE      
      IF(IFLAGS(1) .GE. 10)THEN
       WRITE(6,*) ' @EVCAO2-I, IBAS ',IBAS
      ENDIF
cmn
c Write the record ANGTYP to JOBARC. This is required for MRCC Heff calculations
cmn
      call putrec(20, 'JOBARC', 'ANGTYP', nbasx, ANGTYP)

C
C     Set labels for atoms to which each basis function belongs (just
C     the chemical symbol).
C
      DO 70 I=1,NBASX
      IATOM = ICENT(I)
      ATMLBL(I) = ATSYMB(ATMCHG(IATOM))
      IF(IFLAGS(1).GE.10)THEN
cYAU       write(6,*) i,iatom,atmchg(iatom)
      ENDIF
   70 CONTINUE
C
C     Get eigenvalues.
C
      ITMP = NBAS*IINTFP
      CALL GETREC(20,'JOBARC','SCFEVLA0',ITMP,EVAL)
      IF (IUHF.EQ.0) THEN
         CALL GETREC(20,'JOBARC','SCFEVLA0',ITMP,EVAL(1+NBAS))
      ELSE
         CALL GETREC(20,'JOBARC','SCFEVLB0',ITMP,EVAL(1+NBAS))
      END IF
C
C     Open file AOBASMOS for writing MOs in AO basis.
C
      INQUIRE(FILE='AOBASMOS',EXIST=YESNO)
CSSS
CSSS       OPEN(71,FILE='AOBASMOS',STATUS='OLD',ACCESS='SEQUENTIAL',
CSSS     &      FORM='FORMATTED')
CSSS       CLOSE(UNIT=71,STATUS='DELETE')
CSSS      ENDIF

      IF (.NOT. YESNO) THEN
         OPEN(71,FILE='AOBASMOS',STATUS='NEW',ACCESS='SEQUENTIAL',
     &    FORM='FORMATTED')
      ELSE
         OPEN(71,FILE='AOBASMOS',STATUS='OLD',ACCESS='SEQUENTIAL',
     &    FORM='FORMATTED')
      ENDIF
C
C     Transform the MOs and print eigenvalues and eigenvectors.
C
      ICNT0=0
      DO 100 ISPIN=1,IUHF+1
C
C     Get MOs in SO basis.
C
        IF(ISPIN.EQ.1)THEN
         CALL GETREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,CSO)
#ifdef _DEBUG_LVLM
        write(6,*) 
        write(6, "(a)") "The Alpha SCF natural orbitals"
        call output(cso, 1, nbas, 1, nbas, nbas, nbas, 1)
#endif 
        ELSE
         CALL GETREC(20,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,CSO)
#ifdef _DEBUG_LVLM
        write(6,*) 
        write(6, "(a)") "The beta SCF natural orbitals"
        call output(cso, 1, nbas, 1, nbas, nbas, nbas, 1)
#endif 
        ENDIF
C
C     Get CMP2ZMAT matrix.
C
        CALL GETREC(20,'JOBARC','CMP2ZMAT',NBASX*NBAS*IINTFP,SOZAO)

#ifdef _DEBUG_LVLM
      Write(6,*) "CMP2ZMAT transformation"
      Call output(SOZAO, 1, Nbasx, 1, Nbas, Nbasx, Nbas, 1)
#endif
C
C     Apply CMP2ZMAT to SO eigenvectors.
C
        CALL XGEMM('N','N',NBASX,NBAS,NBAS,1.0D+00,
     &             SOZAO,NBASX,CSO,NBAS,0.0D+00,CZAO,NBASX)

#ifdef _DEBUG_LVLM
      Write(6,*) "The transformed eiegenvectors"
      Call output(CZAO, 1, Nbasx, 1, Nbas, Nbasx, Nbas, 1)
#endif
C
C     Print AO basis eigenvectors to output file.
C
      IF(IFLAGS(1).GE.1)THEN
C
        ICNT1=0
        WRITE(LUOUT,5001)
        ISTART = 0
        DO 105 I=1,NIRREP
          WRITE(LUOUT,5010)I,SPIN(ISPIN)
          NBY4=NBFIRR(I)/4
          DO 110 J=1,NBY4
            ICNT0=4*J+ICNT1
            WRITE(LUOUT,5019)(ICNT0-4+K,K=1,4)
            WRITE(LUOUT,5020)(EVAL((ISPIN-1)*NBAS+IRPOFF(I)+
     &                             (J-1)*4+K),K=1,4)
            WRITE(LUOUT,5025)
c           DO 115 K=1,NBFIRR(I)
            DO 115 K=1,NBASX
              WRITE(LUOUT,5030)K,ATMLBL(K),BASLAB(K),
     &         (CZAO(ISTART+(J-1)*4*NBASX + (L-1)*NBASX + K),L=1,4)
  115       CONTINUE
            WRITE(LUOUT,5035)
  110     CONTINUE
C
C  Now take care of the remaining eigenvectors, if we had a total
C  in this irrep which was not divisible by 4.
C
          IF(NBY4.EQ.0) ICNT0=ICNT1
          IF(4*NBY4.NE.NBFIRR(I)) THEN
            ICNT=NBFIRR(I)-4*NBY4
            IF(ICNT.EQ.1) THEN
              WRITE(LUOUT,5050)ICNT0+1
              WRITE(LUOUT,5021)(EVAL((ISPIN-1)*NBAS+IRPOFF(I)+
     &                               NBY4*4+K),K=1,1)
              WRITE(LUOUT,5026)
              DO 116 K=1,NBASX
                WRITE(LUOUT,5031)K,ATMLBL(K),BASLAB(K),
     &           (CZAO(ISTART+NBY4*4*NBASX + (L-1)*NBASX + K),L=1,1)
  116         CONTINUE
              ICNT1=ICNT0+1
            ELSEIF(ICNT.EQ.2) THEN
              WRITE(LUOUT,5051)(ICNT0+K,K=1,2)
              WRITE(LUOUT,5022)(EVAL((ISPIN-1)*NBAS+IRPOFF(I)+
     &                               NBY4*4+K),K=1,2)
              WRITE(LUOUT,5027)
              DO 117 K=1,NBASX
                WRITE(LUOUT,5032)K,ATMLBL(K),BASLAB(K),
     &           (CZAO(ISTART+NBY4*4*NBASX + (L-1)*NBASX + K),L=1,2)
  117         CONTINUE
              ICNT1=ICNT0+2
            ELSE
              WRITE(LUOUT,5052)(ICNT0+K,K=1,3)
              WRITE(LUOUT,5023)(EVAL((ISPIN-1)*NBAS+IRPOFF(I)+
     &                               NBY4*4+K),K=1,3)
              WRITE(LUOUT,5028)
              DO 118 K=1,NBASX
                WRITE(LUOUT,5033)K,ATMLBL(K),BASLAB(K),
     &           (CZAO(ISTART+NBY4*4*NBASX + (L-1)*NBASX + K),L=1,3)
  118         CONTINUE
              ICNT1=ICNT0+3
            ENDIF
            WRITE(LUOUT,5035)
          ELSE
            ICNT1=ICNT0
          ENDIF
        ISTART = ISTART + NBASX*NBFIRR(I)
  105   CONTINUE
C
      ENDIF
C
c
cmn if we are performing an optimization update AOBASMOS until converged
c
      geom_opt=(iflags2(h_IFLAGS2_geom_opt).ne.0)
      if (geom_opt) then
      istat = 0
      call getrec(20,'JOBARC','JODADONE',1,istat)
      else
         istat = 1
      endif
c
#ifdef _DEBUG_LVLM
      Print*, "Vib calc,  AOBASMOS, GEO?", iFlags(h_IFLAGS_vib),
     $     YESNO, istat
#endif

      If ((iFlags(h_IFLAGS_vib).GT.1) .AND. YESNO
     $     .and. istat .eq. 1) Then
         close(71, status='KEEP')
         RETURN
      endif
C
C     Write eigenvectors to AOBASMOS
C
C     At the beginning write NBAS, NBASX, and IUHF so when we read
C     OLDAOMOS we can check the dimensions are as expected.
C
      IF(ISPIN.EQ.1)THEN
cJDW 11/30/2004
c      WRITE(71,'(3I12)') NBAS,NBASX,IUHF
       WRITE(71,'(11i7)') NBAS,NBASX,IUHF,(NBFIRR(IRREP),IRREP=1,NIRREP)
      ENDIF
C
      NCOLS = 4
      NFULL = NBAS/NCOLS
      NLEFT = NBAS - NCOLS*NFULL
      DO 150 IBLOCK=1,NFULL
      DO 140 IROW  =1,NBASX
      WRITE(71,1020)
     & (CZAO((IBLOCK-1)*NBASX*NCOLS + (J-1)*NBASX + IROW),J=1,NCOLS)
  140 CONTINUE
  150 CONTINUE
C
      IF(NLEFT.NE.0)THEN
       DO 160 IROW  =1,NBASX
       WRITE(71,1020)
     &  (CZAO(     NFULL*NBASX*NCOLS + (J-1)*NBASX + IROW),J=1,NLEFT)
  160  CONTINUE
      ENDIF
C
C     Write occupation numbers and eigenvalues.
C
      DO 170 IBAS=1,NBAS
      WRITE(71,1030) IOCC(IBAS,ISPIN),EVAL((ISPIN-1)*NBAS+IBAS)
  170 CONTINUE
C
C     Write angular type array (length is NBASX).
C
      DO 180 IBAS=1,NBASX
      WRITE(71,*) ANGTYP(IBAS)
  180 CONTINUE
C
  100 CONTINUE
C
C     Write coordinates.
C
      CALL GETREC(20,'JOBARC','COORD',3*NATOMS*IINTFP,COORD)
      DO 250 IATOM=1,NATOMS
      WRITE(71,*) (COORD(I,IATOM),I=1,3)
  250 CONTINUE
C
C     Write computational to ZMAT mapping.
C
      CALL GETREC(20,'JOBARC','MAP2ZMAT',NATOMS,IJUNK)
      DO 260 IATOM=1,NATOMS
      WRITE(71,*) IJUNK(IATOM)
  260 CONTINUE
C
      CLOSE(UNIT=71,STATUS='KEEP')

#ifdef _TRACE_VSCF
      write(*,*) '@VSCF: Leaving ', __FILE__
#endif /* _TRACE_VSCF */

      RETURN
 5001 FORMAT(/,T3,'ORBITAL EIGENVECTORS (ZMAT ordered AO basis)',/)
 5010 FORMAT(T3,'SYMMETRY BLOCK ',I1,' (',A5,')',/)
 5019 FORMAT(T25,2X,'MO #',I3,3(7X,'MO #',I3),/)
 5020 FORMAT(T5,'BASIS/ORB E',9X,F10.5,4X,F10.5,4X,F10.5,
     &       4X,F10.5)
 5025 FORMAT(T25,'----------',4X,'----------',4X,
     &           '----------',4X,'----------')
 5030 FORMAT(T3,I3,';',1X,A4,2X,A4,7X,F10.5,4X,F10.5,4X,
     &        F10.5,4X,F10.5)
 5035 FORMAT(T3,/)
 5021 FORMAT(T5,'BASIS\\ORB E',9X,F10.5,1(4X,F10.5))
 5022 FORMAT(T5,'BASIS\\ORB E',9X,F10.5,2(4X,F10.5))
 5023 FORMAT(T5,'BASIS\\ORB E',9X,F10.5,3(4X,F10.5))
 5026 FORMAT(T25,'----------')
 5027 FORMAT(T25,'----------',4X,'----------')
 5028 FORMAT(T25,'----------',4X,'----------',4X,'----------')
 5031 FORMAT(T3,I3,';',1X,A4,2X,A4,7X,F10.5,1(4X,F10.5))
 5032 FORMAT(T3,I3,';',1X,A4,2X,A4,7X,F10.5,2(4X,F10.5))
 5033 FORMAT(T3,I3,';',1X,A4,2X,A4,7X,F10.5,3(4X,F10.5))
 5050 FORMAT(T25,2X,'MO #',I3,/)
 5051 FORMAT(T25,2X,'MO #',I3,7X,'MO #',I3,/)
 5052 FORMAT(T25,2X,'MO #',I3,7X,'MO #',I3,7X,'MO #',I3,/)
 1020 FORMAT(4F20.10)
 1030 FORMAT(I3,F20.10)
      END
