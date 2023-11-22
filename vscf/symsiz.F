      SUBROUTINE SYMSIZ(NBAS,IRPSZA,IRPSZB,IPKSIZ,MAXMEM,NAOBK,
     &                  MXLSSZ1,MXLSSZ2,IRPPMEM,INOUT,NDIMRPA,
     &                  MAXCOR,IUHF)
C
C  This routine determines the offsets for the full in-core symmetry
C  algorithm, the maximum amount of memory for the full in-core algorithm
C  and partial out-of-core algorithm, and decides which of the three
C  different algorithms to use for the transformation.
C
C      IDXVEC - Lookup vector for determining irreps of orbitals
C
C        NBAS - Number of basis functions before dropped MOs
C               the (aa|aa) and (ab|ab) symmetry integrals in core.
C
C      IRPSZB - Amount of memory (in integer words) required to hold
C               the (aa|bb) symmetry integrals in core.
C
C      TOTSIZ - Total amount of memory of memory required to hold all
C               the AO integrals in core.
C
C      MAXMEM - Maximum amount of memory to be allocated for the
C               transformation (value will be for the algorithm which
C               is selected)
C
C       NAOBK - Number of symmetry combinations of the (aa|aa), (ab|ab),
C               and (aa|bb) symmetry integrals.  This is determined
C               from the number of irreps and is calculated in this
C               routine.  Needed for out-of-core PK list formation.
C
C     MXLSSZ1 - Length of the longest AO (aa|aa) or (aa|bb) symmetry
C               integral list.  Needed for the out-of-core PK list
C               formation.
C
C     MXLSSZ2 - Length of the longest AO (ab|ab) symmetry integral list.
C               Needed for the out-of-core PK list formation.
C
C       INOUT - Variable which tells main routine which algorithm is
C               selected (full in-core, partial out-of-core, or C1)
C
C
C   Definition of SYMM2 common block variables (most are determined in
C   this routine):
C
C       NIRREP - Number of irreps for this molecule.
C
C       NBFIRR - Number of basis functions in each irrep before
C                dropping MOs; dimensioned to max number of irreps
C                in an Abelian group (8).
C
C       IRPSZ1 - This vector contains the sizes of each of the (aa|aa)
C                and (ab|ab) symmetry combinations.  The element of
C                the vector is indexed by INDX(a,b), where "a" and "b"
C                are the irrep numbers.
C
C       IRPSZ2 - Same as IRPSZ1, except for (aa|bb) symmetry combinations.
C                The indexing is INDX(a-1,b).
C
C       IRPDS1 - Vector which contains the distribution size of the
C                (aa|aa) and (ab|ab) symmetry combinations.  Indexed
C                as INDX(a,b), where "a" and "b" are irrep numbers.
C
C       IRPDS2 - Vector which contains the distribution sizes (both
C                bra and ket) of the (aa|bb) symmetry combinations.
C                IRPDS2(odd) contain the bra distribution sizes,
C                IRPDS2(even) contain the ket distribution sizes.  These
C                entries are accessed by 2*INDX(a-1,b)-1 and 2*INDX(a-1,b),
C                or an equivalent counting.
C
C       IRPOFF - This vector contains the orbital number offset for
C                each of the irreps.  For example, if there are 13
C                orbitals in irrep 1, 4 in irrep 2, 7 in irrep 3,
C                and 2 in irrep 4, then IRPOFF(1)=0, IRPOFF(2)=13,
C                IRPOFF(3)=17, IRPOFF(4)=24.
C
C        IREPS - This vector contains the starting orbital number for
C                each irrep, before dropping MOs.  For the example given
C                for IRPOFF, IREPS(1)=1, IREPS(2)=14, IREPS(3)=18,
C                IREPS(4)=25.
C
C       DIRPRD - This is the direct product table; it is read in from
C                JOBARC.
C
C       IWOFF1 - Absolute offset into the integral vector W for the
C                starting of each of the (aa|aa) and (ab|ab) symmetry
C                combinations.  Indexed as INDX(a,b).
C
C       IWOFF2 - Absolute offset into the integral vector W for the
C                starting of each of the (aa|bb) symmetry combinations.
C                Indexed as INDX(a-1,b).
C
C       INEWVC - Vector which contains the relative orbital number within
C                an irrep for each of the orbitals (no dropped MOs).
C
C   Other variables of note:
C
C       IBKOFF - This vector contains cumulative number of the different
C                symmetry combinations.  IBKOFF(1)=0.  IBKOFF(2) is the
C                number of (aa|aa) and (ab|ab) symmetry combinations.
C                IBKOFF(3) is the number of (aa|bb) symmetry combinations
C                plus IBKOFF(2).
C
C       DIRPRD - The direct product table.
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL AOFIL, REACTION_FIELD
      LOGICAL FIRSTSCF, YESNO3
C
      PARAMETER(MINMEM=50 000)
C
      DIMENSION IRPSZA(37),IRPSZB(30)
      DIMENSION IBKOFF(4)
      DIMENSION IADERR(200),IADFCK(200,2)
      DIMENSION IPKOFF(73)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /FOCK/ AOFIL
#include "symm2.com"
      COMMON /OUTCOR/ IBKOFF
      COMMON /PKOFF/ IPKOFF
      COMMON /RPPINF/ ISTART,IORDER,IADERR,IADFCK
C
      DATA I255 /MAXBASFN/
      DATA ILNBUF /600/
      DATA IONE /1/
C
      INDX(I,J)=J+I*(I-1)/2
      IPACK(I,J,K,L)=OR(OR(OR(I,ISHFT(J,IBITWD)),ISHFT(K,2*IBITWD)),
     &                  ISHFT(L,3*IBITWD))
      IUPKI(INT)=AND(INT,IALONE)
      IUPKJ(INT)=AND(ISHFT(INT,-IBITWD),IALONE)
      IUPKK(INT)=AND(ISHFT(INT,-2*IBITWD),IALONE)
      IUPKL(INT)=AND(ISHFT(INT,-3*IBITWD),IALONE)
C
C
CJDW 3/28/97. Check whether this is the first entry to xvscf.
C             FIRSTSCF is used to control whether we print certain
C             things or not.
C
      CALL GETREC(1,'JOBARC','FIRSTRUN',1,I)
      FIRSTSCF = (i.ne.0).or.(iFlags(1).ge.10)
c
      CALL IZERO(INEWVC,I255)
      IRPDX=NIRREP*(NIRREP+1)/2
      IAASZ=0
      IABSZ=0
      ISIZE=0
      JSIZE=0
      IPKSIZ=0
      IMAX=0
      IRPOFF(1)=0
      IPKOFF(1)=0
      IREPS(1)=1
      MXBKSZ=0
C
C  Assign the order for RPP
C
      IORDER=IFLAGS(20)
C
C  Determine the symmetry packed size of each irrep and also place
C  the total in the last location.
C
      ITOT=0
      JTOT=0
      NBAST=0
      DO 325 I=1,NIRREP
        NBAST=NBAST+NBFIRR(I)
        ITRILN(I)=NBFIRR(I)*(NBFIRR(I)+1)/2
        ITOT=ITOT+ITRILN(I)
        ISQRLN(I)=NBFIRR(I)*NBFIRR(I)
        JTOT=JTOT+ISQRLN(I)
  325 CONTINUE
      ITRILN(NIRREP+1)=ITOT
      ISQRLN(NIRREP+1)=JTOT
C
C  Now store the offset values for each of the symmetry blocks for
C  two dimensional arrays using ITRILN.
C
      ITRIOF(1)=1
      ISQROF(1)=1
      DO 326 I=1,NIRREP-1
        ITRIOF(I+1)=ITRIOF(I)+ITRILN(I)
        ISQROF(I+1)=ISQROF(I)+ISQRLN(I)
  326 CONTINUE
C
C  Determine the largest dimensioned irrep and place the square of its
C  value into MXIRR2.  Needed for various dimensioning of scratch arrays.
C
      DO 327 I=1,NIRREP
        IMAX=MAX(IMAX,NBFIRR(I))
  327 CONTINUE
      MXIRR2=IMAX*IMAX
C
C  These loops account for the distribution, total sizes, and offsets of
C  the (aa|aa) and (ab|ab) type symmetry integrals.
C
      DO 10 I=1,NIRREP
        IREPS(I+1)=IREPS(I)+NBFIRR(I)
        IRPOFF(I+1)=IRPOFF(I)+NBFIRR(I)
        IRPDS1(INDX(I,I))=NBFIRR(I)*(NBFIRR(I)+1)/2
        IRPSZ1(INDX(I,I))=IRPDS1(INDX(I,I))*IRPDS1(INDX(I,I))
        ISIZE=ISIZE+IRPSZ1(INDX(I,I))
        IPKSIZ=IPKSIZ+IRPSZ1(INDX(I,I))
        MXLSSZ1=MAX(IRPSZ1(INDX(I,I)),MXLSSZ1)
        IAASZ=IAASZ+IRPSZ1(INDX(I,I))
        DO 11 J=1,I-1
          IRPDS1(INDX(I,J))=NBFIRR(I)*NBFIRR(J)
          IRPSZ1(INDX(I,J))=IRPDS1(INDX(I,J))*IRPDS1(INDX(I,J))
          ISIZE=ISIZE+IRPSZ1(INDX(I,J))
          MXLSSZ2=MAX(IRPSZ1(INDX(I,J)),MXLSSZ2)
          IABSZ=IABSZ+IRPSZ1(INDX(I,J))
   11   CONTINUE
   10 CONTINUE
      IRPSZA(1)=0
      DO 15 I=1,IRPDX
        IRPSZA(I+1)=IRPSZA(I)+IRPSZ1(I)*IINTFP
   15 CONTINUE
      IWOFF1(1)=0
      DO 17 I=2,IRPDX+1
        IWOFF1(I)=IWOFF1(I-1)+IRPSZ1(I-1)
   17 CONTINUE
C
C  This section accounts for the distribution and total size of
C  the (aa|bb) type symmetry integrals.  We need to account for
C  two different distribution sizes because the two "halves" aren't
C  necessarily equal in size (distribution "aa" and distribution "bb").
C
      ICNT=0
      IRPSZB(2)=IRPSZA(IRPDX+1)
      IWOFF2(1)=IWOFF1(IRPDX+1)
      DO 20 I=2,NIRREP
        DO 21 J=1,I-1
          ICNT=ICNT+2
          IRPDS2(ICNT-1)=NBFIRR(I)*(NBFIRR(I)+1)/2
          IRPDS2(ICNT)=NBFIRR(J)*(NBFIRR(J)+1)/2
          IRPSZ2(ICNT/2)=IRPDS2(ICNT-1)*IRPDS2(ICNT)
          MXLSSZ1=MAX(IRPSZ2(ICNT/2),MXLSSZ1)
          IRPSZB((ICNT/2)+2)=IRPSZB((ICNT/2)+1)+IRPSZ2(ICNT/2)*IINTFP
          JSIZE=JSIZE+IRPSZ2(ICNT/2)
          IPKSIZ=IPKSIZ+IRPSZ2(ICNT/2)
   21   CONTINUE
   20 CONTINUE
      IRPDX=(NIRREP*(NIRREP+1)/2)-NIRREP
      DO 25 I=2,IRPDX+1
        IWOFF2(I)=IWOFF2(I-1)+IRPSZ2(I-1)
   25 CONTINUE
C
      ICNT=0
      CALL IZERO(IDXVEC,I255)
      DO 70 I=1,NIRREP
        DO 71 J=1,NBFIRR(I)
          ICNT=ICNT+1
          IDXVEC(ICNT)=I
   71   CONTINUE
   70 CONTINUE
      DO 72 I=1,NBAS
        INEWVC(I)=I-IRPOFF(IDXVEC(I))
   72 CONTINUE
C
C  Set up the offsets for the in-core PK lists.
C
      DO 73 I=1,NIRREP
        DO 74 J=1,I
          IF(I.EQ.J) THEN
            IPKOFF(INDX(I,J)+1)=IPKOFF(INDX(I,J))+IRPSZ1(INDX(I,J))
          ELSE
            IPKOFF(INDX(I,J)+1)=IPKOFF(INDX(I,J))+IRPSZ2(INDX(I-1,J))
          ENDIF
   74   CONTINUE
   73 CONTINUE
      IF(IUHF.EQ.1) THEN
        IPKSIZ=2*IPKSIZ
        IJUNK=INDX(NIRREP,NIRREP)+1
        IPKOFF(37)=IPKOFF(IJUNK)
        DO 75 I=1,NIRREP
          DO 76 J=1,I
            IF(I.EQ.J) THEN
              IPKOFF(INDX(I,J)+37)=IPKOFF(INDX(I,J)+36)+
     &                             IRPSZ1(INDX(I,J))
            ELSE
              IPKOFF(INDX(I,J)+37)=IPKOFF(INDX(I,J)+36)+
     &                             IRPSZ2(INDX(I-1,J))
            ENDIF
   76     CONTINUE
   75   CONTINUE
      ENDIF
C
      TOTSIZ=ISIZE+JSIZE
      NAAAA=NIRREP
      NABAB=(NIRREP*(NIRREP+1)/2)-NIRREP
      NAABB=NABAB
C
C  The following variables are needed for the partial out-of-core
C  transformation.
C
      NAOBK=NAAAA+2*NABAB
      IBKOFF(1)=0
      IBKOFF(2)=IBKOFF(1)+NAAAA+NABAB
      IBKOFF(3)=IBKOFF(2)+NAABB
C
      CALL GETREC(20,'JOBARC','NATOMS  ',IONE,NATOMS)
      CALL GETREC(20,'JOBARC','FULLORDR',IONE,IORDGRP)
C
CJDW 11/13/96.
      CALL GETREC(20,'JOBARC','NAOBASFN',1,NBASX)
C
C  Now determine how much memory is required for in-core and out-of-core
C  algorithms and decide which to do.
C
C getaoso
C
c      IMEM1=7*NBAS+7*MOD(NBAS,2)+8*NBAS+8*NBAS*IINTFP+NBAS*NBAS*IINTFP
       IMEM1=0
C
C in-core
C
C mkpk
C
      IF(.NOT.AOFIL) THEN
       IMEMCOR=IPKSIZ*IINTFP
       IMEMX1=IMEMCOR+ILNBUF*IINTFP+ILNBUF
       IMAXMEM=MAX(IMEM1,IMEMX1)
      ELSE
       IMEMCOR=0
       IMAXMEM=IMEM1
      ENDIF
C
C get1e
C
      IMEMCOR=IMEMCOR+3*ITRILN(NIRREP+1)*IINTFP
      IMEMX1=IMEMCOR+ILNBUF*IINTFP+ILNBUF+NBAS
      IMAXMEM=MAX(IMAXMEM,IMEMX1)
c
c refiled2
c
      REACTION_FIELD = IFLAGS2(118) .EQ. 1
C
      IF (REACTION_FIELD) THEN
         IMEMCOR = 4*NBAS*NBAS*IINTFP + 2*ILNBUF*IINTFP
     &             + MXIRR2*IINTFP + IMEMCOR + 3*MXIRR2
         IMEMX1  = IMEMCOR + IMEMX1
         IMAXMEM = MAX(IMAXMEM,IMEMX1)
      END IF
C
C dovrlp
C
      IMEMX1=IMEMCOR+4*MXIRR2*IINTFP+NBAS*NBAS*IINTFP+NBAS*IINTFP+NBAS
      IMAXMEM=MAX(IMAXMEM,IMEMX1)
C
C initges
C
      IMEMCOR=IMEMCOR+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP
      IMEMX1=IMEMCOR+2*NBAS*NBAS*IINTFP+4*MXIRR2*IINTFP+
     &       (IUHF+1)*NBAS*IINTFP
CJDW
     &      +(IUHF+2)*NBAS *NBAS *IINTFP
     &      +(IUHF+1)*NBASX*NBAS *IINTFP
     &      +         NBASX*NBASX*IINTFP * 4
     &      + MAX(2000,NATOMS+MOD(NATOMS,2),
     &            NBASX+MOD(NBASX,2),
     &            NBASX*NBASX*IINTFP)
     &      +         NBAS *NBASX*IINTFP * 2
     &      +         NBAS *NBASX*IINTFP * 1
     &      +         NBASX      *IINTFP * 21
     &      + MAX(1080,NATOMS+MOD(NATOMS,2),IINTFP*NBASX*5)
     &      +         3 * NATOMS *IINTFP * 3
     &      +         NATOMS     *IINTFP * 2
     &      +         IORDGRP    *IINTFP
      IMAXMEM=MAX(IMAXMEM,IMEMX1)
C
C mkfock
C
      IF(.NOT.AOFIL) THEN
       IMEMCOR=IMEMCOR+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP
       IMEMX1=IMEMCOR+MXIRR2*IINTFP+NBAS*NBAS*IINTFP
       IMAXMEM=MAX(IMAXMEM,IMEMX1)
      ELSE
C
C mkrhff and/or mkuhff
C
       IMEMCOR=IMEMCOR+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP
       IMEMX1=IMEMCOR+ILNBUF*(1+IINTFP)+NBAST*NBAST
       IMAXMEM=MAX(IMAXMEM,IMEMX1)
      ENDIF
C
C mkrohff
C
       IF(IFLAGS(11).EQ.2) THEN
        IMEMX1=IMEMCOR+IINTFP*(2*ITRILN(NIRREP+1)+3*
     &                         ISQRLN(NIRREP+1))
        IMAXMEM=MAX(IMAXMEM,IMEMX1)
       ENDIF
C
C mkener
C
      IMEMX1=IMEMCOR+ITRILN(NIRREP+1)*IINTFP+2*MXIRR2*IINTFP
      IMAXMEM=MAX(IMAXMEM,IMEMX1)
C
C dorpp
C
      IF(IFLAGS(10).NE.0) THEN
        IMEMX1=IMEMCOR+NBAS*NBAS*IINTFP+
     &         (IUHF+1)*ISQRLN(NIRREP+1)*IINTFP+4*MXIRR2*IINTFP+
     &         MAX(MXIRR2,(IORDER+1)*(IORDER+2))*IINTFP+
     &         MAX(MXIRR2,(IORDER+1)*(IORDER+1))*IINTFP
        IMAXMEM=MAX(IMAXMEM,IMEMX1)
      ENDIF
C
C scfit
C
      IMEMCOR=IMEMCOR+(IUHF+1)*NBAS*IINTFP+
     &        (IUHF+1)*ISQRLN(NIRREP+1)*IINTFP + 2*NBAS*NBAS*IINTFP
      IMEMX1=IMEMCOR+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP+
     &       4*MAX(2*NBAS,MXIRR2)*IINTFP + 2*NBAS*NBAS*IINTFP
      if (ndimrpa.ne.0) then
         IMEMX1 = IMEMX1 + 2*NBAS*NBAS*IINTFP
     &          + 7*NDIMRPA*IINTFP + 3*ISQRLN(NIRREP+1)*(IUHF+1)*IINTFP
      end if
      IMAXMEM=MAX(IMAXMEM,IMEMX1)
C
C final mkfock
C
      IF(.NOT.AOFIL) THEN
       IMEMX1=IMEMCOR+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP+MXIRR2*IINTFP+
     &        NBAS*NBAS*IINTFP
       IMAXMEM=MAX(IMAXMEM,IMEMX1)
      ENDIF
C
C final mkener
C
      IMEMX1=IMEMCOR+ITRILN(NIRREP+1)*IINTFP+2*MXIRR2*IINTFP
      IMAXMEM=MAX(IMAXMEM,IMEMX1)
C
C outevc
C
      IMEMX1=IMEMCOR+2*NBAS+2*MOD(NBAS,2)
      IMAXMEM=MAX(IMAXMEM,IMEMX1)
C
C s2calc
C
      IF(IUHF.EQ.1) THEN
        IMEMX1=IMEMCOR+NBAS*NBAS*IINTFP+2*MXIRR2*IINTFP
        IMAXMEM=MAX(IMAXMEM,IMEMX1)
      ENDIF
C
C dmpjob
C
c     IMEMX1=IMEMCOR+NBAS*NBAS*IINTFP+MXIRR2*IINTFP+NBAS+MOD(NBAS,2)+
c    &       4*ISQRLN(NIRREP+1)*IINTFP
      IMEMX1=IMEMCOR+2*NBAS*NBAS*IINTFP+MXIRR2*IINTFP+MAX(2*NIRREP,
     &       NBAS+MOD(NBAS,2))+4*ISQRLN(NIRREP+1)*IINTFP + 3*NBASX*
     &       NBASX*IINTFP
      IMAXMEM=MAX(IMAXMEM,IMEMX1)
C
C evcao
C
      IMEMX1=IMEMCOR+3*NBAS*NBAS*IINTFP
      IMAXMEM=MAX(IMAXMEM,IMEMX1)
C
C  evlout
C
      IMEM2=(IUHF+1)*NBAS*IINTFP+2*NBAS+2*MOD(NBAS,2)+2*NBAS*IINTFP
      IMAXMEM=MAX(IMAXMEM,IMEM2)
C
C  prjden
C
      CALL GETREC(20,'JOBARC','NAOBASFN',IONE,NBASX)
      IMEM2=2*NBASX+7*NATOMS+IORDGRP*NATOMS+1+9*IORDGRP*IINTFP+
     &      2*NBASX*NBASX*IINTFP+9*IINTFP+NBASX+8*NBASX*NBASX*IINTFP+
     &      4*IINTFP*NATOMS
C
      MEMIN=MAX(IMAXMEM,IMEM2)
C
C out-of-core
C
C get1e
C
      IMEMCOR=3*ITRILN(NIRREP+1)*IINTFP
      IMEMX1=IMEMCOR+ILNBUF*IINTFP+ILNBUF+NBAS
      IMAXMEM=MAX(IMEM1,IMEMX1)
C
C dovrlp
C
      IMEMX1=IMEMCOR+4*MXIRR2*IINTFP+NBAS*NBAS*IINTFP+NBAS*IINTFP+NBAS
      IMAXMEM=MAX(IMAXMEM,IMEMX1)
C
C initges
C
      IMEMCOR=IMEMCOR+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP
      IMEMX1=IMEMCOR+2*NBAS*NBAS*IINTFP+4*MXIRR2*IINTFP+
     &       (IUHF+1)*NBAS*IINTFP
CJDW
     &      +(IUHF+2)*NBAS *NBAS *IINTFP
     &      +(IUHF+1)*NBASX*NBAS *IINTFP
     &      +         NBASX*NBASX*IINTFP * 4
     &      + MAX(2000,NATOMS+MOD(NATOMS,2),
     &            NBASX+MOD(NBASX,2),
     &            NBASX*NBASX*IINTFP)
     &      +         NBAS *NBASX*IINTFP * 2
     &      +         NBAS *NBASX*IINTFP * 1
     &      +         NBASX      *IINTFP * 21
     &      + MAX(1080,NATOMS+MOD(NATOMS,2),IINTFP*NBASX*5)
     &      +         3 * NATOMS *IINTFP * 3
     &      +         NATOMS     *IINTFP * 2
     &      +         IORDGRP    *IINTFP
      IMAXMEM=MAX(IMAXMEM,IMEMX1)
C
C mkfock2
C
      IF(.NOT.AOFIL) THEN
      IMEMCOR=IMEMCOR+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP
      IMEMX1=IMEMCOR+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP+
     &       MXIRR2*IINTFP+NBAS*NBAS*IINTFP
      IMAXMEM=MAX(IMAXMEM,IMEMX1)
      ENDIF
C
C mkener
C
      IMEMX1=IMEMCOR+ITRILN(NIRREP+1)*IINTFP+2*MXIRR2*IINTFP
      IMAXMEM=MAX(IMAXMEM,IMEMX1)
C
C mkrohff
C
      IF(IFLAGS(11).EQ.2) THEN
       IMEMX1=IMEMCOR+IINTFP*(2*ITRILN(NIRREP+1)+3*
     &                        ISQRLN(NIRREP+1))
       IMAXMEM=MAX(IMAXMEM,IMEMX1)
      ENDIF
C
C dorpp
C
      IF(IFLAGS(10).NE.0) THEN
        IMEMX1=IMEMCOR+NBAS*NBAS*IINTFP+
     &         (IUHF+1)*ISQRLN(NIRREP+1)*IINTFP+4*MXIRR2*IINTFP+
     &         MAX(MXIRR2,(IORDER+2)*(IORDER+1))*IINTFP+
     &         MAX(MXIRR2,(IORDER+1)*(IORDER+1))*IINTFP
        IMAXMEM=MAX(IMAXMEM,IMEMX1)
      ENDIF
C
C scfit
C
      IMEMCOR=IMEMCOR+(IUHF+1)*NBAS*IINTFP+
     &        (IUHF+1)*ISQRLN(NIRREP+1)*IINTFP
      IMEMX1=IMEMCOR+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP+
     &       (3*MAX(2*NBAS,NBAS*NBAS,MXIRR2)+MAX(2*NBAS,MXIRR2))*IINTFP
      if (ndimrpa.ne.0) then
         IMEMX1 = IMEMX1 + 2*NBAS*NBAS*IINTFP
     &          + 7*NDIMRPA*IINTFP + 3*ISQRLN(NIRREP+1)*(IUHF+1)*IINTFP
      end if
      IMAXMEM=MAX(IMAXMEM,IMEMX1)
C
C final mkfock2
C
      IF(.NOT.AOFIL) THEN
      IMEMX1=IMEMCOR+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP+
     &       MXIRR2*IINTFP+NBAS*NBAS*IINTFP
      IMAXMEM=MAX(IMAXMEM,IMEMX1)
      ENDIF
C
C final mkener
C
      IMEMX1=IMEMCOR+ITRILN(NIRREP+1)*IINTFP+2*MXIRR2*IINTFP
      IMAXMEM=MAX(IMAXMEM,IMEMX1)
C
C outevc
C
      IMEMX1=IMEMCOR+2*NBAS+2*MOD(NBAS,2)
      IMAXMEM=MAX(IMAXMEM,IMEMX1)
C
C s2calc
C
      IF(IUHF.EQ.1) THEN
        IMEMX1=IMEMCOR+NBAS*NBAS*IINTFP+2*MXIRR2*IINTFP
        IMAXMEM=MAX(IMAXMEM,IMEMX1)
      ENDIF
c
C UNO-REF
c
       if (iflags2(148) .ne. 0) then
          imemx1 = imemcor + 8*jtot*iintfp
          IMAXMEM=MAX(IMAXMEM,IMEMX1)
       endif
C
C dmpjob
C
c     IMEMX1=IMEMCOR+NBAS*NBAS*IINTFP+MXIRR2*IINTFP+NBAS+MOD(NBAS,2)+
c    &       4*ISQRLN(NIRREP+1)*IINTFP
      IMEMX1=IMEMCOR+2*NBAS*NBAS*IINTFP+MXIRR2*IINTFP+MAX(2*NIRREP,
     &       NBAS+MOD(NBAS,2))+4*ISQRLN(NIRREP+1)*IINTFP +3*NBASX*
     &       NBASX
      IMAXMEM=MAX(IMAXMEM,IMEMX1)
C
C evcao
C
      IMEMX1=IMEMCOR+3*NBAS*NBAS*IINTFP
      IMAXMEM=MAX(IMAXMEM,IMEMX1)
C
C evlout
C
      IMAXMEM=MAX(IMAXMEM,IMEM2)
C
C  prjden
C
      IMEM2=2*NBAS+5*NATOMS+IORDGRP*NATOMS+9*IORDGRP*IINTFP+
     &      2*NBAS*NBAS*IINTFP+5*NBAS*NBAS*IINTFP+4*NATOMS*IINTFP
C
C memory for out-of-core
C
      MEMOUT=MAX(IMAXMEM,IMEM2)
C
C  Now set up the memory locations for the DIIS error matrices
C
      IF(IFLAGS(10).NE.0) THEN
        ISTART=1
        IRPP1=IORDER*ISQRLN(NIRREP+1)*IINTFP
        IRPPMEM=(IUHF+2)*IRPP1
        IADERR(1)=1
        IADFCK(1,1)=IADERR(1)+IRPP1
        IF(IUHF.EQ.1) THEN
          IADFCK(1,2)=IADFCK(1,1)+IRPP1
        ENDIF
        DO 300 I=2,IORDER
          IADERR(I)=ISQRLN(NIRREP+1)*IINTFP+IADERR(I-1)
          IADFCK(I,1)=ISQRLN(NIRREP+1)*IINTFP+IADFCK(I-1,1)
          IF(IUHF.EQ.1) THEN
            IADFCK(I,2)=ISQRLN(NIRREP+1)*IINTFP+IADFCK(I-1,2)
          ENDIF
  300   CONTINUE
      ELSE
        IRPPMEM=0
      ENDIF
C
C  Now determine whether the SCF will be done using in-core or out-of-core
C  PK lists (i.e., full in-core or out-of-core fock matrix construction).
C  Additionally, if we're going out-of-core, kick up the file record-size
C  for the PK list file so we can get better caching.
C
      IRECSIZ=IFLAGS(37)
      INUMREC=IFLAGS(44)
      ICACHSZ=IRECSIZ*INUMREC
      IF(MEMIN+IRPPMEM+ICACHSZ.GT.MAXCOR) THEN
       IF(AOFIL) CALL ERREX
        IF(MEMOUT+IRPPMEM+ICACHSZ.GT.MAXCOR) THEN
          INOUT=3
        ELSE
          IF(IRECSIZ.LT.65536) THEN
            IRECSIZ=65536
  441       ICACHSZ=IRECSIZ*INUMREC
            IF(MEMOUT+IRPPMEM+ICACHSZ.GT.MAXCOR) THEN
              IRECSIZ=IRECSIZ/2
              GOTO 441
            ENDIF
            WRITE(LUOUT,5001)IFLAGS(37),IRECSIZ
 5001       FORMAT(T2,'@SYMSIZ-I, Out-of-core record size changed ',
     &                'from ',I6,' to ',I6,' words.',/)
            IFLAGS(37)=IRECSIZ
          ENDIF
          MAXMEM=MAX(MAXCOR-ICACHSZ-IRPPMEM,MINMEM)
          INOUT=2
        ENDIF
      ELSE
        MAXMEM=MAX(MEMIN,MINMEM)
        INOUT=1
      ENDIF
C
      IF(IRECSIZ.LT.128) THEN
        WRITE(LUOUT,9101)
 9101   FORMAT(T3,'@SYMSIZ-F, Record size for direct access ',
     &            'files smaller than 128 words.',/,
     &         T14,'This may cause VERY slow I/O.  If this is ',
     &             'to be allowed, then modify routine SYMSIZ',/,
     &         T14,'in this program.',/)
        CALL ERREX
      ENDIF
C
C  The total number of possibilities for each of the symmetry
C  combinations is given below.  The variables are self-explanatory.
C  For the (ab|cd) symmetry integrals, there are 210 possible ways
C  to uniquely combine irrep indices for D2h (largest Abelian subgroup),
C  which arises from taking 8 things 4 at a time, and then multiplying
C  by 3.  Of the 210 possibilities, only 42 are actually allowed by
C  symmetry.  Of course, which ones these are depends on the direct
C  product table (which is why it's read in).
C
      IF (.NOT.AOFIL) THEN
      IF (FIRSTSCF) THEN
         WRITE(LUOUT,5000)NAAAA,NABAB,NAABB,IAASZ,IABSZ,JSIZE,TOTSIZ
 5000 FORMAT(T3,'@SYMSIZ-I, Symmetry integral information:',/,
     &       T8,'There are ',I3,' (aa|aa) symmetry combinations.',/,
     &       T8,'There are ',I3,' (ab|ab) symmetry combinations.',/,
     &       T8,'There are ',I3,' (aa|bb) symmetry combinations.',/,
     &          ' ',/,
     &       T8,'The (aa|aa) integrals require ',I10,' words.',/,
     &       T8,'The (ab|ab) integrals require ',I10,' words.',/,
     &       T8,'The (aa|bb) integrals require ',I10,' words.',/,
     &       T8,'The integrals require ',I10,' words.',/,
     &          ' ')
      END IF
      END IF
C
      IF(AOFIL) THEN
       IF(FIRSTSCF) WRITE(LUOUT,5041) MEMIN+IRPPMEM+ICACHSZ
 5041  FORMAT(T3,'@SYMSIZ-I, Memory information: ',I10,
     &        ' words required.')
       IF(FIRSTSCF) WRITE(LUOUT,5042)
 5042  FORMAT(T3,'@SYMSIZ-I, Fock matrices are constructed',
     &           ' from AO integral file.')
      ELSE
      IF(FIRSTSCF) WRITE(LUOUT,5051)MEMIN+IRPPMEM+ICACHSZ,MAXCOR
 5051 FORMAT(T3,'@SYMSIZ-I, Memory information: ',/,
     &       T8,'Full in-core PK list formation: ',I10,/,
     &       T8,'      Maximum memory available: ',I10,/)
      IF(INOUT.EQ.1) THEN
        IF(FIRSTSCF) WRITE(LUOUT,5061)
 5061   FORMAT(T3,'@SYMSIZ-I, Selecting full in-core PK list ',
     &            'formation.',/)
      ELSEIF(INOUT.EQ.2) THEN
        IF(FIRSTSCF) WRITE(LUOUT,5071)
 5071   FORMAT(T3,'@SYMSIZ-I, Selecting out-of-core PK list ',
     &            'formation.',/)
      ELSEIF(INOUT.EQ.3) THEN
        WRITE(LUOUT,9000)MEMOUT+IRPPMEM+ICACHSZ,MAXCOR
 9000   FORMAT(T3,'@SYMSIZ-F, Insufficient memory to do SCF.  ',
     &            'Sorry.',/,
     &         T8,'Minimum memory for out-of-core algorithm: ',I10,/,
     &         T8,'                Maximum memory available: ',I10)
        CALL ERREX
      ENDIF
      ENDIF
C
      RETURN
      END
