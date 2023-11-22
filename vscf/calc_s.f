      SUBROUTINE CALC_S(SMAT,COORD,IZ,NBASEHT,NBASXEHT,NATOMS)
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
      DOUBLE PRECISION SMAT,COORD
      INTEGER IZ,NBASEHT,NBASXEHT,NATOMS
C-----------------------------------------------------------------------
      DOUBLE PRECISION D,Y,Z,ZETAI,ZETAJ,ZETAA,ZETAB,PAIRS,R,C1,C2,E,T
      DOUBLE PRECISION TEMP
      INTEGER NORB,N,L,M,IATOM,IOFFATOM,LEN,NSHELL,MAPSHELL,
     &        JATOM,NORBI,NORBJ,NSHELLI,NSHELLJ,I,J,K,
     &        ISHELL,JSHELL,NA,NB,LA,LB,MVAL,IORB,JORB,IVAL,JVAL,MAXL
C-----------------------------------------------------------------------
      INTEGER IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      INTEGER IFLAGS
C-----------------------------------------------------------------------
      DOUBLE PRECISION FACT,SS
C-----------------------------------------------------------------------
      DOUBLE PRECISION TWO
C-----------------------------------------------------------------------
      DIMENSION SMAT(NBASXEHT,NBASXEHT),COORD(3,NATOMS),IZ(NATOMS)
      DIMENSION D(0:3,0:3,0:3),Y(0:8,0:8,203),Z(0:16,0:8,0:8)
      DIMENSION NORB(3),N(3),L(4,3),M(4,3),IOFFATOM(100),LEN(100),
     &          MAPSHELL(4,3,3),PAIRS(9,9),C1(3),C2(3),E(3),T(9,9),
     &          TEMP(9,9)
C-----------------------------------------------------------------------
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/  IFLAGS(100)
C-----------------------------------------------------------------------
      DATA TWO /2.0D+00/
C-----------------------------------------------------------------------

      CALL GETREC(20,'JOBARC','COORD',NATOMS*3*IINTFP,COORD)
      CALL GETREC(20,'JOBARC','ATOMCHRG',NATOMS,IZ)
C-----------------------------------------------------------------------
C     Calculate tables needed in SS.
C-----------------------------------------------------------------------
      CALL DTABLE(D)
      CALL YTABLE(Y)
      CALL ZTABLE(Z)
C-----------------------------------------------------------------------
C     Number of orbitals in each shell : our model is minimal basis set
C     where s and p functions with same pqn have same exponent --- ie we
C     have s and sp shells.
C-----------------------------------------------------------------------
      NORB(1)=1
      NORB(2)=4
      NORB(3)=4
C-----------------------------------------------------------------------
C     Principal quantum number of each shell (trivial in our MBS model).
C-----------------------------------------------------------------------
      N(1)=1
      N(2)=2
      N(3)=3
C-----------------------------------------------------------------------
C     Angular momentum array. Angular momentum of each orbital in a
C     shell.
C-----------------------------------------------------------------------
      L(1,1)=0
C
      L(1,2)=0
      L(2,2)=1
      L(3,2)=1
      L(4,2)=1
C
      L(1,3)=0
      L(2,3)=1
      L(3,3)=1
      L(4,3)=1
C-----------------------------------------------------------------------
C     M array.
C-----------------------------------------------------------------------
      M(1,1)= 0
C
      M(1,2)= 0
      M(2,2)= 1
      M(3,2)=-1
      M(4,2)= 0
C
      M(1,3)= 0
      M(2,3)= 1
      M(3,3)=-1
      M(4,3)= 0
C-----------------------------------------------------------------------
C     Calculate atomic offsets in the minimal basis set.
C-----------------------------------------------------------------------
      CALL IZERO(IOFFATOM,NATOMS)
      CALL IZERO(LEN     ,NATOMS)
      DO 30 IATOM =1,NATOMS
C
      NSHELLI = 0
      IF(IZ(IATOM) .EQ.  1 .OR. IZ(IATOM) .EQ.  2)THEN
       NSHELLI = 1
      ELSEIF(IZ(IATOM) .GE.  3 .AND. IZ(IATOM) .LE. 10)THEN
       NSHELLI = 2
      ELSEIF(IZ(IATOM) .GE. 11 .AND. IZ(IATOM) .LE. 18)THEN
       NSHELLI = 3
      ENDIF
C
      IF(NSHELLI .GT. 0)THEN
       DO 20 ISHELL=1,NSHELLI
       LEN(IATOM) = LEN(IATOM) + NORB(ISHELL)
   20  CONTINUE
      ENDIF
C
      IF(IATOM .EQ. 1)THEN
       IOFFATOM(IATOM) = 0
      ELSE
       IOFFATOM(IATOM) = IOFFATOM(IATOM-1) + LEN(IATOM-1)
      ENDIF
C
      IF(IFLAGS(1).GE.10)THEN
       write(6,*) ' @calc_s-i, ioffatom, len ',ioffatom(1),ioffatom(2)
       write(6,*) len(1),len(2)
      ENDIF
C
   30 CONTINUE
C
C-----------------------------------------------------------------------
C     Calculate shell orbital mapping to the usually ordered basis set.
C-----------------------------------------------------------------------
      CALL IZERO(MAPSHELL,4*3)
      NSHELL = 1
      MAPSHELL(1,1,NSHELL) = 1
      NSHELL = 2
      MAPSHELL(1,1,NSHELL) = 1
      MAPSHELL(1,2,NSHELL) = 2
      MAPSHELL(2,2,NSHELL) = 3
      MAPSHELL(3,2,NSHELL) = 4
      MAPSHELL(4,2,NSHELL) = 5
      NSHELL = 3
      MAPSHELL(1,1,NSHELL) = 1
      MAPSHELL(1,2,NSHELL) = 2
c
      MAPSHELL(2,2,NSHELL) = 4
      MAPSHELL(3,2,NSHELL) = 5
      MAPSHELL(4,2,NSHELL) = 6
c orig
      MAPSHELL(1,3,NSHELL) = 2
      MAPSHELL(2,3,NSHELL) = 4
      MAPSHELL(3,3,NSHELL) = 5
      MAPSHELL(4,3,NSHELL) = 6
c correct?
      MAPSHELL(1,3,NSHELL) = 3
      MAPSHELL(2,3,NSHELL) = 7
      MAPSHELL(3,3,NSHELL) = 8
      MAPSHELL(4,3,NSHELL) = 9
C-----------------------------------------------------------------------
C
      CALL ZERO(SMAT,NBASXEHT*NBASXEHT)
C
      DO 500 JATOM=1,NATOMS
      DO 490 IATOM=1,NATOMS
C
      NORBJ=0
      IF(IZ(JATOM) .GE. 1 .AND. IZ(JATOM) .LE. 2)THEN
       NORBJ   = 1
       NSHELLJ = 1
      ELSEIF(IZ(JATOM) .GE. 3 .AND. IZ(JATOM) .LE. 10)THEN
       NORBJ   = 5
       NSHELLJ = 2
      ELSEIF(IZ(JATOM) .GE.11 .AND. IZ(JATOM) .LE. 18)THEN
       NORBJ   = 9
       NSHELLJ = 3
      ENDIF
C
      NORBI=0
      IF(IZ(IATOM) .GE. 1 .AND. IZ(IATOM) .LE. 2)THEN
       NORBI   = 1
       NSHELLI = 1
      ELSEIF(IZ(IATOM) .GE. 3 .AND. IZ(IATOM) .LE. 10)THEN
       NORBI   = 5
       NSHELLI = 2
      ELSEIF(IZ(IATOM) .GE.11 .AND. IZ(IATOM) .LE. 18)THEN
       NORBI   = 9
       NSHELLI = 3
      ENDIF
C
C     Skip dummies. For now this does not allow ghosts, elements with
C     Z>18. But we want to get this working for "simple" elements first !
C
      IF(NSHELLI.EQ.0 .OR. NSHELLJ.EQ.0) GOTO 490
C
      DO 10 I=1,3
      C1(I) = COORD(I,IATOM)
      C2(I) = COORD(I,JATOM)
   10 CONTINUE
C
      CALL RELVEC(R,E,C1,C2)
C
      IF(IFLAGS(1).GE.10)THEN
       write(6,*) ' c1,c2,e,r '
       write(6,*) c1
       write(6,*) c2
       write(6,*) r
       write(6,*) e
      ENDIF
C
C     Exponents are set according to Slater's rules, as described by
C     Pople and Beveridge, pages 27-29. H is an exception.
C
      DO 480 JSHELL=1,NSHELLJ
      DO 470 ISHELL=1,NSHELLI
C
C     n=1
      IF(ISHELL.EQ.1)THEN
       IF(IZ(IATOM) .EQ. 1)THEN
        ZETAI = 1.2D+00
       ELSE
        ZETAI = DFLOAT(IZ(IATOM)) - 0.3D+00
       ENDIF
      ENDIF
C
      IF(JSHELL.EQ.1)THEN
       IF(IZ(JATOM) .EQ. 1)THEN
        ZETAJ = 1.2D+00
       ELSE
        ZETAJ = DFLOAT(IZ(JATOM)) - 0.3D+00
       ENDIF
      ENDIF
C
C     n=2
      IF(ISHELL.EQ.2)THEN
C Li - Ne
       IF(IZ(IATOM).GE. 3 .AND. IZ(IATOM).LE.10)THEN
        ZETAI = 0.325D+00 * DFLOAT(IZ(IATOM)-1)
       ENDIF
C Na - Ar
       IF(IZ(IATOM).GE.11 .AND. IZ(IATOM).LE.18)THEN
        ZETAI = ( DFLOAT(IZ(IATOM)) - 
     &            (7.0D+00 * 0.35D+00 + 2.0D+00 * 0.85D+00) ) / 2.0D+00
       ENDIF
      ENDIF
C
      IF(JSHELL.EQ.2)THEN
       IF(IZ(JATOM).GE. 3 .AND. IZ(JATOM).LE.10)THEN
        ZETAJ = 0.325D+00 * DFLOAT(IZ(JATOM)-1)
       ENDIF
C
       IF(IZ(JATOM).GE.11 .AND. IZ(JATOM).LE.18)THEN
        ZETAJ = ( DFLOAT(IZ(JATOM)) - 
     &            (7.0D+00 * 0.35D+00 + 2.0D+00 * 0.85D+00) ) / 2.0D+00
       ENDIF
      ENDIF
C
C     n=3
      IF(ISHELL.EQ.3)THEN
       ZETAI = 0.65D+00 * DFLOAT(IZ(IATOM)) - 4.95D+00
C Na - Ar
       IF(IZ(IATOM).GE.11 .AND. IZ(IATOM).LE.18)THEN
        ZETAI = ( DFLOAT(IZ(IATOM)) - 
     &           (DFLOAT(IZ(IATOM) - 11) * 0.35D+00 +
     &                           8.0D+00 * 0.85D+00 +
     &                           2.0D+00 * 1.00D+00)) / 3.0D+00
       ENDIF
      ENDIF
C
      IF(JSHELL.EQ.3)THEN
       ZETAJ = 0.65D+00 * DFLOAT(IZ(JATOM)) - 4.95D+00
       IF(IZ(JATOM).GE.11 .AND. IZ(JATOM).LE.18)THEN
        ZETAJ = ( DFLOAT(IZ(JATOM)) - 
     &           (DFLOAT(IZ(JATOM) - 11) * 0.35D+00 +
     &                           8.0D+00 * 0.85D+00 +
     &                           2.0D+00 * 1.00D+00)) / 3.0D+00
       ENDIF
      ENDIF
C
      CALL ZERO(PAIRS,9*9)
C
      IF(IATOM .NE. JATOM)THEN
C
       DO 460 JORB=1,NORB(JSHELL)
       DO 450 IORB=1,NORB(ISHELL)
C
       IF( M(IORB,ISHELL) .EQ. M(JORB,JSHELL) .AND. 
     &     M(IORB,ISHELL) .GE. 0                    )THEN
C
        NA    = N(ISHELL)
        NB    = N(JSHELL)
        LA    = L(IORB,ISHELL)
        LB    = L(JORB,JSHELL)
        MVAL  = M(IORB,ISHELL)
        ZETAA = ZETAI
        ZETAB = ZETAJ
C
        PAIRS(IORB,JORB) = 
     &        DSQRT( (ZETAA*R)**(2*NA+1) * (ZETAB*R)**(2*NB+1) /
     &                 ( FACT(2*NA)      *     FACT(2*NB) )       ) *
     &                   SS(NA,LA,MVAL,NB,LB,ZETAA*R,ZETAB*R,D,Y,Z)
C
        IF(IFLAGS(1).GE.10)THEN
         write(6,*) ' m, ss,iorb,jorb,pairs(iorb,jorb) ',mval
         write(6,*) ss(na,la,mval,nb,lb,zetaa*r,zetab*r,d,y,z)
         write(6,*) pairs(iorb,jorb)
        ENDIF
C
       ENDIF
C
C     Deal with py-py case (M=-1).
C
       IF( M(IORB,ISHELL) .EQ. M(JORB,JSHELL) .AND.
     &     M(IORB,ISHELL) .LT. 0                    )THEN
        PAIRS(IORB,JORB) = PAIRS(IORB-1,JORB-1)
       ENDIF
C
  450  CONTINUE
  460  CONTINUE
C
C     Rotate the shell overlaps.
C
       IF(ISHELL .GE. 2 .OR. JSHELL .GE. 2)THEN
        MAXL = 1
       ELSE
        MAXL = 0
       ENDIF
       CALL HARMTR(T,MAXL,E)
       DO 300 J=1,NORB(JSHELL)
       DO 290 I=1,NORB(ISHELL)
       TEMP(I,J)=0.0D+00
       DO 280 K=1,NORB(JSHELL)
       TEMP(I,J) = TEMP(I,J) + PAIRS(I,K)*T(J,K)
  280  CONTINUE
  290  CONTINUE
  300  CONTINUE
C
       DO 270 J=1,NORB(JSHELL)
       DO 260 I=1,NORB(ISHELL)
       PAIRS(I,J)=0.0D+00
       DO 250 K=1,NORB(ISHELL)
       PAIRS(I,J) = PAIRS(I,J) + T(I,K)*TEMP(K,J)
  250  CONTINUE
  260  CONTINUE
  270  CONTINUE
       
C
      ELSE
C
C
       NA    = N(ISHELL)
       NB    = N(JSHELL)
       ZETAA = ZETAI
       ZETAB = ZETAJ
C
       DO 420 JORB=1,NORB(JSHELL)
       DO 410 IORB=1,NORB(ISHELL)
C
       IF(IORB.EQ.JORB)THEN
        PAIRS(IORB,IORB) =
     &        DSQRT( (TWO*ZETAA)**(2*NA+1) * (TWO*ZETAB)**(2*NB+1) /
     &               ( FACT(2*NA)          *  FACT(2*NB) )          ) *
     &        FACT(NA+NB) / (ZETAA + ZETAB)**(NA+NB+1)
       ENDIF
C
  410  CONTINUE
  420  CONTINUE
C
      ENDIF
C
      IF(IFLAGS(1).GE.10)THEN
       write(6,*) ' @CALC_S-I, PAIRS. ISHELL, JSHELL ',ISHELL,JSHELL
       call output(pairs,1,9,1,9,9,9,1)
      ENDIF
C
C-----------------------------------------------------------------------
C     Fill out appropriate elements of full overlap matrix.
C-----------------------------------------------------------------------
      DO 440 JORB=1,NORB(JSHELL)
      DO 430 IORB=1,NORB(ISHELL)
      IVAL = IOFFATOM(IATOM) + MAPSHELL(IORB,ISHELL,NSHELLI)
      JVAL = IOFFATOM(JATOM) + MAPSHELL(JORB,JSHELL,NSHELLJ)
      SMAT(IVAL,JVAL) = PAIRS(IORB,JORB)
  430 CONTINUE
  440 CONTINUE
C
  470 CONTINUE
  480 CONTINUE
C
  490 CONTINUE
  500 CONTINUE
C
      write(6,*) ' @CALC_S-I, The overlap matrix (AO basis) '
      call output(smat,1,nbasxeht,1,nbasxeht,nbasxeht,nbasxeht,1)
      RETURN
      END
