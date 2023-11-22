      SUBROUTINE LIBINT_READBASIS_DATA(FNAME,IUATMS,NATOMS,NPOP,
     &                                 IREORDER,CARTESIAN,ITFCT, 
     &                                 LNP1,LNPO,NFCT,NBASIS,ALPHA, 
     &                                 IXALPHA,PCOEFF,IXPCOEF,
     &                                 MAX_ATOMS,ATMNAM,COORD,
     &                                 CORD,NAOATM)
C
C     ----INPUT ARGUMENTS----
C
C FNAME     = Unix file name of the MOL file.
C IUATMS    = The number of symmetry unique atoms.
C NATOMS    = The total number of atoms.
C IREORDER  = If there is any reordering of centers, this should
C             give the correspondence. Most often this is a unit vector.
C CARTESIAN = True for Cartesian basis choice.
C NPOP      = number of symmetry equivalent atom for given
C             symmetry unique center.
C NAOATM    = The total number of contracted functions on each atom.
C ITFCT     = Total number of primitive functions.
C LNP1      = Largest value for number of contracted functions for all
C             shells.
C LNPO      = Largest possible value for the product of number of primitives
C             of contracted functions for all shells.
C NFCT      = The total number of  primitive functions on each atom.
C
C     ----OUTPUT ARGUMENTS------
C
C ALPHA     = Exponents for basis functions
C IXALPHA   = Index array, each entry is the beginning of alpha's for the 
C             shell.
C PCOEFF    = The contraction coefficients.
C IXPCOEF   = Index array for PCOEF, each entry is the beginning of
C             pcoeff's for the shell.
C ATMNAM    = Array of Atom names
C COORD     = (x,y,z) coordinate of each shell.
C
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     
      CHARACTER*(*) FNAME 
      CHARACTER*80 XLINE
      CHARACTER*4 ATMNAM(NATOMS)
      LOGICAL CARTESIAN
C     
      DIMENSION NPOP(IUATMS),IREORDER(NATOMS), NFCT(NATOMS), 
     &          NAOATM(NATOMS), PCOEFF(ITFCT*NBASIS), 
     &          ALPHA(ITFCT), SALPHA(LNP1), SPCOEF(LNPO),
     &          COORD(3,*), CORD(3, MAX_ATOMS)

      INTEGER IXALPHA(*), IXPCOEF(*)
 
      DOUBLE PRECISION X, Y, Z 
      INTEGER ISHL(6)

      integer nxt_alpha, nxt_pcoef
C     
      PI=DATAN(1.0D+00)*4.D+00
      PICNST=(0.5D+00/PI)**0.75D+00
C
C Open MOL file for basis set information. The MOL file is created by 
C joda by processing user input file and basis set library.
C
      OPEN(UNIT=10, FILE=FNAME, FORM='FORMATTED', STATUS='OLD')
      REWIND(10)
C
C For the prsent puupose we can ignore the first 5 line of the MOL file.
C
      READ(10,'(A)') XLINE
      READ(10,'(A)') XLINE
      READ(10,'(A)') XLINE
      READ(10,'(A)') XLINE
      READ(10,'(A)') XLINE
C
      INI    = 1
      MAXANG = 0

      ishell = 0
      nxt_alpha = 1
      nxt_pcoef = 1 
      DO 10 IATM = 1, IUATMS
C
         READ(10,1000) ZNUC, IJUNK, NSHL,(ISHL(I),I=1,NSHL)
 1000    FORMAT(F20.1,8I5)
C
         READ(10,1010) ATMNAM(IATM), X, Y, Z
 1010    FORMAT(A4,3F20.12)
C     
         IDGOFF = 0
         IOFF   = 0
         JOFF   = 0
         NUNQSHL= 0

C-The ECP arrays added, Ajith Perera, 04/2012----
C
         CORD(1,IATM) = X
         CORD(2,IATM) = Y
         CORD(3,IATM) = Z
C---
C
         DO 28 LL = 1, NSHL
C
            NPT=0
            NAOT=0
C
            DO 30 II1 = 1, ISHL(LL)
               READ(10,1120) NP1, NAO
C
               ishell = ishell + 1
               coord(1, ishell) = x
               coord(2, ishell) = y
               coord(3, ishell) = z
               NPT  = NPT + NP1
               NAOT = NAOT + NAO
C
               IF (CARTESIAN) THEN
                  NP2 = LL*(LL + 1)/2
               ELSE
                  NP2 = 2*LL - 1
               ENDIF
C
               DO 32 I = 1, NP1
                  READ(10,*) SALPHA(I),(SPCOEF((J-1)*NP1+I),J=1,NAO)
 32            CONTINUE
C     
C Renormalize the atomic orbitals. Multiply the renormalized
C coefficients by the appropriate normalization constants.
C     
               DO 34 INAO = 1, NAO
C
                  SUM = 0.D+00
C
                  DO 36 I = 1, NP1
                     DO 37 J = 1, I
C
                        AI=SALPHA(I)
                        AJ=SALPHA(J)
C
                        TMP=SPCOEF((INAO-1)*NP1+I)*SPCOEF((INAO-1)*
     &                      NP1+J)*(2.0D+00*DSQRT(AI*AJ)/
     &                     (AI+AJ))**(REAL(LL)+0.5D+00)
C
                        SUM = SUM + TMP
                        IF(I .NE. J) SUM = SUM + TMP
C
 37                  CONTINUE
 36               CONTINUE
C
                  XNORM=1.D+00/DSQRT(SUM)
C
                  DO 38 I = 1, NP1
C                   SPCOEF((INAO-1)*NP1+I)= SPCOEF((INAO-1)*NP1+I)*
C     &                                       XNORM*PICNST*
C     &                                       (4.D+00*SALPHA(I))**
C     &                                       (0.5D+00*REAL(LL)+
C     &                                        0.25D+00)
                     SPCOEF((INAO-1)*NP1+I)= SPCOEF((INAO-1)*NP1+I)*
     &                                       XNORM
 38               CONTINUE
 34            CONTINUE
C     
C Place the alpha's and coefficients in their appropriate place in 
C their respective matrices, ALPHA and PCOEFF.
C
CCCCC               DO 40 IPOP = INI, INI + NPOP(IATM) - 1
C
CCCCC                  IATMOFF = 0
CCCCC                  JATMOFF = 0
C
CCCCC                  DO 43 ITMP = 1, IREORDER(IPOP)-1
CCCCC                     IATMOFF = IATMOFF + NFCT(ITMP)
CCCCC                     JATMOFF = JATMOFF + NAOATM(ITMP)*ITFCT
CCCCC 43               CONTINUE
C
CCCCC                  DO 45 I = 1, NP1
CCCCC                     JSHOFF=0
CCCCC                     DO 46 I1 = 1, NP2
CCCCC                        ALPHA(I+IOFF+IATMOFF+(I1-1)*NP1) = SALPHA(I)
C
CCCCC                        DO 57 J=1,NAO
CCCCC                           PCOEFF(JATMOFF+(JOFF+JSHOFF)*ITFCT+I+IOFF+
CCCCC     $                     IATMOFF+(I1-1)*NP1 )= SPCOEF((J-1)*NP1+I)
CCCCC                           JSHOFF=JSHOFF+1
C
CCCCC 57                     CONTINUE
CCCCC 46                  CONTINUE
CCCCC 45               CONTINUE
CCCCC 40            CONTINUE
C
CCCCC               IOFF=IOFF+NP1*NP2
CCCCC               JOFF=JOFF+NP2*NAO
C

c---------------------------------------------------------------------------
c   Move alpha's and pcoef's into their packed return arrays.
c   Save pointers in the ixalpha, ixpcoef arrays so the data may be
c   more easily referenced later.
c---------------------------------------------------------------------------

               ixalpha(ishell) = nxt_alpha
               do i = 1, np1
                  alpha(nxt_alpha+i-1) = salpha(i)
               enddo
       
               nxt_alpha = nxt_alpha + np1

               ixpcoef(ishell) = nxt_pcoef
               do j = 1, nao
               do i = 1, np1
                  pcoeff(nxt_pcoef) = spcoef((j-1)*np1+i)
                  nxt_pcoef        = nxt_pcoef + 1
               enddo
               enddo

 30         CONTINUE
 28      CONTINUE
C
         INI = INI + NPOP(IATM)
         IF(NSHL .GT. MAXANG) MAXANG = NSHL
 10   CONTINUE
C
 1120 FORMAT(2I5)
 1060 FORMAT(4F18.10)
C
      close(10)
      RETURN
      END


