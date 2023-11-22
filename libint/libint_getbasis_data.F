      SUBROUTINE LIBINT_GETBASIS_DATA(fname,MAX_ATOMS,max_shells,
     &                                natoms,nshells,nspc,cartesian,
     &                                Max_prims,Itfct,LNP1,lnpo,nfct,
     &                                Nufct,nbasis,NAOBASIS,nCFpS,
     &                                Npfps,NAOATM,angmom,atomic_label,
     &                                vnn,Maxang,Maxjco,Iecp,Nucstr,
     &                                NRCSTR,IQMSTR,JCOSTR,JSTRT,
     &                                JRS,NAMAT,CHARGE,CENTERS,SCR)
     

C   Simple version of Ajith Perera's INSPECT_MOL routine.  This version 
C   scans the MOL file to determine the number of atoms.  It also assumes
C   no reordering of centers and all atoms are symmetry-unique.
c
c   Mark Ponton 10/2/03
C
C     ----INPUT ARGUMENTS----
c FNAME       = Unix file name of the MOL file.
C MAX_ATOMS   = Maximum number of atoms allowed.
C MAX_SHELLS  = Maximum number of shells
C CARTESIAN   = True for Cartesian basis choice 
C 
C    ----OUTPUT ARGUMENTS----- 
C
C NATOMS    = The total number of atoms. 
c nshells   = total number of shells
c nspc      = Array containing the number of shells per center.
C ITFCT   = Total number of primitive functions.
C LNP1    = Largest possible value for number of contracted functions for
C           all shells.
C LNPO    = Largest possible value for the product of number of primitives
C           of contracted functions for all shells.
C NFCT    = The total number of  primitive functions on each atom.
C NAOATM  = The total number of contracted functions on each atom.
C NUFCT   = The total number of primitives on each sym. unique atom
C NAOUATM = The total number of contracted functions on each sym. 
C           unique atom.
C NBASIS  = Total number of basis functions (contracted)
C NAOBASIS = Total number of basis functions in Cartesian coordinates.
c nCFpS   = Number of contracted functions per shell.
c nPFpS   = Number of primitive functions per shell.
c angmom  = Anuglar momentum of each (sub)shell.
c atomic_label = Array of unique integer ids for the atom of each shell.
c vnn     = Nuclear-nuclear repulsion energy.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER*(*) FNAME
      CHARACTER*4 ATMNAM
      CHARACTER*4 NAMAT(MAX_ATOMS)
      CHARACTER*80 XLINE
      LOGICAL CARTESIAN,ECP,DBG
C     
      DIMENSION NFCT(MAX_ATOMS), NAOATM(MAX_ATOMS),
     &          NUFCT(MAX_ATOMS), NAOUATM(MAX_ATOMS),
     &          COORD(3),jcostr(MAX_ATOMS,Maxjco),
     &          nucstr(MAX_ATOMS,Maxang,Maxjco),
     &          nrcstr(MAX_ATOMS,Maxang,Maxjco),
     &          Jstrt(max_shells),Jrs(max_shells),
     &          iqmstr(MAX_ATOMS),Charge(Max_atoms),
     &          CENTERS(3,MAX_ATOMS),SCR(MAX_ATOMS)

      INTEGER NAOUATM2(MAX_ATOMS)
      INTEGER NAOBASIS, NAOTMP2, NP22
      INTEGER ISHL(6)
      integer idosph, idum(8)
      integer nspc(max_atoms), nCFpS(max_shells), nPFpS(max_shells)
      integer angmom(max_shells)
      integer atomic_label(max_shells)
      double precision vnn, r3, r5 
      double precision nuclear_nuclear_repulsion_energy 
      character*52  dumstring

      integer ihess, jhess, iatom, jatom 
C
C Open MOL file for basis set information. The MOL file is created by
C joda by processing user input file and basis set library.
C
      OPEN(UNIT=10, FILE=FNAME, FORM='FORMATTED', STATUS='OLD')
      REWIND(10)
      ECP = .FALSE. 
      DBG = .FALSE.
      If (Iecp .GE. 1)  ECP = .TRUE.
C     
C Read the first five lines. Information on those five lines 
C are not relevent in the presnt context. 
C
      read (10,'(a6,4x,9i5)') dumstring,(idum(i),i=1,8),icart
      if (icart .eq. 1) then
         cartesian = .true.
         if (dbg)
     *      print *,'*** MOL file contains Cartesian coordinates ***'
      else
         if (dbg)
     *      print *,'*** MOL file contains spherical coordinates ***'
         cartesian = .false.
      endif

      READ(10,'(A)') XLINE
      READ(10,'(A)') XLINE
      READ(10,'(A)') XLINE
      READ(10,'(A)') XLINE
C     
      LNP1   = 0
      LNPO   = 0
      LTNP1  = 0
      LTNPO  = 0
C
      natoms = 0
      nshells = 0
      indx_shell = 0
      DO 10 IATM = 1, MAX_ATOMS
         READ (10, '(A80)') XLINE
         if (xline(1:6) .eq. 'FINISH') go to 2000
C
         natoms = natoms + 1
         READ(xline, 1110) ZNUC, IJUNK, NSHL,
     *       (ISHL(I),I=1,NSHL)
 1110    FORMAT(F20.1,8I5)
         nspc(iatm) = 0

C ECP needs the charges as double precison.

         CHARGE(IATM) = ZNUC
C
         READ(10,1115) ATMNAM,(COORD(I),I=1,3)
 1115    FORMAT(A4,3F20.12)

c----------------------------------------------------------------------------
c   Save the geometry data 
c----------------------------------------------------------------------------

         CALL LIBINT_SET_GEOMETRY(CHARGE,CENTERS,IATM,COORD,ZNUC,
     &                             MAX_ATOMS)
C
         NUFCT(IATM) = 0
         NAOTMP      = 0
         NAOTMP2     = 0

C The ECP arrays added, Ajith Perera, 04/2012----

         IF (ECP) THEN
            IQMSTR(IATM) = NSHL
            NAMAT(IATM) = ATMNAM(1:2)
            DO IMOM = 1, NSHL
               JCOSTR(IATM, IMOM) = ISHL(IMOM)
            ENDDO
         ENDIF 
C----
C
         DO 20 I = 1, NSHL
            NPT  = 0
            NAOT = 0
            iangmom = i - 1
C
            DO 21 I1 = 1, ISHL(I)
               nshells        = nshells + 1
               nspc(iatm)     = nspc(iatm) + 1
               angmom(nshells) = iangmom   ! same a. m. value for each subshell
               atomic_label(nshells) = iatm

               READ(10,1120) NP1, NAO
 1120          FORMAT(2I5)

               IF (NAO .GT. MAX_PRIMS) THEN
               WRITE(6,"(3a,I4)")" The number of contracted functions",
     &                          " per shell exceeded the maximum",
     &                          " allowed", MAX_PRIMS 
               CALL ERREX
               ENDIF
C
               NPT  = NPT  + NP1
               NAOT = NAOT + NAO

C The ECP arrays added, Ajith Perera, 04/2012----

               IF (ECP) THEN
                  NUCSTR(IATM,I,I1) = NP1
                  NRCSTR(IATM,I,I1) = NAO
               ENDIF
C-------
C
               IF (CARTESIAN) THEN 
                  NP2 = I*(I + 1)/2 
                  NP22 = NP2
               ELSE
                  NP2 = 2*I - 1
                  NP22 = I*(I + 1)/2
               ENDIF
C
               NAOTMP = NAOTMP + NP2*NAO
               NAOTMP2 = NAOTMP2 + NP22 * NAO
               NUFCT(IATM) = NUFCT(IATM) + NP2*NP1
               nPFpS(nshells) = np1
               nCFpS(nshells) = nao
C
               NLN = (NAO-3)/4
               IF ((NAO-3) .GT. (NLN*4)) NLN = NLN + 1
               NLN = (NLN + 1)*NP1
C
               DO 30 J=1,NP1
                  READ(10,*) A, (SCR(K),K=1,NAO) 
 30            CONTINUE
C
               IF(NPT .GT. LNP1) THEN
                  IF(NPT .GT. LTNP1) LTNP1 = NPT
               ENDIF
C
               ITMP = NPT*NAOT
               IF(ITMP .GT. LNPO) THEN
                  IF(ITMP .GT. LTNPO) LTNPO = ITMP
               ENDIF
C
 21         CONTINUE
 20      CONTINUE
C
C
         IF (LTNP1 .GT. LNP1) LNP1 = LTNP1
         IF (LTNPO .GT. LNPO) LNPO = LTNPO
C
         NAOUATM(IATM) = NAOTMP
         NAOUATM2(IATM) = NAOTMP2
C
 10   CONTINUE
C     
 2000 continue
      iuatms = natoms
      ITFCT = 0
      DO 110 IATM = 1, IUATMS
            ITFCT = ITFCT + NUFCT(IATM)
 110  CONTINUE
C     
C Fill out NFCT, NAOATM and NMOMFCT for all atoms.
C
      ICNT   = 0
      NBASIS = 0
      NAOBASIS = 0
      DO 1011 II = 1, IUATMS
C
C
            ICNT = ICNT + 1
            NFCT(ICNT) = NUFCT(II)
            NAOATM(ICNT)= NAOUATM(II)
            NBASIS = NBASIS + NAOUATM(II)
            NAOBASIS = NAOBASIS + NAOUATM2(II)
C
 1011 CONTINUE
C
C The ECP arrays added, Ajith Perera, 04/2012----

      JSTRT(1) = 0
      JRS(1)   = 0

      IF (ECP) THEN
         DO ISHELL = 1, NSHELLS-1
            JSTRT(ISHELL+1) = JSTRT(ISHELL) + NPFPS(ISHELL)
            JRS(ISHELL+1)   = JRS(ISHELL) +   NPFPS(ISHELL)*
     &                        NCFPS(ISHELL)
         ENDDO 
      ENDIF
C-------
C
C Compute the Nuclear-Nuclear repulsion energy
C --------------------------------------------
C
CSSS      vnn = nuclear_nuclear_repulsion_energy(natoms)

      close(10)
      RETURN
      END


