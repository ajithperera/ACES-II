#include "flags.h"
#include "program.h"
      subroutine Read_File_Main(Icore, Icrsiz)
c
      implicit double precision (a-h,o-z)
      parameter (mxcoef=30)
c
#include "machsp.com"
#include "jodaflags.com"
#include "mxatms.par"
#include "baslims.par"
C
      character*32 szFile
      logical bExist,NORMCO_EXSIST,POLAR_EXSIST,SHIFT_EXSIST
      Logical Spncc_exsist 
      Character*4 Comp_pgrp, Full_pgrp
      Character*80 Fname 
      Dimension Nocc(16), Atommass(Mxatms), Iatmchrg(Mxatms),  
     &          Coord(3*Mxatms), Norbits_fullG(Mxatms), 
     &          NOrbits_compG(Mxatms),Nbsfns_4irrep(8)
      Dimension Icore(Icrsiz), Orient(3,3)
C
      Data Ione, Ieight /1, 8/
C
      Iuhf = 1
      If (iflags(h_IFLAGS_ref).eq.0) iuhf = 0 
C
      Maxcor   = Icrsiz
      Mxangmom = Nht
      Length   = 0
      Spncc_exsist = .False.
      Shift_exsist = .False. 
c
c Read the JOBARC file for basic data of the molecule. 
c
      INQUIRE(FILE='NORMCO',EXIST=NORMCO_EXSIST)
      INQUIRE(FILE='POLAR',EXIST=POLAR_EXSIST)
      Call GFNAME("CSHIFT",FNAME,ILENGTH)
      INQUIRE(FILE=FNAME(1:ILENGTH),EXIST=SHIFT_EXSIST)
      CALL B_GETREC(0,"JOBARC","JJTENSOR",Length,Ijunk)

      IF (Length .Gt. 0) Spncc_exsist = .True. 
  
      IF (NORMCO_EXSIST) THEN
      
         Call Getrec(-1, 'JOBARC', 'NUMDROPA', Length, Ijunk)       
         If (Length .GT. 0) Then 
             Print*, "Frozen-core is not allowed density plots"
             Call Errex
          Endif
c 
          Call B_GETREC(-1, 'JOBARC', 'NREALATM', Ione, Nreal_atoms)
          Call B_GETREC(-1, 'JOBARC', 'NATOMS  ', Ione, Natoms)
          Call B_GETREC(-1, 'JOBARC', 'FULLNORB', Ione, Iful_unq_atoms)
          Call B_GETREC(-1, 'JOBARC', 'COMPNORB', Ione, Icmp_unq_atoms)
          Call B_GETREC(-1, 'JOBARC', 'COORD   ', 3*Natoms*Iintfp, Coord)
          Call B_GETREC(-1, 'JOBARC', 'ATOMMASS', Natoms*Iintfp, Atommass)
          Call B_GETREC(-1, 'JOBARC', 'ATOMCHRG', Natoms, Iatmchrg)
          Call B_GETREC(-1, 'JOBARC', 'COMPNIRR', Ione, Nirrep)
          Call B_GETREC(-1, 'JOBARC', 'OCCUPYA ', Nirrep, Nocc(1))
          Call B_GETREC(-1, 'JOBARC', 'NBASTOT ', Ione, Nbfns)
          Call B_GETREC(-1, 'JOBARC', 'NAOBASFN', Ione, Naobfns)
          Call B_GETREC(-1, 'JOBARC', 'NUMBASIR', Nirrep, Nbsfns_4irrep)
          Call B_GETREC(-1, 'JOBARC', 'FULLPOPV', Iful_unq_atoms, 
     &                 Norbits_fullG)
          Call B_GETREC(-1, 'JOBARC', 'COMPPOPV', Icmp_unq_atoms, 
     &                 Norbits_compG)
C
          Nvibs  = (3*Nreal_atoms- 6) 
          Iloc   = 1
          Iscr1  = Iloc   + 9*Nvibs*Nvibs*IIntfp
          INext  = IScr1  + 9*Nvibs*Nvibs*IIntfp 
          Ilef   = Maxcor - Inext

          If (INext .ge. Maxcor) Call Insmem("Write_primt_main", Inext, 
     &                                        Maxcor)
C
          Call Read_normco_filE(Nreal_atoms, Nvibs, Icore(Iloc), 
     &                          ICore(IScr1))

      ELSE IF (POLAR_EXSIST) Then

          Call B_GETREC(-1, 'JOBARC', 'NREALATM', Ione, Nreal_atoms)
C
C The ORIENT2 record keep track of the transformation matrix from input
C coordiante to principle axis orientation of the molecule.
C
          CALL B_GETREC(20,'JOBARC','ORIENT2 ',9*IINTFP,Orient)
     
          Call Read_polar_file(Nreal_atoms,Orient)

      ELSE IF (SPNCC_EXSIST) THEN

          Call B_GETREC(-1, 'JOBARC', 'NREALATM', Ione, Nreal_atoms)
          Ndim  = 100
          Ndim2 = Ndim*Ndim

          Iloc   = 1
          Iscr1  = Iloc  + Nreal_atoms2
          Call B_GETREC(20,"JOBARC","JJTENSOR",Ndim2*Iintfp,
     &                  Icore(Iloc))
          Write(6,*) 
          Write(6,"(2a)") "     The total spin-spin coupling tensor",
     &                    " (in Hz)" 
    
CSSS          Call Neatprint(6,Icore(Iloc),Nreal_atoms,Nreal_atoms,100,100)
          Call Prntlo(6,Icore(Iloc),Nreal_atoms,Nreal_atoms,100,100)
          Write(6,*) 
      ELSEIF (SHIFT_EXSIST) THEN

          Icshift = 1
          Ishift  = Icshift + 9*Nreal_atmos 
          Ianiso  = Ishift  + Nreal_atoms 
          Irange  = Ianiso  + Nreal_atoms 
          Iend    = Irange  + Nreal_atoms 

          Call Process_c4shifts(Icore(Icshift),Icore(Ishift),
     &                          Icore(Ianiso),Icore(Irange),
     &                          Nreal_atoms)
      ENDIF 

      Return
      End

