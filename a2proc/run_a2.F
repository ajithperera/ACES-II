      Subroutine Run_a2(Deltaq,Coords,Btmp,Norm_coords,Iatmchrg,AtmMass,
     +                  Nreals,Natoms,Ipoints,Imode,Prop,Plusminus,
     +                  Atmlabel,Label,Imap)

      Implicit Double Precision (A-H,O-z)
C
#include "machsp.com"
#include "jodaflags.com"
#include "const.h"
C
      Dimension Deltaq(3*Nreals), Coords(3*Nreals)
      Dimension Btmp(3*Nreals)
      Dimension Nocc_a(8)
      Dimension Pop_a(8),Pop_b(8)
      Double Precision Norm_coords(3*Nreals)
      Double Precision Atmmass(Nreals)
      Integer Iatmchrg(Nreals)
      Integer Imap(Natoms)

      Character*10 Prop
      Character*5 Plusminus
      Character*5 Atmlabel(Natoms)
      Character*5 Label(Natoms)

      Data Ione /1/

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "@-run_a2: The current normal mode"
      Write(6, "(3F17.13)") (Norm_Coords(i),i=1,3*nreals)
#endif 
     
      Dscale = Deltaq(Imode)*Ipoints 
      Call Dcopy(3*Nreals,Btmp,1,Coords,1)
      Call Daxpy(3*Nreals,Dscale,Norm_coords,1,Coords,1)

      Call Putrec(20, "JOBARC", "COORD  ", 3*Nreals*IINTFP,
     &            Coords)

      Write(6,"(a)") "The coordinates of the current point (in Bohr)"
      Write(6, "(3F17.13)") (Coords(i),i=1,3*Nreals)
      Call Get_nucrep(Coords,Iatmchrg,Nreals,Repulsion)
      Write(6,"(2a,F16.10)") " The nuclear repulsion",
     &                       " energy: ", Repulsion
#ifdef _NOSKIP
      Do Iatm = 1, Natoms 
         If (Imap(Iatm) .Ne. 0) Then
            Label(Iatm) = Atmlabel(Imap(Iatm))
         Endif 
      Enddo 
#else
      Natom = 1
      Do Iatm = 1, Natoms 
         If (Atmlabel(Iatm) .NE. "X") Then
            Label(Natom) = Atmlabel(Iatm)
            Natom = Natom + 1
         Endif 
      Enddo 
#endif 

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "@-run_a2: The atom labels" 
      Write(6,"(6(1x,5a))") (Label(i),i=1,Nreals)
      Write(6,"(a,2(1x,I2))") "@-run_a2: Imode,Ipoints: ",
     +                         Imode,Ipoints 
      Write(6,"(a)") "@-run_a2: The coordinates of the current point"
      Write(6, "(3F17.13)") (Coords(i),i=1,3*Nreals)
      Write(6, "(a)") Plusminus 
#endif 

      Call Dscal(3*Nreals,CONV_a0_a,Coords,1)
      Call Buildzmat(Plusminus,Coords,Label,Nreals,Imode,Ipoints)

      Return

C Down below works fine if there is no symmetry. This work entirly
C from one JOBARC. I find it hard to manage this when there is 
C symmetry. One day I might be able to pull this off but at the
C moment I will prepare ACES II inputs and retun. Ajith Perera
C Oct. 31/2020. 

      If (Prop .EQ. "NMR_SPN_CC") Then 
          Iflags(18)    = 13 
          Iflags(2)     = 10
          Iflags(3)     = 1
          Iflags(4)     = 7
          Iflags(5)     = 7
          Iflags(11)    = 2
          Iflags(12)    = 5
          Iflags(21)    = 1
          Iflags(43)    = 1
          Iflags(97)    = 50
          Iflags(98)    = 5
          Iflags2(106)  = 2
          Iflags2(109)  = 1
          Iflags2(138)  = 1
          Iflags2(171)  = 1
          Iflags2(105)  = 50
          Iflags2(107)  = 1
          Iflags2(122)  = 2
          Iflags2(123)  = 2
          If (Iflags(54) .NE. 0) Iflags(54)   = 0
          Call aces_ja_truncate('JODAOUT ', 1)
          Call Getrec(20,"JOBARC","NIRREP  ",Ione,Nirrep)
          Call Getrec(20,"JOBARC","OCCUPYA0",Nirrep,Nocc_a)
          Call Putrec(20,"JOBARC","OCCUPYB0",Nirrep,Nocc_a)
          Call Putrec(20,"JOBARC","OCCUPYB ",Nirrep,Nocc_a)
          Call Getrec(20,"JOBARC","NUMDROPA",Ione,Numdropa)
         Call Putrec(20,"JOBARC","NUMDROPB",Ione,Numdropa)

      Elseif (Prop .EQ. "SP_ENERGY ") Then

          If (Iflags(54) .NE. 0) Iflags(54)   = 0
          Iflags2(138) = 1
          Iflags2(171) = 1
          If (Iflags(3) .EQ. 2)
     &                       Iflags(3) = 1

      Elseif (Prop .EQ. "NMR_SHIFTS") Then

C This needs to work with CFOUR to do NMR shifts with sphericals
C and beyond MBPT(2). Also, here some work need to be done since the
C final printing is done in joda. It is difficult to set up all the
C flags correctly so that joda can go to the right place. Since the

          If (Iflags(54) .NE. 0) Iflags(54)   = 0
          Iflags(18)    = 3
          Iflags(2)     = 1
          Iflags(3)     = 2
          Iflags(62)    = 0
          Iflags2(138)  = 1
      Endif 

      Call Putrec(1,'JOBARC','IFLAGS  ',   100, iflags)
      Call Putrec(1,'JOBARC','IFLAGS2 ',   500, iflags2)
 
      Call aces_ja_fin
C
#ifdef _DEBUG_LVL0
      Write(6,"(a)") " @-Run_a2: call to runaces2b"
      Write(6,"(a)") " The coordinates of the current point"
      Write(6, "(3F17.13)") (Coords(i),i=1,3*Nreals)
#endif
       Call Runit("runaces2b")
C
       Call aces_ja_init
C
      Call Getrec(20,'JOBARC','TOTENERG', Ione, E_current)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,F15.10)")  " @-Run_a2: Energy: ",E_current 
#endif

      Return
      End

