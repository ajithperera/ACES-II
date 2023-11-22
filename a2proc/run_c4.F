      Subroutine Run_c4(Deltaq,Coords,Norm_coords,Nreals,Ipoints,
     +                  Imode,Prop)

      Implicit Double Precision (A-H,O-z)
C
#include "machsp.com"
#include "jodaflags.com"
C
      Dimension Deltaq(3*Nreals), Coords(3*Nreals)
      Dimension Nocc_a(8)
      Double Precision Norm_coords(3*Nreals)
      Character*10 Prop

      Data Ione /1/
     
      Dscale = Deltaq(Imode)
      Call Daxpy(3*Nreals,Dscale,Norm_coords,1,Coords,1)

      Call Putrec_c4(20, "JOBARC", "COORD  ", Nreals3*IINTFP,
     &               Coords)

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
          Call Getrec_c4(20,"JOBARC","NIRREP  ",Ione,Nirrep)
          Call Getrec_c4(20,"JOBARC","OCCUPYA0",Nirrep,Nocc_a)
          Call Putrec_c4(20,"JOBARC","OCCUPYB0",Nirrep,Nocc_a)
          Call Putrec_c4(20,"JOBARC","OCCUPYB ",Nirrep,Nocc_a)
          Call Getrec_c4(20,"JOBARC","NUMDROPA",Ione,Numdropa)
          Call Putrec_c4(20,"JOBARC","NUMDROPB",Ione,Numdropa)

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

      Call Putrec_c4(1,'JOBARC','IFLAGS  ',   100, iflags)
      Call Putrec_c4(1,'JOBARC','IFLAGS2 ',   500, iflags2)
C 
      Call Dumpja_c4("O")
C
#ifdef _DEBUG_LVL0
      Write(6,"(a)") " @-Run_a2: call to runaces2b"
      Write(6,"(a)") " The coordinates of the current point"
      Write(6, "(3F17.13)") (Coords(i),i=1,3*Nreals)
#endif
       Call Runit("runcfourb")
C
       Call Dumpja_c4("O")
C
      Call Getrec_c4(20,'JOBARC','TOTENERG', Ione, E_current)
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,F15.10)")  " @-Run_a2: Energy: ",E_current 
      Stop
#endif

      Return
      End

