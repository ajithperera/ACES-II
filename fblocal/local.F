#include "flags.h"
      Program Local
C
C Description
C
C   Perform Foster-Boys and Edmiston-Ruedenberg
C   Localizations. Works with both UHF and RHF 
C   reference functions. Frozen core and virtual
C   localiztions are available in near future.
C   Current version was coded by Ajith 05/17/93
C
C Local Variables
C   Itemp  : Temporary storge area
C   Nocupy : Keep No. of occupied orbitals for RHF and UHF runs
C   Spntyp : Corresponds to ALPHA  and BETA spin cases
C   Nbasis : Number of basis functions
C   Nalpa  : No. of ALPHA occupied orbitals
C   
C Declarations
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Nocupy(2), Itemp(12)
      Character*5 Spntyp(2)
      Character*8 Dump(2), Pscf1(2), Pscf2(2)
      Logical Abort
      Parameter(Lirvpr=600, Oned10 = 1.0D-10, Oned8 = 1.0D-8, 
     &          MXNiter = 150, Mxtry = 3)
C
C Common Block informations
C
      Common//Icore(1)
      Common/Istart/I0
      Common/Machsp/Iintln,Ifltln,Iintfp,Ialone,Ibitwd
      Common/Files/Luout,Moints
      Common/Info/Nocco(2),Nvrto(2)
      Common/Flags/ Iflags(100)
      Common/Syminf/Nstart,Nirrep,Irreps(255,2),Dirprd(8,8)
C
      Data Spntyp /'ALPHA', 'BETA'/
      Data Dump /'LOCALMOA', 'LOCALMOB'/
      Data Pscf1 /'SCFEVCA0', 'SCFEVCBO'/
      Data Pscf2 /'SCFEVECA', 'SCFEVECB'/

C Initialize the Memory Allocation
C
      Call Crapsi(Icore, Iuhf, 0)
C
      Ione = 1
      Ithree = 3
      Nbasis = Nocco(1) + Nvrto(1)
C  
      Call Izero(Nocupy, 2)
      Call Getrec(20, 'JOBARC', 'SYMPOPOA', Nirrep, Itemp)
C
      Do 10 I = 1, Nirrep
         Nocupy(1) = Nocupy(1) + Itemp(I)
 10   Continue
      
      If (Iuhf .ne. 0) then 
         Call Getrec(20, 'JOBARC', 'SYMPOPOB', Nirrep, Itemp)
         Do 20 I = 1, Nirrep
            Nocupy(2) = Nocupy(2) + Itemp(I)
 20      Continue
      Endif
C
      Nalpa = Nocupy(1)
C
C Write out some important stuff to the user
C
      Write(Luout,*)
      Write(Luout, 5000)'FOSTER-BOYS ORBITAL LOCALIZATION'
      Write(Luout, 5000)'--------------------------------'
 5000 Format(25X, A)
      Write(Luout,*) 
      Write(Luout,"(a,i2)")'  Maximum number of localization retrys = ',
     &                        Mxtry
      Write(Luout,"(a,i3)")'  Maximum number of iteration per try   = ',
     &                         MXNiter
      Write(Luout,5010)'Threshold for neglect a rotation = ', Oned8
      Write(Luout,5010)'Threshold for convergence        = ', Oned8
      Write(Luout,*)
 5010 Format(2X, A, E7.1)
C      
C Allocate the memory
C      
      I000 = I0
      I010 = I000 + Ithree*Nbasis*Nbasis*Iintfp
      I020 = I010 + Ithree*Nbasis*Nbasis*Iintfp
      I030 = I020 + Nbasis*Nbasis*Iintfp
      I040 = I030 + Nbasis*Nbasis*Iintfp
      I050 = I040 + Nbasis*Nbasis*Iintfp
      I060 = I050 + Lirvpr*Iintfp
      I070 = I060 + Lirvpr
      I080 = I070 + Norbs
      I090 = I080 + Norbs
C
      If (Iuhf .ne. 0) then
         
         Do 30 I = 1, 2
            Write(Luout, 5020)'-UHF WAVE FUNCTION : SPIN TYPE',
     &         Spntyp(I),'-'
 5020       Format(22X, A, 1X, A, A)
            Write(Luout, *)
             Call FBlocal(Iuhf, Spntyp(I), Nalpa, Nocupy(I), Nbasis, 
     &                   Lirvpr, Abort, Icore(I000), Icore(I010),
     &                   Icore(I020), Icore(I030), Icore(I040), 
     &                   Icore(I050), Icore(I060), Icore(I070),
     &                   Icore(I080))
 30      continue

C Write localized MO's into JOBARC file
C
      Call Putrec(20, 'JOBARC', Dump(I), Nbasis*Nbasis*Iintfp,
     &            Icore(I040))

C Overwrite the "SCFEVC{A0,B0} records for post-SCF calculations.

      Call Putrec(20, 'JOBARC', Pscf1(I), Nbasis*Nbasis*Iintfp,
     &            Icore(I040))

      Call Putrec(20, 'JOBARC', Pscf2(I), Nbasis*Nbasis*Iintfp,
     &            Icore(I040))
C      
      Else
C
         Write(Luout, 5030)'-RHF WAVE FUNCTION-' 
         Write(Luout, *) 
 5030    Format(30X, A) 
         Call FBlocal(Iuhf, Spntyp(1), Nalpa, Nocupy(1), Nbasis, 
     &                   Lirvpr, Abort, Icore(I000), Icore(I010),
     &                   Icore(I020), Icore(I030), Icore(I040), 
     &                   Icore(I050), Icore(I060), Icore(I070),
     &                   Icore(I080))
C
C Write localized MO's into JOBARC file for RHF runs
C
      Call Putrec(20, 'JOBARC', Dump(1), Nbasis*Nbasis*Iintfp,
     &            Icore(I040))

C Overwrite the "SCFEVCAO record for post-SCF calculations.

       Call Putrec(20, 'JOBARC', Pscf1(1), Nbasis*Nbasis*Iintfp,
     &             Icore(I040))

       Call Putrec(20, 'JOBARC', Pscf2(1), Nbasis*Nbasis*Iintfp,
     &             Icore(I040))

      Endif

C Check the orthonormality of new vectors

      Call Getrec(20, 'JOBARC', "AOOVRLAP", Nbasis*Nbasis*Iintfp,
     &             Icore(I020))

      Call Xgemm("T", "N", Nbasis, Nbasis, Nbasis, 1.0D0, Icore(I040),
     &            Nbasis, Icore(I020), Nbasis, 0.0D0,  Icore(I030),
     &            Nbasis)
C
      Call Xgemm("N", "N", Nbasis, Nbasis, Nbasis, 1.0D0, Icore(I030),
     &            Nbasis, Icore(I040), Nbasis, 0.0D0,  Icore(I020),
     &            Nbasis)
      Write(6,"(2x,a)") "The overlap (C^tSC) check"
      Call output(Icore(I020), 1, Nbasis, 1, Nbasis, Nbasis, Nbasis, 1)

C Write out the successful completion message and also 
C the time taken
C 
      IFlags(h_IFLAGS_scf_maxcyc) = 0
CSSS      IFlags(h_IFLAGS_orbitals) = 1
CSSS      IFlags(h_IFLAGS_reference)  = 2
      Call Putrec(20, 'JOBARC', "IFLAGS  ", 100, Iflags)

      If ( Abort) Then
         Call Errex
         Call Aces_fin

      Else
         Write(Luout,*)
         Call Timer(Tcpu, Tjunk, 1)
         Write(Luout,5040) Tcpu
 5040    FORMAT(T3,'@LOCAL-I, Orbital localization completed in'
     &   ,F8.3,' cpu seconds.')
C
         Write(Luout, *)
         Call Aces_fin
      Endif
C
      Stop
      End





