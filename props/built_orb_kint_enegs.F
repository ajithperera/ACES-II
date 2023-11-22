      Subroutine Built_orb_kint_enegs(Work,Maxcor,Nbfns,Naobfns,Nirrep,
     &                                Occ,Vrt,Nocc,Nvrt,Iuhf)

      Implicit Double Precision(A-H,O-Z)
     
      Dimension Work(Maxcor)
      Integer Occ(8,2),Vrt(8,2),Nocc(2),Nvrt(2)
      Character*8 Vecs(4) 
      Logical Symmetry 

      Data Vecs /"RELDENSA","RELDENSB","SCFEVECA","SCFEVECB"/
      Data One,Dnull,Two /1.0D0,0.0D0,2.0D0/

#include "machsp.com"
#include "flags.com"
#include "info.com"

      Nbfns2    = Nbfns*Nbfns 
      Naomobfns2= Naobfns*Nbfns 
      Naobfns2  = Naobfns*Naobfns 
      Norbs     = Max(Nocco(1)+Nvrto(1),Nocco(2)+Nvrto(2))
      I000      = 1
      Iorbdens_a = I000 
      Iorbdens_b = Iorbdens_a + Naobfns2*Norbs
      Istart     = Iorbdens_b + Naobfns2*Norbs

      Icorden_a = Istart 
      Iscfvec_a = Icorden_a + Nbfns2
      Imonatorb_a  = Iscfvec_a + Naomobfns2
      Iaonatorb_a  = Imonatorb_a + Nbfns2 
      Iend         = Iaonatorb_a + Naomobfns2
      If (Iend .GT. Maxcor) Call Insmem("built_orb_kint_enegs",
     &                                   Iend,Maxcor)
      If (Iuhf .Ne. 0) Then 
          Icorden_b = Iend
          Iscfvec_b = Icorden_b + Nbfns2
          Imonatorb_b  = Iscfvec_b + Naomobfns2 
          Iaonatorb_b  = Imonatorb_b + Nbfns2 
          Iend         = Iaonatorb_b + Naomobfns2
          If (Iend .GT. Maxcor) Call Insmem("built_orb_kint_enegs",
     &                                       Iend,Maxcor)
      Endif 

C Read in alpha and beta correlated densities from JOBARC. The RELDENSA
C and RELDENSB 

      Call Getrec(20,"JOBARC",Vecs(1),Nbfns2*Iintfp,Work(Icorden_a))
      Call Getrec(20,"JOBARC",Vecs(3),Naomobfns2*Iintfp,Work(Iscfvec_a))
      If (Iuhf .Ne. 0) Then
         Call Getrec(20,"JOBARC",Vecs(2),Nbfns2*Iintfp,Work(Icorden_b))
         Call Getrec(20,"JOBARC",Vecs(4),Naomobfns2*Iintfp,
     &               Work(Iscfvec_b))
      Endif 

      Call Eig(Work(Icorden_a),Work(Imonatorb_a),0,Nbfns,2)
      Call Dgemm("N","N",Naobfns,Nbfns,Nbfns,One,Work(Iscfvec_a),
     &            Naobfns,Work(Imonatorb_a),Nbfns,Dnull,
     &            Work(Iaonatorb_a),Naobfns)
      If (Iuhf .Ne. 0) Then 
          Call Eig(Work(Icorden_b),Work(Imonatorb_b),0,Nbfns,2)
          Call Dgemm("N","N",Naobfns,Nbfns,Nbfns,One,Work(Iscfvec_b),
     &                Naobfns,Work(Imonatorb_a),Nbfns,Dnull,
     &                Work(Iaonatorb_b),Naobfns)
      Endif 

#ifdef _DEBUG_LVLM
      Write(6,"(a)") " The alpha natural eigenvectors"
      call output(Work(Iaonatorb_a),1,Naobfns,1,Nbfns,Nbfns,
     &                Nbfns,1)
      Write(6,"(a)") " The alpha natural occupation numbers"
      call output(Work(Icorden_a),1,Naobfns,1,Nbfns,Nbfns,
     &                Nbfns,1)
      If (Iuhf .Ne. 0) Then
          Write(6,"(a)") " The beta natural eigenvectors"
          call output(Work(Iaonatorb_b),1,Naobfns,1,Nbfns,Nbfns,
     &                Nbfns,1)
      Write(6,"(a)") " The beta natural occupation numbers"
      call output(Work(Icorden_a),1,Naobfns,1,Nbfns,Nbfns,
     &                Nbfns,1)
      Endif
#endif 
      Call Cor_orbdens(Work(Iaonatorb_a),Work(Icorden_a),
     &                 Work(Iorbdens_a),Norbs,Nbfns,
     &                 Naobfns,1)
      If (Iuhf .Ne. 0) Then
         Call Cor_orbdens(Work(Iaonatorb_b),Work(Icorden_b),
     &                    Work(Iorbdens_b),Norbs,Nbfns,
     &                    Naobfns,2)
      Endif 

      If (Iuhf .Ne. 0) Then
         Call Daxpy(Norbs*Naobfns2,One,Work(Iorbdens_b),1,
     &              Work(Iorbdens_a),1)
         Call Daxpy(Norbs*Naobfns2,One,Work(Icorden_b),1,
     &              Work(Icorden_a),1)
      Endif 
 
C Calculate the orbital kinetic energies.

      Memleft = Maxcor - Iend
      Write(6,*)
      Write(6,"(2x,2a)") " The orbital kinetic energies from the ",
     &                   " correlated density matrix"
      Call Orb_Kint_enegs(Work(Iorbdens_a),Work(Icorden_a),
     &                    Work(Iend),Maxcor,Naobfns,Nbfns,
     &                    Norbs) 
     
      Return
      End 
