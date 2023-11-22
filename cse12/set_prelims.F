#include "flags.h"
      Subroutine Set_prelims(Work,Maxcor,Iuhf,Nbasis,Nocc_a,Nocc_b,
     +                       Nvrt_a,Nvrt_b,Frac_occ)

      Implicit Double Precision (A-H,O-Z)
       
      Dimension Work(Maxcor)

      Dimension Nocc_alpha(8),Nocc_beta(8)
      Logical Frac_occ,Occnums_file,Save_ints 
      Character*80 Fname 

#include "machsp.com"
#include "jodaflags.com"
#include "cse12_vars.com"

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(2a,4(1x,I3))") "@-set_prelims, varibales:",
     +                           " nocc_a,nocc_b,nvrt_a,nvrt_b",
     +                             NOcc_a,Nocc_b,Nvrt_a,Nvrt_b
      Write(6,*)
#endif 
      Icontl     = IFLAGS(h_IFLAGS_cc_conv)
      cc_conv    = 10.D0**(-Icontl) 
      cc_maxcyc  = IFLAGS(h_IFLAGS_cc_maxcyc)
      Ring_cc    = IFLAGS(h_IFLAGS_excite)  .EQ. 2 
      Icontl     = IFLAGS(h_IFLAGS_lindep_tol)
      Denom_tol  = 10.D0**(-Icontl) 
      Save_ints  = IFLAGS(h_IFLAGS_save_ints) .GT. 0
      Frac_occ   = .False. 
      Brueck     = IFLAGS(h_IFLAGS_brueckner) .GT. 0
      Icontl     = IFLAGS(h_IFLAGS_bruck_conv)
      Brueck_tol = 10.D0**(-Icontl)

      If (.Not. Save_ints) Then
         Write(6,"(a,a)") " This program requires AO integrals. Since",
     +                    " this program is executed"
         Write(6,"(a,a)") " post-transform, use the SAVE_INTS=ON",
     +                    " key-word to facilitate corret"
         Write(6,"(a)")   " execution."
         Write(6,*)
         Call Errex 
      Endif 
      If (cc_maxcyc .Gt. 500) Then
         Write(6,"(2a)") " The total number of iteration must be",
     +                   " less than or equal to 500"
         Call Errex
      Endif
      Call Gfname("OCCNUMS",Fname,Ilength)
      Inquire(File=Fname(1:7),Exist=Occnums_file)

      If (Occnums_file) Call Read_occnum(Work,Maxcor,Act_min_a,
     +                                   Act_max_a,Act_min_b,
     +                                   Act_max_b,Nbasis,Iuhf)

      Call Getrec(0,"JOBARC","ORB_OCCA",Length,Junk)

      If (Length .Ne. Nbasis) Then
         Write(6,"(2a)") " Inconsistency in the values for number",
     +                  " of basis function and the record length" 
         Write(6,"(a)") " of the doubele precision occupation vector"
         Write(6,*)
         Call Errex 
      Endif 

      Iocn_a = 1
      Iocn_b = Iocn_a + Nbasis
      Iend   = Iocn_b + Nbasis

      If (Occnums_file) Then
          Call Getrec(20,"JOBARC","ORB_OCCA",Nbasis,Work(Iocn_a))
          If (Iuhf .Ne. 0) Call Getrec(20,"JOBARC","ORB_OCCB",
     +                                 Nbasis,Work(Iocn_b))
          Nfrac_a = 0
          Nfrac_b = 0
          Do I =1, Nocc_a
             If (Dabs(Work(Iocn_a+I-1)) .Ne. 1.0D0) Then
                 Nfrac_a = Nfrac_a + 1
             Endif 
          Enddo 
          Do I =1, Nocc_b
             If (Dabs(Work(Iocn_b+I-1)) .Ne. 1.0D0) Then
                 Nfrac_b = Nfrac_b + 1
             Endif 
          Enddo

          If (Nfrac_a .NE. 0 .OR. Nfrac_b .NE. 0) Then
             Frac_occ = .True.
             Write(6,"(2a)") " Fractional occupation CCSD calculation",
     +                       " is performed. The fractionnal occupation" 
             Write(6,"(a)")  " numbers are read from the OCCNUMS file." 
          Else

C Consider the following scenario: There are no fractional occupation
C numbers present but the user has specified occupation numbers using 
C the OCCNUMS file. This is fractional orbital calculation with the
C canonical occupations. 

             Frac_occ = .True.
             Write(6,"(3a)") " Fractional occupation CCSD calculation",
     +                       " with the canonical occupation numbers"
             Write(6,"(2a)")  " (read from the OCCNUMS file) is",
     +                        " performed."
          Endif 
             
          Do I = 1, Nocc_a
             Ocn_oa(i) = Work(Iocn_a+I-1)
          Enddo 
          Do I = 1, Nvrt_a
             Ocn_va(i) = Work(Nocc_a+Iocn_a+I-1)
          Enddo 
          If (Iuhf .GT. 0) Then
             Do I = 1, Nocc_b
                Ocn_ob(i) = Work(Iocn_b+I-1)
             Enddo 
             Do I = 1, Nvrt_b
                Ocn_vb(i) = Work(Iocn_b+Nocc_b+I-1)
             Enddo
          Endif 

C If the fractional occuptation calculation is requested, then the
C occupied and virtual occupation distributions are identical. The
C following logic assuares that.
  
          If (Frac_occ) Then 
             Do I = 1, Nvrt_a
                Ocn_oa(Nocc_a+I) = Ocn_va(I)
             Enddo 
             Do I = 1, Nvrt_b
                Ocn_ob(Nocc_b+I) = Ocn_vb(I)
             Enddo 
             Call Dcopy(Nbasis,Ocn_oa,1,Ocn_va,1)
             Call Dcopy(Nbasis,Ocn_ob,1,Ocn_vb,1)
          Endif 
           
      Else 
          Call Getrec(0,"JOBARC","OCCUPYA ",Length,Junk)
          Call Getrec(20,"JOBARC","OCCUPYA ",Length,Nocc_alpha(1))
          If (Iuhf .Ne. 0) Call Getrec(20,"JOBARC","OCCUPYB ",
     +                                 Length,Nocc_beta(1))

         If (Nocc_alpha(2) .Gt. 0 .OR. Nocc_beta(2) .Gt. 0) Then
             Write(6,"(2a)") " Abelian symmetry is not allowed.",
     +                      " use SYMMETRY=OFF instead."
             Call Errex
          Endif

          Call Dzero(Ocn_va,Nvrt_a)
          Do I = 1, Nocc_alpha(1)
             Ocn_oa(I) = 1.0D0
          Enddo 
          If (Iuhf .Gt. 0) Then
             Call Dzero(Ocn_vb,Nvrt_b)
             Do I = 1, Nocc_beta(1)
                Ocn_ob(I) = 1.0D0
             Enddo 
          Endif 

          Write(6,"(2a)") " Canonical CCSD calculation with the",
     +                    " canonical occupation numbers is performed."
      Endif 

      Active_space = (Act_min_a .Ne. 0 .Or. Act_max_a .Ne. 0 .Or.
     +                Act_min_b .Ne. 0 .Or. Act_max_b .Ne. 0)

      If (Frac_occ) Then
      Write(6,*)
      Write(6,"(a)") " The Alpha and Beta occupied occupation vectors"
      Write(6,"(6(1x,F10.5))") (Ocn_oa(i), i=1,Nbasis)
      Write(6,"(6(1x,F10.5))") (Ocn_ob(i), i=1,Nbasis)
      Write(6,"(a)") " The Alpha and Beta virtual occupation vectors"
      Write(6,"(6(1x,F10.5))") (Ocn_va(i), i=1,Nbasis)
      Write(6,"(6(1x,F10.5))") (Ocn_vb(i), i=1,Nbasis)
      Else
      Write(6,*)
      Write(6,"(a)") " The Alpha and Beta occupied occupation vectors"
      Write(6,"(6(1x,F10.5))") (Ocn_oa(i), i=1,Nocc_a)
      Write(6,"(6(1x,F10.5))") (Ocn_ob(i), i=1,Nocc_b)
      Write(6,"(a)") " The Alpha and Beta virtual occupation vectors"
      Write(6,"(6(1x,F10.5))") (Ocn_va(i), i=1,Nvrt_a)
      Write(6,"(6(1x,F10.5))") (Ocn_vb(i), i=1,Nvrt_b)
      Endif 

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(2a,ES8.2E2,1x,I4,1x,l)") " @-set_prelims, varibales:",
     +                              " cc_conv, cc_maxcyc,Ring_cc= ", 
     +                                cc_conv,cc_maxcyc,Ring_cc
      Write(6,"(a,ES8.2E2)") " Denom_tol= ", Denom_tol
      Write(6,"(2a,1x,4I3,l)")" @-set_prelims, the active space info:",
     +                       " act_min_a-b,act_max_a-b,active_space= "
      Write(6,"(4i3,l)")      act_min_a,act_max_a,act_min_b,act_max_b,
     +                        Active_space
      Write(6,*)
#endif 
      Return
      End 

