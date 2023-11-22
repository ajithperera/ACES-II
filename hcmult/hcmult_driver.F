      Subroutine Hcmult_driver(Work,Maxcor,Iuhf)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Logical SS,SD,DS,DD
      Logical Cis,Rpa,Eomcc,Cisd,Fuldiag,Incore,Readgues 

#include "sympop.com"
#include "syminf.com"

      Common /Drvhbar/SS,SD,DS,DD
      Common /Meth/Cis,Rpa,Eomcc,Cisd,Fuldiag,Incore,Readgues 
      Data Ione /1/
      
C This simply a debugging tool and does a Hbar*C multiplication of
C C1 and C2 (arbitrary). Think of this as mimicing any itertation 
C other than the fist iteration EOM root search. All the EOM lists
C must be present and Hbar must be on the disk. Input would be 
C exactly the same as EOM-CCSD and instead of vee one can run
C this test module. 

      Ss       = .True.
      Sd       = .True.
      Ds       = .True.
      Dd       = .True.
      Cis      = .False.
      Rpa      = .False.
      Eomcc    = .True.
      Cisd     = .False.
      Fuldiag  = .False.
      Incore   = .False. 
      Readgues = .False. 

      Call Make_eom_lists(Work,Maxcor,Iuhf)
      Call Modf(Work,Maxcor,Iuhf,Ione)
      
      Call checkhbar(Work,Maxcor,Iuhf)
      Call check_t2(Work,Maxcor,Iuhf)

      Irrepx = Ione
      Iside  = Ione
      Call Make_newlst(Irrepx,Work,Maxcor,Iuhf)
     
C Get C2 and put it in the appropriate list (here I am using  
C ground state T2 but it can be any type of doubles. 

      Nsize_44 = Idsymsz(Irrepx,Isytyp(1,44),Isytyp(2,44))
      If (Iuhf .Ne. 0) Nsize_45 = Idsymsz(Irrepx,Isytyp(1,45),
     +                                           Isytyp(2,45))
      Nsize_46 = Idsymsz(Irrepx,Isytyp(1,46),Isytyp(2,46))

      Maxdim = Max(Nsize_44,Nsize_46,Nsize_45)
      If (Maxdim .Ge. Maxcor) Call Insmem("hcmult_driver",Maxdim,
     +                                     Maxcor)

      Call Getall(Work,Nsize_44,Irrepx,44)
      Call Dzero(Work,Nsize_44)
      Call Putall(Work,Nsize_44,Irrepx,444)
      If (Iuhf .Ne. 0) Then
         Call Getall(Work,Nsize_45,Irrepx,45)
         Call Dzero(Work,Nsize_45)
         Call Putall(Work,Nsize_45,Irrepx,445)
      Endif 
      Call Getall(Work,Nsize_46,Irrepx,46)
      Call Dzero(Work,Nsize_46)
      Call Putall(Work,Nsize_46,Irrepx,446)

C Get C2 and put it in the appropriate list (here I am using  
C ground state T2 but it can be any type of doubles. 

      Nsize_90_aa =Irpdpd(Irrepx,9)
      Call Getlst(Work,1,1,1,3,90)
      Call Dzero(Work,Nsize_90_aa)
      Work(2) = 1.0D0
      Call Putlst(Work,1,1,1,1,490)

      If (Iuhf .Ne. 0) Then
       Nsize_90_bb =Irpdpd(Irrepx,10)
       Call getlst(Work,1,1,1,4,90)
       Call Dzero(Work,Nsize_90_bb)
       Call Putlst(Work,1,1,1,2,490)
      Endif

      Call Hbardiag(Irrepx,Work,Maxcor,Iuhf) 

      Call Hc_mult(Work,Maxcor,Iuhf,Iside,Irrepx)
  
      Return
      End

      
       


