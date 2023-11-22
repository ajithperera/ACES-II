      Subroutine Make_lccsd_hbar(Icore, Maxcor, Iuhf)

      Implicit double precision (A-H,O-Z)

      Dimension Icore(Maxcor)

      Logical Term1,Term2,Term3,Term4,Term5,Term6

      I0 = 1
 
      Call Stllst(Icore(I0),Maxcor,Iuhf)

      Term1 = .True.
      Term2 = .True.
      Term3 = .True.
      Term4 = .True.
      Term5 = .True.
      Term6 = .True.
C 
C Form Hbariakjk
C
      Call Formw5(Icore(I0),Maxor,Iuhf,Term1,Term2,Term3,Term4,
     &            Term5,Term6,.True.,.True.)
C
C Form Hbarabci
C
      Call Formw5(Icore(I0),Maxor,Iuhf,Term1,Term2,Term3,Term4,
     &            Term5,Term6,.True.,.True.)
   
C 
C  Form Hbarijka (no difference between CCSD and LCCSD Hbar(ijka)
C
      Call Hbrijka0(Icore(I0),Maxcor,Iuhf)
C
C Form Hbaraibc (no difference between CCSD and LCCSD Hbar(aibc)
C

      Call Modaibc(Icore,Maxcor,Iuhf)
C
C Now we can call genint from CC code to form the one-paticle 
C intermediates, Hbarijkl and W(mb,ej). Then call the routine 
C apply fixes so at the end we have Hbar elements. 
 
      Call genint(Icore(i0),Maxcor,Iuhf,2)
C
C This will fix  Hbar(ij) and Hbar(ea)

      Call Fixfbar(Icore(I0),Maxcor,Iuhf)
C 
C  This will fix the W(mb,ej) so that it would become Hbar(mb,ej)
C
      
      Call Formwl(Icore(I0),Maxcor,Iuhf)
      Call Wtwtw(Icore(i0),Maxcor,Iuhf)
   
      Write(6,"(a)") "The Hbar calcualation for LCCSD is completed"
      Call aces_fin
      Stop

      Return
      End
