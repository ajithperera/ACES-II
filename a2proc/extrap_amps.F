      Subroutine Extrap_amps(Work,Maxcor,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Maxcor)

      Parameter(Max_xp=100)
      Dimension Eta_val(Max_xp)

      Data Ione,Inull,I30,Dnull,Two,Twl /1,0,30,0.0D0,2.0D0,12.0D0/

      N_xp            = I30
      Eta_val(1)      = Dnull
      Eta_val(2)      = Two
      Eta_val(n_xp+2) = Twl

      Do I = 3, N_xp+1
         Eta_val(I) = Eta_val(I-1) + 0.125D0
      Enddo

      Call aces_ja_fin

      Do I_xp = 1, N_xp+1
CSSS      Do I_xp = 1, 1
         Eta_o = - 10.0d0**(-Eta_val(I_xp))
         Call aces_ja_init
         Call Putrec(20,"JOBARC","ETA_4EXP",Ione,Eta_o)
         Call Putrec(20,"JOBARC","ETA_ITER",Ione,I_xp)
         Call aces_ja_fin

         Call Runit("xlsccresp")
CSSS         Call Runit("xhbar") 
CSSS         Call Runit("xvcceh")

      Enddo 
      
      Return
      End 
