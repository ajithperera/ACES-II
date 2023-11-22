      Subroutine Finalize(Work,Iatmchrg,Isotopes,Atmlabel,Temp,Maxcor,
     +                    Nreals,Natoms,Ndata,NMR_spincc,NMR_shifts)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Dimension Iatmchrg(Nreals)
      Dimension Isotopes(Nreals)
      Dimension Temp(Ndata)
      Character*80 Fname 
      Character*6 Atmlabel(Nreals)
      Logical NMR_spincc 
      Logical NMR_shifts
      Logical eJJ_SPNCC_present
      Logical eNMR_SHIFTS_present
      Data Ione /1/

      Write(6,*)
      Write(6,"(6x,2a)") "-------Finalizing the vibrational",
     +                   " corrections-------"
#ifdef _DEBUG_LVL0
      Print*, Temp(1)
#endif
C Read the equilibrium tensors for 

      Iunit = 99 
      If (NMR_spincc) Then
         Call Gfname("eJJ_SPNCC",Fname,Length)
         Inquire(File=Fname(1:Length),Exist=eJJ_SPNCC_present)
         If (eJJ_SPNCC_present) Then
             Open(Unit=Iunit,File=Fname(1:Length),Form="Formatted")
         Else
             Write(6,"(2a)") " The equilibrium tensor stored in",
     +                       " eJJ_SPNCC is needed to proceed"
             Call Errex
         Endif 

         Ndim = Nreals*Nreals 
         I000 = Ione
         I010 = I000 + Nreals*Nreals*Ndata
         Iend = I010 + Nreals*Nreals
         If (Iend .Ge. Maxcor) Call Insmem("finalize",Iend,Maxcor)

         Call Dzero(Work(I010),Ndim)
         Call Get_prop_values(Work(I010),NMR_spincc,NMR_shits,Nreals,
     +                             Iunit)
         Call Get_sum_spincc(Work(I000),Work(I010),Nreals,Ndata)
          
         Ioff = 0
         Do K = 1, Ndata 
            Call Print_spncc_tensor(Work(I000+Ioff),Iatmchrg,Isotopes,
     +                              Atmlabel,Nreals,Natoms) 
            Ioff = Ioff + Nreals*Nreals 
         Enddo 
         Close(Iunit)
      Elseif(NMR_shifts) Then

         Call Gfname("eNMR_SHIFT",Fname,Length)
         Inquire(File=Fname(1:Length),Exist=eNMR_SHIFTS_present)
         If (eNMR_SHIFTS_present) Then
             Open(Unit=Iunit,File=Fname(1:Length),Form="Formatted")
         Else
             Write(6,"(2a)") " The equilibrium tensor stored in",
     +                       " eNMR_SHIFT is needed to proceed"
             Call Errex
         Endif

         I000 = Ione
         I010 = I000 + Nreals*Ndata
         Iend = I010 + Nreals
         If (Iend .Ge. Maxcor) Call Insmem("finalize",Iend,Maxcor)

         Call Get_prop_values(Work(I010),.False.,NMR_shifts,Nreals,
     +                            Iunit)
         Call Get_sum_shifts(Work(I000),Work(I010),Nreals,Ndata)
         Ioff = 0
         Do K = 1, Ndata 
            Call Print_shifts_tensor(Work(I000+Ioff),Iatmchrg,Isotopes,
     +                               Atmlabel,Nreals) 
            Ioff = Ioff + Nreals 
         Enddo 
         Close(Iunit)

       Endif      

       Return
       End 
