      Subroutine Vibavg_lvl0(Prop_derv1,Prop_derv2,Cfc,Vibfrq,Work,
     +                       Imap,Temp,Maxcor,Ndata,Nmodes,Ndim,
     +                       Nreals,B2ang,Au2invcm,NMR_spincc,
     +                       NMR_shifts,Iamlinear)
     
      Implicit Double Precision(A-H,O-Z)

      Dimension Prop_derv1(Ndim,Nmodes)
      Dimension Prop_derv2(Ndim,Nmodes)
      Dimension Cfc(Nmodes*Nmodes*Nmodes)
      Dimension Temp(Ndata)
      Dimension Work(Maxcor)
      Dimension Imap(Nreals)
      Dimension Vibfrq(3*Nreals)
      Logical NMR_spincc
      Logical NMR_shifts

      Data Ione /1/

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,a)")" @-vibavg_lvl0: The cubic force constant",
     +                " matrices (per mode)"
      call print_cfc(Cfc,Nmodes)
      Write(6,*) 
      If (NMR_spincc) Then
      Write(6,"(a)") " The fist-derivative of NMR spin-spin coupling"
      do i=1,Nmodes
      call output(Prop_derv1(1,i),1,Nreals,1,Nreals,Nreals,Nreals,1)
      enddo 
      Write(6,"(a)") " The second-derivative of NMR spin-spin coupling"
      do i=1,Nmodes
      call output(Prop_derv2(1,i),1,Nreals,1,Nreals,Nreals,Nreals,1)
      enddo 
      elseif (NMR_shifts) then
      Write(6,"(a)") " The fist-derivative of NMR chemical shifst"
      do j=1,Nmodes
      Write(6,"(6(1x,F15.6))") (Prop_derv1(i,j),i=1,Nreals)
      enddo 
      Write(6,"(a)") " The second-derivative of NMR chemical shifts"
      do j=1,Nmodes
      Write(6,"(6(1x,F15.6))") (Prop_derv2(i,j),i=1,Nreals)
      enddo
      endif 
      Write(6,*)
      Write(6,"(a)") " The harmonic frequencies"
      Write(6,"(4(1x,F15.6))") (Vibfrq(i),i=1,3*Nreals)
#endif

      Icol = 6
      if (Iamlinear .Gt. 0) Icol = 5

      Ndim1 = Nreals*Nreals*Ndata
      Ndim2 = Nreals*Ndata

      Call Trans_2nodim_coords(Prop_derv1,Prop_derv2,Vibfrq,Ndim,
     +                         Nreals,Nmodes,Icol,NMR_spincc,
     +                         NMR_shifts)

      I000 = Ione
      If (NMR_spincc) I010 = I000 + Ndim1
      If (NMR_shifts) I010 = I000 + Ndim2
      I020 = I010 + Ndata*Nmodes 
      Iend = I020 + Ndata*Nmodes

      If (Iend .Ge. Maxcor) Call Insmem("vibavg_lvl0",Iend,Maxcor)

      Call Dzero(Work(I010),Nmodes*Ndata)
      Call Get_cfc_contr(Cfc,Work(I010),Work(I020),Vibfrq,
     +                   Temp,Ndata,Nmodes,Nreals,Icol)

      If (NMR_spincc) Then
         Call Dzero(Work(I000),Ndim1)
         Call Vibavg_lvl0_spincc_1der(Prop_derv1,Work(I010),Vibfrq,
     +                                Nreals,Nmodes,Ndata,Icol,
     +                                Work(I000))
#ifdef _SEPARATE_12DER_CONTRS 
      If (NMR_spincc) Call Dzero(Work(I000),Ndim1)
#endif
         Call Vibavg_lvl0_spincc_2der(Prop_derv2,Work(I020),Nreals,
     +                               Nmodes,Ndata,Work(I000))
      Elseif (NMR_shifts) Then 
         Call Dzero(Work(I000),Ndim2)
         Call Vibavg_lvl0_shifts_1der(Prop_derv1,Work(I010),Vibfrq,
     +                                Nreals,Nmodes,Ndata,Icol,
     +                                Work(I000))
         Call Vibavg_lvl0_shifts_2der(Prop_derv2,Work(I020),Nreals,
     +                                Nmodes,Ndata,Work(I000))
      Endif 

      Return
      End 
