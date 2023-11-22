      Subroutine Tdee_form_restart_lists(Work,Memleft,Iuhf,Irrepx,Nsize)

      Implicit Double Precision (A-H,O-Z)
      Dimension Work(Memleft)

C We need place holder on the disk to store the current Vector. The
C Vector is 2*Nsize long (both real and imaginary components). It
C actually consitute the Dipole function Mu(1) AA (BB for RHF) and
C Mu(2) (for UHF AAAA BBBB) ABAB ordered as AA-BB-AAAA-BBBB-ABAB.
C Lets not try to split this into pices but write as a single 
C Coulmn Vector.

      Iright = 1
      Ileft  = 2*Nsize

      Write(6,*) Iright,Ileft, Irrepx, Nsize 
      Call Updmoi(Iright,Ileft,Irrepx,495,0,0)
      Call Updmoi(Iright,Ileft,Irrepx,496,0,0)

      Return
      End
