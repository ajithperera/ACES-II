      Subroutine Osc_bd_fin(Osc_qm,Osc_om,Osc_am,Osc_dm,Osc_mq,
     +                      Qm_2,Om_2,Am_2,Dm_2,Mq_2,Jroot_1,
     +                      Irrepy_1,Iroot_1,Irrepx_1,E,Length)
    
      Implicit Double Precision(A-H,O-Z)
#include "machsp.com"

      Double precision Mq_2,Mq2_t,M_iso
    
      Dimension Qm_2(4,8,4,6,4,8,4,6,9)
      Dimension Om_2(4,8,4,6,4,8,4,6,10)
      Dimension Am_2(4,8,4,6,4,8,4,6,3)
      Dimension Mq_2(4,8,4,6,4,8,4,6,27)
      Dimension Dm_2(4,8,4,6,4,8,4,6,3)

      Dimension Qm2_t(9)
      Dimension Om2_t(9)
      Dimension Am2_t(3)
      Dimension Mq2_t(27)
      Dimension Dm2_t(3)
    
      Dimension Osc_Qm(1024)
      Dimension Osc_Om(1024)
      Dimension Osc_Am(1024)
      Dimension Osc_Mq(1024)
      Dimension Osc_Dm(1024)
      
      Qm2_t(1) = Qm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,1) 
      Qm2_t(2) = Qm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,2) 
      Qm2_t(3) = Qm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +    Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,3) 
      Qm2_t(4) = Qm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,4) 
      Qm2_t(5) = Qm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,5) 
      Qm2_t(6) = Qm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,6) 
      Qm2_t(7) = Qm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,7) 
      Qm2_t(8) = Qm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,8) 
      Qm2_t(9) = Qm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,9) 

      Om2_t(1) = Om_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,1) 
      Om2_t(2) = Om_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,2)
      Om2_t(3) = Om_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,3) 
      Om2_t(4) = Om_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,4)
      Om2_t(5) = Om_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,5) 
      Om2_t(6) = Om_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,6) 
      Om2_t(7) = Om_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,7) 
      Om2_t(8) = Om_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,8) 
      Om2_t(9) = Om_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,9)

      Am2_t(1) = Am_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,1) 
      Am2_t(2) = Am_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,2) 
      Am2_t(3) = Am_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,3) 

      Dm2_t(1) = Dm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,1) 
      Dm2_t(2) = Dm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,2) 
      Dm2_t(3) = Dm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,3) 

      Mq2_t(1) =  Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,1) 
      Mq2_t(2) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,2)
      Mq2_t(3) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,3) 
      Mq2_t(4) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,4) 
      Mq2_t(5) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,5) 
      Mq2_t(6) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,6) 
      Mq2_t(7) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,7) 
      Mq2_t(8) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,8) 
      Mq2_t(9) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,9) 
      Mq2_t(10) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,10) 
      Mq2_t(11) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,11)
      Mq2_t(12) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,12) 
      Mq2_t(13) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,13) 
      Mq2_t(14) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,14) 
      Mq2_t(15) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,15) 
      Mq2_t(16) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,16) 
      Mq2_t(17) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,17) 
      Mq2_t(18) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,18) 
      Mq2_t(19) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,19) 
      Mq2_t(2) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,20)
      Mq2_t(21) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,21) 
      Mq2_t(22) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,22) 
      Mq2_t(23) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,23) 
      Mq2_t(24) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,24)
      Mq2_t(25) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,25) 
      Mq2_t(26) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,26) 
      Mq2_t(27) = Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1,27) 

      Call Genrlizd_multipole_fin(Qm2_t,Om2_t,Am2_t,Dm2_t,Mq2_t,E,
     +                            D_iso,A_iso,O_iso,Q_iso,M_iso)

      Osc_Qm(Length) = Q_iso
      Osc_Om(Length) = O_iso
      Osc_Am(Length) = A_iso
      Osc_Mq(Length) = M_iso
      Osc_Dm(Length) = D_iso

      Return
      End
