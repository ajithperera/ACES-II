      Subroutine Osc_bd_init(Qm,Om,Am,Dm,Mq,Qm_2,Om_2,Am_2,Dm_2,Mq_2,
     +                       Work,Maxcor,Irrepx_1,Iroot_1,Irrepy_1,
     +                       jroot_1,Irrepx_2,Iroot_2,Irrepy_2,
     +                       jroot_2)
    
      Implicit Double Precision(A-H,O-Z)

      Double precision Mq,Mq_2,Mq2_t,Mq_t
      Dimension Work(Maxcor),Qm(4,8,4,8,9),Om(4,8,4,8,10),
     +                       Am(4,8,4,8,3),Mq(4,8,4,8,9),
     +                       Dm(4,8,4,8,3)

      Dimension Qm_tr(9),Om_tr(10),Am_tr(3),Dm_tr(3),Mq_tr(9)
      Dimension Qm_tl(9),Om_tl(10),Am_tl(3),Dm_tl(3),Mq_tl(9)
      Dimension Qm2_t(9),Om2_t(10),Am2_t(3),Dm2_t(3),Mq2_t(27)

      Dimension Qm_2(4,8,4,6,4,8,4,6,9)
      Dimension Om_2(4,8,4,6,4,8,4,6,10)
      Dimension Am_2(4,8,4,6,4,8,4,6,3)
      Dimension Mq_2(4,8,4,6,4,8,4,6,27)
      Dimension Dm_2(4,8,4,6,4,8,4,6,3)

#ifdef _DEBUG_LVL0
      Write(6,"(a,1x4i3)") "Iroot_1,Irrepx_1,jroot_1,Irrepy_1: ",
     +                      Iroot_1,Irrepx_1,jroot_1,Irrepy_1
      Write(6,"(a,1x,4i3)") "Jroot_2,Irrepy_2,Iroot_2,Irrepx_2: ",
     +                      Jroot_2,Irrepy_2,Iroot_2,Irrepx_2
C      Write(6,"(3(1x,F10.5))") 
C     +   Dm(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,1),
C     +   Dm(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,2),
C     +   Dm(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,3)
#endif 
      Qm_tl(1) = Qm(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,1)
      Qm_tl(2) = Qm(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,2)
      Qm_tl(3) = Qm(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,3)
      Qm_tl(4) = Qm(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,4)
      Qm_tl(5) = Qm(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,5)
      Qm_tl(6) = Qm(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,6)
      Qm_tl(7) = Qm(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,7)
      Qm_tl(8) = Qm(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,8)
      Qm_tl(9) = Qm(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,9)

      Om_tl(1)  = Om(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,1)
      Om_tl(2)  = Om(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,2)
      Om_tl(3)  = Om(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,3)
      Om_tl(4)  = Om(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,4)
      Om_tl(5)  = Om(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,5)
      Om_tl(6)  = Om(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,6)
      Om_tl(7)  = Om(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,7)
      Om_tl(8)  = Om(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,8)
      Om_tl(9)  = Om(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,9)
      Om_tl(10) = Om(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,10)
 
      Dm_tl(1) = Dm(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,1)
      Dm_tl(2) = Dm(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,2)
      Dm_tl(3) = Dm(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,3)

      Am_tl(1) = Am(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,1)
      Am_tl(2) = Am(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,2)
      Am_tl(3) = Am(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,3)

      Mq_tl(1) = Mq(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,1)
      Mq_tl(2) = Mq(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,2)
      Mq_tl(3) = Mq(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,3)
      Mq_tl(4) = Mq(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,4)
      Mq_tl(5) = Mq(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,5)
      Mq_tl(6) = Mq(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,6)
      Mq_tl(7) = Mq(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,7)
      Mq_tl(8) = Mq(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,8)
      Mq_tl(9) = Mq(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,9)

      Qm_tr(1) = Qm(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,1)
      Qm_tr(2) = Qm(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,2)
      Qm_tr(3) = Qm(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,3)
      Qm_tr(4) = Qm(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,4)
      Qm_tr(5) = Qm(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,5)
      Qm_tr(6) = Qm(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,6)
      Qm_tr(7) = Qm(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,7)
      Qm_tr(8) = Qm(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,8)
      Qm_tr(9) = Qm(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,9)

      Om_tr(1)  = Om(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,1)
      Om_tr(2)  = Om(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,2)
      Om_tr(3)  = Om(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,3)
      Om_tr(4)  = Om(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,4)
      Om_tr(5)  = Om(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,5)
      Om_tr(6)  = Om(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,6)
      Om_tr(7)  = Om(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,7)
      Om_tr(8)  = Om(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,8)
      Om_tr(9)  = Om(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,9)
      Om_tr(10) = Om(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,10)

      Dm_tr(1) = Dm(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,1)
      Dm_tr(2) = Dm(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,2)
      Dm_tr(3) = Dm(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,3)

      Am_tr(1) = Am(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,1)
      Am_tr(2) = Am(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,2)
      Am_tr(3) = Am(tJroot_1,Irrepy_1,Iroot_1,Irrepx_1,3)

      Mq_tr(1) = Mq(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,1)
      Mq_tr(2) = Mq(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,2)
      Mq_tr(3) = Mq(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,3)
      Mq_tr(4) = Mq(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,4)
      Mq_tr(5) = Mq(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,5)
      Mq_tr(6) = Mq(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,6)
      Mq_tr(7) = Mq(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,7)
      Mq_tr(8) = Mq(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,8)
      Mq_tr(9) = Mq(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,9)

       
      Call Genrlizd_multipole_init(Qm_tl,Om_tl,Am_tl,Mq_tl,Dm_tl,
     +                             Qm_tr,Om_tr,Am_tr,Mq_tr,Dm_tr,
     +                             Qm2_t,Om2_t,Am2_t,Mq2_t,Dm2_t)

      Qm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,1) = Qm2_t(1)
      Qm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,2) = Qm2_t(2)
      Qm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,3) = Qm2_t(3)
      Qm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,4) = Qm2_t(4)
      Qm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,5) = Qm2_t(5)
      Qm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,6) = Qm2_t(6)
      Qm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,7) = Qm2_t(7)
      Qm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,8) = Qm2_t(8)
      Qm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,9) = Qm2_t(9)

      Om_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,1) = Om2_t(1)
      Om_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,2) = Om2_t(2)
      Om_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,3) = Om2_t(3)
      Om_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,4) = Om2_t(4)
      Om_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,5) = Om2_t(5)
      Om_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,6) = Om2_t(6)
      Om_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,7) = Om2_t(7)
      Om_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,8) = Om2_t(8)
      Om_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,9) = Om2_t(9)

      Am_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,1) = Am2_t(1)
      Am_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,2) = Am2_t(2)
      Am_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,2) = Am2_t(3)

      Dm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,1) = Dm2_t(1)
      Dm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,2) = Dm2_t(2)
      Dm_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,3) = Dm2_t(3)

      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,1) = Mq2_t(1)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,2) = Mq2_t(2)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,3) = Mq2_t(3)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,4) = Mq2_t(4)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,5) = Mq2_t(5)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,6) = Mq2_t(6)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,7) = Mq2_t(7)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,8) = Mq2_t(8)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,9) = Mq2_t(9)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,10) = Mq2_t(10)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,11) = Mq2_t(11)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,12) = Mq2_t(12)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,13) = Mq2_t(13)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,14) = Mq2_t(14)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,15) = Mq2_t(15)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,16) = Mq2_t(16)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,17) = Mq2_t(17)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,18) = Mq2_t(18)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,19) = Mq2_t(19)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,20) = Mq2_t(20)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,21) = Mq2_t(21)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,22) = Mq2_t(22)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,23) = Mq2_t(23)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,24) = Mq2_t(24)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,25) = Mq2_t(25)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,26) = Mq2_t(26)
      Mq_2(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +     Jroot_2,Irrepy_2,Iroot_2,Irrepx_2,27) = Mq2_t(27)

      Return
      End
