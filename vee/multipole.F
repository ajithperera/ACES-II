      Subroutine multipole(Root)
   
      Implicit Double Precision(A-H,O-Z)
      Double Precision Mq_r,Mq_l,M_strength,M_iso,M_fac
      Dimension Qm_r(6),Qm_l(6),Om_r(10),Om_l(10)
      Dimension Am_r(3),Am_l(3),Dm_r(3),Dm_l(3)
      Dimension Mq_r(3),Mq_l(3)
      Dimension A_strength(3),Q_strength(9),O_strength(10)
      Dimension D_strength(3),M_strength(27)
      Character*8 Q_string(2),O_string(2),A_string(2)
      Character*8 D_string(2),M_string(2)

      Data D_string/"TMRIGHT","TMLEFT"/
      Data Q_string/"QMRIGHT","QMLEFT"/
      Data O_string/"OMRIGHT","OMLEFT"/
      Data A_string/"AMRIGHT","AMLEFT"/
      Data M_string/"MQRIGHT","MQLEFT"/
      Data Isix,Iten,Ithr,Inin /6,10,3,9/
      Data Dzero,Two,Three,Twozero,Fourfive/0.0D0,2.0D0,3.0D0,
     +                                      20.0D0,45.0D0/
      Data Done,Dnine,Six/1.0D0,9.0D0,6.0D0/
      Data Recp_fsc/137.03599/

      Call Getrec(20,"JOBARC",D_string(1),Ithr,Dm_R)
      Call Getrec(20,"JOBARC",D_string(2),Ithr,Dm_L)

      Call Getrec(20,"JOBARC",Q_string(1),Isix,Qm_R)
      Call Getrec(20,"JOBARC",Q_string(2),Isix,Qm_L)

      Call Getrec(20,"JOBARC",O_string(1),Iten,Om_R)
      Call Getrec(20,"JOBARC",O_string(2),Iten,Om_L)

      Call Getrec(20,"JOBARC",A_string(1),Ithr,Am_R)
      Call Getrec(20,"JOBARC",A_string(2),Ithr,Am_L)

      Call Getrec(20,"JOBARC",M_string(1),Inin,Mq_R)
      Call Getrec(20,"JOBARC",M_string(2),Inin,Mq_L)

C Dipole (rx,ry,rz)

      Do Icomp = 1, 3
        D_strength(Icomp) = Dm_r(Icomp)*Dm_l(Icomp)
      Enddo 

C Orbital angular momentum (px,py,pz)

      Do Icomp = 1, 3
        A_strength(Icomp) = Am_r(Icomp)*Am_l(Icomp)
      Enddo 

C Quadrupole (Qxx,Qyy,Qzz,Qxy,Qxz,Qyz)

      Do Icomp = 1, 6
        Q_strength(Icomp) = Qm_r(Icomp)*Qm_l(Icomp)
      Enddo 
      Q_strength(7) =  Qm_r(1)*Qm_l(2)+Qm_r(2)*Qm_l(1)
      Q_strength(8) =  Qm_r(1)*Qm_l(3)+Qm_r(3)*Qm_l(1)
      Q_strength(9) =  Qm_r(2)*Qm_l(3)+Qm_r(3)*Qm_l(2)

C Dipole/Octupole (Oxxx;rx,Oyyy;ry,Ozzz;rz,Oxxy;ry,
C Oxxz;rz,Oyyx;rx,Oyyz;rz,Ozzx;rx,Ozzy;ry)

      O_strength(1) = Om_r(1)*Dm_l(1)+Dm_r(1)*Om_l(1)
      O_strength(2) = Om_r(2)*Dm_l(1)+Dm_r(2)*Om_l(2)
      O_strength(3) = Om_r(3)*Dm_l(1)+Dm_r(3)*Om_l(3)
      O_strength(4) = Om_r(4)*Dm_l(2)+Dm_r(2)*Om_l(4)
      O_strength(5) = Om_r(5)*Dm_l(3)+Dm_r(3)*Om_l(5)
      O_strength(6) = Om_r(6)*Dm_l(1)+Dm_r(1)*Om_l(6)
      O_strength(7) = Om_r(7)*Dm_l(3)+Dm_r(3)*Om_l(7)
      O_strength(8) = Om_r(8)*Dm_l(1)+Dm_r(1)*Om_l(8)
      O_strength(9) = Om_r(9)*Dm_l(2)+Dm_r(2)*Om_l(9)

      M_strength(1)  = Mq_r(1)*Dm_l(1)+Mq_l(1)*Dm_r(1)
      M_strength(2)  = Mq_r(1)*Dm_l(2)+Mq_l(1)*Dm_r(2)
      M_strength(3)  = Mq_r(1)*Dm_l(3)+Mq_l(1)*Dm_r(3)
      M_strength(4)  = Mq_r(2)*Dm_l(1)+Mq_l(1)*Dm_r(1)
      M_strength(5)  = Mq_r(2)*Dm_l(2)+Mq_l(2)*Dm_r(2)
      M_strength(6)  = Mq_r(2)*Dm_l(3)+Mq_l(2)*Dm_r(3)
      M_strength(7)  = Mq_r(3)*Dm_l(1)+Mq_l(3)*Dm_r(1)
      M_strength(8)  = Mq_r(3)*Dm_l(2)+Mq_l(3)*Dm_r(2)
      M_strength(9)  = Mq_r(3)*Dm_l(3)+Mq_l(3)*Dm_r(3)
      M_strength(10) = Mq_r(4)*Dm_l(1)+Mq_l(4)*Dm_r(1)
      M_strength(11) = Mq_r(4)*Dm_l(2)+Mq_l(4)*Dm_r(2)
      M_strength(12) = Mq_r(4)*Dm_l(3)+Mq_l(4)*Dm_r(3)
      M_strength(13) = Mq_r(5)*Dm_l(1)+Mq_l(5)*Dm_r(1)
      M_strength(14) = Mq_r(5)*Dm_l(2)+Mq_l(5)*Dm_r(2)
      M_strength(15) = Mq_r(5)*Dm_l(3)+Mq_l(5)*Dm_r(3)
      M_strength(16) = Mq_r(6)*Dm_l(1)+Mq_l(6)*Dm_r(1)
      M_strength(17) = Mq_r(6)*Dm_l(2)+Mq_l(6)*Dm_r(2)
      M_strength(18) = Mq_r(6)*Dm_l(3)+Mq_l(6)*Dm_r(3)
      M_strength(19) = Mq_r(7)*Dm_l(1)+Mq_l(7)*Dm_r(1)
      M_strength(20) = Mq_r(7)*Dm_l(2)+Mq_l(7)*Dm_r(2)
      M_strength(21) = Mq_r(7)*Dm_l(3)+Mq_l(7)*Dm_r(3)
      M_strength(22) = Mq_r(8)*Dm_l(1)+Mq_l(8)*Dm_r(1)
      M_strength(23) = Mq_r(8)*Dm_l(2)+Mq_l(8)*Dm_r(2)
      M_strength(24) = Mq_r(8)*Dm_l(3)+Mq_l(8)*Dm_r(3)
      M_strength(25) = Mq_r(9)*Dm_l(1)+Mq_l(9)*Dm_r(1)
      M_strength(26) = Mq_r(9)*Dm_l(2)+Mq_l(9)*Dm_r(2)
      M_strength(27) = Mq_r(9)*Dm_l(3)+Mq_l(9)*Dm_r(3)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " The tensor components" 
      Write(6,"(a)") " The <(L_0|mu|n><n|mu|0> "
      Write(6,"(3(1x,F15.10))") (D_strength(i),i=1,3)
      Write(6,"(a)") " The <L_0|r^p|n><n|r^p>"
      Write(6,"(3(1x,F15.10))") (A_strength(i),i=1,3)
      Write(6,"(a)") " The <L_0|Q|n><n|Q|0>  "
      Write(6,"(3(1x,F15.10))") (Q_strength(i),i=1,9)
      Write(6,"(a)") " The <L_0|mu|n><n|O|0> "
      Write(6,"(3(1x,F15.10))") (O_strength(i),i=1,9)
      Write(6,"(a)") " The <L_0|mu|n><n|Mq|0> "
      Write(6,"(3(1x,F15.10))") (M_strength(i),i=1,27)
#endif 
C Istropic averaging 

      Two_ov_three = Two/Three   
      One_ov_three = Done/Three
      One_ov_tzero = Done/Twozero 
      One_ov_nine  = Done/Dnine 
      Half         = Done/Two 
      Fsc   = Done/Recp_fsc

      D_iso = Dzero
      Do Icomp = 1, 3
         D_iso = D_iso + D_strength(Icomp)
      Enddo
      D_iso = Two_ov_three*Root*D_iso

      A_iso = Dzero
      Do Icomp = 1, 3
         A_iso = A_iso + A_strength(Icomp)
      Enddo
      A_iso = Two_ov_three*Root*A_iso

      Q_fac = One_ov_three*One_ov_tzero*Fsc**2
      Q_iso_1 = Two*(Q_strength(1)+Q_strength(2)+Q_strength(3))
      Q_iso_2 = Six*(Q_strength(4)+Q_strength(5)+Q_strength(6))
C The elements 7, 8 and 9 are build with a factor of two.
      Q_iso_3 = -(Q_strength(7)+Q_strength(8)+Q_strength(9))

      Q_iso = Q_iso_1 + Q_iso_2 + Q_iso_3
      Q_iso = Q_fac*Q_iso*Root**3

      O_iso = Dzero
      O_fac = Done/Fourfive
      O_fac = O_fac*Fsc**2
      Do Icomp = 1, 9
         O_iso = O_iso + O_Strength(Icomp)
      Enddo
      O_iso = O_fac*O_iso*Root**3

C The factor 1/3c associated with the integral is also absorbed 
C in here. The factor half is to account for the factor two 
C introduced in the symmetrization (see above)

      M_iso = Dzero 
      M_fac = Half*One_ov_nine*Fsc**2 
      Do Icomp = 1, 27
         M_iso = M_iso + M_strength(Icomp)
      Enddo
      M_iso = M_fac*M_iso*Root**2
       
      Write(6,*)
      Write(6,100) "Electric dipole contribution                   :",
     +              D_iso
      Write(6,100) "Magnetic dipole contribution                   :",
     +              A_iso
      Write(6,100) "Electric quadrupole contribution               :",
     +              Q_iso 
      Write(6,100) "Electric dipole/electric octupole contribution :",
     +             O_iso 
      Write(6,100) "Electric dipole/Magnetic octupole contribution :",
     +             M_iso 
 100  Format(1x,a,1x,es12.4e2)

      Return
      End 
