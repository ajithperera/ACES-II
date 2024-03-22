      Subroutine Genrlizd_multipole_fin(Q_strength,O_strength,
     +                                  A_strength,D_strength,
     +                                  M_strength,Root,D_iso,A_iso,
     +                                  O_iso,Q_iso,M_iso)
   
      Implicit Double Precision(A-H,O-Z)
      Double Precision Mq_r,Mq_l,M_strength,M_iso,M_fac

      Data Isix,Iten,Ithr,Inin /6,10,3,9/
      Data Dzero,Two,Three,Twozero,Fourfive/0.0D0,2.0D0,3.0D0,
     +                                      20.0D0,45.0D0/
      Data Done,Dnine,Six/1.0D0,9.0D0,6.0D0/
      Data Recp_fsc/137.03599/

      Dimension A_strength(3),Q_strength(9),O_strength(10)
      Dimension D_strength(3),M_strength(27)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " @-genrlizd_multipole_fin The tensor components" 
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

