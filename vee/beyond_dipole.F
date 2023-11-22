      Subroutine Beyond_dipole(Dens,Scr1,Scr2,Tm,Work,Maxcor,Nao, 
     +                         Iside)
    
      Implicit Double Precision(A-H,O-Z)

      Double precision Mq
      Dimension Dens(Nao,Nao),Scr1(Nao,Nao),Scr2(Nao,Nao)
      Dimension Work(Maxcor),Tm(3),Qm(6),Om(10),Am(3),Mq(9),Dm(3)

      Character*8 Label_Angm(3)
      Character*8 Label_dmom(3)
      Character*8 Label_octp(10)
      Character*8 Label_quad(6)
      Character*8 Label_magd(9)
      Character*8 Q_string(2),O_string(2),A_string(2),M_string(2)

      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),IOFFV(8,2),IOFFO(8,2),
     +             IRPDPDAO(8),IRPDPDAOS(8),
     +             ISTART(8,8),ISTARTMO(8,3)
 
      Data Label_dmom/'DIPOLE_X','DIPOLE_Y','DIPOLE_Z'/
      Data Label_quad/'QUAD_XX ','QUAD_YY ','QUAD_ZZ ','QUAD_XY ',
     +                'QUAD_XZ ','QUAD_YZ '/
      Data Label_octp/'OCTUPXXX','OCTUPYYY','OCTUPZZZ','OCTUPXXY',
     +                'OCTUPXXZ','OCTUPXYY','OCTUPYYZ','OCTUPXZZ',
     +                'OCTUPYZZ','OCTUPXYZ'/
      Data Label_angm/'   AMX  ', '   AMY  ','   AMZ  '/
      Data Label_magd/'   MQXX ','   MQYX ','   MQZX ','   MQXY ',
     +                '   MQYY ','   MQZY ','   MQXZ ','   MQYZ ',
     +                '   MQZZ '/
      Data Q_string/"QMRIGHT","QMLEFT"/
      Data O_string/"OMRIGHT","OMLEFT"/
      Data A_string/"AMRIGHT","AMLEFT"/
      Data M_string/"MQRIGHT","MQLEFT"/
      Data Isix,Iten,Ithr,Inin /6,10,3,9/

      Call Dzero(Dm,3)
      Call Dzero(Qm,6)
      Call Dzero(Om,10)
      Call Dzero(Am,3)
      Call Dzero(Mq,9)
     
C Lets do the quadrapole contriubtion 

      Length = Nao*(Nao+1)/2
      Do Icomp = 1, 6
         Call Getrec(20,"JOBARC",Label_quad(Icomp),Length,Scr1)
         Call Expnd2(Scr1,Scr2,Nao)
         Qm(Icomp) = Qm(Icomp) + Ddot(Nao*Nao,Dens,1,Scr2,1)
      Enddo

      Call Putrec(20,"JOBARC",Q_string(Iside),Isix,Qm)

C Similarly Octapole contribution

      Do Icomp = 1, 10
         Call Getrec(20,"JOBARC",Label_Octp(Icomp),Length,Scr1)
         Call Expnd2(Scr1,Scr2,Nao)
         Om(Icomp) = Om(Icomp) + Ddot(Nao*Nao,Dens,1,Scr2,1)
      Enddo

      Call Putrec(20,"JOBARC",O_string(Iside),Iten,Om)

C Angular momentum contribution (rxp + s)

      Do Icomp = 1,3 
         Call Getrec(20,"JOBARC",Label_Angm(Icomp),Length,Scr1)
         Call Expnd2(Scr1,Scr2,Nao)
         Am(Icomp) = Am(Icomp) + Ddot(Nao*Nao,Dens,1,Scr2,1)
      Enddo

C For excited states <0|S|Psi_k> is zero. This term contribute
C only to ioniized or attached states. 

      Call Putrec(20,"JOBARC",A_string(Iside),Ithr,Am)

C Lets do the orbital magnetic quadrupole contribution 
C r(rxp) - (rxp)r

      Do Icomp = 1,3
         Call Getrec(20,"JOBARC",Label_Dmom(Icomp),Length,Scr1)
         Call Expnd2(Scr1,Scr2,Nao)
         Dm(Icomp) = Dm(Icomp) + Ddot(Nao*Nao,Dens,1,Scr2,1)
      Enddo

      Do Icomp = 1, 9
         Call Getrec(20,"JOBARC",Label_magd(Icomp),Length,Scr1)
         Call Expnd2(Scr1,Scr2,Nao)
         Mq(Icomp) = Mq(Icomp) + Ddot(Nao*Nao,Dens,1,Scr2,1)
      Enddo 

C These are the contributions from the second part of the 
C magtentic quadrupole operator (L*r) piece. This translate to 
C using dipole moment integrals. 
C This is (r^p) r is evaluated as 
C  (yd/dz-zd/dy)(x,y,z) =>  (0,-z,y)
C -(xd/dz-zd/dx)(x,y,z) => -(-z,0,x)
C  (yd/dx-xd/dy)(x,y,z) =>  (y,-x,0)

      Mq(1) = Mq(1)
      Mq(2) = Mq(2) - Dm(3)
      Mq(3) = Mq(3) + Dm(2)

      Mq(4) = Mq(4) + Dm(3)
      Mq(5) = Mq(5) 
      Mq(6) = Mq(6) - Dm(1)

      Mq(7) = Mq(7) + Dm(2)
      Mq(8) = Mq(8) - Dm(1)
      Mq(9) = Mq(9) 

      Call Putrec(20,"JOBARC",M_string(Iside),Inin,Mq)

      Return 
      End 
