      Subroutine Genrlizd_multipole_init(Qml,Oml,Aml,Mql,Dml,Qmr,Omr,
     +                                   Amr,Mqr,Dmr,Qm_2,Om_2,Am_2, 
     +                                   Mq_2,Dm_2)
   
      Implicit Double Precision(A-H,O-Z)
      Double Precision Mq,Mq_2
      Dimension Qmr(6),Omr(10),Amr(3),Dmr(3),Mqr(9)
      Dimension Qml(6),Oml(10),Aml(3),Dml(3),Mql(9)
      Dimension Qm_2(9),Om_2(10),Am_2(3),Dm_2(3),Mq_2(27)

      Data Isix,Iten,Ithr,Inin /6,10,3,9/
      Data Dzero,Two,Three,Twozero,Fourfive/0.0D0,2.0D0,3.0D0,
     +                                      20.0D0,45.0D0/
      Data Done,Dnine,Six/1.0D0,9.0D0,6.0D0/
      Data Recp_fsc/137.03599/
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "@-Enrty Genrlizd_multipole_init"
      Write(6,"(a)") " The tensor components" 
      Write(6,"(a)") " The <(L_0|mu|n><n|mu|0> "
      Write(6,"(3(1x,F15.10))") (Dmr(i),i=1,3)
      Write(6,"(a)") " The <L_0|r^p|n><n|r^p>"
      Write(6,"(3(1x,F15.10))") (Amr(i),i=1,3)
      Write(6,"(a)") " The <L_0|Q|n><n|Q|0>  "
      Write(6,"(3(1x,F15.10))") (Qmr(i),i=1,9)
      Write(6,"(a)") " The <L_0|mu|n><n|O|0> "
      Write(6,"(3(1x,F15.10))") (Omr(i),i=1,10)
      Write(6,"(a)") " The <L_0|mu|n><n|Mq|0> "
      Write(6,"(3(1x,F15.10))") (Mqr(i),i=1,9)
#endif 

C Dipole (rx,ry,rz)

      Do Icomp = 1, 3
        Dm_2(Icomp) = Dmr(Icomp)*Dml(Icomp)
      Enddo 

C Orbital angular momentum (px,py,pz)

      Do Icomp = 1, 3
        Am_2(Icomp) = Amr(Icomp)*Aml(Icomp)
      Enddo 

C Quadrupole (Qxx,Qyy,Qzz,Qxy,Qxz,Qyz)

      Do Icomp = 1, 6
        Qm_2(Icomp) = Qmr(Icomp)*Qml(Icomp)
      Enddo 
      Qm_2(7) =  Qmr(1)*Qml(2)+Qmr(2)*Qml(1)
      Qm_2(8) =  Qmr(1)*Qml(3)+Qmr(3)*Qml(1)
      Qm_2(9) =  Qmr(2)*Qml(3)+Qmr(3)*Qml(2)

C Dipole/Octupole (Oxxx;rx,Oyyy;ry,Ozzz;rz,Oxxy;ry,
C Oxxz;rz,Oyyx;rx,Oyyz;rz,Ozzx;rx,Ozzy;ry)

      Om_2(1) = Omr(1)*Dmr(1)+Dmr(1)*Oml(1)
      Om_2(2) = Omr(2)*Dmr(1)+Dmr(2)*Oml(2)
      Om_2(3) = Omr(3)*Dmr(1)+Dmr(3)*Oml(3)
      Om_2(4) = Omr(4)*Dmr(2)+Dmr(2)*Oml(4)
      Om_2(5) = Omr(5)*Dmr(3)+Dmr(3)*Oml(5)
      Om_2(6) = Omr(6)*Dmr(1)+Dmr(1)*Oml(6)
      Om_2(7) = Omr(7)*Dmr(3)+Dmr(3)*Oml(7)
      Om_2(8) = Omr(8)*Dmr(1)+Dmr(1)*Oml(8)
      Om_2(9) = Omr(9)*Dmr(2)+Dmr(2)*Oml(9)

      Mq_2(1)  = Mqr(1)*Dml(1)+Mql(1)*Dmr(1)
      Mq_2(2)  = Mqr(1)*Dml(2)+Mql(1)*Dmr(2)
      Mq_2(3)  = Mqr(1)*Dml(3)+Mql(1)*Dmr(3)
      Mq_2(4)  = Mqr(2)*Dml(1)+Mql(1)*Dmr(1)
      Mq_2(5)  = Mqr(2)*Dml(2)+Mql(2)*Dmr(2)
      Mq_2(6)  = Mqr(2)*Dml(3)+Mql(2)*Dmr(3)
      Mq_2(7)  = Mqr(3)*Dml(1)+Mql(3)*Dmr(1)
      Mq_2(8)  = Mqr(3)*Dml(2)+Mql(3)*Dmr(2)
      Mq_2(9)  = Mqr(3)*Dml(3)+Mql(3)*Dmr(3)
      Mq_2(10) = Mqr(4)*Dml(1)+Mql(4)*Dmr(1)
      Mq_2(11) = Mqr(4)*Dml(2)+Mql(4)*Dmr(2)
      Mq_2(12) = Mqr(4)*Dml(3)+Mql(4)*Dmr(3)
      Mq_2(13) = Mqr(5)*Dml(1)+Mql(5)*Dmr(1)
      Mq_2(14) = Mqr(5)*Dml(2)+Mql(5)*Dmr(2)
      Mq_2(15) = Mqr(5)*Dml(3)+Mql(5)*Dmr(3)
      Mq_2(16) = Mqr(6)*Dml(1)+Mql(6)*Dmr(1)
      Mq_2(17) = Mqr(6)*Dml(2)+Mql(6)*Dmr(2)
      Mq_2(18) = Mqr(6)*Dml(3)+Mql(6)*Dmr(3)
      Mq_2(19) = Mqr(7)*Dml(1)+Mql(7)*Dmr(1)
      Mq_2(20) = Mqr(7)*Dml(2)+Mql(7)*Dmr(2)
      Mq_2(21) = Mqr(7)*Dml(3)+Mql(7)*Dmr(3)
      Mq_2(22) = Mqr(8)*Dml(1)+Mql(8)*Dmr(1)
      Mq_2(23) = Mqr(8)*Dml(2)+Mql(8)*Dmr(2)
      Mq_2(24) = Mqr(8)*Dml(3)+Mql(8)*Dmr(3)
      Mq_2(25) = Mqr(9)*Dml(1)+Mql(9)*Dmr(1)
      Mq_2(26) = Mqr(9)*Dml(2)+Mql(9)*Dmr(2)
      Mq_2(27) = Mqr(9)*Dml(3)+Mql(9)*Dmr(3)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "@-Exit Genrlizd_multipole_init"
      Write(6,"(a)") " The tensor components" 
      Write(6,"(a)") " The <(L_0|mu|n><n|mu|0> "
      Write(6,"(3(1x,F15.10))") (Dm_2(i),i=1,3)
      Write(6,"(a)") " The <L_0|r^p|n><n|r^p>"
      Write(6,"(3(1x,F15.10))") (Am_2(i),i=1,3)
      Write(6,"(a)") " The <L_0|Q|n><n|Q|0>  "
      Write(6,"(3(1x,F15.10))") (Qm_2(i),i=1,9)
      Write(6,"(a)") " The <L_0|mu|n><n|O|0> "
      Write(6,"(3(1x,F15.10))") (Om_2(i),i=1,9)
      Write(6,"(a)") " The <L_0|mu|n><n|Mq|0> "
      Write(6,"(3(1x,F15.10))") (Mq_2(i),i=1,27)
#endif 

      Return
      End 
