      Subroutine Genrlizd_Beyond_dipole(Dens,Scr1,Scr2,Qm,Om,Am,
     +                                  Dm,Mq,Work,Maxcor,Nao,
     +                                  Jroot,Irrepy,Iroot,Irrepx)
    
      Implicit Double Precision(A-H,O-Z)

      Double precision Mq
      Dimension Dens(Nao,Nao),Scr1(Nao,Nao),Scr2(Nao,Nao)
      Dimension Work(Maxcor),Qm(4,8,4,6,6),Om(4,8,4,8,10),
     +                       Am(4,8,4,8,3),Mq(4,8,4,8,9),
     +                       Dm(4,8,4,8,3)

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
      Data Isix,Iten,Ithr,Inin /6,10,3,9/

      Length = Nao*(Nao+1)/2

C See below for the dipole contribution 

C Lets do the quadrapole contriubtion 
      Do Icomp = 1, 6
         Call Getrec(20,"JOBARC",Label_quad(Icomp),Length,Scr1)
         Call Expnd2(Scr1,Scr2,Nao)
         Qm(Jroot,Irrepy,Iroot,Irrepx,Icomp) = 
     +   Qm(Jroot,Irrepy,Iroot,Irrepx,Icomp) + 
     +   Ddot(Nao*Nao,Dens,1,Scr2,1)
      Enddo

C Similarly Octapole contribution

      Do Icomp = 1, 10
         Call Getrec(20,"JOBARC",Label_Octp(Icomp),Length,Scr1)
         Call Expnd2(Scr1,Scr2,Nao)
         Om(Jroot,Irrepy,Iroot,Irrepx,Icomp) = 
     +   Om(Jroot,Irrepy,Iroot,Irrepx,Icomp) + 
     +   Ddot(Nao*Nao,Dens,1,Scr2,1)
      Enddo

C Angular momentum contribution (rxp + s)

      Do Icomp = 1,3 
         Call Getrec(20,"JOBARC",Label_Angm(Icomp),Length,Scr1)
         Call Expnd2(Scr1,Scr2,Nao)
         Am(Jroot,Irrepy,Iroot,Irrepx,Icomp) = 
     +   Am(Jroot,Irrepy,Iroot,Irrepx,Icomp) + 
     +   Ddot(Nao*Nao,Dens,1,Scr2,1)
      Enddo

C For excited states <0|S|Psi_k> is zero. This term contribute
C only to ioniized or attached states. 

      Call Putrec(20,"JOBARC",A_string(Iside),Ithr,Am)

C Lets do the orbital magnetic quadrupole contribution 
C r(rxp) - (rxp)r

      Do Icomp = 1,3
         Call Getrec(20,"JOBARC",Label_Dmom(Icomp),Length,Scr1)
         Call Expnd2(Scr1,Scr2,Nao)
         Dm(Jroot,Irrepy,Iroot,Irrepx,Icomp) = 
     +   Dm(Jroot,Irrepy,Iroot,Irrepx,Icomp) + 
     +   Ddot(Nao*Nao,Dens,1,Scr2,1)
      Enddo

      Do Icomp = 1, 9
         Call Getrec(20,"JOBARC",Label_magd(Icomp),Length,Scr1)
         Call Expnd2(Scr1,Scr2,Nao)
         Mq(Jroot,Irrepy,Iroot,Irrepx,Icomp) = 
     +   Mq(Jroot,Irrepy,Iroot,Irrepx,Icomp) + 
     +   Ddot(Nao*Nao,Dens,1,Scr2,1)
      Enddo 

C These are the contributions from the second part of the 
C magtentic quadrupole operator (L*r) piece. This translate to 
C using dipole moment integrals. 
C This is (r^p) r is evaluated as 
C  (yd/dz-zd/dy)(x,y,z) =>  (0,-z,y)
C -(xd/dz-zd/dx)(x,y,z) => -(-z,0,x)
C  (yd/dx-xd/dy)(x,y,z) =>  (y,-x,0)

      Mq(Jroot,Irrepy,Iroot,Irrepx,1) = Mq(Jroot,Irrepy,Iroot,Irrepx,1)
      Mq(Jroot,Irrepy,Iroot,Irrepx,2) = Mq(Jroot,Irrepy,Iroot,Irrepx,2)
     +                                - Dm(Jroot,Irrepy,Iroot,Irrepx,3)
      Mq(Jroot,Irrepy,Iroot,Irrepx,3) = Mq(Jroot,Irrepy,Iroot,Irrepx,3) 
     +                                + Dm(Jroot,Irrepy,Iroot,Irrepx,2)

      Mq(Jroot,Irrepy,Iroot,Irrepx,4) = Mq(Jroot,Irrepy,Iroot,Irrepx,4) 
     +                                + Dm(Jroot,Irrepy,Iroot,Irrepx,3)
      Mq(Jroot,Irrepy,Iroot,Irrepx,5) = Mq(Jroot,Irrepy,Iroot,Irrepx,5) 
      Mq(Jroot,Irrepy,Iroot,Irrepx,6) = Mq(Jroot,Irrepy,Iroot,Irrepx,6) 
     +                                + Dm(Jroot,Irrepy,Iroot,Irrepx,1)

      Mq(Jroot,Irrepy,Iroot,Irrepx,7) = Mq(Jroot,Irrepy,Iroot,Irrepx,7)
     +                                + Dm(Jroot,Irrepy,Iroot,Irrepx,2)
      Mq(Jroot,Irrepy,Iroot,Irrepx,8) = Mq(Jroot,Irrepy,Iroot,Irrepx,8)
     +                                - Dm(Jroot,Irrepy,Iroot,Irrepx,8)
      Mq(Jroot,Irrepy,Iroot,Irrepx,9) = Mq(Jroot,Irrepy,Iroot,Irrepx,9) 

      Return 
      End 
