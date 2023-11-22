      Subroutine Form_genrlizd_oscstrength(Work,Maxcor,Iuhf,Eigval)

      Implicit Double Precision(A-H,O-Z)

      Character*8 Label_xyz(3),Label
      Character*2 Suffix(32)

      Dimension Tm(3,4,8,4,8), Eigval(100,8)
      Dimension Tm_x(4,8,4,8), Tm_y(4,8,4,8), Tm_z(4,8,4,8)
      Dimension Osc(1024)
      Dimension Work(Maxcor)

#include "flags.com"
#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"
#include "machsp.com"
       
      Common /Calcinfo/Nroot(8)

      Data Label_xyz/'DIPOLE_X','DIPOLE_Y','DIPOLE_Z'/
      Data Suffix /"01","02","03","04","05","06","07","08","09","10",
     +             "11","12","13","14","15","16","17","18","19","20",
     +             "21","22","23","24","25","26","27","28","29","30",
     +             "31","32"/
      Data Ione, Ithree, Fact, Done, Dnull /1,3,0.66666666666666D0,
     +                                      1.0D0,0.0D0/

      Call Getrec(20,"JOBARC","NBASTOT ",Ione,Nao)
      Nmo  = Nocco(1) + Nvrto(1)
      Nmo2 = Nmo*Nmo
      Nao2 = Nao*Nao

      Itmp1 = Ione
      Itmpx = Itmp1 + Nao2
      Itmpy = Itmpx + Nao2 
      Itmpz = Itmpy + Nao2
      Itmp3 = Itmpz + Nao2
      Iend  = Itmp3 + Nao2
      
      Len_drec = Nmo*(Nmo+1)/2
      Do Ixyz = 1, 3
         Call Getrec(20,"JOBARC",Label_xyz(Ixyz),Len_drec*IIntfp,
     +               Work(Itmp1))
         If (Ixyz.Eq.1) Call Expnd2(Work(Itmp1),Work(Itmpx),Nao)
         If (Ixyz.Eq.2) Call Expnd2(Work(Itmp1),Work(Itmpy),Nao)
         If (Ixyz.Eq.3) Call Expnd2(Work(Itmp1),Work(Itmpz),Nao)
      Enddo
      
      Call Dzero(Tm_x,1024)
      Call Dzero(Tm_y,1024)
      Call Dzero(Tm_z,1024)

      Ndone = 0
      Do Irrepx = 1, Nirrep
         Do Iroot = 1, Nroot(Irrepx)
            Do Irrepy = 1, Nirrep
               Do Jroot = 1, Nroot(Irrepy)
                   Do Ispin = 1, Iuhf+1
                      Ndone = Ndone + 1
                      Label = "GTDENS"//Suffix(Ndone)
                      Call Getrec(20,"JOBARC",LABEL,Nao2*Iintfp, 
     +                            Work(Itmp3))
#ifdef _DEBUG_LVL0
      Write(6,"(a,2(1x,i2))") "Irrepx,Irrepy : ",
     +                         Irrepx,Irrepy
      Write(6,"(a,2(1x,i2))") "Rootx,Rooty : ",
     +                         Iroot,Jroot 
      Call Checksum("TDENS-READ ",Work(Itmp3),Nao2,s)
#endif
                      Do Ixyz = 1, 3
                         IF (Ixyz .Eq.1) 
     +                                Tm_x(Jroot,Irrepy,Iroot,Irrepx) = 
     +                                Tm_x(Jroot,Irrepy,Iroot,Irrepx) +
     +                                Ddot(Nao2,Work(Itmp3),1,
     +                                Work(Itmpx),1)
                         If (Ixyz .Eq.2) 
     +                                Tm_y(Jroot,Irrepy,Iroot,Irrepx) = 
     +                                Tm_y(Jroot,Irrepy,Iroot,Irrepx) + 
     +                                Ddot(Nao2,Work(Itmp3),1,
     +                                Work(Itmpy),1)
                         If (Ixyz .Eq.3) 
     +                               Tm_z(Jroot,Irrepy,Iroot,Irrepx) = 
     +                               Tm_z(Jroot,Irrepy,Iroot,Irrepx) + 
     +                                Ddot(Nao2,Work(Itmp3),1,
     +                                Work(Itmpz),1)
                      Enddo 
                   Enddo
               Enddo
            Enddo 
         Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,"(20x,a)") "Ixyz,root(irrepy),irrepy,root(irrepx), irrepx"
      Write(6,"(20x,a)") "---------------------------------------------"
      Do Irrepx = 1, Nirrep
      Do Iroot  = 1, Nroot(Irrepx)
      Do Irrepy = 1, Nirrep
      Do Jroot  = 1, Nroot(Irrepy)
         write(6,"(a,(4(2x,i2)),3(1x,f10.5))") "mu(pq)<L|(p+q)|R>: ", 
     +            Jroot,Irrepy,Iroot,Irrepx,Tm_x(Jroot,Irrepy,
     +            Iroot,Irrepx),Tm_y(Jroot,Irrepy,Iroot,Irrepx),
     +            Tm_z(Jroot,Irrepy,Iroot,Irrepx)
      Enddo  
      Enddo  
      Enddo  
      Enddo  
      Write(6,"(20x,a)") "---------------------------------------------"
#endif
      I000 = Ione
      I00y = I00x + 1048576
      I00z = I00y + 1048576
      Iend = I00z + 1048576
      If (Iend .Ge. Maxcor) Call Insmem("form_genrlizd_oscstrength",
     +                                   Iend,Maxor)

      Length = 0 
      Do Irrepx_1 = 1, Nirrep
         Do Iroot_1 = 1, Nroot(Irrepx_1)
            E_i = Eigval(Iroot_1,Irrepx_1)
            Do Irrepy_1 = 1, Nirrep
               Do Jroot_1 = 1, Nroot(Irrepy_1)
                  E_j = Eigval(jroot_1,Irrepy_1)

#ifdef _DEBUG_LVL0
          write(6,*)
          Write(6,"(2a,(4(1x,i2)),2(1x,F15.8))") "Irrepx_1,Iroot_1,",
     +    "Irrepy_1,Jroot_1,E_i,E_j:",Irrepx_1,Iroot_1,Irrepy_1,
     +    Jroot_1,E_i,E_j
#endif 
                  If (E_i .eq. E_j) Then
                      E = Dnull 
                  Else
                      E = Dabs(E_i - E_j)
                  Endif 

                  Do Irrepx_2 = 1, Nirrep
                     Do Iroot_2 = 1, Nroot(Irrepx_2)
                        Do Irrepy_2 = 1, Nirrep
                           Do Jroot_2 = 1, Nroot(Irrepy_2)
                           
#ifdef _DEBUG_LVL0
          Write(6,"(2a,(4(1x,i2),1x,F15.8))")"Irrepx_2,Iroot_2",
     +                                       "Irrepy_2,Jroot_2,E:",
     +                                        Irrepx_2,Iroot_2,
     +                                        Irrepy_2,Jroot_2,E
#endif 
                           Call osc_init(Work(I00X),Work(I00y),
     +                                   Work(I00z),TM_x,TM_y,TM_z,
     +                                   E,Fact,Jroot_1,Irrepy_1,
     +                                   Iroot_1,Irrepx_1,Jroot_2,
     +                                   Irrepy_2,Iroot_2,Irrepx_2)
                           Enddo 
                        Enddo
                     Enddo
                  Enddo

                  Length = Length + 1 
                  Call osc_fin(Osc,Work(I00X),Work(I00y),Work(I00z),
     +                                  Jroot_1,Irrepy_1,Iroot_1,
     +                                  Irrepx_1,Length)

               Enddo
            Enddo
         Enddo
      Enddo

      Call Putrec(20,"JOBARC","GEN_OSCS",Length,Osc)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(20x,2a)") "Ixyz,root(irrepy),irrepy,root(irrepx)",
     +                     " irrepx", " Osc. strength"
      Write(6,"(20x,a)") "---------------------------------------------"
      Ndone = 0
      Do Irrepx = 1, Nirrep
      Do Iroot  = 1, Nroot(Irrepx)
      Do Irrepy = 1, Nirrep
      Do Jroot  = 1, Nroot(Irrepy)
         Ndone = Ndone + 1
         write(6,"(a,(4(7x,i2)),f10.5)") "Osc. Strength: ",
     +            Jroot,Irrepy,Iroot,Irrepx,Osc(Ndone)
      Enddo
      Enddo
      Enddo
      Enddo
      Write(6,"(20x,a)") "---------------------------------------------"
#endif
      Return
      End
