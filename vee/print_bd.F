      Subroutine Print_bd1(Array1,Ndim,Nloop,Nirrep,Nroot,Nc)

      Implicit Double Precision(A-H,O-Z)
      Dimension Array1(4,8,4,8,Ndim),Nroot(8),Array2(4,8,4,8,4,8,4,8,
     +          Ndim)

#ifdef _DEBUG_LVL0
      If (Nloop .Eq. 4) Then
      Do Irrepx = 1, Nirrep
         Do Iroot = 1, Nroot(Irrepx)
            Do Irrepy = 1, Nirrep
               Do Jroot = 1, Nroot(Irrepy)
               Write(6,"(2a,(4(1x,i2)))") "Irrepx,Iroot,",
     +        "Irrepy,Jroot:",Irrepx,Iroot,Irrepy,Jroot
               If (Nc.eq.6)Write(6,"(6(1x,F15.9))")(Array1(Jroot,Irrepy,
     +                                    Iroot,Irrepx,Ic),Ic=1,Nc)
               If (Nc.eq.3)Write(6,"(3(1x,F15.9))")(Array1(Jroot,Irrepy,
     +                                    Iroot,Irrepx,Ic),Ic=1,Nc)
               If (Nc.eq.10)Write(6,"(6(1x,F15.9))")(Array1(Jroot,
     +                                Irrepy,Iroot,Irrepx,Ic),Ic=1,Nc)
               If (Nc.eq.9)Write(6,"(6(1x,F15.9))") (Array1(Jroot,
     +                                Irrepy,Iroot,Irrepx,Ic),Ic=1,Nc)
               Enddo
            Enddo
         Enddo
      Enddo
      Elseif(Nloop .Eq. 8) Then

      Do Irrepx_1 = 1, Nirrep
         Do Iroot_1 = 1, Nroot(Irrepx_1)
            Do Irrepy_1 = 1, Nirrep
               Do Jroot_1 = 1, Nroot(Irrepy_1)
#ifdef _DEBUG_LVL0
          write(6,*)
          Write(6,"(2a,(4(1x,i2)),2(1x,F15.8))") "Irrepx_1,Iroot_1,",
     +    "Irrepy_1,Jroot_1,E_i,E_j:",Irrepx_1,Iroot_1,Irrepy_1,
     +    Jroot_1,E_i,E_j
#endif
                  Do Irrepx_2 = 1, Nirrep
                     Do Iroot_2 = 1, Nroot(Irrepx_2)
                        Do Irrepy_2 = 1, Nirrep
                           Do Jroot_2 = 1, Nroot(Irrepy_2)
#ifdef _DEBUG_LVL0
          Write(6,"(2a,(4(1x,i2),1x,F15.8))")"Irrepx_2,Iroot_2",
     +                                       "Irrepy_2,Jroot_2:",
     +                                        Irrepx_2,Iroot_2,
     +                                        Irrepy_2,Jroot_2,E
               If (Nc.eq.6)Write(6,"(6(1x,F15.9))")
     +            (Array2(Irrepx_1,Iroot_1,Irrepy_1,Jroot_1,Irrepx_2,
     +                    Iroot_2,Irrepy_2,jroot_2,Ic),Ic=1,Nc)
               If (Nc.eq.3)Write(6,"(3(1x,F15.9))") 
     +            (Array2(Irrepx_1,Iroot_1,Irrepy_1,Jroot_1,Irrepx_2,
     +                    Iroot_2,Irrepy_2,jroot_2,Ic),Ic=1,Nc)
               If (Nc.eq.10)Write(6,"(6(1x,F15.9))") 
     +            (Array2(Irrepx_1,Iroot_1,Irrepy_1,Jroot_1,Irrepx_2,
     +                   Iroot_2,Irrepy_2,jroot_2,Ic),Ic=1,Nc)
               If (Nc.eq.9)Write(6,"(6(1x,F15.9))") 
     +            (Array2(Irrepx_1,Iroot_1,Irrepy_1,Jroot_1,Irrepx_2,
     +                   Iroot_2,Irrepy_2,jroot_2,Ic),Ic=1,Nc)


#endif
                           Enddo
                        Enddo
                     Enddo
                  Enddo

             Enddo
           Enddo 
         Enddo
       Enddo
      Endif 

#endif

      Return
      End
