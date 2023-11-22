      Subroutine Osc_init(Osc_x,Osc_y,Osc_z,Tm_x,Tm_y,Tm_z,E,Fact,
     +                    Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,Jroot_2,
     +                    Irrepy_2,Iroot_2,Irrepx_2)

      Implicit Double Precision(A-H,O-z)

      Dimension Tm_x(4,8,4,8) 
      Dimension Tm_y(4,8,4,8) 
      Dimension Tm_z(4,8,4,8) 
      Dimension Osc_x(4,8,4,8,4,8,4,8)
      Dimension Osc_y(4,8,4,8,4,8,4,8)
      Dimension Osc_z(4,8,4,8,4,8,4,8)
 
      Osc_x(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +      Jroot_2,Irrepy_2,Iroot_2,Irrepx_2)=
     +                     Tm_x(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1)*
     +                     Tm_x(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2)*
     +                     E*Fact


      Osc_y(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +      Jroot_2,Irrepy_2,Iroot_2,Irrepx_2)=
     +                     Tm_y(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1)*
     +                     Tm_y(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2)*
     +                     E*Fact

      Osc_z(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +      Jroot_2,Irrepy_2,Iroot_2,Irrepx_2)=
     +                     Tm_z(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1)*
     +                     Tm_z(Jroot_2,Irrepy_2,Iroot_2,Irrepx_2)*
     +                     E*Fact

      Return
      End


