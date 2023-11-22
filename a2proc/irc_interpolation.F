      Subroutine Irc_interpolation(Cpi,Vec1,Cpr,Coords_K1CM,Coords_K1C,
     +                             GTcu,GTpr,Gcu,Gpr,Grad_K1C,AtmMass,
     +                             Nreals)

      Implicit Double Precision(A-H,O-Z)

#include <mxatms.par>

      Dimension Vec1(3*Nreals),Cpi(3*Nreals),Cpr(3*Nreals)
      Dimension GTcu(3*Nreals),GTpr(3*Nreals)
      Dimension Gcu(3*Nreals),Gpr(3*Nreals)
      Dimension Vec2(3*Mxatms),Vec3(3*Mxatms)
      Dimension AtmMass(Nreals)
      Dimension Coords_K1CM(3*Nreals),Coords_K1C(3*Nreals)
      Dimension Grad_K1C(3*Nreals)
      Data One,Pt9/1.0D0, 1.0D-09/

C Vec1 - Cpi-Cur (Cur is the current point)
C Cpr  - previous point
c Cpi  - pivot point 
C Gcu  - current gradient 
C Gpr  - previous gradient 
C Gtpr - Gradient along previous vector Cpi-Cpr
C Gtcu - Gradient along current vector Cpi-Cur
C all these quantities are mass-weigthed 

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "@-Irc_interpolate: At the entry"
      Write(6,*) "The M.W. current tangent gradient"
      Write(6, "(3F17.13)") (GTcu(j),j=1,3*Nreals)
      Write(6,*) "The M. W. previous tangent gradient"
      Write(6, "(3F17.13)") (GTpr(j),j=1,3*Nreals)
      Write(6,*) "The M. W. pivot point"
      Write(6, "(3F17.13)") (Cpi(j),j=1,3*Nreals)
      Write(6,*) "The M. W. previous point"
      Write(6, "(3F17.13)") (Cpr(j),j=1,3*Nreals)
      Write(6,*) "The M. W. Cpi_Ccu point"
      Write(6, "(3F17.13)") (Vec1(j),j=1,3*Nreals)
      Write(6,*) "The M.W. current gradient"
      Write(6, "(3F17.13)") (Gcu(j),j=1,3*Nreals)
      Write(6,*) "The M. W. previous gradient"
      Write(6, "(3F17.13)") (Gpr(j),j=1,3*Nreals)
#endif
C Eqns. 12-15 in Gonzales and Schlegel, JCP,90,2154,1989.

      Acut  = 5.0D0
      Todeg = 180.0D0/Acos(-One)
      Acut  = Acut/Todeg 

      Do Iatom = 1, Nreals
         Ioff = 1 + 3*(Iatom - 1)
         Call Vec(Cpi(Ioff),Cpr(Ioff),Vec2(Ioff),0)
      Enddo
   
      Dlen_1  = Ddot(3*Nreals,Vec1,1,Vec1,1)
      Dlen_2  = Ddot(3*Nreals,Vec2,1,Vec2,1)
      Dlen_12 = Dsqrt(Dlen_1)*DSqrt(Dlen_2)
      Dlen    = Ddot(3*Nreals,Vec1,1,Vec2,1)
      
      Ang1 = Dlen/Dlen_12
      If (Abs(Ang1) .Gt. One) Ang1 = Ang1/Dabs(Ang1)

      Ang = Acos(Ang1)

      Do Iatom = 1, Nreals
         Ioff = 1 + 3*(Iatom - 1)
         Call Vec(Vec2(Ioff),Vec1(Ioff),Vec3(Ioff),0)
      Enddo

      G2prim = Ddot(3*Nreals,Vec3,1,GTpr,1)
      G1prim = Ddot(3*Nreals,Vec3,1,GTcu,1) 
      Gdiff  = (G2prim-G1prim)
    
      If (Dabs(Gdiff) .Gt. Pt9) Then
         Theta = Ang*G2prim/Gdiff 
      Else
         Write(6,"(3a)") " Division by zero occurs in",
     &                   " evaluating theta and hence", 
     &                   " interpolation is skipped."
     &                   
          Return
      Endif 

      If ((Dabs(Theta) .Gt. Acut) .AND. 
     &   (Dabs(Ang1-Theta)) .GT. Acut) Then
          Error = Theta*Todeg
          Write(6,"(a,F10.5,a,F10.5,a)") " The angle ", Error,
     &                                   " exceed the allowed",
     &                                   " limit: ",Acut*Todeg,
     &                                   " Interpolation is skipped"
          Return
      Endif

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(2a,2F10.5)") "@-IRC_interpolate: The angles theta",
     &                       " and theta': ",Theta,Ang
#endif 

      T1 = Cos(Theta) - Sin(Theta)*Cos(Ang)/Sin(Ang)
      T2 = Sin(Theta)/Sin(Ang)
      T3 = Theta/Ang

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,3F10.5)") "@-IRC_interpolate: The T1-3: ",
     &                       T1,T2,T3
#endif 

      Do Ixyz = 1, 3*Nreals
         Vec3(Ixyz) = T1*Vec2(Ixyz) + T2*Vec1(Ixyz)
         Gcu(Ixyz)  = Gpr(Ixyz)*(ONE-T3) + Gcu(Ixyz)*T3
      Enddo

      Call Dcopy(3*Nreals,Vec3,1,Vec1,1)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(2a)") "@-IRC_interpolate: The M. W. interpolted",
     &                " Cpi_Ccu point point"
      Write(6, "(3F17.13)") (Vec1(j),j=1,3*Nreals)
#endif

      Call Daxpy(3*Nreals,One,Cpi,1,Vec3,1)
      Call Dcopy(3*Nreals,Vec3,1,Coords_K1CM,1)

      Ioff = 0
      Do Iatom = 1, Nreals
         Sqrtmass = One/Dsqrt(AtmMass(Iatom))
            Do Ixyz = 1, 3
               Ioff = Ioff + 1
               Vec3(Ioff) = Vec3(Ioff)*Sqrtmass
               Grad_K1C(Ioff) = Gcu(Ioff)*Dsqrt(AtmMass(Iatom))
            Enddo
      Enddo

      Call Dcopy(3*Nreals,Vec3,1,Coords_K1C,1)

#ifdef _DEBUG_LVL0
      Write(6,"(2a)") "@-IRC_interpolate: Interpolated current points",
     &               " gradients"
      Write(6,"(a)")  "The cur. point (mass weighted)"
      Write(6, "(3F17.13)") (Coords_K1CM(i),i=1,3*Nreals)
      Write(6,"(a)")  "The cur. point (no. mass weight)"
      Write(6, "(3F17.13)") (Coords_K1C(i),i=1,3*Nreals)
      Write(6,"(a)")  "The cur. gradient (mass weighted)"
      Write(6, "(3F17.13)") (Gcu(i),i=1,3*Nreals)
      Write(6,"(a)")  "The cur. gradient (no mass weighing)"
      Write(6, "(3F17.13)") (Grad_K1C(i),i=1,3*Nreals)
#endif 

      Return
      End 

     
