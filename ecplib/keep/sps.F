      Subroutine Sps(Nlq, Lprj, La, Lb, Ltot, Lamalo, Lamahi, Lamblo, 
     &                 Lambhi, Lamau, Lambu, Alpha2, Rc, Rc2, Beta1, 
     &                 Beta2, Rad2)

#include "ecp.par"
#include "ecp.com"
      
      Dimension Rad2(0:2*Maxang,0:2*Maxang,0:2*maxang)
C
      Talpha = 2.0D0*Alpha2

      If (Beta1 .GT. Beta2) Then
C
         Alpha   = 1.0D0
         Sqalpha = Dsqrt(Alpha2)
         Beta    = Beta1/Sqalpha
         Xval    = Beta**2/4.0D0
         Tka     = Beta2*Beta2/Talpha

         Do Lama = Lamalo, Lamahi

            Na_min = Iabs(Lprj-Lama)

            Do Lamb = Lamblo, Lambhi

               Nb_min = Iabs(Lprj-Lamb) 

               Nlo = Na_min + Nb_min
               Nhi = (Ltot - Mod(La-Na_min,2)) - Mod(Lb-Nb_min,2)

               If (Nlo .EQ. Nhi .AND. Tka .EQ. 0) Then

                   Nstart = Nlq+NLo+Lama
                   Qj0 = Qcomp(Nstart, Lama)/Fac2(2*Lamb+1)
                   Fctr(Nlo) = 1.0D0
                   Term(Nlo) = Qj0
                   Sum(Nlo)  = Qj0

                   Do N = Nlo, Nhi, 2

                      If (Lamb .Eq. 0) Then
                         Prefac = Prd/Sqalpha**(Nlq+N+La+1)
                      Else
                         Prefac = Prd*Betaa**Lamb/Sqalpha**(Nlq+N+Lb+1)
                      Endif
C
                      Rad2(N, Lama, Lamb) = Rad2(N, Lama, Lamb) +
     &                                      Prefac*Sum(N)
                    Enddo

               Else
                    Nprime = Nlq + Nlo + Lama + 1
                    Qj1 = Qcomp(Nstart, Lama)/Fac2(2*Lamb+1)
                   
                    IF (NLo .Eq. Nhi) Then

                        Nlo2      = Nlo + 2
                        Sum(Nlo2) = Qj0

                        If (Nlo2 .Eq. Nhi .AND. Tka .EQ. 0 Then
                           
                           If (Lamb .Eq. 0) Then
                              Prefac = Prd/Sqalpha**(Nlq+N+La+1)
                           Else
                              Prefac = Prd*Betaa**Lamb/Sqalpha**(Nlq+N+Lb+1)
                           Endif
C
                           Rad2(N, Lama, Lamb) = Rad2(N, Lama, Lamb) +
     &                                      Prefac*Sum(N)
                        Enddo

                    Else

                        Fctr(Nlo+2) = Fctr(Nlo)
                        F1          = (2*Lama + 3) 
                        Fctr(Nlo)   = Tk/F1
                        Term(Nlo)   = Fctr(Nlo)*Qj1
                        Sum(Nlo)    = Sum(Nlo) + Term(Nlo)
                        Qj0         = Fctr(Nlo)*Qj0
                        Qj1         = Term(Nlo)

               Endif
            Enddo
         Enddo

      Else

         Alpha   = 1.0D0
         Sqalpha = Dsqrt(Alpha2)
         Beta    = Beta2/Sqalpha
         Xval    = Beta**2/4.0D0
         Tka     = Beta1*Beta1/Talpha
      Do Lama = Lamalo, Lamahi
        
         Do lamb = Lamblo, Lambhi
     
            If (Tka .EQ. 0) Then

                Nstart = Nlq+N+Lama
                Qj0 = Qcomp(Nstart, Lamb)/Fac2(2*Lama+1)
                Fctr(N) = 1.0D0
                Term(N) = Qj0
                Sum(N)  = Qj0

                Do N = 0, Ltot 

                   If (Lama .Eq. 0) Then
                      Prefac = xxx*/Sqalpha**(Nlq+N+La+1)
                   Else
                      Prefac = xxx*Betaa**Lama/Sqalpha**(Nlq+N+La+1)
                   Endif
C
                   Rad2(N, Lama, Lamb) = Rad2(N, Lama, Lamb) + 
     &                                   Prefac*Sum(N)
                Enddo

            Else

                Do  N = 0, Ltot 
                    J = 0
                    Do while (Term .Lt. ) 
 
C Evaluate J=0 term Eqn. 51 of MD paper.
                
                    Nstart = Nlq+N+Lama 
                    Qj0 = Qcomp_sps(Nstart, Lamb)/Fac2(2*Lama+1)
                    Fctr(J) = 1.0D0
                    Term(J) = Qj0
                    Sum(J)  = Qj0

C Evaluate J=1 term
                     Nstart = Nstart + 2 
                     Qj1 = Qcom_sps(Nstart, Lamb)/Fac2(2*Lama+3)
                     Fctr(J+1) = Tka
                     Term(J+1) = Tka*Qj1
                     Sum (J+1) = Sum(N) + Term(N)


                Nstart = Nstart + 2
                    F1 = (2*Nstart - 5)/2.0D0
                    F2 = (Lamb - Nstart -4)*(Lamb+Nstart-3)/4.0D0
                    Qjn= ((Xval + F1)*Qj1 + F2*Qj0))
    
                Do while (Term .Lt. )

  
            
                Enddo
            Enddo
          Enddo
      Enddo

      Return
      End

