      Subroutine Set_atm_4shell(Atm_4shell,Ncenters,Nshells,
     &                          Iqmstr,Jcostr)

      Implicit Integer (A-Z)

      Dimension Atm_4Shell(Nshells)
      Dimension Iqmstr(Ncenters)
      Dimension Jcostr(Ncenters,Nshells)

      Nshll = 0
      Do Iat=1,Ncenters
         Do Ilq=1,Iqmstr(iat)
           Do Ijco=1,Jcostr(iat,ilq)
              Nshll=Nshll+1
              Atm_4shell(Nshll)=Iat
           Enddo
         Enddo
      Enddo

      Return
      End 


     
