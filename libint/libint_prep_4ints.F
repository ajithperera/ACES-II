      Subroutine Libint_prep_4ints(Angmom,Ncfs,Npfs,Pcoef,Ioff_pcoef,
     &                             Nfps,Index_ccf,End_nfps,Cc_begin,
     &                             CC_end,Nshells,Npcoef,Ntfct,
     &                             Cartesian,Spherical)

      Implicit Integer (A-Z)

      Dimension Npfs(Nshells),Ncfs(Nshells),Nfps(Nshells)
      Dimension Cc_begin(Ntfct),Cc_end(Ntfct)
      Dimension Index_ccf(Nshells),End_nfps(Nshells)
      Dimension Ioff_Pcoef(Nshells)
      Dimension Angmom(Nshells)

      Double Precision Pcoef(Npcoef)
    
      Logical Cartesian, Spherical 

      Data Dzero /0.0D0/

CSSS      Print*, Angmom
CSSS      Print*, Ncfs 
CSSS      Print*, Ioff_pcoef
CSSS      Print*, Nshells, Npcoef,Cartesian,Spherical
CSSS      Print*, Pcoef
   
      If (Cartesian) Then
         Do Ishell = 1, Nshells 
            Degenaracy = (Angmom(Ishell)+1)*(Angmom(Ishell)+2)/2 
            Nfps(Ishell) = Degenarcy*Ncfs(Ishell)
         Enddo
      Endif 

      If (Spherical) Then
         Do Ishell = 1, Nshells
            Degenarcy = (2*Angmom(Ishell)+1)
            Nfps(Ishell) = Degenarcy*Ncfs(Ishell)
         Enddo
      Endif
      
      Icpount = 0
      Do Ishell = 1, Nshells
         End_Nfps(Ishell) = Icpount + Nfps(Ishell)
         Icpount           = End_Nfps(Ishell)
      Enddo 

      Icc = 1
      Do Ishell = 1, Nshells
         Index_ccf(Ishell) = Icc
         Do K = 1, Ncfs(Ishell)
            l = (k-1)*Npfs(Ishell) + 1
            Do I = 1, Npfs(Ishell) 
               If (Pcoef(Ioff_pcoef(Ishell)+l-1) .Ne. Dzero) Then
                  CC_begin(Icc) = I
                  Go to 10
               Endif
               l = l + 1
            Enddo
 10   Continue 
            l = k*Npfs(Ishell) 
            Do I = Npfs(Ishell), 1, -1 
               If (Pcoef(Ioff_pcoef(Ishell)+l-1) .Ne. Dzero) Then
                  CC_end(Icc) = I
                  Go to 20
               Endif 
               l = l - 1
            Enddo
 20   Continue 
            Icc = Icc + 1
         Enddo
      Enddo 

#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,"(a)") " The Nfps array"
      Write(6,"(6(1x,I4))") (Nfps(i),    i=1, Nshells)
      Write(6,"(a)") " The End_Nfps array"
      Write(6,"(6(1x,I4))") (End_Nfps(i),i=1, Nshells)
      Write(6,"(a)") " The Index_ccf array"
      Write(6,"(6(1x,I4))") (Index_ccf(i),i=1, Nshells)
      Write(6,"(a)") " The Cc_begin array"
      Write(6,"(16(1x,I4))") (Cc_begin(i),i=1, Ntfct)
      Write(6,"(a)") " The Cc_end array"
      Write(6,"(16(1x,I4))") (Cc_end(i),  i=1, Ntfct)
#endif 

      Return
      End
