      Subroutine Track_occ_ntos(CNtos,Ntos,Cmos,Doo,Numi,Nbas,
     +                          ICOLL,IMAP_A,IMAP_B,Lunitn)

      Implicit Double Precision (A-H, O-Z)

#include "maxbasfn.par"    

      Double Precision Ntos
      Integer Track

      Dimension CNtos(Nbas,Numi), Cmos(Nbas,Numi),Doo(Numi,Numi)
      Dimension Ntos(Numi,Numi)
      Dimension IMAP_A(*),IMAP_B(*)
C     Dimension Ovlp(Maxbasfn,Maxbasfn), Track(Maxbasfn)
      Dimension Ovlp(Numi,Numi), Track(Numi)
#ifdef _DEBUG_LVLM
      Write(6,"(a)") "The NTOs"
      Call output(Ntos,1,Numi,1,Numi,Numi,Numi,1)
      Write(6,"(a)") "The NTO transformed HF eigenvectors"
      Call output(CNtos,1,Nbas,1,Numi,Nbas,Numi,1)
      Write(6,"(a)") "The HF eigenvectors"
      Call output(Cmos,1,Nbas,1,Numi,Nbas,Numi,1)
      Write(6,"(a)") "Lambda"
      Write(6,"(6(1X,F10.5))") (Doo(i,i), i=1,Numi)
#endif 

      Do i = 1, Numi
         Do j = 1, Numi
C           Ovlp(j,i) = Ddot(Nbas, CNtos(1,i), 1, Cmos(1,j), 1)
            Ovlp(j,i) = Ntos(j,i)
         Enddo 
      Enddo 
      
#ifdef _DEBUG_LVLM
      Write(6, "(a)") "The maximum overlap"
      Write(6,"(5(1X,F10.5))") ((Ovlp(j,i), j=1,Numi),i=1,Numi)
#endif 
C
C For a given i (NTO) find the MO that has the largest overlap
C
      Do i = 1, Numi
         Test = 0.0D0
         Do j = 1, Numi 
             If (Dabs(Ovlp(j,i)) .Gt. Test) Then
                 Test = Dabs(Ovlp(j,i))
C                Track(i) = j
                 Track(i) = IMAP_A(ICOLL+j-1)
             Endif 
         Enddo
      Enddo 
#ifdef _DEBUG_LVLM
      Write(6, "(a)") "The NTO -> SCF MO map"
      Write(6,"(5(1X,I3))") (Track(I), I=1,Numi)
#endif 
      Do i = 1, Numi
         WRITE(LUNITN,'(A,I7,1X,F15.10)') "OCC ",Track(i),
     +                                     2.0D0*Doo(i,i)
      Enddo 
C
      Return
      End
