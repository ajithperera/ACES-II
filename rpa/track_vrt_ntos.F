      Subroutine Track_vrt_ntos(CNtos,Ntos,Cmos,Dvv,Numb,Nbas,
     +                          ICOLL,IMAP_A,IMAP_B,Lunitn)

      Implicit Double Precision (A-H, O-Z)

#include "maxbasfn.par"    
#include "info.com"    

      Double Precision Ntos
      Integer Track

      Dimension CNtos(Nbas,Numb), Cmos(Nbas,Numb),Dvv(Numb,Numb)
      Dimension Ntos(Numb,Numb)
      Dimension IMAP_A(*),IMAP_B(*)
C     Dimension Ovlp(Maxbasfn,Maxbasfn), Track(Maxbasfn)
      Dimension Ovlp(Numb,Numb), Track(Numb)
#ifdef _DEBUG_LVLM
      Write(6,"(a)") "The NTOs"
      Call output(Ntos,1,Numb,1,Numb,Numb,Numb,1)
      Write(6,"(a)") "The NTO transformed HF eigenvectors"
      Call output(CNtos,1,Nbas,1,Numb,Nbas,Numb,1)
      Write(6,"(a)") "The HF eigenvectors"
      Call output(Cmos,1,Nbas,1,Numb,Nbas,Numb,1)
      Write(6,"(a)") "Lambda"
      Write(6,"(6(1X,F10.5))") (Dvv(i,i), i=1,Numb)
#endif 

      Do i = 1, Numb
         Do j = 1, Numb
C           Ovlp(j,i) = Ddot(Nbas, CNtos(1,i), 1, Cmos(1,j), 1)
            Ovlp(j,i) = Ntos(j,i)
         Enddo 
      Enddo 
      
#ifdef _DEBUG_LVLM
      Write(6, "(a)") "The maximum overlap"
      Write(6,"(5(1X,F10.5))") ((Ovlp(j,i), j=1,Numb),i=1,Numb)
#endif 
C
C For a given i (NTO) find the MO that has the largest overlap
C
      Do i = 1, Numb
         Test = 0.0D0
         Do j = 1, Numb 
            If (Dabs(Ovlp(j,i)) .Gt. Test) Then
                Test = Dabs(Ovlp(j,i))
C               Track(i) = j + Nocco(1)
                Track(i) = IMAP_A(ICOLL+j-1)
            Endif 
         Enddo
      Enddo 
#ifdef _DEBUG_LVLM
      Write(6, "(a)") "The NTO -> SCF MO map (overlap)"
      Write(6,"(5(1X,I3))") (Track(I), I=1,Numb)
#endif 
C     Do i = 1, Numb
C        Test = 0.0D0
C        Do j = 1, Numb
C            If (Dabs(Ntos(j,i)) .Gt. Test) Then
C                Test = Dabs(Ntos(j,i))
C                Track(i) = j + Nocco(1)
C                Track(i) = IMAP_A(ICOLL+j-1) 
C            Endif
C        Enddo
C     Enddo

#ifdef _DEBUG_LVLM
      Write(6, "(a)") "The NTO -> SCF MO map"
      Write(6,"(5(1X,I3))") (Track(I), I=1,Numb)
#endif 
      Do i = 1, Numb
         WRITE(LUNITN,'(A,I7,1X,F15.10)') "VRT ",Track(i),
     +                                     2.0D0*Dvv(i,i)
      Enddo
C
      Return
      End
