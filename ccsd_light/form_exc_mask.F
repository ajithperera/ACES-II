      Subroutine Form_exc_mask(P,Length,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension P(Length)
      Integer I,J,A,B,Inill 
      Double Precision Dnill,One 
      Data Ione,One,INill,DNill/1,1.0D0,0,0.0D0/

#include "ccsdlight_vars.com"

C This is a simple minded way to implement an active space scheme
C within the fractional occupation CCD/CCSD codes. This simply 
C project out (set to zero) non-active amplitudes. One could argues
C that this fairly inefficient since all the sums are still intact.
C this is sufficient to conduct preliminary tests.

      Do I = 1, Length
         P(I) = One
      Enddo 
     
C AA singles block       

      Index = INill
      Do I = 1, Nocc_a
      Do A = 1, Nvrt_a
        Index = Index + Ione 
        If ((I .Ge. Act_min_a .And. I .Le. Act_max_a) .And.
     +      (A .Gt. Act_max_a)) Then
           P(Index) = Dnill
        Endif
        If ((I .Lt. Act_min_a) .And. 
     +      (A .Ge. Act_min_a .And. A .Le. Act_max_a)) Then
           P(Index) = Dnill
        Endif
      Enddo 
      Enddo 
       
C BB singles block       

      Do i = 1, Nocc_b
      Do a = 1, Nvrt_b
        Index = Index + Ione 
        If ((i .Ge. Act_min_b .And. i .Le. Act_max_b) .And.
     +      (a .Gt. Act_max_b)) Then
           P(Index) = Dnill
        Endif
        If ((i .Lt. Act_min_b) .And.
     +      (a .Ge. Act_min_b  .And. a .Le. Act_max_b)) Then
           P(Index) = Dnill
        Endif
      Enddo 
      Enddo 

C AAAA doubles block 
      
      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         Index = Index + Ione 
         If ((I .Lt. Act_min_a .Or.  J .Lt. Act_min_a)   .And. 
     +      ((A .Ge. Act_min_a .And. A .Le. Act_max_a)   .Or.
     +       (B .Ge. Act_min_a .And. B .Le. Act_max_a))) Then
              P(Index) = Dnill
         Endif
         If ((I .Ge. Act_min_a  .And.  I .Le. Act_max_a   .Or.  
     +        J .Ge. Act_min_a  .And.  J .Le. Act_max_a)  .And. 
     +       (A .Gt. Act_max_a  .Or.   B .Gt. Act_max_a)) Then
              P(Index) = Dnill
         Endif
      Enddo
      Enddo 
      Enddo
      Enddo 

C BBBB doubles block 

      Do j = 1, Nocc_b
      Do i = 1, Nocc_b
      Do b = 1, Nvrt_b
      Do a = 1, Nvrt_b
         Index = Index + Ione 
         If ((i .Lt. Act_min_b  .Or. j .Lt. Act_min_b)   .And.
     +      ((a .Ge. Act_min_b .And. a .Le. Act_max_b)   .Or.
     +       (b .Ge. Act_min_b .And. b .Le. Act_max_b))) Then
              P(Index) = Dnill
         Endif
         If ((i .Ge. Act_min_b  .And. i .Le. Act_max_b  .And.
     +        J .Ge. Act_min_b  .And. J .Le. Act_max_b) .And. 
     +       (a .Gt. Act_max_b  .Or.  b .Gt. Act_max_b)) Then
             P(Index) = Dnill
        Endif
      Enddo
      Enddo 
      Enddo
      Enddo 

C ABAB doubles block 

      Do j = 1, Nocc_b
      Do I = 1, Nocc_a
      Do b = 1, Nvrt_b
      Do A = 1, Nvrt_a
         Index = Index + Ione 
         If ((I .Lt. Act_min_a  .Or. j .Lt. Act_min_b)   .And.
     +      ((A .Ge. Act_min_a .And. A .Le. Act_max_a)   .Or.
     +       (b .Ge. Act_min_b .And. b .Le. Act_max_b))) Then
              P(Index) = Dnill
         Endif
         If ((I .Gt. Act_min_a  .And. I .Lt. Act_max_a   .And.
     +        j .Ge. Act_min_b  .And. j .Le. Act_max_b)  .And. 
     +       (A .Gt. Act_max_a  .Or.  b .Gt. Act_max_b)) Then
             P(Index) = Dnill
         Endif
      Enddo
      Enddo 
      Enddo
      Enddo 

#ifdef _DEBUG_LVL0
      Sum = 0.0D0
      Do i=1,Length
         Sum = Sum + P(i)*P(i)
      Enddo 
      Write(6,"(a,1x,2(1x,ES8.2E2))") " @-Form_exc_mask,length,P*P = ",
     +                                 Dble(Length),Sum
#endif 
      Return
      End
