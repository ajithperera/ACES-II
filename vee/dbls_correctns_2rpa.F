      Subroutine Dbls_correctns_2rpa(Work,Length,Iuhf,Rcc_vecl,
     &                               Len_ph_aa,Len_ph_bb,Len_ph_pq,
     &                               E2_rcc_aa,E2_rcc_bb,Listz)
   
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Length),Rcc_vecl(Len_ph_pq)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"
        
C Do the dot product with the left EOM(SF)-CCD vectors to
C compute the second-order correcion.

#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,"(a)") "--------Enter dbls_correctns_2rpa-------"
      Write(6,*) "The incoming vectors" 
      Write(6,"(4(1x,F15.10))") (Rcc_vecl(i),i=1,Len_ph_pq)
#endif 
      I000 = 1
      Ioff = 1
      Do Ispin = 1, Iuhf+1
         If (Ispin .Eq. 1) Then
            I010 = I000 + Len_ph_aa
         Elseif (Ispin .Eq. 2) Then 
            I010 = I000 + Len_ph_bb
         Endif 

         Call Getlst(Work,1,1,1,2+Ispin,Listz)
         Ioff = Ioff + Len_ph_aa * (Ispin - 1)
#ifdef _DEBUG_LVL0
         If (Ispin .eq. 1) call checksum("CR(a,i)",Work,Len_ph_aa,s)
         If (Ispin .eq. 2) call checksum("CR(b,j)",Work,Len_ph_bb,s)
#endif
         If (Ispin .Eq. 1) Then

C This was done to compare with the UHF for debugging purposes.
C            If (Iuhf .EQ. 0) Call Dscal(Len_ph_aa,1.00/Dsqrt(2.0D0),
C     +                                  Rcc_vecl(Ioff),1)

#ifdef _DEBUG_LVL0
            call checksum("CL(a,i)",Rcc_vecl(Ioff),Len_ph_aa,s)
#endif
            E2_rcc_aa = Ddot(Len_ph_aa,Work(I000),1,Rcc_vecl(Ioff),1)

         Else if (Ispin .Eq. 2) Then
#ifdef _DEBUG_LVL0
       call checksum("CL(b,j)",Rcc_vecl(Ioff),Len_ph_bb,s)
#endif
            E2_rcc_bb = Ddot(Len_ph_bb,Work(I000),1,Rcc_vecl(Ioff),1)

         Endif
      Enddo 

      Return
      End
