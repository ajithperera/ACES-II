      Subroutine Assign_states_uhf(Coefs,Nature,Bgn,End,Bgn_irp,
     &                             End_irp,Iuhf,Irrepx,Iroot)

      Implicit Double Precision (A-H, O-Z) 

      Dimension Coefs(*)
      Integer Bgn(100,8),End(100,8),NBfirr(8),Bgn_irp(100,8),
     &        End_irp(100,8)
      
      Character*1 Nature(100,8)
      Logical Singlet, Triplet

#include "machsp.com"
#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Nalpha = 0
      Call Getrec(20,"JOBARC",'NUMBASIR',Nirrep, Nbfirr)

      Do Irepr = 1, Nirrep
         irepl = Dirprd(Irepr,Irrepx)
         Do I = 1, Pop(irepr,1)
            Do J = 1, vrt(irepl,1)
               Nalpha =  Nalpha + 1
            Enddo
         Enddo
      Enddo

      Triplet = .False.
      Singlet = .False.

      Thres1  = 1.0D-05
      Thres2  = 0.05D0
      Ncount = 1 
      Call Getlst(Coefs, 1, 1, 1, 1, 490)
      Call Getlst(Coefs(Nalpha+1), 1, 1, 1, 2, 490)

      Do Irepr = 1, Nirrep
         irepl = Dirprd(Irepr,Irrepx)
         Do I = 1, Pop(irepr,1)
            Do J = 1, vrt(irepl,1)
               If (Coefs(Ncount)        .ne. 0.0D0 .and. 
     &             Coefs(Ncount+Nalpha) .ne. 0.0D0) then

               Diff = Coefs(Ncount)/Coefs(Ncount+Nalpha)
#ifdef _DEBUG_LVLM
               Write(6,"(a,2(1x,F10.6))") "The A and B vecs: ", 
     &               Coefs(Ncount),Coefs(Ncount+Nalpha)
               Write(6,"(a,F10.6)") "The A/B ratio:", Dabs(Diff+1.0d0)
#endif 
               If (Dabs(1.0D0+Diff) .lt. Thres1) Then
                   Triplet = .True. 
               Else if (Dabs(1.0D0-Diff) .lt. Thres1) then 
                   Singlet = .True. 
               Endif 

               Endif 
               If (Singlet .or. Triplet)  Go To 10
               Ncount = Ncount + 1 

            Enddo 
         Enddo
      Enddo

 10   Continue
      
      If (.NOT. (Singlet .or. Triplet)) Then
          Nature(Iroot,Irrepx) = "U"
      ElseIf (Triplet) Then
          Nature(Iroot,Irrepx) = "T"
      ElseIf (Singlet) Then
          Nature(Iroot,Irrepx) = "S"
      Endif 
#ifdef _DEBUG_LVLM
      Write(6,"(a,2(1x,I3))") "The root number and Irrep: ", Iroot,
     &                         Irrepx
      Write(6,"(a,l,1x,l)") "Singlet or Triplet state? :", Singlet,
     &                       Triplet
#endif       
      Ncount = 1
      Weight_max = Thres2
      Do Irepr = 1, Nirrep
         irepl = Dirprd(Irepr,Irrepx)
         Do I = 1, Pop(irepr,1)
            Do J = 1, vrt(irepl,1)
               Weight1    = Dabs(Coefs(Ncount))
               Save       = Weight1
               Weight_max = max(Weight1,Weight_max) 

               If (Dabs(Save-Weight_max) .eq. 0.0D0) Then
                   Bgn(Iroot,Irrepx)       = i
                   Bgn_irp(Iroot,Irrepx)   = irepr
                   End(Iroot,Irrepx)       = j
                   End_irp(Iroot,Irrepx)   = irepl
               Endif  
               Ncount = Ncount + 1 
            Enddo
         Enddo
      Enddo

#ifdef _DEBUG_LVLM
      Write(6,"(a,1x,a)") "Singlet or triplet?:", Nature(Iroot,Irrepx)
      Write(6,"(a,i3)") "The origin      irrep:", Bgn_irp(Iroot,Irrepx)
      Write(6,"(a,i3)") "The destination irrep:", end_irp(Iroot,Irrepx)
      Write(6,"(a,i3)") "The origin        orb:", Bgn(Iroot,Irrepx)
      Write(6,"(a,i3)") "The destination   orb:", end_irp(Iroot,Irrepx)
#endif       
     
      Return
      End
