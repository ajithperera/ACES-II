      Subroutine Process_natorbs(Natorb_f,Scr1,Scr2,Iscr1,Iaopop,Nmo,
     +                           Ndim,Ispin)

      Implicit Double Precision(A-H,O-Z)
      
      Double Precision Natorb_f

      Dimension Natorb_f(Nmo*Nmo)
      Dimension Scr1(Nmo,Nmo),Scr2(Nmo,Nmo)
      Dimension Iscr1(0:Nmo),Iaopop(Nmo)

      Character*1 Type(2)

      Data Done,DZero,Thres /1.0D0,0.0D0,1.0D-05/
      Data Type /"A","B"/

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"

      Call Dcopy(Nmo*Nmo,Natorb_f,1,Scr2,1)

#ifdef _DEBUG_LVLM
      Write(6,"(2a,1x,i2)") "The MO basis natural orbitals of",
     +                      " spin",Ispin
      call output(Scr2,1,Nmo,1,Nmo,Nmo,Nmo,1)
#endif 

      Call Getrec(20,"JOBARC","NUMBASIR",Nirrep,Iscr1(1))

      Ithru = 0
      Imin  = 0
      Imax  = 0
      Iscr1(0) = 1
      Do J = 1, Nirrep 
         Imin = Imin + Iscr1(j-1)
         Imax = Imax + Iscr1(j) 
         Do K = 1, Nmo
            Do I = Imin, Imax
               If (Dabs(Scr2(I,K)) .Gt. Thres)  Then 
                   Iaopop(K) = J
               Endif
            Enddo
         Enddo
      Enddo
       
#ifdef _DEBUG_LVLM
      Write(6,"(a)") "The sym reorder array for natural orbs: "
      Write(6,"(15(1x,I3))") (Iaopop(I),I=1,Nmo)
      Write(6,"(a,10(1x,I3))") "The AO population array   : ",
     + (Iscr1(I),I=1,Nirrep)
#endif 

      I1 = 0
      I2 = Iscr1(1)
      I3 = Iscr1(2)+I2
      I4 = Iscr1(3)+I3
      I5 = Iscr1(4)+I4
      I6 = Iscr1(5)+I5
      I7 = Iscr1(6)+I6
      I8 = Iscr1(7)+I7

      Do I = 1, Nmo
         If (Iaopop(I) .Eq. 1) Then
            I1 = I1 + 1
            Iscr1(I) =I1
         Elseif (Iaopop(I) .Eq. 2) Then
            I2 = I2 + 1
            Iscr1(I) =I2
         Elseif (Iaopop(I) .Eq. 3) Then
            I3 = I3 + 1
            Iscr1(I) =I3
         Elseif (Iaopop(I) .Eq. 4) Then
            I4 = I4 + 1
            Iscr1(I) =I4
         Elseif (Iaopop(I) .Eq. 5) Then
            I5 = I5 + 1
            Iscr1(I) =I5
         Elseif (Iaopop(I) .Eq. 6) Then
            I6 = I6 + 1
            Iscr1(I) =I6
         Elseif (Iaopop(I) .Eq. 7) Then
            I7 = I7 + 1
            Iscr1(I) =I7
         Elseif (Iaopop(I) .Eq. 8) Then
            I8 = I8 + 1
            Iscr1(I) =I8
         Endif 
      Enddo

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "The sym reorder array for natural orbs: "
      Write(6,"(15(1x,I3))") (Iscr1(I),I=1,Nmo)
#endif 
      Call Getrec(20,"JOBARC","NUMBASIR",Nirrep,Iaopop(1))
      Do Imo = 1, Nmo
         Idestn = Iscr1(Imo)
         Call Dcopy(Nmo,Scr2(1,Imo),1,Scr1(1,Idestn),1)
      Enddo

#ifdef _DEBUG_LVLM
      Write(6,"(2a,1x,i2)") "The ZMAT symmetry SCF eigevalue ordered",
     +                " natural orbitals of spin",Ispin
      call output(Scr1,1,Nmo,1,Nmo,Nmo,Nmo,1)
#endif 
      Iaopop(0) = 1
      Ioff      = 1 
      Joff      = 1
      Ibgn      = 0
      Iend      = 0
      Ndim      = 0
      Do I = 1, Nirrep
         Length = Iaopop(I)
         Ibgn   = Ibgn + Iaopop(I-1)
         Iend   = Iend + Iaopop(I)
         Do Imo = Ibgn, Iend
            Call Dcopy(Length,Scr1(Ioff,Imo),1,Natorb_f(Joff),1)
            Joff = Joff + Length
         Enddo
         Ioff = Ioff + Length 
         Ndim = Ndim + Length*Length 
      Enddo 

#ifdef _DEBUG_LVLM
      Ioff = 1 
      Do I = 1, Nirrep
         Write(6,*) 
         Write(6,"(2a,1x,i2,a,i2)") "The symmetry packed",
     +                " natural orbitals of spin",Ispin," irrep ",i
         Nrow = Iaopop(I)
         Ncol = Iaopop(I)
         call output(Natorb_f(Ioff),1,Nrow,1,Ncol,Nrow,Ncol,1)
         Ioff = Ioff + Nrow*Ncol
      Enddo 
#endif 
       Return
       End
