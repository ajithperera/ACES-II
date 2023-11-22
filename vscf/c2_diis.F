      Subroutine C2_Diis(B_aug, B_bare, BVec, B_keep, Error, IOrder,
     &                   Ikeep)

      Implicit Double Precision(A-H, O-Z)

      Dimension B_aug(Iorder+1,Iorder+2)
      Dimension B_bare(Iorder,Iorder), Bvec(Iorder,Iorder)
      Dimension B_keep(Iorder,Iorder), Error(Iorder)
  
      Double Precision Normj
   
#ifdef _DEBUG_LVLM
      Write(*,*) "The agumented RPP-Matrix"
      Call output(B_aug,1,Iorder+1,1,Iorder+2,Iorder+1,Iorder+2,1)
#endif

      Do J = 2, IORDER+1
          Do I = 2, IORDER+1
         
            B_bare(I-1,J-1) = B_aug(I,J)
          
          Enddo
      Enddo

#ifdef _DEBUG_LVLM
      Write(*,*) "The Bare RPP-Matrix"
      Call output(B_bare,1,Iorder,1,Iorder,Iorder,Iorder,1)
#endif
      Call Dcopy(Iorder*Iorder, B_bare, 1, B_keep, 1)

      Call Eig(B_bare, Bvec, 0, Iorder, 1)

#ifdef _DEBUG_LVLM
      Write(*,*) "The eigen-vectors of the RPP-Matrix"
      Call output(Bvec, 1, Iorder, 1, Iorder, Iorder, Iorder,1)
      Write(*,*)
      Write(*,*) "The eigen-values of the RPP-Matrix"
      Call output(B_bare, 1, Iorder, 1, Iorder, Iorder, Iorder,1)
      Write(*,*)
#endif
      Do J = 1, Iorder
         Sumj = 0.0D0

         Do I = 1, Iorder
            Sumj = Sumj + Bvec(i,j)
         Enddo

         Normj = 1.0D0/(Sumj)

         Do I = 1, Iorder
            Bvec(I,j) = Normj * Bvec(i,j)
         Enddo 
      Enddo

#ifdef _DEBUG_LVLM
      Write(*,*) "The Normalized eigen-vectors of the RPP-Matrix"
      Call output(Bvec, 1, Iorder, 1, Iorder, Iorder, Iorder,1)
      Write(*,*)
#endif

      Do J = 1, Iorder
         Error(j) = 0.0D0
         DO I= 1, Iorder
            Error(J) = Error(J) + Dabs(Bvec(I,J) * B_keep(I,J))
          Enddo 
      Enddo
#ifdef _DEBUG_LVLM
      Write(*,*) "The Error vectors"
      Write(*,"(6(1x,F12.8))") (Error(J), j= 1, Iorder) 
      Write(*,*)
#endif
      
      Call Simple_sort(Error, Bvec, Iorder)
      Ikeep = 1
        
#ifdef _DEBUG_LVLM
      Write(*,*) "The sorted vectors"
      Write(*,"(6(1x,F12.8))") (Error(J), j= 1, Iorder) 
      Write(*,*)
#endif

      Return
      End

  
