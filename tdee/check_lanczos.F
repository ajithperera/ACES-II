      Subroutine check_lanczos(Work,Memleft,Nsize,irrepx,Nlim)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Memleft), Overlap(Nlim,Nlim)

      I000 = 1
      I010 = I000 + Nsize
      I010 = I010 + Nsize

      Do I=1, Nlim
      Do J=1, Nlim

         Call Getlst(Work(I000),I,1,0,Irrepx,497)
         Call Getlst(Work(I010),J,1,0,Irrepx,498)
         Overlap(I,J) = Ddot(Nsize,Work(I000),1,Work(I010),1)

      Enddo
      Enddo

      Write(6,"(a)") "The <R|L> of Lanczos vectors" 
      call output(Overlap,1,Nlim,1,Nlim,Nlim,Nlim,1)

      Return
      End 
