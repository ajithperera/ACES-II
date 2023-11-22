      subroutine initGuess(C,row,cols)
      integer, intent(in):: row,cols
      double precision,intent(inout):: C(row,cols)
      integer:: i
      Data Done /1.0D0/

!! Zero out entire matrix
!! Create Columns of identiy matrix according to the number of 
!! Initial Guess Vectors

      call dzero(C,row*cols)
      do i=1,cols
         C(i,i)=Done
      enddo
      
      return
      end 
