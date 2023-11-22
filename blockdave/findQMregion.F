        subroutine findQMregion(indx,indxSize,QM1,QM1Size,
     &                             QM2,QM2Size,QMregOut)!,output)
           integer,intent(in)::indxSize,QM1Size,QM2Size
           integer,intent(in)::indx(indxSize),QM1(QM1Size),
     &                                        QM2(QM2Size)
           integer, intent(inout)::QMregOut(indxSize)
 
           integer:: i,j
           logical:: finalOut(indxSize) 
           logical:: test(indxSize)
!******************************************************************
!******************************************************************
!  Purpose: receives the <aj||ib> indices in 'indx'. Iterates thru
!           each element of indx(4) at a time. For each element
!           ask:
!               1) is index in QM1? if no, then index in QM2
!               2) if yes, then ask is index also in QM2? 
!               3a)if no, then index is only in QM1
!               3b) if yes, then index is in both QM1&QM2
!
!
!           Returns this logic in 'QMregOut(4)'. Each element of 
!           array holds either:
!                       1 - index only in QM1
!                       2 - index only in QM2
!                       3 - index in both QM1&QM2
!
!******************************************************************
!******************************************************************
          QMregOut=0
!          print*, 'QM1', QM1
!          print*
!          print*, 'QM2', QM2
!          print*
          do i=1,indxSize
                if (any(indx(i).eq.QM1)) then
                   if (any(indx(i).eq.QM2)) then
                      QMregOut(i)=3
                   else
                      QMregOut(i)=1
                   endif
                else
                   QMregOut(i)=2
                endif
          enddo
!          print*,'**************************************'
!          print*,'** INside findQM  **'
!          print*,'**************************************'
!          print*,'**************************************'
        end subroutine 
