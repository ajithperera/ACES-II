program CMOL
implicit double precision (a-h,o-z)
!
!I am going to comment out all the lines which are
!in the nddo output for the initial guess routines. 
!
!
!
print*,'-----------------------------------------------'
print*,'|          Projected NDDO Guess               |'
print*,'|       WRITTEN BY DeCARLOS TAYLOR (2003)     |'
print*,'-----------------------------------------------'
print*,''
!write(*,*)' '
!write(*,*)' '
!write(*,*)'@Projected NDDO Initial Guess'
call aces_init_rte
call aces_ja_init
call cpusec(time1)
call readin
call cpusec(time2)
write(*,200)time2-time1
call aces_ja_fin 
!write(*,*)''
!write(*,*)''


200 format(' Total job time (sec) =  ',F10.3)
!write(*,*)''
!write(*,*)'Calculation completed. Program will stop'



end program cmol
