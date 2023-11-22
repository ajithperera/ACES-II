       subroutine traSGC(sgrad,cgrad,tra,ncoor)
       implicit none
       integer ncoor,i,j
       double precision sgrad(ncoor),cgrad(ncoor),tra(ncoor,ncoor)
       double precision one,azero
       data azero,one /0.d0,1.d0/
c       DO 299 I = 1, NCOOR
c        WRITE (*,'(9F6.2)')(TRA(I,J),J=1,NCOOR)
c 299   CONTINUE  

      call xgemm('n','n',ncoor,1,ncoor,one,
     *           TRa,NCOOR,sgrad,ncoor,azero,
     *           cgrad,ncoor)
      RETURN
      END  
