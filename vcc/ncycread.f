      subroutine ncycread
c
c the motivation for this change is to avoid situation that a CC job crashes
c because of too lw max number of iterations was specified
c if the user observes this is going to happen, (s)he can signal the xvcc 
c process to change the value
c
cjp this is called from a signal handler and it reads
cjp the ncycle value from a file
c
cjp another possibility is to use unix IPC message queues or socket communication
cjp for socket version, note, that in contrast to typical daemons, we must
cjp allow more xaces2 processes running, so they cannot all listen() at the same
cjp port, some dymical portmapping should be used or another port allocation scheme
c
cjp this is just an easy way, how to get this functionality
cjp - allow the user to change max number of CC iterations when the program
cjp is already running
cjp the user prepares a file in the working directory of the job and  issus the command
cjp kill -USR1 <PID of xvcc>
c
      implicit none
      integer iflags,ncycle
      COMMON /FLAGS/  IFLAGS(100)
      common /ncycle/ncycle
cjp
      open(unit=1,file='ncycle',form='formatted',
     + status='old',err=9)
      read(1,*) ncycle
      close(1)
      if(ncycle.gt.0) then
        write(6,*) '@VCC-I CC_MAXCYC changed to ',ncycle
        IFLAGS(7)=ncycle
        CALL PUTREC(20,'JOBARC','IFLAGS  ',100,IFLAGS)
      else
        NCYCLE=IFLAGS(7)
      endif
9     continue
      return
      end
