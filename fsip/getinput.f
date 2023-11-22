   
      subroutine getinput(nactive,iactive,mxbas,IUHF)
      implicit integer (a-z)
      character*3 string
      dimension nactive(2),iactive(mxbas,2)
c
      open(unit=70,file='ZMAT',status='old',form='formatted')
1     read(70,'(A3)',end=900)string
      if(string.ne.'*fs')goto1
      read(70,*)nactive(1)
      if(iuhf.eq.0)then
        nactive(2)=nactive(1)
      else
        read(70,*)nactive(2)
      endif
      read(70,*)(iactive(i,1),i=1,nactive(1))
      do 11 i=1,nactive(1)
       ix=iactive(i,1)
       iactive(i,1)=map2norm(ix,1,'OCC',ijunk,ijunk2)
11    continue
      if(iuhf.eq.0)then
        call icopy(mxbas, iactive(1,1),1,iactive(1,2),1)
      elseif(iuhf.ne.0.and.nactive(2).ne.0)then
       read(70,*)(iactive(i,2),i=1,nactive(2))
       do 12 i=1,nactive(2)
        ix=iactive(i,2)
        iactive(i,2)=map2norm(ix,2,'OCC',ijunk,ijunk2)
12     continue
      endif
      close(unit=70,status='keep')
      return
900   write(6,1000)
1000  format(T3,'@GETINPUT-F, FS input card not found.')
      call errex
      end
