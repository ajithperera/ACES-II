      SUBROUTINE mksyv(nu,SYVEC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SYVEC(1)
      integer syvec,dirprd,a,b,c,d,pop,vrt,spop,svrt
      character*8 lb
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),N(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      common/sym/pop(8,2),vrt(8,2),ntaa(6)
      common/sumpop/spop(8),svrt(8)
      spop(1)=1
      svrt(1)=1
c      write(6,*)'numery kolejne'
      do 98 i=2,nirrep
         spop(i)=spop(i-1)+pop(i-1,1)
         svrt(i)=svrt(i-1)+vrt(i-1,1)
 98   continue
c      do 99 i=1,nirrep
c         write(6,*)spop(i),svrt(i)
c 99   continue
      ksyv=0
      do 101 ir=1,nirrep
c         write(6,*)'repr nr.',ir
         do 102 jr=1,nirrep
            kr=dirprd(ir,jr)            
c            write(6,*)'skladowe reprezentacje',kr,jr
c            write(6,*)'wymiary:',vrt(jr,1),vrt(kr,1)
            if(vrt(jr,1).eq.0.or.vrt(kr,1).eq.0)goto 102
c        write(6,*)'przedzialy wierszy :',svrt(kr),svrt(kr)+vrt(kr,1)-1
c         write(6,*)'przedzialy kolumn :',svrt(jr),svrt(jr)+vrt(jr,1)-1
         ilow =svrt(jr)
         ihigh=ilow+vrt(jr,1)-1
         jlow =svrt(kr)
         jhigh=jlow+vrt(kr,1)-1
c         write(6,*)'ilow,ihigh,jlow,jhigh:',ilow,ihigh,jlow,jhigh
         do 50 i=ilow,ihigh
            jlim=i
            if(jhigh.lt.i)jlim=jhigh
            do 51 j=jlow,jlim
               ksyv=ksyv+1
               syvec(ksyv)=(i-1)*nu+j
 51         continue
 50      continue
 102     continue
 101  continue
c      do 70 i=1,nu2
c         write(6,*)syvec(i)
 70   continue
      RETURN
      END
