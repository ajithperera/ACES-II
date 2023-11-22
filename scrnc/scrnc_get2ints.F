C*************************************************************
      subroutine Scrnc_get2ints(G2,buf,ibuf,Nbfns)
C     
C     Get two electron integrals
C         
C*************************************************************
      implicit double precision (a-h,o-z)
c
      logical FileExist
      integer and,or,nut,Nbfns,ifile,ldim2
      character*4 filename(4)

#include "maxbasfn.par"
c
      Dimension Nbfirr(8),Iao_sym(Maxbasfn) 
      common /machsp/ iintln,ifltln,iintfp,ialone,ibitwd

      dimension G2(Nbfns*Nbfns,Nbfns*Nbfns),buf(600),ibuf(600)
      Data half /0.50D0/
c
      iupki(int)=and(int,ialone)
      iupkj(int)=and(ishft(int,-ibitwd),ialone)
      iupkk(int)=and(ishft(int,-2*ibitwd),ialone)
      iupkl(int)=and(ishft(int,-3*ibitwd),ialone)
      ipack(i,j,k,l)=or(or(or(i,ishft(j,ibitwd)),ishft(k,2*ibitwd)),
     &                  ishft(l,3*ibitwd))
      indx(i,j)=j+(i*(i-1))/2
      indx2(i,j)=i+(j-1)*Nbfns 

      Call Getrec(20,"JOBARC","NIRREP ",1,Nirrep)
      Call Getrec(20,"JOBARC","NUMBASIR",Nirrep, NbfIrr)

      Index = 0
      Do Irrep=1,Nirrep
         Do Ibf = 1, Nbfirr(Irrep)
            Index = Index + 1
            Iao_sym(index) = Irrep
         Enddo 
      Enddo 
c

      filename(1)='IIII'
      filename(2)='IIJJ'
      filename(3)='IJIJ'
      filename(4)='IJKL'

      Call Dzero(G2,Nbfns**4)

      do 10 ifile=1,4
         FileExist=.false.
         inquire(file=filename(ifile),exist=FileExist)
         if (FileExist) then
            open(unit=ifile,file=filename(ifile),form='UNFORMATTED',
     &           access='SEQUENTIAL')
            rewind ifile
c     
            call locate(ifile,'TWOELSUP')
            nut = 1
            do while (nut.ge.1)
               read(ifile) buf, ibuf, nut
               do int = 1, nut
                  ix=iupki(ibuf(int))
                  jx=iupkj(ibuf(int))
                  kx=iupkk(ibuf(int))
                  lx=iupkl(ibuf(int))

                  i=max(ix,jx)
                  j=min(ix,jx)
                  k=max(kx,lx)
                  l=min(kx,lx)

                  ind1=indx(i,j)
                  ind2=indx(k,l)

                  if(ind1.lt.ind2) then
                     ix=i
                     jx=j
                     i=k
                     j=l
                     k=ix
                     l=jx
                  endif

                  x = buf(int)
#ifdef _DEBUG_LVL0
      write(6,"(a,1x,9i4,1x,F15.9)") "ifile,I,J,K,L,4irreps:", 
     +                                ifile,i,j,k,l,
     +                                Iao_sym(i),Iao_sym(j),Iao_sym(k),
     +                                Iao_sym(l),X
#endif 
                  G2(indx2(i,j),indx2(k,l))=x
                  G2(indx2(i,j),indx2(l,k))=x
                  G2(indx2(j,i),indx2(k,l))=x
                  G2(indx2(j,i),indx2(l,k))=x
                  G2(indx2(k,l),indx2(i,j))=x
                  G2(indx2(k,l),indx2(j,i))=x
                  G2(indx2(l,k),indx2(i,j))=x
                  G2(indx2(l,k),indx2(j,i))=x
               end do
            end do
            close(ifile)
         else
         end if   
 10   continue  

#ifdef _DEBUG
      Indj = 0
      Do k = 1, Nbfns
         Do l = 1, Nbfns
            Indi = 0
            Indj = Indj + 1
            Do j = 1, Nbfns
               Do i = 1, Nbfns
                  Indi = Indi + 1
#ifdef _DEBUG_LVL0
      write(6,"(a,1x,4i4,1x,F15.9)") "I,J,L,K:", i,j,l,k,G4(i,j,l,k)
#endif 
                  G2(Indi,Indj) = G4(i,j,l,k)
               Enddo
            Enddo 
         Enddo  
      Enddo
#endif 
c     
      return
      end
