      SUBROUTINE buildKSpot(V,Z1,screxc,scr,
     &                        valao,valgradao,totwt,
     &                        intnumradpts,max_angpts,ncount,
     &                        Ksa,F,
     &                        naobasfn,iuhf,ntotal)

      implicit none

      integer i,naobasfn,iuhf,ntotal
      integer intnumradpts,max_angpts,ncount
      integer Z1(naobasfn,2)

      double precision V(naobasfn,naobasfn,iuhf+1),
     &                 ksa(ntotal),F(ntotal),
     &                 screxc(naobasfn,naobasfn),
     &                 scr(naobasfn,naobasfn,iuhf+1),
     & valao(naobasfn,max_angpts,intnumradpts,ncount),
     & valgradao(naobasfn,max_angpts,intnumradpts,ncount,3),
     & totwt(ncount,intnumradpts,max_angpts)
       call dzero(ksa,ntotal) 
       call dzero(V,naobasfn*naobasfn*(1+iuhf))
       call integxc(V,z1,screxc,scr,.False.,
     &              valao,valgradao,totwt,
     &              intnumradpts,max_angpts,ncount,.false.)
c       write(*,*) ' let us prin V naosbaf=',naobasfn
c       call kkk(naobasfn*naobasfn,V)
       call mat_trans(0,2,v(1,1,1),ksa,0)
       call DAXPY(ntotal,1.0d0,ksa,1,F,1)
       do i=1,ntotal
c        write(*,*) i,ksa(i),F(i)
       end do
       return
       end 
