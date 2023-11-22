      SUBROUTINE CONVCRIT(FOCK,DENS,SOVRLP,SCR1,SCR2,SCRTMP,LDIM,
     &                    NBAS,IUHF,ERRMAX)
      implicit none
C Common blocks
#include "symm2.com"
#include "machsp.com"
      integer iflags(100)
      common/flags/iflags
C Input variables
      integer iuhf,nbas,ldim
      double precision fock((iuhf+1)*ldim),dens((iuhf+1)*ldim)
C Output variables
      double precision errmax
C Pre-allocated local variables
      double precision sovrlp(nbas*nbas),scr1(isqrln(nirrep+1)),
     &                 scr2(mxirr2),scrtmp(mxirr2)
C Local variables
      logical rohf
      integer i,ispin,len,jj,kk
      double precision dev,devmax(2),one,onem,zilch
      data one,onem,zilch /1.0d0,-1.0d0,0.0d0/
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C This actually needs to use the same convergence criterion as in RPP 
      CALL GETREC(20,'JOBARC','AOOVRLAP',NBAS*NBAS*IINTFP,SOVRLP)
      rohf=iflags(11).eq.2
      if (rohf) then
        CALL GETREC(20,'JOBARC','MOVECTOR',ISQRLN(NIRREP+1)*IINTFP,
     &             SCR1)
        call FAO2MOROHF(FOCK,SOVRLP,SCR1,SCR2,SCRTMP,NBAS,.TRUE.)
      endif
      devmax(1)=zilch
      devmax(2)=zilch
      do 100 i=1,nirrep
        if (nbfirr(i).gt.0) then
          len=nbfirr(i)
          do 101 ispin=1,iuhf+1
              if (rohf) then
                CALL EXPND2(FOCK(ITRIOF(I)),SCR1,len)
              else
                CALL EXPND2(FOCK((ISPIN-1)*LDIM+ITRIOF(I)),SCR1,len)
              endif
              CALL EXPND2(DENS((ISPIN-1)*LDIM+ITRIOF(I)),SCR2,len)
C Form FOCK * DENS {F*D}
              call xgemm('n','n',len,len,len,one,scr1,len,scr2,len,
     &                   zilch,scrtmp,len)
C Form FOCK * DENS * SOVRLP {F*D*S}
              CALL GETBLK(SOVRLP,SCR1,len,NBAS,IREPS(I))
              call xgemm('n','n',len,len,len,one,scrtmp,len,scr1,len,
     &                   zilch,scr2,len)
C S*D*F = (F*D*S)^T, use that to calculate Err=F*D*S-S*D*F
              devmax(ispin)=zilch
              do jj=1,len
                do kk=1,jj-1
                  dev=abs(scr2((jj-1)*len+kk)-scr2((kk-1)*len+jj))
                  if (dev.gt.devmax(ispin)) devmax(ispin)=dev
                end do
              end do
101       continue
          dev=devmax(1)+(iuhf-1)*devmax(2)
          if (i.eq.1) then
            errmax=dev
          else if (dev.gt.errmax) then
            errmax=dev
          endif
        endif
100   continue
      if (rohf) then
        CALL GETREC(20,'JOBARC','MOVECTOR',ISQRLN(NIRREP+1)*IINTFP,
     &             SCR1)
        call FAO2MOROHF(FOCK,SOVRLP,SCR1,SCR2,SCRTMP,NBAS,.FALSE.)
      endif
      return
      end
