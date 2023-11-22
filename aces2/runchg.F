
      subroutine runchg
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "aces2.h"
#include "maxbasfn.par"
      DIMENSION ITOTREC(5),ITOTWRD(5)
      DIMENSION MOIO(10,500),MOIOWD(10,500),MOIOSZ(10,500),
     &          MOIODS(10,500),MOIOFL(10,500)
      DIMENSION NOIO(10,500),NOIOWD(10,500),NOIOSZ(10,500),
     &          NOIODS(10,500),NOIOFL(10,500)
      DIMENSION ISYTYP(2,500),isyty2(2,500)

#ifdef _TRACE_ACES2
      write(*,*) '@ACES2: Entered ', __FILE__
#endif

      call a2getrec(20,'JOBARC','TOTRECMO',5,ITOTREC)
      call a2getrec(20,'JOBARC','TOTWRDMO',5,ITOTWRD)
      call a2getrec(20,'JOBARC','NDROTVRT',1,ndrvrtt)

      if (fno) call resetfno(eval,nbas,'DM')

      call runit('mv JOBARC JOBARC_DM')
      call runit('mv JAINDX JAINDX_DM')
      call runit('mv JOBARC_AM JOBARC')
      call runit('mv JAINDX_AM JAINDX')

      ndrgeo=2

      call a2putrec(20,'JOBARC','TOTRECMO',5,ITOTREC)
      call a2putrec(20,'JOBARC','TOTWRDMO',5,ITOTWRD)
      call a2putrec(20,'JOBARC','NDROPGEO',1,NDRGEO)

      call runit('xvtran')
      call a2putrec(20,'JOBARC','NDROTVRT',1,ndrvrtt)
      call runit('xintprc')
    
      if (fno) call resetfno(eval,nbas,'AM')

      MOIOSIZ=5000

      CALL A2GETREC(20,'JOBARC','MOIOVEC',MOIOSIZ,MOIO)
      CALL A2GETREC(20,'JOBARC','MOIOWRD',MOIOSIZ,MOIOWD)
      CALL A2GETREC(20,'JOBARC','MOIOSIZ',MOIOSIZ,MOIOSZ)
      CALL A2GETREC(20,'JOBARC','MOIODIS',MOIOSIZ,MOIODS)
      CALL A2GETREC(20,'JOBARC','MOIOFIL',MOIOSIZ,MOIOFL)
      CALL A2GETREC(20,'JOBARC','ISYMTYP',1000,ISYTYP)

      call runit('mv JOBARC JOBARC_AM')
      call runit('mv JAINDX JAINDX_AM')
      call runit('mv JOBARC_DM JOBARC')
      call runit('mv JAINDX_DM JAINDX')

      CALL A2GETREC(20,'JOBARC','MOIOVEC',MOIOSIZ,NOIO)
      CALL A2GETREC(20,'JOBARC','MOIOWRD',MOIOSIZ,NOIOWD)
      CALL A2GETREC(20,'JOBARC','MOIOSIZ',MOIOSIZ,NOIOSZ)
      CALL A2GETREC(20,'JOBARC','MOIODIS',MOIOSIZ,NOIODS)
      CALL A2GETREC(20,'JOBARC','MOIOFIL',MOIOSIZ,NOIOFL)
      CALL A2GETREC(20,'JOBARC','ISYMTYP',1000,ISYTY2)

      if (ndrvrtt.eq.0) then
         istart=291
      else
         istart=201
      end if
      do i=istart,400
         isyty2(1,i)=isytyp(1,i)
         isyty2(2,i)=isytyp(2,i)
         do k=1,10
            noio(k,i)  =moio(k,i)
            noiowd(k,i)=moiowd(k,i)
            noiosz(k,i)=moiosz(k,i)
            noiods(k,i)=moiods(k,i)
            noiofl(k,i)=moiofl(k,i)
         end do
      end do
      do i=istart,400
         isytyp(1,i)=0
         isytyp(2,i)=0
         do k=1,10
            moio(k,i)  =0
            moiowd(k,i)=0
            moiosz(k,i)=0
            moiods(k,i)=0
            moiofl(k,i)=0
         end do
      end do

      CALL A2PUTREC(20,'JOBARC','MOIOVEC',MOIOSIZ,NOIO)
      CALL A2PUTREC(20,'JOBARC','MOIOWRD',MOIOSIZ,NOIOWD)
      CALL A2PUTREC(20,'JOBARC','MOIOSIZ',MOIOSIZ,NOIOSZ)
      CALL A2PUTREC(20,'JOBARC','MOIODIS',MOIOSIZ,NOIODS)
      CALL A2PUTREC(20,'JOBARC','MOIOFIL',MOIOSIZ,NOIOFL)
      CALL A2PUTREC(20,'JOBARC','ISYMTYP',1000,ISYTY2)

      call runit('mv JOBARC JOBARC_DM')
      call runit('mv JAINDX JAINDX_DM')
      call runit('mv JOBARC_AM JOBARC')
      call runit('mv JAINDX_AM JAINDX')

      CALL A2PUTREC(20,'JOBARC','MOIOVEC',MOIOSIZ,MOIO)
      CALL A2PUTREC(20,'JOBARC','MOIOWRD',MOIOSIZ,MOIOWD)
      CALL A2PUTREC(20,'JOBARC','MOIOSIZ',MOIOSIZ,MOIOSZ)
      CALL A2PUTREC(20,'JOBARC','MOIODIS',MOIOSIZ,MOIODS)
      CALL A2PUTREC(20,'JOBARC','MOIOFIL',MOIOSIZ,MOIOFL)
      CALL A2PUTREC(20,'JOBARC','ISYMTYP',1000,ISYTYP)

      call runit('mv JOBARC JOBARC_AM')
      call runit('mv JAINDX JAINDX_AM')
      call runit('mv JOBARC_DM JOBARC')
      call runit('mv JAINDX_DM JAINDX')

#ifdef _TRACE_ACES2
      write(*,*) '@ACES2: ===== end of runchg ====='
#endif

      return
      end

