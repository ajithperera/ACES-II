      subroutine pccsd_modification(icore, maxcor, iuhf)
c     
c PCCSD modifications. 
c     
c     This subroutine is modelled after GENINT, and subsequent multiplication calls (see also vcc.f)
c     
c     Here we calculate most of the modifications due to parameterized CCSD
c     The one-particle intermediates are processed here: F(IM), F(AE), F(ME)
c     Also the paralpha 2-body contribution is calculated here.
c     Both the contributions to singles and doubles residuals are considered.
c     
c     The two-body parbeta contributions is included through a modification to 
c     dwmbej.F
c     
c     
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION TCPU,TSYS
      DOUBLE PRECISION C1,C2,C3,C4,FACT
C     
      DIMENSION ICORE(MAXCOR)
      COMMON/FLAGS/IFLAGS(100)
      COMMON/FLAGS2/IFLAGS2(500)
      COMMON/SYM/POP(8,2),VRT(8,2),NT1AA,NT1BB,NF1AA,NF1BB,
     &     NF2AA,NF2BB
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &     DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      EQUIVALENCE(METHOD,IFLAGS(2))
c     
#include "parcc.com"
c     
c     
c     these are the contributions to the doubles residuals as coded by Ondrej Demel
c     
      if (abs(paralpha -1.0d0) .gt. 1.0d-6) then
         write(6,*) ' Evaluate 2-body paralpha contributions '
         iii=5   
         jjj=21
c     
c     The subroutine quad1 creates the (paralpha-1) contribution to the wiiii * taaii term
c     
c     set the intermediate list used in quad1 explicitly to zero
c     
         if (iuhf .eq. 1) then
            ilow = 1
         else
            ilow = 3
         endif
         do ispin = ilow, 3
            listw = 50 + ispin
            NSIZE=0
            DO IRREP=1,NIRREP
               NUMSYM=IRPDPD(IRREP,ISYTYP(2,listw))
               DSZSYM=IRPDPD(IRREP,ISYTYP(1,listw))
               NSIZE=NSIZE+NUMSYM*DSZSYM
            ENDDO
            CALL ZERO(ICORE,NSIZE)
            CALL PUTALL(ICORE,NSIZE,1,listw)
         enddo
c
         call QUAD1(ICORE,MAXCOR,iii,IUHF, paralpha-1.0d0)
         call drlad(icore,maxcor,iuhf,jjj)
c     
c     set the intermediate list created in quad1 to zero
c     
         if (iuhf .eq. 1) then
            ilow = 1
         else
            ilow = 3
         endif
         do ispin = ilow, 3
            listw = 50 + ispin
            NSIZE=0
            DO IRREP=1,NIRREP
               NUMSYM=IRPDPD(IRREP,ISYTYP(2,listw))
               DSZSYM=IRPDPD(IRREP,ISYTYP(1,listw))
               nSZE=NSIZE+NUMSYM*DSZSYM
            ENDDO
            CALL ZERO(ICORE,NSIZE)
            CALL PUTALL(ICORE,NSIZE,1,listw)
         enddo
      endif
c     
c     calculate the remaining one-particle terms to the t2 residual
c     
      CALL ZERLST(ICORE,NF1AA,1,1,1,91)
      CALL ZERLST(ICORE,NF2AA,1,1,1,92)
      IF(IUHF.NE.0)THEN
         CALL ZERLST(ICORE,NF1BB,1,1,2,91)
         CALL ZERLST(ICORE,NF2BB,1,1,2,92)
      ENDIF
c     
C      if (abs(paralpha - 1.0d0) .gt. 1.0d-6) then
         write(6,*) 'Evaluate 1-body paralpha contributions '
         iii = 5                ! no tau formation
         call quad2(icore,maxcor,iuhf,iii,(paralpha-1.0d0)/2.0d0)
         call fmicont2(icore,maxcor,iuhf)
C      endif
     
C      if (abs(parbeta - 1.0d0) .gt. 1.0d-6) then
         write(6,*) 'Evaluate 1-body parbeta terms '
         iii = 5                ! no tau formation
         call quad3(icore,maxcor,iuhf,iii, parbeta - 1.0d0)
         call feacont2(icore,maxcor,iuhf)
C      endif
c     
      CALL ZERLST(ICORE,NF1AA,1,1,1,91)
      CALL ZERLST(ICORE,NF2AA,1,1,1,92)
      IF(IUHF.NE.0)THEN
         CALL ZERLST(ICORE,NF1BB,1,1,2,91)
         CALL ZERLST(ICORE,NF2BB,1,1,2,92)
      endif
c     
c     modification for singles equations (parepsilon parameter)
c     f(ea) contribution
c     
C      if (abs(parepsilon - 1.0d0) .gt. 1.0d-6) then
         write(6,*) 'Evaluate pccsd epsilon contributions  '
         CALL ZERLST(ICORE,NF2AA,1,1,1,92)
         IF(IUHF.NE.0)THEN
            CALL ZERLST(ICORE,NF2BB,1,1,2,92)
         endif
c     
         iii = 5                ! no tau
         call quad3(icore,maxcor,iuhf,iii, parepsilon - 1.0d0)
         call feacont1(icore, maxcor, iuhf)
c     
         CALL ZERLST(ICORE,NF2AA,1,1,1,92)
         IF(IUHF.NE.0)THEN
            CALL ZERLST(ICORE,NF2BB,1,1,2,92)
         endif
c     
         CALL ZERLST(ICORE,NT1AA,1,1,1,93)
         IF(IUHF.NE.0)THEN
            CALL ZERLST(ICORE,NT1BB,1,1,2,93)
         endif
         call makfme(icore,maxcor,iuhf,1)
         if (iuhf .ne. 0) call makfme(icore,maxcor,iuhf,2)
c     
c     scale contributions by (parepsilon - 1.0)
c     
         listf = 93
         do ispin = 1, iuhf+1
            nsizf =irpdpd(1,8+ispin)
            call getlst(icore,1,1,1,ispin,listf)
            call sscal(nsizf, parepsilon-1.0d0, icore, 1)
            call putlst(icore,1,1,1,ispin,listf)
         enddo
         call fmecont(icore, maxcor, iuhf, 1)
         if (iuhf .ne. 0) call fmecont(icore, maxcor, iuhf, 2)
c     
         CALL ZERLST(ICORE,NT1AA,1,1,1,93)
         IF(IUHF.NE.0)THEN
            CALL ZERLST(ICORE,NT1BB,1,1,2,93)
         endif
C      endif
c     
      return
      end
