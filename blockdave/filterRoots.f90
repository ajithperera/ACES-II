#include "flags.h"
subroutine filterRoots(cVals,cEVECS,Nsize,Nblocks,NUMSOL,&
           CISfilterFlag,filterIndx)
  IMPLICIT INTEGER (A-Z)
  logical,intent(in)::CISfilterFlag
  integer,intent(in)::Nsize,Nblocks
  integer,intent(inout)::NUMSOL,filterIndx
  double precision,intent(inout)::cVals(Nblocks),cEVECS(Nsize,Nblocks)
  INTEGER::terminalIndex,qm1switch,qm2switch,counti,countj
  double precision::THRESH
  integer::jj,zz,SSPIN,ICOUNT,IIRREP,AIRREP,I,A
CHARACTER*12 STRINGI(2), STRINGA(2)
LOGICAL incQM1, FIRST
  double precision::SCR(Nsize),tmpEVAL
  double precision::bkupeVecs(NSize,Nblocks),bkupeVals(Nblocks)
  integer::QM1count
  logical::purge(Nblocks)


      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/ROOTS/EIGVAL(100,8),EIGVAL_T(100,8),OSCSTR(100,8),&
                  BGN(100,8),BGN_IRP(100,8),END(100,8),&
                  END_IRP(100,8),NATURE

      INTEGER::LUNITN
      integer*8, dimension(2) :: scratch
      integer ierr,natm,nbas,nocc,nvirt
      INTEGER QM1num,NLMOnum,NLMOnum2,NLMOorigin,NLMOorigin2,indx
      INTEGER,allocatable::QM1atoms(:),QM2atoms(:),NLMO(:,:)
      integer,allocatable:: NLMOQM1(:),NLMOQM2(:)
      integer::compareIA(2),QMregIA(2),compareIJ(2),QMregIJ(2)
#include "info.com"
#include "flags.com"
#include "flags2.com"
!! First, read in QM1 information
      call getrec(1,'JOBARC','NBASTOT',1,nbas)
      call getrec(1,'JOBARC','NREALATM',1,natm)
        call getrec(1,'JOBARC','NOCCORB',2,scratch)
      nocc=scratch(1)
      nvirt=nbas-nocc
      allocate(NLMOQM1(nbas),NLMOQM2(nbas))
      NLMOQM1=0
      NLMOQM2=0

      LUNITN = 30
      INQUIRE(FILE='QMcenter',EXIST=NLS_EXIST)
      if (NLS_EXIST) THEN
        open(unit=LUNITN,file='QMcenter',status='old',iostat=ierr)
        if (ierr.eq.0) then
          rewind LUNITN
          read(LUNITN,*) QM1num
          read(LUNITN,*)
          if (natm-QM1num .gt. 0) then
            allocate(QM1atoms(QM1num),QM2atoms(natm-QM1num))
          else 
            allocate(QM1atoms(QM1num))
          endif
          do i=1,QM1num
            read(LUNITN,*) site
            QM1atoms(i)=site
          enddo
          read(LUNITN,*)
          read(LUNITN,*) QM2num
          if (QM2num.gt.0) then
            do i=1,QM2num
              read(LUNITN,*) site
              QM2atoms(i)=site
            enddo
          endif
        else
          print*, "error opening QMcenter"
        endif
        close(LUNITN)  
      else
        print*,'FATAL ERROR: QMcenter-NLS input file-can not be found'
        call exit(1)
      endif


!******************************************************************
! * Read output of N. Flocke NLMO code 'nbocenters'
! * 'NLMO(x,y)' stores linked list of NLMOs associated with atom site
!               where x is the NLMO number and y are the atoms sharing
!               the NLMO
!******************************************************************
      INQUIRE(FILE='nbocenters',EXIST=NBO_EXIST)
      if (NBO_EXIST) THEN
        NBOfile=300
        allocate(NLMO(nbas,5))
        open(unit=NBOfile,file='nbocenters',status='old',iostat=ierr)
        rewind NBOfile
        read(NBOfile,30,IOSTAT=ierr) NLMOnum,NLMOorigin
        print*,'line1',NLMOnum,NLMOorigin 
        NLMO=100
        NLMO(NLMOnum,1)=NLMOorigin
        indx=1
        do while (ierr.eq.0)
           read(NBOfile,30,IOSTAT=ierr) NLMOnum2,NLMOorigin2
           print*,NLMOnum2,NLMOorigin2
           if (ierr.ne.0) exit
           if (NLMOnum.eq.NLMOnum2) then
              indx=indx+1
           else
              indx=1
           endif
           NLMO(NLMOnum2,indx)=NLMOorigin2
           NLMOnum=NLMOnum2
           NLMOorigin2=NLMOorigin
        enddo
        close(NBOfile)
      ELSE
        print*,'FATAL ERROR: nbocenters-output of N. Flocke NLMO code-&
                    can not be found'
        call exit(1)
      endif     






!******************************************************************
! * Determine if the NLMO read from nbocenters is entirely within
! * QM1, QM2, or both QM1&QM2
!******************************************************************
        counti=1
        countj=1
        do i=1,nbas
          terminalIndex=findloc(NLMO(i,1:5),value=100,dim=1)
          qm1switch=0
          qm2switch=0
          do j=1,terminalIndex-1
                if (any(NLMO(i,j).eq.QM1atoms)) then
                   qm1switch=1
                else if (any(NLMO(i,j).eq.QM2atoms)) then
                   qm2switch=1
                endif
          enddo
!         * Handles if NLMO is shared between QM1 and QM2
          if ((qm1switch.eq.1) .and. (qm2switch.eq.1)) then
                NLMOQM1(counti)=i
                counti=counti+1

                NLMOQM2(countj)=i
                countj=countj+1
!         * Only if NLMO is in QM1
          else if ((qm1switch.eq.1) .and. (qm2switch.eq.0)) then
                NLMOQM1(counti)=i
                counti=counti+1
!         * Only if NLMO is in QM2
          else if ((qm1switch.eq.0) .and. (qm2switch.eq.1)) then
                NLMOQM2(countj)=i
                countj=countj+1
          endif
        enddo
!#ifdef _DEBUG_LVL0
        print*
        print*,'NLMOQM1', NLMOQM1
        print*,'NLMOQM2',NLMOQM2

!! Then determine origin of largest magnitude of excitation;
!! if 'i' is not in QM1, omit that root
STRINGI(1)='I  [I_SYM]  '
STRINGI(2)='i  [i_SYM]  '
STRINGA(1)='A  [A_SYM]  '
STRINGA(2)='a  [a_SYM]  '
      THRESH = 0.05D0
bkupeVecs=0.0d0
bkupeVals=0.0d0
purge=.False.
   do zz=1,Nblocks
      print*,'Starting root', zz,'out of ',Nblocks
      ICOUNT=1
      DO SSPIN = 1, 1!1+IUHF
!        WRITE(6,1300I)
        PRINT*,'STARTING ICOUNT',ICOUNT
        IF (SSPIN .EQ. 1) THEN
          WRITE(6,*)'   SINGLE EXCITATION COEFFICIENTS AA'
        ELSE
          WRITE(6,*)'   SINGLE EXCITATION COEFFICIENTS BB'
        ENDIF
        FIRST =.TRUE.
!        CALL GETLST(SCR, zz, 1, 1, SSPIN, 94)
!        SCR=cEVECS(:,zz)!cEVECS(1+(SSPIN-1)*(NSize/2):(NSize/2)+(NSize/2)*(SSPIN-1),zz)
        if (.not.CISfilterFlag) then
          SCR=0.0d0
          call getlst(SCR,zz,1,1,1,497)
          call getlst(tmpEVAL,zz,1,1,1,95)
        endif
!        ICOUNT = 1
        PRINT*,SCR(1)
        DO IIRREP = 1, 1!NIRREP
          AIRREP=DIRPRD(IIRREP,1)
          print*,'airrep',AIRREP,POP(IIRREP, SSPIN)
          DO 1 I = 1, POP(IIRREP, SSPIN)
            DO 2 A = 1, VRT(AIRREP,SSPIN)
              IF (ABS(SCR(ICOUNT)).GT.THRESH) THEN
                print*,'significant coefficient found'
                IF (FIRST) THEN
                  WRITE(6,999) STRINGI(SSPIN), STRINGA(SSPIN)
  999             FORMAT(6X,A12,3X,A12)
                  WRITE(6,*)
                ENDIF
                FIRST =.FALSE.
!                IF (ABS(SCR(ICOUNT)).GT.0.10D0) THEN
                  print*,'significant coefficient found'
                  compareIA=(/ I,A+nocc /)
                  call findQMregion(compareIA,size(compareIA),NLMOQM1,&
                  size(NLMOQM1),NLMOQM2,size(NLMOQM2),QMregIA)

                  if (QMregIA(1).eq.2) then ! *i is outside QM1
                    if (CISfilterFlag) then
                      filterIndx=zz
                      RETURN
                    endif
                    purge(zz)=.True.
                    print*,'Purging root',zz
                  endif
!                ENDIF
                !print*,QMregIA
                WRITE(6, 1001) I, IIRREP, A, AIRREP, SCR(ICOUNT)
              ENDIF
              ICOUNT = ICOUNT + 1
    2       CONTINUE
    1     CONTINUE
        ENDDO
! 1001   FORMAT(3X,I4,3X, ' [',I1,'],  ',1X, I4,3X, ' [',I1,']   ;',&
!            10X, F12.6)
       enddo
       if (purge(zz)) then
         print*,'** Purging root:',zz
       endif
       print*,purge
       print*,'Ending root',zz,purge(zz),.not.(purge(zz))
! Both alpha and beta portion of vector inside QM1; so
! save this vector and value
       if (.not.(purge(zz))) then ! vector in QM1:: save
          print*,'Saving root',zz,QM1count+1
          PRINT*,'SAVING EVAL',tmpEVAL
          QM1count=QM1count+1
          bkupeVecs(:,QM1count)=SCR!cEVECS(:,zz)
          bkupeVals(QM1count)=tmpEVAL !cVals(zz)
          print*,'leaving if'
       endif
        print*,'outside of if'
  enddo

  print*,'** Number of vectors inside QM1:',QM1count
! Dump bkupeVecs and bkupeVals to disk
  do jj=1,QM1count
    call putlst(bkupeVals(jj),jj,1,1,1,95)
  enddo

  do jj=1,QM1count
    call putlst(bkupeVecs(:,jj),jj,1,1,1,497)
  enddo
!  cEVECS=bkupeVecs
!  cVals=bkupeVals
  NUMSOL=QM1count


30    FORMAT (I4,2X,I4)
 1001   FORMAT(3X,I4,3X, ' [',I1,'],  ',1X, I4,3X, ' [',I1,']   ;',&
            10X, F12.6)

end subroutine
