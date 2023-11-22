subroutine NLMOcorrect(cVals,cEVECS,Nsize,NUMSOL)
  integer,intent(in)::Nsize
  integer,intent(inout)::NUMSOL
  double precision,intent(inout)::cVals(NUMSOL),cEVECs(Nsize,NUMSOL)

  double precision::THRESH
  integer::zz,SSPIN,ICOUNT,IIRREP,AIRREP,I,A
CHARACTER*12 STRINGI(2), STRINGA(2)
LOGICAL incQM1, FIRST
  double precision::SCR(Nsize)
  double precision::bkupeVecs(NSize,NUMSOL),bkupeVals(NUMSOL)
  integer::QM1count


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
!! First, read in QM1 information
      call getrec(1,'JOBARC','NBASTOT',1,nbas)
      call getrec(1,'JOBARC','NREALATM',1,natm)
        call getrec(1,'JOBARC','NOCCORB',2,scratch)
      nocc=scr(1)
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
!! Then determine origin of largest magnitude of excitation;
!! if 'i' is not in QM1, omit that root
STRINGI(1)='I  [I_SYM]  '
STRINGI(2)='i  [i_SYM]  '
STRINGA(1)='A  [A_SYM]  '
STRINGA(2)='a  [a_SYM]  '
      THRESH = 0.05D0
bkupeVecs=0.0d0
bkupeVals=0.0d0
   do zz=1,NUMSOL


      DO SSPIN = 1, 1+IUHF
!        WRITE(6,1300)
        IF (SSPIN .EQ. 1) THEN
          WRITE(6,*)'   SINGLE EXCITATION COEFFICIENTS AA'
        ELSE
          WRITE(6,*)'   SINGLE EXCITATION COEFFICIENTS BB'
        ENDIF
        FIRST =.TRUE.
!        CALL GETLST(SCR, zz, 1, 1, SSPIN, 94)
        SCR=cEVECS(1+(SSPIN-1)*(NSize/2):(NSize/2)+(NSize/2)*(SSPIN-1),zz)
        ICOUNT = 1
        DO IIRREP = 1, NIRREP
          AIRREP=DIRPRD(IIRREP,IRREPX)
          DO 1 I = 1, POP(IIRREP, SSPIN)
            DO 2 A = 1, VRT(AIRREP,SSPIN)
              IF (ABS(SCR(ICOUNT)).GT.THRESH) THEN
                IF (FIRST) THEN
                  WRITE(6,999) STRINGI(SSPIN), STRINGA(SSPIN)
  999             FORMAT(6X,A12,3X,A12)
                  WRITE(6,*)
                ENDIF
                FIRST =.FALSE.
                WRITE(6, 1001) I, IIRREP, A, AIRREP, SCR(ICOUNT)
              ENDIF
              ICOUNT = ICOUNT + 1
    2       CONTINUE
    1     CONTINUE
        ENDDO
 1001   FORMAT(3X,I4,3X, ' [',I1,'],  ',1X, I4,3X, ' [',I1,']   ;',&
            10X, F12.6)
       enddo

! Both alpha and beta portion of vector inside QM1; so
! save this vector and value
       if (incQM1) then ! vector in QM1:: save
          QM1count=QM1count+1
          bkupeVecs(:,QM1count)=cEVECS(:,zz)
          bkupeVals(QM1count)=cVals(zz)
       endif
  enddo

  print*,'** Number of vectors inside QM1:',QM1count
  Call Dzero(cVals,NUMSOL)
  Call Dzero(cEVECS,Nsize)
  cEVECS=bkupeVecs
  cVals=bkupeVals
  NUMSOL=QM1count


30    FORMAT (I4,2X,I4)

end subroutine
