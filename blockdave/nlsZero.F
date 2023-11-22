#include "flags.h"
        subroutine NLSZERO(SCR,MAXCOR,IRREPX,IUHF,IROOT,NSIZE)
C
C The core EEs are obtained by setting the R(ij,ab) = 0 when i and
C j do  not belong to the core region. The present version works
C state specifc fashion.
C
C                           all; each
C LIST 444:    C(IJ,AB )       A<B ; I<J     AA AA
C      445:    C(ij,ab )       a<b ; i<j     BB BB
C      446:    C(Ij,Ab )       A,b ; I,j     AB A
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION SCR(MAXCOR)
      DOUBLE PRECISION VALUE,O_SCALE,V_SCALE
      DOUBLE PRECISION VAL1,VAL2,DIFF,THRES2
      INTEGER INTOO(2),INTOV(2)
C
#include "maxbasfn.par"
C
C
      CHARACTER*12 STRINGI(2), STRINGJ(2), STRINGA(2),STRINGB(2)
      DIMENSION LS2OUT(2,2)
      INTEGER LS1OUT
      LOGICAL NLS_EXIST,NBO_EXIST,PROJECT_SINGLES
      integer*8, dimension(2) :: scratch
      integer ierr,natm,nbas,nocc,nvirt
      INTEGER QM1num,NLMOnum,NLMOnum2,NLMOorigin,NLMOorigin2,indx
      INTEGER,allocatable::QM1atoms(:),QM2atoms(:),NLMO(:,:)
      integer,allocatable:: NLMOQM1(:),NLMOQM2(:)
      integer::compareIA(2),QMregIA(2),compareIJ(2),QMregIJ(2)
        double precision, allocatable :: Waa(:),Wab(:) 
      DOUBLE PRECISION::OLDSCR(NSIZE)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /PROJECT/ IPROJECT, IPATTERN, NCALC, ICALC, IWINDOW(8)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
C
      DATA BLANK /" "/

#include "info.com"
#include "flags.com"
#include "flags2.com"
C
      PROJECT_SINGLES = (Iflags2(h_IFLAGS2_prjct_single) .NE. 0)
      print*, 'value of proj singles: ', PROJECT_SINGLES
      call getrec(1,'JOBARC','NBASTOT',1,nbas)
      call getrec(1,'JOBARC','NREALATM',1,natm)
        call getrec(1,'JOBARC','NOCCORB',2,scratch)
      nocc=scr(1)
      nvirt=nbas-nocc
      allocate(NLMOQM1(nbas),NLMOQM2(nbas))
      NLMOQM1=0
      NLMOQM2=0
!******************************************************************
! * Read NLS file detailing QM regions 'QMcenter'
!******************************************************************
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
        print*,'FATAL ERROR: nbocenters-output of N. Flocke NLMO code-
     &               can not be found'
        call exit(1)
      endif     

!******************************************************************
        print*, 'INSIDE DRIVE_NLS.F'
#ifdef _DEBUG_LVL0
        print*, 'INSIDE DRIVE_NLS.F'
        print*, 'Contents of NLMO()'
        
        print*, NLMO
#endif

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
!#endif
!******************************************************************
!******************************************************************
!       *       DONE OBTAINING AUXILLARY NLS INFO       *
!       *       NOW MAKE MODIFCATIONS TO GUESS VECTOR   *
!******************************************************************
!******************************************************************

      LS1OUT      = 490
      LS2OUT(1,1) = 444
      LS2OUT(1,2) = 446
      LS2OUT(2,1) = 446
      LS2OUT(2,2) = 445
C
      STRINGI(2)='i  [i_SYM]  '
      STRINGJ(2)='j  [j_SYM]  '
      STRINGA(2)='a  [a_SYM]  '
      STRINGB(2)='b  [b_SYM]  '
      STRINGI(1)='I  [I_SYM]  '
      STRINGJ(1)='J  [J_SYM]  '
      STRINGA(1)='A  [A_SYM]  '
      STRINGB(1)='B  [B_SYM]  '
!      CALL ZERO(SCR,NSIZEC)
!     pre-dumps the R vector from block Davidson so that existing
!     code should run
      call LANCZOS_DUMP_VEC(irrepx,scr,maxcor,490,0,0,
     &                      443,0,IUHF,.False.)

!******************************************************************
!******************************************************************
!******************************************************************
! * Modified version of 'drive_nto_ee_state_prjct.F' lines 165-233
! * Works on the single excitation coefficient C(i,a) by zeroing out
! * the i not in QM1
! *** FIRST TEST FOR JUST THE UHF CASE, THEN TRY RHF WHICH HAS NEVER 
! *** BEEN DONE IN THE NLS SCHEME 12/1/2020


      IF (PROJECT_SINGLES) THEN
        DO 500 SSPIN = 1, 1+IUHF
          CALL GETLST(SCR, 1, 1, 1, SSPIN, LS1OUT)
          ICOUNT = 1
          DO IIRREP = 1, NIRREP
            AIRREP=DIRPRD(IIRREP,IRREPX)

          DO 1 I = 1, POP(IIRREP, SSPIN)
            DO 2 A = 1, VRT(AIRREP,SSPIN)

            compareIA=0
            QMregIA=0
            compareIA=(/ i,a+nocc /)
            call findQMregion(compareIA,size(compareIA),NLMOQM1,
     &           size(NLMOQM1),NLMOQM2,size(NLMOQM2),QMregIA)
!            print*,'ia',I,A
!            print*,'QMredIA', QMregIA
!            print*, 'value of Eig?', SCR(ICOUNT)
            if (QMregIA(1).eq.2) then
              SCR(ICOUNT) = 0.0D0
            else
!              SCR(ICOUNT) =1.0d0
            endif
            print*, SCR(ICOUNT)
!                IF (IIRREP.EQ.INTOO(1) .AND. AIRREP.EQ.INTOV(1)) THEN
!                   IF (I.EQ.INTOO(2) .AND. A.EQ.INTOV(2)) THEN
!                      WRITE(6,'(2(A,I3,A,I1,A))')'R1:',I,'[',IIRREP,']',
!     &                                           ' ->',A,'[',AIRREP,']'
!                      SCR(ICOUNT) = 1.0D0
!                   ELSE
!                      SCR(ICOUNT) = 0.0D0
!                   ENDIF
!                ELSE
!                   SCR(ICOUNT) = 0.0D0
!                ENDIF
              ICOUNT = ICOUNT + 1
    2       CONTINUE
    1     CONTINUE
        print*,'Evecs modified'
        print*,SCR(1:ICOUNT)
        ENDDO
        CALL PUTLST(SCR, 1, 1, 1, SSPIN, LS1OUT)
  500 CONTINUE
      ENDIF



30    FORMAT (I4,2X,I4)
!******************************************************************
!******************************************************************
!       * BEGIN MODIFYING R(IJ,AB)      *
!******************************************************************
!******************************************************************

      DO ICASE = 1, 1 + 2 * IUHF
        IF (ICASE .EQ. 1) THEN
          ASPIN = 1
          BSPIN = 2
        ELSEIF (ICASE .EQ. 2) THEN
          ASPIN = 1
          BSPIN = 1
        ELSEIF (ICASE .EQ. 3) THEN
          ASPIN = 2
          BSPIN = 2
        ENDIF

      LISTS2EX = LS2OUT(ASPIN, BSPIN)

      DO 10 RIRREP = 1, NIRREP
        LIRREP = DIRPRD(RIRREP,IRREPX)
        IF (ICASE .EQ. 1) THEN
          DISSYS = IRPDPD(LIRREP, ISYTYP(1,LISTS2EX))
          NUMDSS = IRPDPD(RIRREP, ISYTYP(2,LISTS2EX))
          DISSYA = DISSYS
          NUMDSA = NUMDSS
        ELSE
          DISSYA = IRPDPD(LIRREP, ISYTYP(1,LISTS2EX))
          NUMDSA = IRPDPD(RIRREP, ISYTYP(2,LISTS2EX))
          DISSYS = IRPDPD(LIRREP, 18 + ASPIN)
          NUMDSS = IRPDPD(RIRREP, 20 + ASPIN)
        ENDIF



       IF (DISSYA * NUMDSA .GT. 0) THEN
          I000 = 1
          I010 = I000 + DISSYS*NUMDSS *IINTFP
          IF (ICASE .EQ. 1) THEN
            CALL GETLST(SCR(I000), 1, NUMDSS, 1, RIRREP,
     $       LISTS2EX)
            CALL ZERO(SCR(I000),DISSYS*NUMDSS)
          ELSE ! expands I<J, A<B into IJ;AB
            CALL GETLST(SCR(I000), 1, NUMDSA, 1, RIRREP,
     $       LISTS2EX)
            CALL SYMEXP(RIRREP, POP(1,ASPIN),DISSYA,SCR(I000))
            CALL SYMEXP2(LIRREP,VRT(1,ASPIN),DISSYS, DISSYA,
     &         NUMDSS, SCR(I000), SCR(I000))
            CALL ZERO(SCR(I000),DISSYS*NUMDSS)
          ENDIF
          ICOUNT = I000
          DO 20 JIRREP = 1, NIRREP
            IIRREP = DIRPRD(JIRREP, RIRREP)
!       ******** Begin NLS logic *******
!       If i,j outside of user-defined QM1, then that
!       element in the vector is 0; 1 otherwise.

            DO 300 J= 1, POP(JIRREP,BSPIN)
              DO 40 I = 1, POP(IIRREP, ASPIN)
                DO 50 BIRREP = 1, NIRREP
                  AIRREP = DIRPRD(LIRREP, BIRREP)
                  DO 60 B = 1, VRT(BIRREP, BSPIN)
                    DO 70 A = 1, VRT(AIRREP, ASPIN)

! **NOTE** : Will only implement normal NLS scheme
!            Once tested, will implement charge transfer
!            modification. ZW 12/3/2020

            compareIJ=0
            QMregIJ=0
            compareIJ=(/ i,j /)
            call findQMregion(compareIJ,size(compareIJ),NLMOQM1,
     &           size(NLMOQM1),NLMOQM2,size(NLMOQM2),QMregIJ)
!            print*,'ia',I,J
!            print*,'QMredIA', QMregIJ
            if (any(QMregIJ.eq.2)) then
              SCR(ICOUNT) = 0.0D0
            else
!              SCR(ICOUNT) =1.0d0
            endif
!                print*,'SCR val:', SCR(ICOUNT)
!                       IF (ICASE .EQ. 1) THEN
!
!                          SCR(ICOUNT) = 1.0D0
!
!                       ELSE IF (ICASE .EQ. 2) THEN
!
!                          SCR(ICOUNT) = 1.0D0
!
!
!                       ELSE IF (ICASE .EQ. 3) THEN
!
!                          SCR(ICOUNT) = 1.0D0
!
!                       ENDIF

                    ICOUNT = ICOUNT + 1
   70             CONTINUE
   60           CONTINUE
   50         CONTINUE
   40       CONTINUE
  300     CONTINUE
   20   CONTINUE

      ENDIF
! ***********************************************************
! ***********************************************************!
! ***********************************************************
! * ICASE=1 is only case in RHF; other ICASE available in UHF
      IF (ICASE .EQ. 1) THEN
         CALL PUTLST(SCR(I000), 1, NUMDSS, 1, RIRREP,
     $               LISTS2EX)
      ELSE
C
C SQUEEZE ARRAY IN PROPER FORM A<B, I<J
C
         CALL SQSYM(LIRREP, VRT(1,ASPIN),DISSYA,DISSYS,
     &              NUMDSS,SCR(I010), SCR(I000))
         CALL TRANSP(SCR(I010), SCR(I000),NUMDSS,DISSYA)
         CALL SQSYM(RIRREP,POP(1,ASPIN), NUMDSA, NUMDSS,
     &              DISSYA, SCR(I010), SCR(I000))
         CALL TRANSP(SCR(I010), SCR(I000),DISSYA,NUMDSA)
         CALL PUTLST(SCR(I000), 1, NUMDSA, 1, RIRREP,
     $               LISTS2EX)
      ENDIF

   10 CONTINUE
      ENDDO
C
      CALL LOADVEC1(IRREPX,SCR,MAXCOR,IUHF,490,0,443,NSIZEC,
     &              .FALSE.)
C
C  PUT NORMALIZED PROJECTION VECTORS BACK ON LIST
C
#ifdef _DEBUG_LVLM
      Write(6,"(a)")" The checksum of the state spe. core projection"
      Call checksum(" @-DRIVE_CORE_EE",SCR,NSIZEC,S1)
#endif
      CALL UPDATES(IRREPX,SCR,444,0,490,IUHF)
C
! **PERFORM THE NLS-EOM PART NOW THAT SCR IS LOADED/DUMPEED TO DISK
      OLDSCR=SCR(1:NSIZEC)
      Call Hbarxc(SCR,Maxcor,Iuhf,1,Irrepx)
      Call Loadvec1(Irrepx,SCR,Maxcor,Iuhf,490,2,460,
     +                      NSIZEC,.False.)

!     oldscr^t * scr === R^t H R === E
      print*,'NLS-EOM-CCSD energy:', dot_product(OLDSCR,SCR(1:NSIZE))

      RETURN

          end subroutine


