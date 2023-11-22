#include "flags.h"
C
C THIS ROUTINE CREATES THE VARIOUS FOCK MATRIX LISTS AND
C  WRITES THEM TO MOINTS.
C
cjp for BWCC it adds the contributions to fock and vacuum energies
cjp for the reference configurations
cjp it also sets indices for Heff elements corresponding to internal
cjp monoexcitations (internfrom1 etc.)
cjp (for double excitations, this is done in dijab)
c
      SUBROUTINE FOCKLIST(F,BUF,Work,Isize,NBAS,NBAST,IUHF)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION F(NBAS,NBAS),BUF(4*NBAST*NBAST),X,ONEM,DENOM
      DOUBLE PRECISION Work(Isize)
      CHARACTER*1 SP(2)
      LOGICAL NONHF
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFEA(2),NFMI(2)
      COMMON /FLAGS/ IFLAGS(100)
CJDW 10/4/95. KKB stuff
      COMMON /SHIFT/ ISHIFT,NDRGEO
CJDW END
      LOGICAL Response, X_field, Y_field, Z_field
cjp
      character*80 fname
       real*8 ruhf,enuc,escf
#include "bwcc.com"
       real*8 enermo
       dimension enermo(maxorb)
cjp
      real*8 denmax
      integer imax,amax

      integer findex
      findex(i,j)=(i*(i-1))/2+j
cjp
      DATA ONEM /-1.0/
      DATA SP /'A','B'/
C
C CREATE SINGLE LISTS FOR METHODS INCLUDING SINGLE EXCITATIONS
C
      CALL UPDMOI(1,NT(1),1,90,0,0)
cYAU - initialize T1(alpha) to zero
      call aces_list_memset(1,90,0)
      IF (IUHF.NE.0) THEN
         CALL UPDMOI(1,NT(2),2,90,0,0)
cYAU - initialize T1(beta) to zero
         call aces_list_memset(2,90,0)
      END IF
C
C  SKIP THIS PART FOR HF REFERENCE FUNCTIONS (check whether we 
C  are doing a finite differenc response property calculations. 
C
C Useful in doing finite difference response property
C calculations. Ajith Perera, 23/2014 (Redid something that I had done
C in 1992-94, but loss track of it).

      Response = .FALSE.
      X_field  = .FALSE.
      Y_field  = .FALSE.
      Z_field  = .FALSE.

      Response = (Iflags(h_IFLAGS_density) .EQ. 1)
      X_field  = (Iflags(h_IFLAGS_xfield) .NE. 0)
      Y_field  = (Iflags(h_IFLAGS_yfield) .NE. 0)
      Z_field  = (Iflags(h_IFLAGS_zfield) .NE. 0)

      If (Response .AND. (X_field .OR. Y_field .OR. Z_field)) Then

          If (Iflags(h_IFLAGS_nonhf) .EQ. 0) Then
              Write(*,*)
              Write(6,"(a,a,a)") " Finite field response property",
     &                           " calculations must turn on",
     &                           " the NONHF flag."
              Call Errex
          Endif  
      Endif 

      IF (.not.isbwcc    .and.
     &    IFLAGS(38).EQ.0.AND.
     &    IFLAGS(77).EQ.0.AND.
     &    IFLAGS(11).NE.2     )  RETURN
C
      NBAS2=NBAS*NBAS
      NBAST2=NBAST*NBAST
      IOFFBF=NBAST*NBAST+1
c
cjp spin loop, some actions like file opening must be outside!
cjp necessary changes with respect to RHF code, when this loop
      call gfname('FOCKCD',fname,ilength)
      open(unit=99,file=fname(1:ilength),form='unformatted',
     &     status='unknown')
      write(99)nbas
cjp
c
      DO 1000 ISPIN=1,1+IUHF
C
C PICK UP AO BASIS FOCK MATRIX FROM JOBARC, TRANSFORM TO
C   MO BASIS AND GET RID OF DIAGONAL ELEMENTS AFTER WRITING THEM
C   TO JOBARC
c
cjp in mr-bwcc the diagonal elements will not be written to jobarc,
cjp but to a separate file in the form of differences wrto the RHF eigenvalues
cjp since jobarc is not replicated for each reference the offdiagonals are
cjp stored incl. the differences in the moints files for each reference
c
         if (ispin.eq.1) then
            CALL GETREC(20,'JOBARC','FOCKA   ',NBAST2*IINTFP,F)
         else
            CALL GETREC(20,'JOBARC','FOCKB   ',NBAST2*IINTFP,F)
         end if
#ifdef _DEBUG_LVL
        If(ispin .eq. 1)  then
          write(6,"(a)") "The Alpha Fock Matrix"
          call output(f, 1, nbas, 1, nbas, nbas, nbas, 1)
        else 
          write(6,"(a)") "The Beta Fock Matrix"
          call output(f, 1, nbas, 1, nbas, nbas, nbas, 1)
        endif 
#endif 
         CALL AO2MO2(F,F,BUF,BUF(IOFFBF),NBAS,NBAST,ispin)

         IF (IFLAGS(1).GT.10) THEN
            WRITE(6,5000) ispin
5000        FORMAT(T3,' MO basis Fock matrix (',i1,' spin)')
            WRITE(6,'((2I5,F10.5,20X,2I5,F10.5))')
     &                                  ((I,J,F(I,J),I=1,NBAS),J=1,NBAS)
        END IF
#ifdef _DEBUG_LVLM
        If(ispin .eq. 1)  then
          write(6,"(a)") "The Alpha Fock Matrix"
          call output(f, 1, nbas, 1, nbas, nbas, nbas, 1)
        else 
          write(6,"(a)") "The Beta Fock Matrix"
          call output(f, 1, nbas, 1, nbas, nbas, nbas, 1)
        endif 
#endif 

cjp store diagonal elements extra
cjp NOTE: in order not to replicate JOBARC, we store only elements of 1st
cjp reference's fock matrix there NOTE: the fock contributions we add are
cjp indexed in energy order, but fock matrix here is in IR,energy order.
cjp we ignore symmetry at first anyway here orbital energies are replaced by
cjp fock diagonal elements
c
         CALL dcopy(NBAS,F,NBAS+1,BUF,1)
         if (ispin.eq.1) then
            CALL PUTREC(20,'JOBARC','SCFEVALA',NBAS*IINTFP,BUF)
         else
            CALL PUTREC(20,'JOBARC','SCFEVALB',NBAS*IINTFP,BUF)
         end if
c
cjp compute one-determinantal energies
         if (isbwcc) then
cjp for HF prepare the hcore by backcomputation from eigenvalues and 2el
cjp contributions
            if (iref.eq.1) then
               do i = 1, nbas
                  hcore(i,ispin) = buf(i)-fock2elcontr0(i,ispin)
               end do
            end if
cjp initialize before alpha spin is added
            if (ispin.eq.1) enerscf(iref)=0.0
            if (ispin.eq.1.and.iref.eq.1) enerscf0=0.0
cjp one-electron integral contributions
            if (iuhf.ne.0) then
               ruhf=1.
            else
               ruhf=2.0
            end if
            do i = 1, nbas
               if (iocc(i,iref,ispin).ne.0) then
                  enerscf(iref)=enerscf(iref)+ruhf*hcore(i,ispin)
               end if
               if (iref.eq.1.and.iocc0(i,ispin).ne.0) then
                  enerscf0=enerscf0+ruhf*hcore(i,ispin)
               end if
            end do
         end if
         IF (IUHF.EQ.0) THEN
            CALL PUTREC(20,'JOBARC','SCFEVALB',NBAS*IINTFP,BUF)
         END IF
cjp add here the contributions; store diagonal contribs;
cjp in the fockdiag files leave diagonal ones
cjp otherwise jobarcs would have to be replicated too
         if (isbwcc) then
            do i = 1, nbas
            do j = 1, i
               if (i.eq.j) then
                  f(i,j)=       fockcontr(findex(i,j),ispin)
               else
                  f(i,j)=f(i,j)+fockcontr(findex(i,j),ispin)
                  f(j,i)=f(j,i)+fockcontr(findex(i,j),ispin)
               end if
            end do
            end do
cjp compute one-determinantal energies of the current  fermi vacuum
            if (iuhf.ne.0) then
               ruhf=0.5
            else
               ruhf=1.0
            end if
            do i = 1, nbas
               if (iocc(i,iref,ispin).ne.0) then
                  enerscf(iref)=enerscf(iref)+ruhf*fock2elcontr(i,ispin)
               end if
               if (iref.eq.1.and.iocc0(i,ispin).ne.0) then
                  enerscf0=enerscf0+ruhf*fock2elcontr0(i,ispin)
               end if
            end do
cjp write diagonal differences and vacuum energies into file for use in other
cjp modules
            do i = 1, nbas
               write(99)f(i,i)
            end do
cjp write also vacuum energy difference to HF vacuum
cjp NOTE: do not rewrite enerscf0 into jobarc: MOs can be DROPMOed and we would
cjp compute nonsence. However, it does not matter if only differences are
cjp computed here and the SCF energy of the 1st reference is computed in XVSCF
cjp this has to be handled in xvscf, if foreign MOs are read in, as well as the
cjp Fock matrix in AO, which is here only read from JOBARC
c
cjp .... postponed after additional write in dijab.f!!!        close(unit=99)
         end if

cjp zero out the diagonal elements again
c YAU : old
c        CALL SAXPY(NBAS,ONEM,F,NBAS+1,F,NBAS+1)
c YAU : new
         do i = 1, nbas
            f(i,i) = 0.d0
         end do
c YAU : end
       

C --- Now the Fock digonal has been set to zero, we can add a external
C perturbation. Useful in doing finite difference response property
C calculations. Ajith Perera, 23/2014 (Redid something that I had done 
C in 1992-94, but loss track of it). See at the begining of this 
C routine for setting of these variables. 

      If (Response .AND. (X_field .OR. Y_field .OR. Z_field)) Then

          Call Add_ext_pert(Work, Isize, X_field, Y_field, Z_field,
     &                       F, Nbas, Nbast, Ispin)
     
      Endif 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c now split the fock offdiagonals to OO, VV and OV parts, form also denominators
c from the diagonal
c
C
C FORM SYMMETRY PACKED OCCUPIED-OCCUPIED PART.
C
         ITHRU = 0
         IOFF  = 0
         DO 10 IRREP=1,NIRREP
         DO 20 J=1,POP(IRREP,ISPIN)
         DO 21 I=1,POP(IRREP,ISPIN)
            ITHRU=ITHRU+1
cjp changed indexing for general reference
cjp this and the following is unclean when symmetry should be introduced!!!
            if (isbwcc) then
               iind=invhnum(i,iref,ispin)
               jind=invhnum(j,iref,ispin)
               buf(ithru)=f(iind+ioff,jind+ioff)
            else
               BUF(ITHRU)=F(I+IOFF,J+IOFF)
            end if
21       CONTINUE
20       CONTINUE
         IOFF=IOFF+POP(IRREP,ISPIN)
10       CONTINUE
C
CJDW 10/4/95. KKB stuff. ISHIFT added.
         CALL UPDMOI(1,ITHRU,ISPIN+2,91+ISHIFT,0,0)
         CALL UPDMOI(1,ITHRU,ISPIN+4,91+ISHIFT,0,0)
         CALL PUTLST(BUF,1,1,1,ISPIN+2,91+ISHIFT)
         CALL PUTLST(BUF,1,1,1,ISPIN+4,91+ISHIFT)
CJDW END
C
C FORM SYMMETRY PACKED VIRTUAL-VIRTUAL PART.
C
         ITHRU = 0
         IOFF  = NOCCO(ISPIN)
         DO 110 IRREP=1,NIRREP
         DO 120 J=1,VRT(IRREP,ISPIN)
         DO 121 I=1,VRT(IRREP,ISPIN)
            ITHRU=ITHRU+1
cjp changed indexing for general reference
            if (isbwcc) then
               iind=invpnum(i,iref,ispin)
               jind=invpnum(j,iref,ispin)
               buf(ithru)=f(iind+ioff-nocco(ispin),
     &                      jind+ioff-nocco(ispin))
            else
               BUF(ITHRU)=F(I+IOFF,J+IOFF)
            end if
121      CONTINUE
120      CONTINUE
         IOFF=IOFF+VRT(IRREP,ISPIN)
110      CONTINUE
C
CJDW 10/4/95. KKB stuff. ISHIFT added.
c         print *, '@FOCKLIST: list ',ISPIN+2,92+ISHIFT
c         call dmat_prt(buf,1,ithru,1)
         CALL UPDMOI(1,ITHRU,ISPIN+2,92+ISHIFT,0,0)
         CALL UPDMOI(1,ITHRU,ISPIN+4,92+ISHIFT,0,0)
         CALL PUTLST(BUF,1,1,1,ISPIN+2,92+ISHIFT)
         CALL PUTLST(BUF,1,1,1,ISPIN+4,92+ISHIFT)
CJDW END
C
C FORM SYMMETRY PACKED VIRTUAL-OCCUPIED PART.
C
       ITHRU=0
       IOFFA=NOCCO(ISPIN)
       IOFFI=0
       I000=1
       I010=I000+NOCCO(ISPIN)*NVRTO(ISPIN)
       I020=I010+NOCCO(ISPIN)+NVRTO(ISPIN)
       IF(IFLAGS(38)+IFLAGS(77).NE.0)THEN
C
CJDW 10/4/95. Extra logic for KKB stuff.
C
         IF(NDRGEO.EQ.2)THEN
           CALL GETREC(20,'JOBARC','SCFEVL'//SP(ISPIN)//'0',IINTFP*NBAS,
     &                 BUF(I010))
         ELSE
           CALL GETREC(20,'JOBARC','SCFEVAL'//SP(ISPIN),IINTFP*NBAS,
     &                 BUF(I010))
         ENDIF
CJDW END
cjp correct for diagonal contributions
        IOFFO=I010-1
        IOFFV=I010+NOCCO(ISPIN)-1
cjp add diagonal corrections for non-hf character of higher references
      if(isbwcc) then
      do i=1,nbas
        buf(ioffo+i)=buf(ioffo+i)+fockcontr(findex(i,i),ispin)
      enddo
      endif
cjp modifications for bwcc - blow denominators and
cjp find indices of Heff nondiagonal elements of single excitations
cjp
        imax=0
        amax=0
        denmax=0.
        DO 210 IRREP=1,NIRREP
         DO 220 I=1,POP(IRREP,ISPIN)
          DO 221 A=1,VRT(IRREP,ISPIN)
           ITHRU=ITHRU+1
cjp changed indexing and special treatment of internal excitations
           if(isbwcc) then
              iind=invhnum(i,iref,ispin)
              aind=invpnum(a,iref,ispin)
              buf(ithru)=f(aind+ioffa-nocco(ispin),iind+ioffi)
              if(isactive(iind,ispin) .and. isactive(aind,ispin)) then
              denom=1.d0/denomblow
              if(abs(buf(ithru)).gt.1e-8) then
                  write(6,*)'@FOCKLIST-I: nonzero internal amplitude
     &             encountered, ref ',
     &           iref,iind,aind,buf(ithru),(buf(ioffo+iind)
     &    -buf(ioffo+aind)),buf(ithru)/(buf(ioffo+iind)-buf(ioffo+aind))
              endif
cjp here also set internfrom1, internto1, internindex1
cjp candidate for Heff matrix element
cjp find for which pair of references it plays the role of Heff(offdiagonal)
cjp and store this info
               if(nref.gt.1 .and. iuhf.eq.1 ) then
                  k=ifindref(nbas,1,iind,aind,ispin)
                  if(k.gt.0) then
cjp                 write the information to array and exit searching
                   internnum1(iref,ispin)=internnum1(iref,ispin)+1
                   internindex1(internnum1(iref,ispin),iref,ispin)=ithru
                   internfrom1(internnum1(iref,ispin),iref,ispin)=iref
                   internto1(internnum1(iref,ispin),iref,ispin)=k
                  endif
                endif
cjp end setting Heff indices
              else
cjp not indices from active space
              denom=1.0/(buf(ioffo+iind)-buf(ioffv+aind-nocco(ispin)))
              if(abs(denom).gt.abs(denmax)
     +           .and.abs(denom).lt.intruder) then
                   denmax=denom
                   imax=iind
                   amax=aind
              endif
              endif
           else
              BUF(ITHRU)=F(A+IOFFA,I+IOFFI)
              DENOM=1.0/(BUF(IOFFO+I)-BUF(IOFFV+A))
           endif
          BUF(I020+ITHRU-1)=BUF(ITHRU)*DENOM
221       CONTINUE
220      CONTINUE
        IOFFI=IOFFI+POP(IRREP,ISPIN)
        IOFFA=IOFFA+VRT(IRREP,ISPIN)
        IOFFO=IOFFO+POP(IRREP,ISPIN)
        IOFFV=IOFFV+VRT(IRREP,ISPIN)
210     CONTINUE
c        CALL UPDMOI(1,ITHRU,ISPIN,90,0,0)
CJDW 10/5/95. Should we or should we not have this PUTLST call when we
C             go through xintprc for the second time in dropped core runs ?
        CALL PUTLST(BUF(I020),1,1,1,ISPIN,90)
       ELSE
        if(isbwcc) stop 'this place not assumed to be entered in BWCC'
        DO 1210 IRREP=1,NIRREP
         DO 1220 I=1,POP(IRREP,ISPIN)
          DO 1221 A=1,VRT(IRREP,ISPIN)
           ITHRU=ITHRU+1
           BUF(ITHRU)=F(A+IOFFA,I+IOFFI)
1221      CONTINUE
1220     CONTINUE
         IOFFI=IOFFI+POP(IRREP,ISPIN)
         IOFFA=IOFFA+VRT(IRREP,ISPIN)
1210    CONTINUE
       ENDIF
CJDW 10/4/95. KKB stuff. ISHIFT added.
       CALL UPDMOI(1,ITHRU,ISPIN+2,93+ISHIFT,0,0)
       CALL PUTLST(BUF,1,1,1,ISPIN+2,93+ISHIFT)
       CALL UPDMOI(1,ITHRU,ISPIN+4,93+ISHIFT,0,0)
       CALL PUTLST(BUF,1,1,1,ISPIN+4,93+ISHIFT)
CJDW END
C
C PUT OUT FIRST-ORDER T1 FOR QRHF
C

cjp intruder diagnostics
      if(isbwcc) then
       write(6,*)'@FOCKLIST-I: MAX denomin. ref, spin, i, a, denom',
     &  iref,ispin,imax,amax,denmax
       if(abs(denmax).gt.abs(totmaxdenom)) totmaxdenom=denmax
      endif

cjp some debug output
      if(isbwcc.and.bwgossip) then
      if(ispin.eq.1) then
          call getrec(20,'JOBARC','SCFEVALA',nbas*iintfp,enermo(1))
      else
          call getrec(20,'JOBARC','SCFEVALB',nbas*iintfp,enermo(1))
      endif
      write(6,*)' reference ',iref,' spin ',ispin,' fock diagonal'
      write(6,*) (enermo(i)+f(i,i),i=1,nbas)
      write(6,*)' reference ',iref,' spin ',ispin,' fock offdiagonal'
      do i=2,nbas
      do j=1,i-1
       write(6,*) 'i,j,f(i,j)',i,j,f(i,j)
      enddo
      enddo
      endif

cjp end of ispin loop
1000  CONTINUE
cjp enerscf output must be postponed outside the spin loop
        write(99) enerscf(iref)-enerscf0
        if(isbwcc.and.iref.eq.1)
     & write(6,*)'@FOCKLIST-I Evacuum(SCF)= ',enerscf0,' proc. ',1
        if(isbwcc) write(6,*)'@FOCKLIST-I Evacuum('
     &          ,iref,')-Evacuum(SCF): ',enerscf(iref)-enerscf0
cjp

c
cjp finish preparation of data for masik program - THIS IS ONLY FOR RHF CASE
cjp if relevant
      if(isbwcc.and.masik) then
cjp store fock diagonal and offdiagonal
cjp store "SCF" energy
      call getrec(20,'JOBARC','SCFEVALA',nbas*iintfp,enermo(1))
      do i=1,nbas
        f(i,i)=enermo(i)+f(i,i)
      enddo
      write(33) (f(i,i),i=1,nbas)
      enuc=0d0
      ehf=0d0
      call getrec(20,'JOBARC','NUCREP', iintfp,enuc)
      call getrec(20,'JOBARC','SCFENEG',iintfp,escf)
      write(33) enuc,escf-enuc,escf
      close(33)
      write(6,*)' masik data part 2 prepared'
cjp write also file with the fock matrix of closed shell first reference
      open(unit=33,file='ftn33',status='unknown', form='unformatted')
      do i=1,nbas
        write(33) (f(j,i),j=1,nbas)
      enddo
      close(33)
      stop 'masik data prepared'
      endif
c
      RETURN
      END
