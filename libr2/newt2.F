
      SUBROUTINE NEWT2(ICORE,MAXCOR,IUHF)
C
C THIS ROUTINE PICKS UP THE FINAL T2 OR T2 INCREMENTS, DENOMINATOR
C  WEIGHTS THEM AND THEN OVERWRITES THE T2 INCREMENT LIST WITH THE NEW VALUES.
C  GOES ONE IRREP AT A TIME OVER ALL SPIN CASES.  FOR RHF CASES,
C  IT ALSO CALLS A ROUTINE WHICH FORMS THE AA AMPLITUDES FROM THE 
C  AB VALUES.
cjp
c This routine is modified for mr-bw-CC in following way:
c
c 1. ecorrbw variable (contains epsilon0-Heff(mm)) is additional
c     parameter of subroutine
c 2. offdiagonal Heff elements are identified as the internal amplitudes
c 3. 1/D(ij,ab) is replaced by 1/(D(ij,ab)+ecorrbw): instead of vecprd
c     routine a new routine bwvecprd is called
c 4. terms from disconnected diagrams arred added to the RHS of bwcc T2 
c
c This task is split in two passes - the routine will be called twice
c In the first pass the offdiagonal Heff elements are gathered
c The second pass must be dine after Heff diagonalization is finished!
c since until then ecorrbw is not known.
c For standard cc, the routine performs everything in one pass
c
cjp

C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ELAST,FACT
      DIMENSION ICORE(MAXCOR)
      LOGICAL LINCC,CICALC
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/ POP1(8),POP2(8),VRT1(8),VRT2(8),NTAA,NTBB,NF1AA,
     &             NF1BB,NF2AA,NF2BB
      COMMON /LINEAR/ LINCC,CICALC
cjp needed common switch
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,QCISD
      LOGICAL NONHF,TRIPIT,TRIPNI,TRIPNI1,T3STOR,PRESNT,INIT,
     &        DORESET,UCC,RESTART,BRUECK
      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC
      COMMON /CORENG/ ELAST
cjp
      COMMON /INFO/ NOCCO(2),NVRTO(2)
cjp
#include "bwcc.com"
        integer intptr
         real*8 getheff
         real*8 scalfactor
cjp help variables to pass uniquelly through all 3 spin cases
        integer ntxx,nocvrxx,ispinxx
        dimension ntxx(3),nocvrxx(3),ispinxx(3)
        
cjp
        nocca=nocco(1)
      noccb=nocco(2)
      nvrta=nvrto(1)
      nvrtb=nvrto(2)
      nocvrxx(1)=nocca+nvrta
      nocvrxx(2)=nocca+nvrta
      nocvrxx(3)=nocca+nvrta
      ntxx(1)=ntaa
      ntxx(2)=ntbb
      ntxx(3)=ntaa
      ispinxx(1)=1
      ispinxx(2)=2
      ispinxx(3)=1

      if(isbwcc.and.bwgossip) 
     +    write(6,*)'entering bwnewt2, ecorrbw= ',ecorrbw
      IBOT=1
      IF(IUHF.EQ.0)IBOT=3
      I000=1
      DO 5 ISPIN=IBOT,3
       LSTDEN=47+ISPIN
       LSTINC=60+ISPIN
       DO 10 IRREP=1,NIRREP
        NDSSYM=IRPDPD(IRREP,ISYTYP(2,LSTINC))
        DISSYM=IRPDPD(IRREP,ISYTYP(1,LSTINC))
        NSIZE=NDSSYM*DISSYM
        I010=I000+NSIZE*IINTFP
        I020=I010+NSIZE*IINTFP
        IF(I020.GT.MAXCOR)CALL INSMEM('NEWT2',I020,MAXCOR)
cjp read amplitude increments
        CALL GETLST(ICORE,1,NDSSYM,1,IRREP,LSTINC)
#ifdef _DEBUG_LVLM
        call checksum("T2-increment",ICORE, NDSSYM*DISSYM)
#endif 
cjp
cjp identify and store Heff offdiagonal elements
cjp skip this if ibwpass ne 1
        if(isbwcc .and. nref.gt.1 .and. ibwpass.eq.1) then
cjp FOR LATER NOTE THAT WE ARE HERE INSIDE THE ISPIN AND IRREP LOOPS!!!
cjp
cjp@@@ NOTE: we should also check if some internal T2 amplitudes
cjp are nonzero, as it is done in e4seng for t1. It might be done
cjp also in bwvecprd - but would probably slow down the code
cjp
        intptr=1
        do i=1,nsize
        if(intptr.le.internnum(iref,ispin) .and. 
     +      internindex(intptr,iref,ispin).eq.i) then
cjp             this is offdiag element, but maybe such one which does not belong to current iref
                if(iref.eq.internfrom(intptr,iref,ispin)) then
cjp                     relevant one
                        heff(internto(intptr,iref,ispin),iref)=
     +                   getheff(icore(I000),i)
                endif
                intptr=intptr+1
cjp             they are sorted acc. to increasing internindex, look simply at next
        endif
        enddo
        endif
cjp end heff gathering
cjp we are done with first pass in bwcc - continue irrep and spin loops
        if(ibwpass.eq.1) goto 10
c
c
c
c
cjp
cjp add unlinked contributions according to Hubac eq. 4.28
cjp in the size-ext correction, these unlinked contributions will be omitted
cjp
        if(isbwcc.and.sing1) then
cjp we need core for temp storage of t1
        i030=i020+ntxx(ispin)*iintfp
        if(ispin.eq.3 .and. iuhf.ne.0) then
        i025=i030
        i030=i025+ntbb*iintfp
        else    
        i025=i020
        endif
        if(i030.gt.maxcor)call insmem('newt2',i030,maxcor)
cjp read old t1 amplitudes (aa case only for rhf)
        call getlst(icore(i020),1,1,1,ispinxx(ispin),90)
        if(ispin.eq.3 .and. iuhf.ne.0)
     +     call getlst(icore(i025),1,1,1,2,90)
c
c
c
cjp the use of useeq429=true is already obsolete
        if(useeq429.and.iuhf.eq.0) then
cjp simplified implementation for comparison with masik's' program
cjp add the t(i,a)t(j,b)-t(j,a)t(i,b) term scaled by ecorrbw
        call ftau(icore(I000),icore(i020),icore(i020),dissym,ndssym,
     +  pop1,pop2,vrt1,vrt2,irrep,ispin,ecorrbw)
        else
c
c
cjp add 2 terms of RHS of 4.28 with single excitation contributions
cjp first some allocation
        i040=i030+ntxx(ispin)*iintfp
        if(ispin.eq.3 .and. iuhf.ne.0) then
        i035=i040
        i040=i035+ntbb*iintfp
        else
        i035=i030
        endif
        i050=i040+ntxx(ispin)*iintfp
        if(ispin.eq.3 .and. iuhf.ne.0) then
        i045=i050
        i050=i045+ntbb*iintfp
        else
        i045=i040
        endif
        i060=i050+nocvrxx(ispin)*iintfp
        if(ispin.eq.3 .and. iuhf.ne.0) then
        i055=i060
        i060=i055+(noccb+nvrtb)*iintfp
        else
        i055=i050
        endif
        if(i060.gt.maxcor)call insmem('newt2',i060,maxcor)
c
cjp read t1 right hand sides and fix them
cjp for this is essential that second pass of e4s is called AFTER second pass of newt2
        call getlst(icore(i030),1,1,1,2+ispinxx(ispin),90)
        if(ispin.eq.3 .and. iuhf.ne.0) 
     +     call getlst(icore(i035),1,1,1,4,90)
c
cjp but these are not the <...> terms of 4.28 yet, denomin(not zeroed for t1)*t(i,a) must be added to them!
c
        if(ispin.eq.1 .or. ispin.eq.3) 
     +      call fixrhs(icore(i030),icore(i020),icore(i040),icore(i050),
     +         ntaa,pop1,vrt1,nocca,nvrta,1,iuhf)
        if(ispin.eq.2) 
     +      call fixrhs(icore(i030),icore(i020),icore(i040),icore(i050),
     +         ntbb,pop2,vrt2,noccb,nvrtb,2,iuhf)
        if(ispin.eq.3 .and. iuhf.ne.0) 
     +  call fixrhs(icore(i035),icore(i025),icore(i045),icore(i055),
     +         ntbb,pop2,vrt2,noccb,nvrtb,2,iuhf)
c
cjp NOTE: at ecorrbw=0, fixrhs result vanishes for all but internal amplitudes
c
c
c
cjp add rhs contributions to t2 equations, permute Pij Pab
cjp therefore use gtau procedure
        if(correctiontype.eq.0) then
          if(ihubaccorr.lt.2) then
             scalfactor=hfakt
          else
             scalfactor=0.0
          endif
        else
          scalfactor=1.0
        endif
       if(scalfactor.ne.0.0) call gtau(icore(I000),icore(i020),
     +    icore(i030),icore(i025),icore(i035),
     +    dissym,ndssym,pop1,pop2,vrt1,vrt2,irrep,ispin,scalfactor)
c
cjp these are monoexcitations on the LHS of 4.28
cjp this is OK, there is only one permutation and both quantities are t1
cjp so we can use simply ftau
c
cjp this is skipped in a posteriori correction and scaled in iterative one
cjp (ecorrbw is already scaled appropriatelly)
        if(ihubaccorr.lt.2) 
     +      call ftau(icore(I000),icore(i020),icore(i025),
     +         dissym,ndssym,pop1,pop2,vrt1,vrt2,irrep,ispin,-ecorrbw)
        endif
cjp
        endif
cjp
cjp denominator-weight the T2 amplitude increments
cjp
cjp read denominators
        CALL GETLST(ICORE(I010),1,NDSSYM,1,IRREP,LSTDEN)
#ifdef _DEBUG_LVL0
        call checksum("@-NEWT2=D2 :",ICORE(I010),NDSSYM*DISSYM)
        call checksum("@-NEWT2-T2 :", ICORE(I000), NSIZE)   
#endif 
        IF(CICALC)CALL CIDENOM(NSIZE,ELAST,ICORE(I010))
cjp
        if(isbwcc) then
cjp ICORE(I000) are t2 increments
        call bwvecprd(ecorrbw,ICORE(I000),ICORE(I010),ICORE(I000),NSIZE)
        else
        CALL VECPRD(ICORE(I000),ICORE(I010),ICORE(I000),NSIZE)
        endif

cjp store new amplitudes
        CALL PUTLST(ICORE(I000),1,NDSSYM,1,IRREP,LSTINC)
#ifdef _DEBUG_LVL0
        call checksum("@-NEWT2-T2:", ICORE(I000), NSIZE)   
        write(*,*)
#endif 
10     CONTINUE
c end of the ispin loop
5     CONTINUE

cjp done with pass one
       if(ibwpass.eq.1) return
C
C FORM T2(AA) AMPLITUDES AND DUMP THEM.
cjp in spin-restricted, the aa and bb amplitudes are redundant 
cjp and can be generated from the ab amplitudes
C
      IF(IUHF.EQ.0)THEN
       DO 200 IRREP=1,NIRREP
        NUMAB=IRPDPD(IRREP,ISYTYP(2,63))
        DSZAB=IRPDPD(IRREP,ISYTYP(1,63))
        NUMAA=IRPDPD(IRREP,ISYTYP(2,61))
        DSZAA=IRPDPD(IRREP,ISYTYP(1,61))
        NSIZ1=DSZAB*NUMAA
        NSIZ2=DSZAB*NUMAB
        I000=1
        I010=I000+IINTFP*NSIZ1
        I020=I010+IINTFP*NSIZ2
cjp new t2 amplitudes
        CALL GETLST(ICORE(I010),1,NUMAB,1,IRREP,63)
        CALL ASSYM(IRREP,POP1,DSZAB,DSZAB,ICORE(I000),ICORE(I010))
        CALL SQSYM(IRREP,VRT1,DSZAA,DSZAB,NUMAA,ICORE(I010),ICORE(I000))
        CALL PUTLST(ICORE(I010),1,NUMAA,1,IRREP,61)
200    CONTINUE 
      ENDIF
      RETURN
      END
