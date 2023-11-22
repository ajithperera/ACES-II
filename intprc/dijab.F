      SUBROUTINE DIJAB(E,D,ISCR,NBAS,IUHF)
C
C THIS ROUTINE WRITES THE 1/D(ijab) LISTS FOR MBPT/CC
C CALCULATIONS.  THE LISTS ARE WRITTEN OUT IN SYMMETRY
C PACKED FORMAT, WITH THE AA EIGENVALUES ON LIST 48,
C THE BB EIGENVALUES (UHF ONLY) ON LIST 49, AND THE AB
C EIGENVALUES ON LIST 50.  D IS A SCRATCH ARRAY WHICH MUST
C BE AT LEAST AS LARGE AS THE TOTAL SIZE OF THE LARGEST IRREP,
C E MUST BE 2*NBAS LONG, WHERE NBAS IS THE NUMBER OF BASIS
C FUNCTIONS, AND ISCR MUST BE (NVRTA*NVRTB+NOCCA*NOCCB+NIRREP)*2
C IN LENGTH.  RECIPROCAL DENOMINATORS ARE WRITTEN OUT FOR RLE STUFF
C IF THIS IS REQUIRED.
C
CEND
cjp
cjp in BWCC this is a very important routine: it blows denominators
cjp corresponding to internal excitations and identifies to which
cjp Heff matrix elements such RHSides will contribute - it sets
cjp index arrays internto, internfrom, internindex
cjp (for single excitations, internto1 etc. are set in focklist)
cjp
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION E(*),D(*)
      LOGICAL RLE
C
      DIMENSION ISCR(*)
C
      COMMON/FLAGS/IFLAGS(100)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
cjp 
#include "bwcc.com"
C
      DATA ONE/1.D0/
cjp for ifind on biexcitations
      integer efrom,eto,espin
      dimension efrom(2),eto(2),espin(2)
cjp
      integer imax,jmax,amax,bmax
      real*8 denmax

      integer findex
      findex(i,j)=(i*(i-1))/2+j


C
      NNM1O2(I)=(I*(I-1))/2
      ILRG(IX)=INT(0.5*(1.0+SQRT(8.0*IX-7))+1.D-8)+1
      RLE=.FALSE.
      IF(IFLAGS(21).NE.3)RLE=.TRUE.
      NOCA=NOCCO(1)
      NOCB=NOCCO(2)
      NVRTA=NVRTO(1)
      NVRTB=NVRTO(2)
      NBAS=NOCA+NVRTA
C
C THIS ROUTINE FORMS THE SYMMETRY-PACKED DENOMINATOR
C  LISTS.  THESE ARE WRITTEN OUT IN A SYMMETRY PACKED
C  WAY (A,B;I,J) FOR AB, (A<B;I<J) FOR AA AND (a<b;i<j) FOR BB.
C
      CALL GETREC(20,'JOBARC','SCFEVALA',NBAS*IINTFP,E)
      CALL GETREC(20,'JOBARC','SCFEVALB',NBAS*IINTFP,E(NBAS+1))
cjp add corrections for non-hf character of higher references
      if(isbwcc) then
      do i=1,nbas
        e(i)=e(i)+fockcontr(findex(i,i),1)
        if(iuhf.ne.0) then
        e(nbas+i)=e(nbas+i)+fockcontr(findex(i,i),2)
        else
cjp for RHF take the same for A and B
        e(nbas+i)=e(nbas+i)+fockcontr(findex(i,i),1)
        endif
      enddo
      endif

C For QRHF lets always print the diagonals 

       WRITE(6,*)' Alpha diagonal element of the reference state'
       WRITE(6,*)' ---------------------------------------------'
       WRITE(6,'((7F10.5))')(E(I),I=1,NBAS) 
       WRITE(6,*)' Beta  diagonal element of the reference state'
       WRITE(6,*)' ---------------------------------------------'
       WRITE(6,'((7F10.5))')(E(NBAS+I),I=1,NBAS)
       IF (IFLAGS(77) .NE. 0) THEN
       WRITE(6,"(A,A,A)") " CAUTION! For QRHF the Fock Matrix of the",
     +                    " parent state is different from the"
       WRITE(6,"(A)")     " reference state."
       ENDIF
C
C AB SPIN CASE.
C
cjp j,b are beta
cjp
      imax=0
      jmax=0
      amax=0
      bmax=0
      denmax=0.
      NDSSIZ=NVRTA*NVRTB
      NUMDIS=NOCA*NOCB
      I000=1
      I010=I000+2*NIRREP
      I020=I010+NUMDIS
      I030=I020+2*NIRREP
      I040=I030+NDSSIZ
      CALL GETREC(20,'JOBARC','SOAOB2X ',NIRREP,ISCR(I000))
      CALL GETREC(20,'JOBARC','SOAOB2  ',NUMDIS,ISCR(I010))
      CALL GETREC(20,'JOBARC','SVAVB2X ',NIRREP,ISCR(I020))
      CALL GETREC(20,'JOBARC','SVAVB2  ',NDSSIZ,ISCR(I030))
      IOFFL=0
      IOFFR=0
      DO 10 IRREP=1,NIRREP
       NDSSYM=ISCR(I000+IRREP-1)
       DISSYM=ISCR(I020+IRREP-1)
       ITHRU=0
       DO 20 IDIS=1,NDSSYM
        IDISSYM=IDIS+IOFFR
        IDISFULL=ISCR(I010+IDISSYM-1)
        J=1+(IDISFULL-1)/NOCA
        I=IDISFULL-(J-1)*NOCA
        DO 30 IELM=1,DISSYM
         ITHRU=ITHRU+1
         IELMSYM=IELM+IOFFL
         IELMFULL=ISCR(I030+IELMSYM-1)
         B=1+NOCB+(IELMFULL-1)/NVRTA
         IB=B-NOCB
         A=NOCA+IELMFULL-(IB-1)*NVRTA
cjp modification for mr-bwcc
cjp store ithru index of the amplitudes which become offdiagonal H elements!!
         if(isbwcc) then
                itrue=invhnum(i,iref,1)
                jtrue=invhnum(j,iref,2)
                atrue=invpnum(a-noca,iref,1)
                btrue=invpnum(b-nocb,iref,2)
                if(isactive(itrue,1).and.isactive(jtrue,2)
     +            .and.isactive(atrue,1).and.isactive(btrue,2)) then
                d(ithru)=1d0/denomblow
cjp find for which pair of references it plays the role of Heff(offdiagonal)
cjp and store this info
cjp we collect information for all such elements, in spite that only
cjp those where iref plays role are relevant
cjp Note that ithru for given excitation is specific for each iref
cjp the tables must thus be iref-replicated
cjp for rhf only i=j, a=b cases are relevant
cjp otherwise more general treatment is necessary
          if(nref.gt.1 .and. (iuhf.eq.1 .or. itrue.eq.jtrue
     +                 .and. atrue.eq.btrue)) then
cjp			candidate for heff element (offdiagonal)
                     efrom(1)=itrue
                     efrom(2)=jtrue
                     eto(1)=atrue
                     eto(2)=btrue
                     espin(1)=1
                     espin(2)=2
                     k=ifindref(nbas,2,efrom,eto,espin)
                     if(k.gt.0) then
        internnum(iref,3)=internnum(iref,3)+1
        internindex(internnum(iref,3),iref,3)=ithru
        internfrom(internnum(iref,3),iref,3)=iref
        internto(internnum(iref,3),iref,3)=k
                     endif
          endif
cjp end heff offdiagonal element
                else
         d(ithru)=1.0/(e(itrue)+e(nbas+jtrue)
     +                       -e(atrue)-e(nbas+btrue))
                 if(abs(d(ithru)).gt.abs(denmax)
     +            .and.abs(d(ithru)).lt.intruder ) then
                 imax=itrue
                 jmax=jtrue
                 amax=atrue
                 bmax=btrue
                 denmax=d(ithru)
                 endif
                endif
         else
         D(ITHRU)=1.0/(E(I)+E(NBAS+J)-E(A)-E(NBAS+B))
         endif
cjp end
30      CONTINUE
20     CONTINUE
       IOFFR=IOFFR+NDSSYM
       IOFFL=IOFFL+DISSYM
       CALL PUTLST(D,1,NDSSYM,1,IRREP,50)
       IF(RLE)THEN 
        DO 15 IINV=1,ITHRU
         D(IINV)=1.0/D(IINV)
15      CONTINUE
        CALL PUTLST(D,1,NDSSYM,1,IRREP,66)
       ENDIF
10    CONTINUE
cjp
cjp intruder diagnostics
cjp
       if(isbwcc) then
       write(6,*)'@DIJAB-I: AB MAX denomin.: ref,spin,i,j,a,b,den',
     & iref,3,imax,jmax,amax,bmax,denmax
       if(abs(denmax).gt.abs(totmaxdenom)) totmaxdenom=denmax
       endif

C
C AA AND BB SPIN CASES.
C
      DO 99 ISPIN=1,1+IUHF
       IOFF=(ISPIN-1)*NBAS
       NDSSIZ=NNM1O2(NVRTO(ISPIN))
       NUMDIS=NNM1O2(NOCCO(ISPIN))
       NONV=NVRTO(ISPIN)*NOCCO(ISPIN)
       I000=1
       I010=I000+NIRREP
       I020=I010+NUMDIS
       I030=I020+NIRREP
       I040=I030+NDSSIZ
       IF(ISPIN.EQ.1)THEN
        CALL GETREC(20,'JOBARC','SOAOA0X ',NIRREP,ISCR(I000))
        CALL GETREC(20,'JOBARC','SOAOA0  ',NUMDIS,ISCR(I010))
        CALL GETREC(20,'JOBARC','SVAVA0X ',NIRREP,ISCR(I020))
        CALL GETREC(20,'JOBARC','SVAVA0  ',NDSSIZ,ISCR(I030))
       ELSEIF(ISPIN.EQ.2)THEN
        CALL GETREC(20,'JOBARC','SOBOB0X ',NIRREP,ISCR(I000))
        CALL GETREC(20,'JOBARC','SOBOB0  ',NUMDIS,ISCR(I010))
        CALL GETREC(20,'JOBARC','SVBVB0X ',NIRREP,ISCR(I020))
        CALL GETREC(20,'JOBARC','SVBVB0  ',NDSSIZ,ISCR(I030))
       ENDIF
       IOFFL=0
       IOFFR=0
cjp
       imax=0
       jmax=0
       amax=0
       bmax=0
       denmax=0.
       DO 110 IRREP=1,NIRREP
        NDSSYM=ISCR(I000+IRREP-1)
        DISSYM=ISCR(I020+IRREP-1)
        ITHRU=0
        DO 120 IDIS=1,NDSSYM
         IDISSYM=IDIS+IOFFR
         IDISFULL=ISCR(I010+IDISSYM-1)
         I=ILRG(IDISFULL)
         J=IDISFULL-NNM1O2(I-1)
         DO 130 IELM=1,DISSYM
          ITHRU=ITHRU+1
          IELMSYM=IELM+IOFFL
          IELMFULL=ISCR(I030+IELMSYM-1)
          A=NOCCO(ISPIN)+ILRG(IELMFULL)
          B=NOCCO(ISPIN)+IELMFULL-NNM1O2(ILRG(IELMFULL)-1)
cjp
cjp modification for mr-bwcc
         if(isbwcc) then
                itrue=invhnum(i,iref,ispin)
                jtrue=invhnum(j,iref,ispin)
                if(ispin.eq.1) then
                 nocx=noca
                else
                 nocx=nocb
                endif
                atrue=invpnum(a-nocx,iref,ispin)
                btrue=invpnum(b-nocx,iref,ispin)
                if(isactive(itrue,ispin).and.isactive(jtrue,ispin)
     +    .and.isactive(atrue,ispin).and.isactive(btrue,ispin)) then
                d(ithru)=one/denomblow
cjp candidate for Heff matrix element
cjp find for which pair of references it plays the role of Heff(offdiagonal)
cjp and store this info
                if(nref.gt.1 .and. iuhf.eq.1 ) then
                     efrom(1)=itrue
                     efrom(2)=jtrue
                     eto(1)=atrue
                     eto(2)=btrue
                     espin(1)=ispin
                     espin(2)=ispin
                     k=ifindref(nbas,2,efrom,eto,espin)
                     if(k.gt.0) then
cjp                    write the information to array and exit searching
                   internnum(iref,ispin)=internnum(iref,ispin)+1
                   internindex(internnum(iref,ispin),iref,ispin)=ithru
                   internfrom(internnum(iref,ispin),iref,ispin)=iref
                   internto(internnum(iref,ispin),iref,ispin)=k
                  endif
                endif
cjp end treatment of Heff matrix element
                else
                d(ithru)=one/(e(ioff+itrue)+e(ioff+jtrue)
     +                   -e(ioff+atrue)-e(ioff+btrue))
                 if(abs(d(ithru)).gt.abs(denmax)
     +            .and.abs(d(ithru)).lt.intruder ) then
                 imax=itrue
                 jmax=jtrue
                 amax=atrue
                 bmax=btrue
                 denmax=d(ithru)
                 endif
                endif
         else
CSSS          Write(6,"(a,4(1x,i4))") "I,J,A,B: ", I,J,A,B
          D(ITHRU)=ONE/(E(IOFF+I)+E(IOFF+J)-E(IOFF+A)-E(IOFF+B))
          endif
cjp end
130      CONTINUE
120     CONTINUE
        CALL PUTLST(D,1,NDSSYM,1,IRREP,47+ISPIN)
        IF(RLE)THEN 
         DO 115 IINV=1,ITHRU
          D(IINV)=ONE/D(IINV)
115       CONTINUE
         CALL PUTLST(D,1,NDSSYM,1,IRREP,63+ISPIN)
        ENDIF
        IOFFR=IOFFR+NDSSYM
        IOFFL=IOFFL+DISSYM
110    CONTINUE
cjp
cjp intruder diagnostics
       if(isbwcc) then
       write(6,*)'@DIJAB-I: AA,BB MAX denomin.: ref,spin,i,j,a,b,den',
     & iref,ispin,imax,jmax,amax,bmax,denmax
       if(abs(denmax).gt.abs(totmaxdenom)) totmaxdenom=denmax
       endif
cjp
C
C NOW MAKE SYMMETRY PACKED EPS(I)-EPS(A) LISTS FOR RLE.
C The following line was commented to get lambda code properly
C working when no convergence accelerations are used. The CC
C codes have no effect from this since it built denominator arrays
C it self. However, the lambda code relies on having denominators
C built here. I have no idea why this was done only only when
C RLE and DISS is on in the original implementation. The list
C (9, 64) is used in l2inl1.f in lambda code. Ajith Perera 12/2003.
C 
C       IF(RLE)THEN
C
        I000=1
        I010=I000+NIRREP
        I020=I010+NONV
        IF(ISPIN.EQ.1)THEN
         CALL GETREC(20,'JOBARC','SVAOA2X ',NIRREP,ISCR(I000))
         CALL GETREC(20,'JOBARC','SVAOA2  ',NONV,ISCR(I010))
        ELSE
         CALL GETREC(20,'JOBARC','SVBOB2X ',NIRREP,ISCR(I000))
         CALL GETREC(20,'JOBARC','SVBOB2  ',NONV,ISCR(I010))
        ENDIF
        NDSSYM=ISCR(I000)
        CALL UPDMOI(1,NDSSYM,9,63+ISPIN,0,0)
        CALL UPDMOI(1,NDSSYM,9,47+ISPIN,0,0)
        ITHRU=0
        DO 300 NUMDIS=1,NDSSYM
         ITHRU=ITHRU+1
         IDISFULL=ISCR(I010+NUMDIS-1)
         I=1+(IDISFULL-1)/NVRTO(ISPIN)
         A=NOCCO(ISPIN)+IDISFULL-(I-1)*NVRTO(ISPIN)
cjp modification for mr-bwcc
cjp here is code only for rle, so internindex1 etc. must be done elsewhere
         if(isbwcc) then
                itrue=invhnum(i,iref,ispin)
                if(ispin.eq.1) then
                 nocx=noca
                else
                 nocx=nocb
                endif
                atrue=invpnum(a-nocx,iref,ispin)
                if(isactive(itrue,ispin).and.isactive(atrue,ispin))then
                d(ithru)=1d0/denomblow
                else
                d(ithru)=1.0/(e(ioff+itrue)-e(ioff+atrue))
                endif
         else
         D(ITHRU)=ONE/(E(IOFF+I)-E(IOFF+A))
         endif
300     CONTINUE
        CALL PUTLST(D,1,1,1,9,47+ISPIN)
        DO 301 IINV=1,ITHRU
         D(IINV)=ONE/D(IINV)
301     CONTINUE
        CALL PUTLST(D,1,1,1,9,63+ISPIN)
C the endif for the IF (RLE) commented on 12/2003, A. Perera. 
C       ENDIF
cjp end ispin loop
99    CONTINUE
cjp
cjp write info into to fockcd file (opened in focklist as unit=99)
cjp close the file here
cjp
cjp only relevant spin cases are written to the file
      if(isbwcc) then
      ibot=1
      if(iuhf.eq.0) ibot=3
      do is=ibot,3
      write(99)internnum(iref,is)
      write(99)(internfrom(j,iref,is),j=1,(maxref*(maxref-1)/2))
      write(99)(internto(j,iref,is),j=1,(maxref*(maxref-1)/2))
      write(99)(internindex(j,iref,is),j=1,(maxref*(maxref-1)/2))
      enddo
      if(iuhf.ne.0) then
      do is=1,2
      write(99)internnum1(iref,is)
      write(99)(internfrom1(j,iref,is),j=1,(maxref*(maxref-1)/2))
      write(99)(internto1(j,iref,is),j=1,(maxref*(maxref-1)/2))
      write(99)(internindex1(j,iref,is),j=1,(maxref*(maxref-1)/2))
      enddo
      endif
      close(unit=99)
      endif
      RETURN
      END
