      SUBROUTINE E4SENG(W,D,DEN,E,NTSIZ,POP,VRT,NOCC,NVRT,SE4,ISPIN)
C
C     THIS ROUTINE CALCULATES THE NEW T1 AMPLITUDES AND
C     FOR FOURTH ORDER PERTURBATION THEORY THE SINGLE
C     CONTRIBUTION TO THE ENERGY
C
C     THIS ROUTINE HAS TO BE CALLED SEPARATELY FOR ALPHA AND BETA
C     CASE
cjp
cjp in bwcc, this routine is called twice. in the first pass it
cjp collects matrix elements of Heff, which correspond to 
cjp single excitations (uhf only), in the second pass, called after
cjp Heff diagonalalization is complete, the amplitudes are finally updated
cjp
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINCC,CICALC
      INTEGER DIRPRD,POP,VRT
cjp
      integer itrue,atrue,intptr
      DIMENSION W(NTSIZ),E(NOCC+NVRT),
     &          D(NTSIZ),DEN(NTSIZ),POP(8),VRT(8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /LINEAR/ LINCC,CICALC
      COMMON /CORENG/ ELAST
C
      EQUIVALENCE(IFLAGS(2),METHOD)
C
C
cjp
#include "bwcc.com"
C
C     GET ORBITAL ENERGIES FROM DISK
C
      IF(ISPIN.EQ.1)THEN
       CALL GETREC(20,'JOBARC','SCFEVALA',IINTFP*(NOCC + NVRT),
     &             E)
      ELSE
       CALL GETREC(20,'JOBARC','SCFEVALB',IINTFP*(NOCC + NVRT),
     &             E)
      ENDIF
cjp fix diagonal fock elements
      if(isbwcc) then
        do i=1,nocc+nvrt
          e(i)=e(i)+fockcd(i,iref,ispin)
        enddo
      endif
C
C     LOOP OVER IRREPS 
C
      IND=0
      INDI=0
      INDA0=0
C
      DO 100 IRREP=1,NIRREP
C
      NOCCI=POP(IRREP)
      NVRTI=VRT(IRREP)
c
c
C
cjp
cjp we do this both independently on ibwpass, although it might be splitted 
cjp in the 2 passes
cjp identify also the internal amplitudes
cjp and check them if they come out zero;
cjp set the denominator to blowdenom
cjp the following is also completely wrong if symmetry is on
cjp NOTE for parallelization ... parts of these arrays are meaningless
      intptr=1
      DO  10 I=1,NOCCI
      INDI=INDI+1
      INDA=INDA0
      DO  10 IA=1,NVRTI
      INDA=INDA+1
      IND=IND+1
      if(isbwcc) then
           itrue=invhnum(indi,iref,ispin)
           atrue=invpnum(inda,iref,ispin)
           if(isactive(itrue,ispin).and.isactive(atrue,ispin))then
             den(ind)=denomblow
cjp
cjp collect matrix elements of effective hamiltonian
cjp
          if(ibwpass.eq.1.and.nref.gt.1) then
            if(intptr.le.internnum1(iref,ispin) .and.
     +      internindex1(intptr,iref,ispin).eq.ind) then
cjp             this is offdiag element, but maybe such one which does not belong to current iref
                if(iref.eq.internfrom1(intptr,iref,ispin)) then
cjp                     relevant one
                        heff(internto1(intptr,iref,ispin),iref)=w(ind)
                endif
            intptr=intptr+1
cjp         they are sorted acc. to increasing internindex, look simply at next
           else
c
cjp if such internal excitation is not included in the model space
cjp notice the user about it
             if(abs(w(ind)).gt.1d-8) then
        write(6,*)'@E4SENG-W: ref. ',iref,
     &  ' nonzero internal T1 amplitude:',
     &               itrue,atrue,w(ind),e(itrue)-e(atrue),w(ind)
     &    /(e(itrue)-e(atrue)),
     & '  ADDITIONAL REFERENCES MIGHT BE IMPORTANT'

             bwwarning(2)=.true.
             bwwarntext(2)='@E4SENG-W: Nomnzero internal T1 amplitudes
     & encountered, additional references might be important'
        if(useeq429) stop 'algorithm to treat  nonzero 
     & internal T1s not activated'
        endif
              endif
cjp             intptr.le.internnum1(iref,ispin) .and. ...
            endif
cjp             ibwpass.eq.1.and.nref.gt.1
cjp end collect matrix elements of effective hamiltonian
cjp
cjp    was not internal excitation
           else
             den(ind)=e(itrue)-e(atrue)
           endif
      else
      DEN(IND) = E(INDI) - E(NOCC + INDA)
      endif
   10 CONTINUE
      INDA0=INDA0+NVRTI
100   CONTINUE
C
C     CALCULATE NEW AMPLITUDES
C
cjp do not do it if ibwpass ==1
        if(ibwpass.eq.1) then
cjp store the intermediates as if they were final amplitudes
cjp they will be needed anyway in the general term of eq 4.28 in Hubac paper
        do i=1,ntsiz
        d(i)=w(I)
        enddo
        return
        endif
      FACT=0.0
      IF(CICALC)FACT=ELAST
cjp BW-CC shift of the denominator
      if(isbwcc) fact=ecorrbw
      DO 20 I=1,NTSIZ
C
      tmp=1.d0/(DEN(I)+FACT)
cjp intruder prevention
      if(isbwcc .and. abs(tmp).gt.intruder) tmp=0.d0
      D(I)=W(I)*tmp
20    CONTINUE
C
C     FOR FINITE ORDER PERTURBATION THEORY 
C     CALCULATE ENERGY
C
      IF(METHOD.EQ.3.OR.METHOD.EQ.4) THEN
       SE4=SE4+SDOT(NTSIZ,D,1,W,1)
      ENDIF
      RETURN
      END
