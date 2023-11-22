      SUBROUTINE RCC_FINISH(ICYCLE)
C
C RUN-DOWN ROUTINE FOR CC CALCULATIONS.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*6 STAR
      CHARACTER*13 ITCALC(50)
      LOGICAL LINCC,CICALC

      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
cjp
      logical o
      COMMON /LINEAR/ LINCC,CICALC
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /ENERGY/ ECORR(500,2),IXTRLE(500)
      COMMON /FLAGS/ IFLAGS(100)
      DATA ITCALC /
     & '             ','             ','             ','             ',
     & '  LCCD       ','  LCCSD      ',' UCCSD(4)    ','   CCD       ',
     & '  UCC(4)     ','  CCSD       ','  CCSD       ','  CCSD       ',
     & '  CCSDT-1a   ','  CCSDT-1b   ','  CCSDT-2    ','  CCSDT-3    ',
     & '  CCSDT-4    ','  CCSDT      ','  LCCSDT     ','   CCD       ',
     & '  QCISD      ','  CCSD       ','  QCISD      ','   CID       ',
     & '   CISD      ','  QCISD      ','  CCSD       ','  CCSD       ',
     & '  CCSDT      ','  CCSDT      ','  CCSD       ','  CCSD       ',
     & '     CC3     ',' CCSDT-T1T2  ','  CCSDTQ-1   ','  CCSDTQF-1  ',
     & '  CCSDTQ-2   ',' CCSDTQ-3    ','  CCSDTQ     ','  ACCSD      ',
     & '             ',' ACCSD(T)    ','  CCSD(TQf)  ','  CCSDT(Qf)  ',
     & ' OO-MP(2)    ',' OO-MBPT(2)  ','  CC2        ','  rCCD       ',
     & ' drCCD       ','             '/
C
cjp
#include "bwcc.com"
cjp
cjp in parallel run, assure that this printout will appear 
cjp at the end of output
c
cjp output only from one processor
      o=.true.

      CALL GETREC(20,'JOBARC','SCFENEG ',IINTFP,ESCF)
      if(o)WRITE(6,99)
99    FORMAT(T3,'     Summary of iterative solution of CC equations ')
      if(o)WRITE(6,100)
100   FORMAT(6X,59('-'))
      if(o)WRITE(6,101)
101   FORMAT(T24,'Correlation',T46,'Total',/,T8,'Iteration',T26,
     &       'Energy',T46,'Energy')
      if(o)WRITE(6,100)
      DO 20 I=1,ICYCLE
       STAR='JACOBI'
       IF(IXTRLE(I).EQ.1)STAR=' RLE  '
       IF(IFLAGS(21).EQ.1)STAR=' DIIS '
cjp write a note about correction of BWCC
       if(isbwcc.and.I.EQ.ICYCLE.and.ihubaccorr.eq.3) STAR=' CORR '
cjp write a note about homotopic correction
       if(isbwcc.and.i-1.gt.ihomotop.and.ihomotop.gt.0) then
        hfakt=lambdahomotop**(i-1-ihomotop)
        if(hfakt .lt. diishonset .and. IFLAGS(21).EQ.1 ) then
           STAR=' HDIIS'
        else
           STAR='HJACOB'
        endif
       endif
       IF(I.EQ.ICYCLE.AND.UCC)THEN
        if(o)WRITE(6,1001)
1001    FORMAT(T15,' Adding <T2^+T2^+WT2> energy contribution.')
       ENDIF
cjp info about homotopicsolution
       if(isbwcc .and. i-1.gt.ihomotop.and.ihomotop.gt.0) then
          if(o)WRITE(6,1002)I-1,ECORR(I,1),ECORR(I,1)+ESCF,STAR,hfakt
       else
          if(o)WRITE(6,1000)I-1,ECORR(I,1),ECORR(I,1)+ESCF,STAR
       endif
1000   FORMAT(T10,I4,T19,F18.12,T39,F19.12,T59,A6)
1002   format(T10,I4,T19,F18.12,T39,F19.12,T59,A6,1x,e9.3)
20    CONTINUE
      if(o)WRITE(6,100)
C
CJDW 10/1/96. Write out energy for iterative part of method.
C
cjp
       if(isbwcc ) then
         if(ihubaccorr.eq.3) then
        if(o)WRITE(6,1010) 'cMR-BW-'//ITCALC(IFLAGS(2)),
     &     ECORR(ICYCLE,1)+ESCF
         else
          if(ihomotop.gt.0) then
        if(o)WRITE(6,1010) 'hMR-BW-'//ITCALC(IFLAGS(2)),
     &    ECORR(ICYCLE,1)+ESCF
          else
         if(o) WRITE(6,1010) 'MR-BW-'//ITCALC(IFLAGS(2)),
     &ECORR(ICYCLE,1)+ESCF
          endif
         endif
        else
       if(o)WRITE(6,1010) ITCALC(IFLAGS(2)),ECORR(ICYCLE,1)+ESCF
        endif
C
      IF(.NOT.CICALC)THEN
cjp
       if(isbwcc) then
        if(o)write(6,210)
        else
       if(o)WRITE(6,202)
        endif
      ELSE
       if(o)WRITE(6,203)
       if(o)WRITE(6,205)
       if(o)WRITE(6,204)
      ENDIF
cjp summarize serious warnings
      do i=1,maxbwwarnings
      if(bwwarning(i).and.o)write(6,*) bwwarntext(i)
      enddo
c write total energy to jobarc
      CALL PUTREC(1,'JOBARC','TOTENERG',IINTFP,ECORR(ICYCLE,1)+ESCF)
cmn again for vibronic purposes:
      CALL PUTREC(1,'JOBARC','TOTENER2',IINTFP,ECORR(ICYCLE,1)+ESCF)
      CALL PUTREC(1,'JOBARC','PARENERG',IINTFP,ECORR(ICYCLE,1)+ESCF)
c
210   FORMAT(T7,'A multi-miracle has come to pass. ',
     &          'The MR-BW-CC iterations have converged.')
202   FORMAT(T7,'A miracle has come to pass. ',
     &          'The CC iterations have converged.')
203   FORMAT(T7,'A miracle has come to pass. ',
     &          'The CI iterations have converged.')
205   FORMAT(/,T12, ' ****************WARNING!!!!************* ')
204   FORMAT(T3,'For the most part, the authors of this program ',
     &         'feel that CI is ',/
     &         ,T3,'a bad method.  Please run CC calculations in the',
     &         ' future!',/)
 1010 FORMAT(/,T14,A13,' energy is ',F20.12,' a.u. ',/)
      RETURN
      END
