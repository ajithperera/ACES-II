      SUBROUTINE active_SMICAN(FOCK,SCR1,SCR2,SCR3,EVEC,EVAL,
     &     POP,VRT,iclosed, iactive, NIRREP)
C
C THIS SUBROUTINE PERFORMS A TRANSFORMATION TO Active space SEMICANONICAL ORBITALS,
C  WHICH ARE DEFINED BY
C
C             F(IJ) = 0 (I.NE.J)
C             F(MN) = 0 (M.NE.N)
C             F(AB) = 0 (A.NE.B)
C
C THE FOCK MATRIX WHICH IS PASSED IN MUST BE IN THE AO BASIS.
C
CEND
      IMPLICIT INTEGER (A-Z)
C
      CHARACTER*1 SP(2)
      CHARACTER*8 WLABEL(2,2)
      DOUBLE PRECISION FOCK,SCR1,SCR2,SCR3,SCR4,SCR5,EVEC,ONE,ZILCH,EVAL
      DOUBLE PRECISION SCR,OSCR2
      DIMENSION FOCK(1),SCR1(1),SCR2(1),SCR3(1),EVEC(1),EVAL(1)
      DIMENSION POP(NIRREP),VRT(NIRREP), iclosed(nirrep),
     $     iactive(nirrep)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/ IFLAGS(100)
C
      DATA ONE /1.0/
      DATA ZILCH /0.0/
      DATA WLABEL /'WMATRXOA','WMATRXVA','WMATRXOB','WMATRXVB'/
      DATA SP /'A','B'/
      logical print
C
      NNP1O2(I)=(I*(I+1))/2
      print = .true.
C
C TRANSFORM AO BASIS FOCK MATRIX TO THE ORTHOGONAL BASIS
C
      write(6,*) ' Enter active_smican '
c
      LENGTH=0
      DO 222 I=1,NIRREP
       LENGTH=LENGTH+NNP1O2(POP(I)+VRT(I))
222   CONTINUE
      CALL PUTREC(20,'JOBARC','FOCK    ',LENGTH*IINTFP,FOCK)
      CALL FAO2MO(FOCK,EVEC,SCR2,SCR3,POP,VRT,NIRREP)
C
      IPOS=1
      IOFFF=1
      IOFFE=1
      IOFFEVL=1
C
C LOOP OVER IRREPS
C
      IOFFWO=1
      IOFFWV=1
      DO 10 IRREP=1,NIRREP
       NOCC=POP(IRREP)
       NVRT=VRT(IRREP)
       NSIZ=NOCC+NVRT
       nclosed = iclosed(irrep)
       nactive = iactive(irrep)
       CALL EXPND2(FOCK(IOFFF),SCR1,NSIZ)
#ifdef _DEBUG_LVLM1
      Write(6,*) "The MO basis FOCK matrix in standard orbitals"
      call output(Scr1, 1, NSIZ, 1, NSIZ, NSIZ, NSIZ, 1)
#endif
C
C DIAGONALIZE closed-closed BLOCK OF FOCK MATRIX
C
       CALL BLKCPY2(SCR1,NSIZ,NSIZ,SCR2,nclosed,nclosed,1,1)
       CALL EIG(SCR2,SCR3,nclosed,nclosed,0)
cmn       CALL SCOPY(nclosed*nclosed,SCR3,1,SCR4(IOFFWO),1)
cmn       IOFFWO=IOFFWO+nclosed*nclosed
       CALL SCOPY(nclosed,SCR2,nclosed+1,EVAL(IOFFEVL),1)
       do i = 1, nclosed
          eval(ioffevl+i-1) = eval(ioffevl+i-1) - 100.0d0
       enddo
       if (print) then
       Write(6,*)
       Write(6,*) "The closed-shell eigenvalues, shifted by -100 a.u. ",
     $      irrep, nclosed
       call output(eval(ioffevl), 1, 1, 1, nclosed, 1, nclosed, 1)
c       Write(6,"(6F12.8)") (EVAL(IOFFEVL - 1+k),k=1,nclosed)
       endif

       IOFFEVL=IOFFEVL+nclosed
C
C NOW USE EIGENVECTORS TO ROTATE MOLECULAR ORBITALS
C
       CALL XGEMM('N','N',NSIZ,Nclosed,Nclosed,ONE,EVEC(IOFFE),NSIZ,
     &            SCR3,nclosed,ZILCH,SCR2,NSIZ)
       CALL SCOPY(NSIZ*nclosed,SCR2,1,EVEC(IOFFE),1)
C
       IOFFE=IOFFE+NSIZ*nclosed


C
C NOW DIAGONALIZE active-active BLOCK OF FOCK MATRIX
C
       CALL BLKCPY2(SCR1,NSIZ,NSIZ,SCR2,nactive,nactive,
     $      nclosed+1,nclosed+1)
       CALL EIG(SCR2,SCR3,nactive,nactive,0)
cmn       CALL SCOPY(nactive*nactive,SCR3,1,SCR5(IOFFWV),1)
cmn       IOFFWV=IOFFWV+nactive*nactive
       CALL SCOPY(nactive,SCR2,nactive+1,EVAL(IOFFEVL),1)
       do i = 1, nactive
          eval(ioffevl+i-1) = eval(ioffevl+i-1) - 10.0d0
       enddo
       if (print) then
       Write(6,*)
       Write(6,*) "The active eigenvalues, shifted by -10 a.u. ",
     $      "in irrep", irrep, nactive
       call output(eval(ioffevl), 1, 1, 1, nactive, 1, nactive, 1)
       endif
       IOFFEVL=IOFFEVL+nactive
       CALL XGEMM('N','N',NSIZ,nactive, nactive,ONE,EVEC(IOFFE),NSIZ,
     &            SCR3,nactive,ZILCH,SCR2,NSIZ)
       CALL SCOPY(NSIZ*nactive,SCR2,1,EVEC(IOFFE),1)
C
       IOFFE=IOFFE+NSIZ*nactive
C

C
C NOW DIAGONALIZE VIRTUAL-VIRTUAL BLOCK OF FOCK MATRIX
C
       CALL BLKCPY2(SCR1,NSIZ,NSIZ,SCR2,NVRT,NVRT,NOCC+1,NOCC+1)
       CALL EIG(SCR2,SCR3,NVRT,NVRT,0)
cmn       CALL SCOPY(NVRT*NVRT,SCR3,1,SCR5(IOFFWV),1)
cmn       IOFFWV=IOFFWV+NVRT*NVRT
       CALL SCOPY(NVRT,SCR2,NVRT+1,EVAL(IOFFEVL),1)
       if (print) then
       Write(6,*)
       Write(6,*) "The virtual eigenvalues in irrep", irrep, nvrt
       call output(eval(ioffevl), 1, 1, 1, nvrt, 1, nvrt, 1)
       endif
       IOFFEVL=IOFFEVL+NVRT
       CALL XGEMM('N','N',NSIZ,NVRT,NVRT,ONE,EVEC(IOFFE),NSIZ,
     &            SCR3,NVRT,ZILCH,SCR2,NSIZ)
       CALL SCOPY(NSIZ*NVRT,SCR2,1,EVEC(IOFFE),1)
C
       IOFFE=IOFFE+NSIZ*NVRT
C
       IOFFF=IOFFF+NNP1O2(NSIZ)
       IPOS=IPOS+NVRT+NOCC
10    CONTINUE
C
C  Now dump the canonical to semi-canonical transformation matrices.
C
cmn      IF(SEMI) THEN
cmn        CALL PUTREC(20,'JOBARC',WLABEL(1,ISPIN),(IOFFWO-1)*IINTFP,SCR4)
cmn        CALL PUTREC(20,'JOBARC',WLABEL(2,ISPIN),(IOFFWV-1)*IINTFP,SCR5)
cmn      ENDIF
C
c restore original fock matrix in AO basis
c
c write eval to jobarc
c
      nsiz = 0
      do irrep = 1, nirrep
         nsiz = nsiz + pop(irrep) + vrt(irrep)
      enddo
      call putrec(20,'JOBARC','NSMICAN ', 1, nsiz)
      call putrec(20,'JOBARC','EVALSEMI', nsiz*IINTFP, eval)
c
      CALL GETREC(20,'JOBARC','FOCK    ',LENGTH*IINTFP,FOCK)
      RETURN
      END
