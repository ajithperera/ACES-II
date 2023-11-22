
      SUBROUTINE SCFIT(ONEH,   XFORM,  DENSO,  FOCK,
     &                 EVAL,   EVEC,   DENSN,  SCR1,
     &                 SCRA,   SCRB,   SCRTMP, COLD,
     &                 SOVRLP, REPULS, LDIM1,  LDIM2,
     &                 LDIM3,  NOCONV, TOL,    ITER,
     &                 DMAX,   IUHF,   NBAS,
c QC-SCF args
     &                 FVEC,   XTRIAL, CVECS,  SOVER,
     &                 RESIDUE,D1,     D2,     D3,
     &                 D4,     RAA,    G1AA,   G2AA,
     &                 NPRVOCC,IPRVOCC,
     &                 NDIM,   NDET,   QC_THRESH, QC_PHASE,
     &                 scrsqr,scfks,screvec,scrsort,screval,
     &                 scraks,naobasfn,scrocc)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER SWAP,READMO,WRITMO
      INTEGER ILINOF
      LOGICAL NOCONV
      LOGICAL GSSOPT,GSSALW,GSSALT,GSSLOK,GSSRED,GSSWRT,GSSUFR
      LOGICAL ROHFMO
      DIMENSION ONEH(LDIM1),XFORM(LDIM1),DENSO((IUHF+1)*LDIM1)
      DIMENSION FOCK((IUHF+1)*LDIM1),DENSN((IUHF+1)*LDIM1)
      DIMENSION EVAL((IUHF+1)*NBAS),EVEC((IUHF+1)*LDIM3)

C prakash ks variable
       logical scfks
       dimension screval(nbas,iuhf+1),scrsqr(naobasfn*naobasfn)
       integer scrocc(nbas,2)
       dimension screvec(nbas,nbas),scrsort(naobasfn,nbas,2)
       dimension scraks(naobasfn*naobasfn)
        character*1 spinc(2)
#include<aces.h>
#include<maxbasfn.par>
C
C     COLD and SOVRLP are used when orbital locking (monitoring of
C     C(OLD)T * S * C) has been requested. In this case they have
C     length NBAS**2. When orbital locking has not been requested, they
C     are equivalenced to SCRB and SCRTMP. Note that the allocation of
C     memory for SCR1,SCRA,SCRB, and SCRTMP in VSCF is not the same as
C     in OOFCOR. The use of LDIM2 below is consistent with OOFCOR
C     allocation, but not with VSCF allocation.

      DIMENSION SCR1(LDIM2),SCRA(LDIM2),SCRB(LDIM2),SCRTMP(LDIM2),
     &          COLD(1),SOVRLP(1)
      DIMENSION NOCC(16), NSUM(16)
c QC-SCF args
      DIMENSION FVEC(NDIM),XTRIAL(NDIM),
     &          CVECS(NBAS*NBAS),SOVER(NBAS*NBAS),
     &          RESIDUE(NDIM),D1(NDIM),
     &          D2(NDIM),D3(NDIM),DOCC(MAXBASFN*2),
     &          D4(NDIM),RAA(LDIM3*(IUHF+1)),
     &          G1AA(LDIM3*(IUHF+1)),G2AA(LDIM3*(IUHF+1)),NDET(16)
      DIMENSION NPRVOCC(16)
      LOGICAL QC_PHASE
C
CJDW 9/18/96. ILINOF is Offset array to point to beginning of each
C             symmetry block in EVAL.
      DIMENSION ILINOF(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
#include "flags.h"
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /POPUL/ NOCC
#include "symm2.com"
C
      COMMON /GSCOMA/ GSSOPT,GSSALW,GSSALT,GSSLOK,GSSRED,GSSWRT,GSSUFR
      COMMON /GSCOMB/ SWAP(4,8,2),LOCK(8,2),IPRTGS(8,2),ISTOP1,ISTOP2,
     1                READMO,WRITMO,IUHFRHF,LUGSS
      COMMON /GSCOMC/ ROHFMO
      COMMON /LSHCOM/ ALPHA1,BETA1
C
      DATA ONEM,ONE,AZERO /-1.D0,1.D0,0.D0/


       INDX2(I,J,N)=I+(J-1)*N
C
        spinc(1)='A'
        spinc(2)='B'
C
C-----------------------------------------------------------------------
CJDW 9/18/96. Compute ILINOF.
C
      DO 10 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
       ILINOF(IRREP) = 0
      ELSE
       ILINOF(IRREP) = ILINOF(IRREP-1) + NBFIRR(IRREP-1)
      ENDIF
   10 CONTINUE
C-----------------------------------------------------------------------
C QC-SCF (will skip the regular orbital-generation part)
      IF (QC_PHASE) THEN
         CALL QCSCF(NBAS,NDIM,NDET,FOCK,LDIM1,LDIM2,LDIM3,IUHF,
     &              EVEC,FVEC,XTRIAL,CVECS,SOVER,
     &              SCR1,SCRA,SCRB,SCRTMP,RESIDUE,
     &              D1,D2,D3,D4,RAA,G1AA,G2AA)
         GOTO 112
      END IF
C-----------------------------------------------------------------------
      IF(GSSLOK)THEN
C
C     Get the overlap matrix in the SO basis.
C
      CALL GETREC(20,'JOBARC','AOOVRLAP',NBAS*NBAS*IINTFP,SOVRLP)
C     WRITE(6,1010)
C1010 FORMAT(' @SCFIT-I, Overlap matrix in SO basis. ')
C     CALL OUTPUT(SOVRLP,1,NBAS,1,NBAS,NBAS,NBAS,1)
C
      ENDIF
C
C  Transform the Fock matrix to the orthogonal basis; follow this
C  by diagonalizing the new Fock matrix and obtain the eigenvectors
C  in the orthogonal basis.  Transform the coefficients to the
C  non-orthogonal basis.
C
      ICNT=0
      CALL ZERO(DENSN,LDIM1)
C
      NFOCK=IUHF+1
      IF(IFLAGS(11).EQ.2) NFOCK=1 

      DO 100 ISPIN=1,NFOCK
C
      IF(GSSLOK)THEN
C
C     Get full eigenvectors from previous iteration.
C
      IF(ISPIN.EQ.1)THEN
CSSS       CALL GETREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,COLD)
       CALL GETREC(20,'JOBARC','GUESSVA0',NBAS*NBAS*IINTFP,COLD)
C     WRITE(6,1020)
C1020 FORMAT(' @SCFIT-I, Previous alpha MOs in SO basis. ')
C     CALL OUTPUT(COLD,1,NBAS,1,NBAS,NBAS,NBAS,1)
      ELSE
CSSS       CALL GETREC(20,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,COLD)
       CALL GETREC(20,'JOBARC','GUESSVB0',NBAS*NBAS*IINTFP,COLD)
C     WRITE(6,1030)
C1030 FORMAT(' @SCFIT-I, Previous beta  MOs in SO basis. ')
C     CALL OUTPUT(COLD,1,NBAS,1,NBAS,NBAS,NBAS,1)
      ENDIF
C
      ENDIF
C
ctest  IF(ITER.NE.1.AND.IFLAGS(11).EQ.2) THEN
       IF(ROHFMO   .AND.IFLAGS(11).EQ.2) THEN
CSSS        Write(6,*) "The MOVECTOR", ISPIN
CSSS        Write(6, "(7F10.5)") (EVEC(I), i=1,ISQRLN(NIRREP+1))
        CALL GETREC(20,'JOBARC','MOVECTOR',IINTFP*ISQRLN(NIRREP+1),EVEC)
       ENDIF
        DO 101 I=1,NIRREP
          IF(NBFIRR(I).EQ.0)GOTO 101
          CALL EXPND2(FOCK(((ISPIN-1)*LDIM1)+ITRIOF(I)),SCR1,
     &                NBFIRR(I))
ctest     IF(ITER.EQ.1.OR.IFLAGS(11).NE.2) THEN
          IF(.NOT.ROHFMO.OR.IFLAGS(11).NE.2) THEN
           CALL EXPND2(XFORM(ITRIOF(I)),SCRA,NBFIRR(I))
           CALL TRANSP(SCRA,SCRB,NBFIRR(I),NBFIRR(I))
           CALL MXM(SCRB,NBFIRR(I),SCR1,NBFIRR(I),SCRTMP,NBFIRR(I))
           CALL MXM(SCRTMP,NBFIRR(I),SCRA,NBFIRR(I),SCR1,NBFIRR(I))
          ENDIF
          CALL ZERO(SCRB,LDIM2)
C
C     ----- ROHF level shifting -----
C
      if(iflags(11).eq.2)then
        nd = min(nocc(i),nocc(i+8))
        ns = abs(nocc(i) - nocc(i+8))
        call eigshf(scr1,nbfirr(i),nd,ns,alpha1,beta1)
      endif
          CALL EIG(SCR1,SCRB,NBFIRR(I),NBFIRR(I),0)
#ifdef _DEBUG_LVLM
          write(6,*) "The Fock matrix @SCFIT:"
          call realprt(scr1, nbfirr(i)*nbfirr(i)) 
          Write(6,*) "The eigenvectors @SCFIT:"
          call realprt(scrb, nbfirr(i)*nbfirr(i)) 
#endif 
C
C  Copy the eigenvalues into the EVAL vector.
C
          ICNT0=ICNT
          DO 110 K=1,NBFIRR(I)
            ICNT=ICNT+1
            EVAL(ICNT)=SCR1(INDX2(K,K,NBFIRR(I)))
  110     CONTINUE
          IF(IFLAGS(11).EQ.2) THEN
           ICNT=ICNT0
           DO 111 K=1,NBFIRR(I)
            ICNT=ICNT+1
            EVAL(ICNT+NBAS)=EVAL(ICNT)
  111      CONTINUE
          ENDIF
C
C  Now get the coefficients for the non-orthogonal basis and form
C  the new density matrix.
C
ctest     IF(ITER.EQ.1.OR.IFLAGS(11).NE.2) THEN
          IF(.NOT.ROHFMO.OR.IFLAGS(11).NE.2) THEN
           CALL MXM(SCRA,NBFIRR(I),SCRB,NBFIRR(I),SCR1,NBFIRR(I))
          ELSE
           CALL XGEMM('N','N',NBFIRR(I),NBFIRR(I),NBFIRR(I),ONE,
     &                EVEC(ISQROF(I)),NBFIRR(I),SCRB,NBFIRR(I),AZERO,
     &                SCR1,NBFIRR(I))
          ENDIF
C
C   check on signs of eigenvectors
C
          call chksgn(scr1,nbfirr(i))
C
c          call expnd2(oneh(itriof(i)),scrtmp,nbfirr(i))
c          one=1.
c          zilch=0.
c          call xgemm('T','N',nbfirr(i),nbfirr(i),nbfirr(i),one,
c     &               scr1,nbfirr(i),scrtmp,nbfirr(i),zilch,
c     &               scrb,nbfirr(i))
c          call xgemm('N','N',nbfirr(i),nbfirr(i),nbfirr(i),one,
c     &               scrb,nbfirr(i),scr1,nbfirr(i),zilch,scrtmp,
c     &               nbfirr(i))
C
C-----------------------------------------------------------------------
CJDW 9/18/96.
C
C  Non-orthogonal basis set eigenvectors are in SCR1. Determine which
C  ones (if any) have zero norm. These eigenvectors are shifted to the
C  end of the symmetry block and the "higher" eigenvectors are shifted
C  down. The corresponding thing is done with the eigenvalues, but the
C  zeros are replaced by 1.0D+06. These steps are intended to hide the
C  unused eigenvectors away so they do not interfere. Note that there
C  is a possible complication here if one wants to drop top virtuals
C  and we have a linear dependency. One may simply have to drop more
C  orbitals. Note also that the algorithm assumes that eigenvectors
C  and eigenvalues are ordered starting with most negative eigenvalue.
C
C  NLINDP --- number of linear dependencies in this symmetry block. This
C             is measured by the number of eigenvectors with essentially
C             zero norms.
C  NNEG   --- number of eigenvectors with negative eigenvalues. Note
C             this is initialized to NBAS and only has meaning if there
C             is at least one zero eigenvector.
C  NSHIFT --- number of eigenvectors/eigenvalues which have to shifted
C             down. Equal to number of eigenvectors with positive eigen-
C             value.
C  NFIRST --- First vector/value to be shifted down.
C
      NLINDP = 0
      NNEG   = NBAS
      CTOL = 1.0D-06
      DO 20 IBAS=1,NBFIRR(I)
      CNORM = 0.0D+00
      CNORM = SNRM2(NBFIRR(I),SCR1(1+(IBAS-1)*NBFIRR(I)),1)
c     WRITE(LUOUT,1050) I,IBAS,CNORM
      IF(CNORM.LT.CTOL)THEN
       NLINDP = NLINDP+1
       IF(NNEG.EQ.NBAS)THEN
        NNEG = IBAS-1
       ENDIF
      ENDIF
   20 CONTINUE
C
C  Do the shifting.
C
      NSHIFT = NBFIRR(I) - NNEG - NLINDP
      NFIRST = NNEG + NLINDP + 1
      IF(NLINDP.GT.0)THEN
       DO 30 IBAS=NFIRST,NBFIRR(I)
       CALL DCOPY(NBFIRR(I),SCR1(1+(IBAS-1       )*NBFIRR(I)),1,
     &                      SCR1(1+(IBAS-1-NLINDP)*NBFIRR(I)),1)
       EVAL((ISPIN-1)*NBAS + ILINOF(I)+IBAS-NLINDP) =
     & EVAL((ISPIN-1)*NBAS + ILINOF(I)+IBAS)
       IF(IFLAGS(11).EQ.2)THEN
        EVAL(NBAS + ILINOF(I)+IBAS-NLINDP) = EVAL(NBAS + ILINOF(I)+IBAS)
       ENDIF
   30  CONTINUE
C
C  Zero out top NLINDP eigenvectors of this block.
C
       DO 40 IBAS=NBFIRR(I)-NLINDP+1,NBFIRR(I)
       CALL ZERO(SCR1(1+(IBAS-1)*NBFIRR(I)),NBFIRR(I))
       EVAL((ISPIN-1)*NBAS + ILINOF(I)+IBAS) = 1.0D+06
       IF(IFLAGS(11).EQ.2)THEN
        EVAL(NBAS + ILINOF(I)+IBAS) = 1.0D+06
       ENDIF
   40  CONTINUE
      ENDIF
 1050 FORMAT(' @SCFIT-I, Symmetry ',I2,' MO ',I3,' Norm is ',F10.8)
C-----------------------------------------------------------------------
C
      IF(GSSLOK)THEN
c     IF(GSSLOK.and.iter.gt.25)THEN
C
C     SCR1 now has SO MOs for this spin and symmetry. Compute
C     COLD(T) * S * C and print. This will approach the identity as
C     convergence is reached and should allow us to detect if the
C     order of vectors has changed.
C
      CALL GETBLK(COLD  ,SCRA,NBFIRR(I),NBAS,IREPS(I))
      CALL GETBLK(SOVRLP,SCRB,NBFIRR(I),NBAS,IREPS(I))
C
C     COLD(T) * S ---> SCRTMP
      CALL XGEMM('T','N',NBFIRR(I),NBFIRR(I),NBFIRR(I),1.0D+00,
     &           SCRA,NBFIRR(I),SCRB,NBFIRR(I),0.0D+00,SCRTMP,NBFIRR(I))
C
C     SCRTMP * SCR1 ---> SCRA
C
      CALL XGEMM('N','N',NBFIRR(I),NBFIRR(I),NBFIRR(I),1.0D+00,
     &           SCRTMP,NBFIRR(I),SCR1,NBFIRR(I),0.0D+00,SCRA,NBFIRR(I))
C
C     Print out COLD(T) * S * C
C
CSSS      WRITE(6,1040) ISPIN,I
CSSS1040  FORMAT(' @SCFIT-I, COLD(T)*S*CNEW for spin ',I2,' symmetry ',I3)
CSSS      CALL OUTPUT(SCRA,1,NBFIRR(I),1,NBFIRR(I),NBFIRR(I),NBFIRR(I),1)
      IF(LOCK(I,ISPIN).GT.0)THEN
      CALL NEWOLD(SCRA,NBFIRR(I),
     &            SCR1,SCRTMP,NOCC((ISPIN-1)*8+I),I)
      ENDIF
C
      ENDIF
C
C         CALL MKDEN(SCR1,SCRB,NOCC(((ISPIN-1)*8)+I),NBFIRR(I),IUHF)
C         CALL SQUEZ2(SCRB,DENSN(((ISPIN-1)*LDIM1)+ITRIOF(I)),
C    &                NBFIRR(I))
C
C
C  Now copy the non-orthogonal basis set eigenvectors into EVEC.
C
          CALL DCOPY(ISQRLN(I),SCR1,1,
     &                         EVEC(((ISPIN-1)*LDIM3)+ISQROF(I)),1)
C This is for ROHF.
          IF(IFLAGS(11).EQ.2) THEN
           CALL DCOPY(ISQRLN(I),SCR1,1,EVEC(LDIM3+ISQROF(I)),1)
          ENDIF        
C
  101   CONTINUE
  100 CONTINUE
C
C  Now do the occupancy check and form the density (IFLAGS(41)=LOCK_ORBOCC is
C  set to zero by default)
C
      IF(IFLAGS(41).EQ.0 .AND. .NOT. GSSLOK) THEN
        CALL DCOPY((IUHF+1)*NBAS,EVAL,1,SCRA,1)
        CALL OCCUPY(NIRREP,NBFIRR,NBAS,SCRA,SCR1,IUHF)
      ENDIF
c  See if the occupancies have changed.
      J = 0
      DO I = 1, 16
         J = J + ABS(NOCC(I)-NPRVOCC(I))
      END DO
      IPRVOCC = IPRVOCC + 1
      IF (J.NE.0) THEN
#ifdef _DEBUG
         if (iter.ne.1) then
            print *, '@SCFIT: occupancy changed'
            if (iFlags(h_IFLAGS_print).ne.0) then
               print *, 'OLD: ',(NPRVOCC(I),I=1,NIRREP),' /',
     &                          (NPRVOCC(J),J=9,8+NIRREP)
               print *, 'NEW: ',(NOCC(I),I=1,NIRREP),' /',
     &                          (NOCC(J),J=9,8+NIRREP)
            end if
         end if
#endif
         DO I = 1, 16
            NPRVOCC(I) = NOCC(I)
         END DO
         IPRVOCC = 0
      END IF
C
C  Jump-to location for QC-SCF
  112 CONTINUE
C
#ifdef _DPOCC
        CALL MKDP_OCCNOS(DOCC, NOCC, NSUM, NBAS, IUHF, .FALSE.)
#endif
      DO 120 ISPIN=1,(IUHF+1)
        IROFF = (ISPIN - 1)*NBAS
        DO 121 I=1,NIRREP
          IF(NBFIRR(I).EQ.0)GOTO 121
#ifdef _DPOCC 
          IOFF = IROFF + 1
          CALL MKDDEN_4IRREP(EVEC((ISPIN-1)*LDIM3+ISQROF(I)),
     &                       SCR1,SCRA,NSUM(((ISPIN-1)*8)+I),
     &                       DOCC(IOFF),NBFIRR(I),ISPIN,IUHF)
           IROFF = IROFF + NBFIRR(I) 
#else 
          CALL MKDEN(EVEC(((ISPIN-1)*LDIM3)+ISQROF(I)),SCR1,
     &               NOCC(((ISPIN-1)*8)+I),NBFIRR(I),IUHF)
#endif
          CALL SQUEZ2(SCR1,DENSN(((ISPIN-1)*LDIM1)+ITRIOF(I)),
     &                NBFIRR(I))
  121   CONTINUE
  120 CONTINUE

C  prakash 
C ks variable
CC
      if (scfks) then
          call DCOPY((iuhf+1)*nbas,eval,1,screval,1)
          call int_occupy(screval,scrocc)
       do 70 ispin=1,iuhf+1
       call mat_trans(MAT_SQR,MAT_FULL,evec((ispin-1)*ldim3+1),
     &                scrsqr,0)
       call putrec(1,'JOBARC','SCFEVL'//spinc(ispin)//'0',
     &              iintfp*nbas,eval((ispin-1)*nbas+1))
       call putrec(1,'JOBARC','SCFEVC'//spinc(ispin)//'0',
     &              iintfp*nbas*nbas,scrsqr)
       call putrec(1,'JOBARC','SCFOCC'//spinc(ispin),
     &              nbas,scrocc(1,ispin))
       call t2fulao(screvec,scrsort,scraks,nbas,nbas,
     &              max(nbas,naobasfn),
     &    'SCFEVC'//spinc(ispin)//'0','EVECAO_'//spinc(ispin),
     &    'SCFEVL'//spinc(ispin)//'0','EVALORDR',1,0)
   70  continue
      end if

#ifdef _DEBUG
      CALL SAXPY((IUHF+1)*ITRILN(NIRREP+1),ONEM,DENSN,1,DENSO,1)
      DTMP=0.D0
      DO 200 I=1,(IUHF+1)*ITRILN(NIRREP+1)
        DTMP=MAX(DTMP,ABS(DENSO(I)))
  200 CONTINUE
      Write(6, "(a, F10.5)") "@-SCFIT Density Criteria: ", DTMP
#endif
      IF (IFLAGS(10).EQ.2.AND..NOT.QC_PHASE) THEN
         QC_PHASE = DMAX.LT.QC_THRESH.AND.NOCONV.AND.ITER.GT.1
     &              .AND.IPRVOCC.GT.3
         IF (QC_PHASE) THEN
            PRINT *, '@SCFIT: entered quadratic convergence phase'
         END IF
      END IF
C
C At the moment I do not know what this suppose to do. Iflags(17)
C corresponds to the occupation. I am not sure about the relation 
C between convergence and occupation. 
C
      IF(IFLAGS(17).EQ.0.AND.ITER.EQ.1)NOCONV=.TRUE.
C
      CALL DCOPY((IUHF+1)*ITRILN(NIRREP+1),DENSN,1,DENSO,1)

#ifdef _PRINTCRAP 
      ntotal = ITRILN(NIRREP+1)
      write(*,*) "The alpha density matrix @SCFIT:"
      call realprt(DENSO,ntotal)
      write(*,*) "The beta density matrix @SCFIT:"
      call realprt(DENSO(ntotal),ntotal)
#endif
C
C  If the calculation is almost converged, lock the occupation.
C
      ZHECTO=100.
      IF(DMAX.LT.TOL*ZHECTO) THEN
        IFLAGS(41)=1
      ENDIF
C
      IF(GSSLOK.OR.IFLAGS(h_IFLAGS_scf_extrapo).EQ.2)THEN
C
C     Write the current MOs to disk.
C
      DO  300 ISPIN=1,IUHF+1
      CALL ZERO(COLD,NBAS*NBAS)
      DO  290 I    =1,NIRREP
      CALL PUTBLK(COLD,EVEC((ISPIN-1)*LDIM3+ISQROF(I)),NBFIRR(I),
     1            NBAS,IREPS(I))
  290 CONTINUE
      IF(ISPIN.EQ.1)THEN
      CALL PUTREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,COLD)
      ELSE
      CALL PUTREC(20,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,COLD)
      ENDIF
  300 CONTINUE
C
      ENDIF
      IF(IFLAGS(11).EQ.2) THEN
       CALL PUTREC(20,'JOBARC','MOVECTOR',ISQRLN(NIRREP+1)*IINTFP,EVEC)
       ROHFMO = .TRUE.
      ENDIF
C
C If the QC-SCF has converged, then turn off QC and do one last iteration.
      IF (QC_PHASE.AND..NOT.NOCONV) THEN
         IFLAGS(h_IFLAGS_scf_extrapo) = 0
         QC_PHASE = .FALSE.
         NOCONV = .TRUE.
         TOL = 0.1
      END IF
      RETURN
      END
