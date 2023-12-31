      SUBROUTINE TRNMOM1(VDIP,VECTOR,ISCR,MXCOR,NROOT,NROOT_ORG,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD
      LOGICAL CIS,RPA,EOMCC,CISD,FULDIAG,INCORE,READGUES
      LOGICAL WARN
      DIMENSION NROOT_ORG(8)
      DIMENSION VDIP(*),VECTOR(*),TRNMOM(3),NROOT(8),ISCR(MXCOR)
      COMMON/METH/CIS,RPA,EOMCC,CISD,FULDIAG,INCORE,READGUES
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/STATSYM/IRREPX
C
      FACT=(2.0D0*DFLOAT(2-IUHF))/3.0D0
      FACTEV = 27.2113957D0
C
      WRITE(6,1003)
      WRITE(6,1000)
      WRITE(6,1001)
      WRITE(6,1003)
      CALL GETREC(20,'JOBARC','SCFENEG ',IINTFP,EGS)
C
      DO 5 IRREPX=1,NIRREP
       IF(CISD)CALL NEWLST(IRREPX,ISCR,MXCOR,IUHF)
       LENA=IRPDPD(IRREPX,9)
       LENB=IRPDPD(IRREPX,10)
       LENGTH=LENA+IUHF*LENB
C
C LOAD UP VIRTUAL-OCCUPIED DIPOLE INTEGRALS
C
       IOFF=1
       DO 10 IXYZ=1,3
        CALL GETLST(VDIP(IOFF),IXYZ,1,1,IRREPX,480)
        IOFF=IOFF+LENA
        IF(IUHF.NE.0)THEN
         CALL GETLST(VDIP(IOFF),IXYZ,1,1,IRREPX,481)
         IOFF=IOFF+LENB
        ENDIF
10     CONTINUE
C
C NOW READ IN EIGENVECTORS
C
       WARN = .FALSE.
       DO 20 IROOT=1,NROOT(IRREPX)
C
C LOOP OVER SPIN CASES
C
        CALL GETLST(VECTOR,IROOT,1,1,IRREPX,94)
        CALL GETLST(ROOT  ,IROOT,1,1,IRREPX,95)
        IF (NROOT(IRREPX) .NE. NROOT_ORG(IRREPX)) THEN
            WARN = .TRUE.
        ENDIF 
C
C NOW LOOP OVER CARTESIAN DIRECTIONS
C
        IOFF=1
        DO 21 IXYZ=1,3
         TRNMOM(IXYZ)=SDOT(LENGTH,VECTOR,1,VDIP(IOFF),1)
         IOFF=IOFF+LENGTH
21      CONTINUE
        X=FACT*ROOT*SNRM2(3,TRNMOM,1)**2
C
        WRITE(6,1002)IROOT,IRREPX,ROOT*FACTEV,(TRNMOM(I),I=1,3),X
        WRITE(6,1004)EGS+ROOT
        CALL TDAANAL(VECTOR,ISCR,IUHF,IRREPX)
20     CONTINUE
5     CONTINUE
      CALL PUTREC(20,'JOBARC','TOTENER2',IINTFP,EGS+ROOT)
      WRITE(6,1003)

      IF (WARN) THEN
          WRITE(6,*)
          WRITE(6,"(2a)") " The number of roots printed per irrep. is",
     &                    " is less than requested. This means"
          WRITE(6,"(a)")  " fewer than requested were located."
      ENDIF 
C 
1000  FORMAT(T5,'Tamm-Dancoff (CIS) Excitation Energies and ',
     &        'Transition Moments',/)
1001  FORMAT(T17,'Excitation',T36,'Transition Dipole',T62,'Oscillator',
     &       /,T2,'Root',T8,'Symmetry',T18,'Energy(eV)',T33,'x',T45,
     &       'y',T57,'z',T63,'Strength') 
1002  FORMAT(T2,I3,T12,I1,T14,F10.5,T25,F11.4,T37,F11.4,T49,F11.4,
     &       T62,D13.8)
1003  FORMAT(74('_'))
1004  FORMAT(T3,'Total TDA electronic energy ',F15.8,' a.u.')
      RETURN
      END
