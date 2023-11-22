      SUBROUTINE DT1IA(ICORE,MAXCOR,ANTI,GRAD1)
C
C THIS ROUTINE COMPUTES THE DERIVATIVES OF THE FIRST-ORDER T1
C AMPLITUDES. IT IS ONLY REQUIRED FOR ROHF OR GENERAL NON-HF
C CASES.
C
C  d T(I,A)    a -1   d fai          d fmi
C  -------- = D     { ------- - SUM  ----- T(M,A)
C   d chi      i       d chi     m   d chi
C
C                                  d fae
C                     + SUM T(I,E) ----- 
C                        e         d chi
C
C                               d T(M,A)
C                     - SUM fim --------
C                        m       d chi
C
C                               d T(I,E)
C                     + SUM fae --------
C                                d chi
C
C SEE J. GAUSS, UNPUBLISHED NOTES
C
CEND
C
C CODED MAY/91 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,POP,VRT
      LOGICAL QRHF,ROHF,SEMI,FIELD,GEOM,ANTI
      CHARACTER*2 SPCASE(2)
C
      DIMENSION ICORE(MAXCOR)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/DSYM/IRREPX,IPERT,NDT(2),NDF1(2),NDF2(2),
     &            IOFFIJ(8,2),IOFFAB(8,2),IOFFAI(8,2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
      COMMON/DTRAN/FIELD,GEOM,ROHF,QRHF,SEMI
C
      DATA ONE,ONEM /1.D0,-1.D0/
      DATA SPCASE /'AA','BB'/
C
C LOOP OVER SPIN CASES
C
      DO 1 ISPIN=1,2
C
      WRITE(6,3000) SPCASE(ISPIN)
3000  FORMAT('  @DT1IA-I, Beginning iterative solution for ', 
     &        A2,' perturbed T1 amplitudes.')
C
C ALLOCATE MEMORY FOR THE DERIVATIVES OF THE FOCK 
C MATRICES
C
      IFOO=1
      IFVV=IFOO+IINTFP*IRPDPD(IRREPX,20+ISPIN)
      IFVO=IFVV+IINTFP*IRPDPD(IRREPX,18+ISPIN)
      IDT=IFVO+IINTFP*IRPDPD(IRREPX,8+ISPIN)
      IT=IDT+IINTFP*IRPDPD(IRREPX,8+ISPIN)
      IDDEN=IT+IINTFP*NT(ISPIN)
      IEND=IDDEN+IINTFP*IRPDPD(IRREPX,8+ISPIN)
      IF(.NOT.SEMI) THEN
       IFIJ=IEND
       IFAB=IFIJ+IINTFP*NF1(ISPIN)
       IDTAR=IFAB+IINTFP*NF2(ISPIN)
       IEND=IDTAR+IINTFP*NDT(ISPIN)
      ENDIF
       
C
C READ THE DERIVATIVES OF THE FOCK MATRICES FROM DISK
C
      CALL GETLST(ICORE(IFOO),IPERT,1,1,IRREPX,175+ISPIN)
      CALL GETLST(ICORE(IFVV),IPERT,1,1,IRREPX,177+ISPIN)
      CALL GETLST(ICORE(IFVO),IPERT,1,1,IRREPX,179+ISPIN)
C
C READ UNPERTURBED T1 AMPLITUDES
C
      CALL GETLST(ICORE(IT),1,1,1,ISPIN,90)
C
C READ ORBITAL ENERGY DENOMINATOR
C
      CALL GETLST(ICORE(IDDEN),1,1,1,9,447+ISPIN)
C
C IF STANDARD ORBITALS ARE USED, READ ALSO OCC.-OCC. 
C AND VIRT.-VIRT. BLOCK OF FOCK MATRICES
C
      IF(.NOT.SEMI) THEN
       CALL GETLST(ICORE(IFIJ),1,1,1,2+ISPIN,91)
       CALL GETLST(ICORE(IFAB),1,1,1,2+ISPIN,92)
      ENDIF
c
c calculate gradient contribution d f(,i)/d x  t(a,i)
c
      if(irrepx.eq.1) then 
       grad=sdot(nt(ispin),icore(it),1,icore(ifvo),1)
       GRAD1=GRAD1+GRAD
C
       write(6,*) 'gradient contribution df*t', grad
C
      ENDIF
C
C LOOP OVER ALL IRREPS
C
      IOFFD=0
      IOFFFD1=0
      IOFFFD2=0
      IOFFT2=0
C
      DO 1000 IRREPR=1,NIRREP
C
       IRREPL=DIRPRD(IRREPR,IRREPX)
C
C  -  T(M,A) FD(IM)
C
       NOCCR=POP(IRREPR,ISPIN)
       NOCCL=POP(IRREPL,ISPIN)
       NVRTR=VRT(IRREPR,ISPIN)
       NVRTL=VRT(IRREPL,ISPIN)
C
       IOFFT1=0
       DO 89 IRREP=1,IRREPL-1
        IOFFT1=IOFFT1+IINTFP*POP(IRREP,ISPIN)*VRT(IRREP,ISPIN)
89     CONTINUE
C
       CALL XGEMM('N','N',NVRTL,NOCCR,NOCCL,ONEM,ICORE(IT
     &            +IOFFT1),NVRTL,ICORE(IFOO+IOFFFD1),NOCCL,ONE,
     &            ICORE(IFVO+IOFFD),NVRTL)
C
       CALL XGEMM('N','N',NVRTL,NOCCR,NVRTR,ONE,ICORE(IFVV
     &            +IOFFFD2),NVRTL,ICORE(IT+IOFFT2),NVRTR,
     &            ONE,ICORE(IFVO+IOFFD),NVRTL)
C
C UPDATE OFFSETS
C
       IOFFFD1=IOFFFD1+IINTFP*NOCCR*NOCCL
       IOFFFD2=IOFFFD2+IINTFP*NVRTR*NVRTL
       IOFFD=IOFFD+IINTFP*NVRTL*NOCCR
       IOFFT2=IOFFT2+IINTFP*NVRTR*NOCCR
C
1000  CONTINUE
C
c YAU : old
c     CALL ICOPY(IINTFP*NDT(ISPIN),ICORE(IFVO),1,ICORE(IDT),1)
c YAU : new
      CALL DCOPY(NDT(ISPIN),ICORE(IFVO),1,ICORE(IDT),1)
c YAU : end
C
      CALL VECPRD(ICORE(IDT),ICORE(IDDEN),ICORE(IDT),
     &            NDT(ISPIN))
C
      IF(.NOT.SEMI) THEN
C
C FOR STANDARD ORBITALS, ITERATE
C
       IFIJ1=IFIJ
       IDT1=IDT
       IDT2=IDTAR
       IFAI1=IFVO
       IDD1=IDDEN
C
       DO 2000 IRREPR=1,NIRREP
C
        IRREPL=DIRPRD(IRREPX,IRREPR)
C
        ITER=0
C 
C CALCULATE D T(EI)/D chi * F(EA)
C
        NVRTL=VRT(IRREPL,ISPIN)
        NOCCR=POP(IRREPR,ISPIN)
        LENGTH=NVRTL*NOCCR
C
        IFAB1=IFAB
        DO 189 IRREP=1,IRREPL-1
         IFAB1=IFAB1+IINTFP*VRT(IRREP,ISPIN)*VRT(IRREP,ISPIN)
189      CONTINUE
C
        IF(LENGTH.EQ.0) GO TO 2002
C
2001    CALL IZERO(ICORE(IDT2),LENGTH)
C
        CALL XGEMM('N','N',NVRTL,NOCCR,NVRTL,ONE,ICORE(IFAB1),
     &             NVRTL,ICORE(IDT1),NVRTL,ONE,ICORE(IDT2),
     &             NVRTL)
C
C CALCULATE D T(AM)/D chi * F(M,I)
C
        CALL XGEMM('N','N',NVRTL,NOCCR,NOCCR,ONEM,ICORE(IDT1),
     &             NVRTL,ICORE(IFIJ1),NOCCR,ONE,ICORE(IDT2),
     &             NVRTL)
C
C ADD IN D F(AI)/D chi ... PIECES
C      
        CALL SAXPY(LENGTH,ONE,ICORE(IFAI1),1,ICORE(IDT2),1)
        CALL VECPRD(ICORE(IDT2),ICORE(IDD1),ICORE(IDT2),LENGTH)
C
C COMPARE WITH PREVIOUS ITERATE
C
        CALL SAXPY(LENGTH,ONEM,ICORE(IDT2),1,ICORE(IDT1),1)
        X=FNDLRGAB(ICORE(IDT1),LENGTH)
        IF(X.GT.1.E-10) THEN
         ITER=ITER+1
c YAU : old
c        CALL ICOPY(LENGTH*IINTFP,ICORE(IDT2),1,ICORE(IDT1),1)
c YAU : new
         CALL DCOPY(LENGTH,ICORE(IDT2),1,ICORE(IDT1),1)
c YAU : end
         GO TO 2001
        ELSE
         ITER=ITER+1
         WRITE(6,1001) IRREP,ITER
1001     FORMAT(T3,'Irrep ',I2,' perturbed amplitudes converged',
     &          ' in ',I4,' iterations.')
c YAU : old
c        CALL ICOPY(LENGTH*IINTFP,ICORE(IDT2),1,ICORE(IDT1),1)
c YAU : new
         CALL DCOPY(LENGTH,ICORE(IDT2),1,ICORE(IDT1),1)
c YAU : end
        ENDIF
2002    CONTINUE
        IFIJ1=IFIJ1+NOCCR*NOCCR*IINTFP
        IDD1=IDD1+NVRTL*NOCCR*IINTFP
        IDT1=IDT1+NVRTL*NOCCR*IINTFP
        IDT2=IDT2+NVRTL*NOCCR*IINTFP
        IFAI1=IFAI1+NVRTL*NOCCR*IINTFP
2000   CONTINUE
C
      ENDIF
C
C SAVE PERTURBED AMPLITUDES ON DERGAM
C
      CALL PUTLST(ICORE(IDT),1,1,1,ISPIN,490)
C
C CALCULATE GRADIENT CONTRIBUTION
C
      IF(IRREPX.EQ.1) THEN
       CALL GETLST(ICORE(IDTAR),1,1,1,ISPIN+2,93)
       GRAD=SDOT(NT(ISPIN),ICORE(IDT),1,ICORE(IDTAR),1)
       GRAD1=GRAD1+GRAD
C
       write(6,*) 'gradient contribution to f*dt', grad
C
      ENDIF
        
1     CONTINUE
C
C ALL DONE
C
      RETURN
      END