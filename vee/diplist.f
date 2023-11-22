      SUBROUTINE DIPLIST(F,BUF,NBAS,NBAST,IUHF,IRREPX)
C
C THIS ROUTINE CREATES THE VARIOUS DIPOLE INTEGRAL LISTS AND
C WRITES THEM OUT.  HACKED UP FROM FOCKLIST.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DIMENSION IOFFO(8),IOFFV(8)
      DOUBLE PRECISION F(NBAS,NBAS),BUF(4*NBAST*NBAST)
      LOGICAL VPROP,PRINT
      CHARACTER*8 LABEL(3)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),IOFFVS(8,2),IOFFOS(8,2),
     &             IRPDPDAO(8),IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFEA(2),NFMI(2)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /VDINTLEN/ LENS(8),LENT(8)
      COMMON /VDINTPRT/NTPERT,NPERT(8),KPERT(8),IDIPSM(3),
     &                 IYZPERT,IXZPERT,IXYPERT,ITRANSX,
     &                 ITRANSY,ITRANSZ,NUCIND
      COMMON /INTPROG/ VPROP
      DATA LABEL /'DIPOLE_X','DIPOLE_Y','DIPOLE_Z'/
C
      NNP1O2(I)=(I*(I+1))/2
C
      PRINT = (IFLAGS(1) .GT. 30)

      NBAST2=NBAST*NBAST
      IOFFBF=NBAST*NBAST+1
C
C LOOP OVER IRREPX
C
      DO 1000 IRREPX=1,NIRREP
C
C LOOP OVER SPIN CASE
C
       DO 1001 ISPIN=1,1+IUHF
C
C CALCULATE OFFSETS
C
        IOFFO(1)=0
        IOFFV(1)=NOCCO(ISPIN)
        DO 1 IRREP=1,NIRREP-1
         IOFFO(IRREP+1)=IOFFO(IRREP)+POP(IRREP,ISPIN)
         IOFFV(IRREP+1)=IOFFV(IRREP)+VRT(IRREP,ISPIN)
1       CONTINUE
C
C LOOP OVER CARTESIAN DIRECTIONS
C
        DO 1002 IXYZ=1,3
C
C PICK UP DIPOLE MOMENT INTEGRALS FROM JOBARC AND TRANSFORM TO THE MO BASIS
C
         IF(VPROP)THEN
          LENGTH=NNP1O2(NBAST)
          CALL GETREC(-1,'JOBARC',LABEL(IXYZ),LENGTH*IINTFP,F)
          CALL EXPND2(F,BUF,NBAST)
          CALL SCOPY (NBAST2,BUF,1,F,1)
         ELSE
          IRRPRT=IDIPSM(IXYZ)
          LENGTH=LENT(IRRPRT)
          CALL GETREC(-1,'JOBARC',LABEL(IXYZ),LENGTH*IINTFP,F)
          CALL VMINUS(F,LENGTH)
          CALL ZERO(BUF,NBAST*NBAST)
          CALL MATEXP(IRRPRT,IAOPOP,F,BUF)
          CALL MATEXP2(IRRPRT,BUF,F,NBAST)
         ENDIF
         CALL AO2MO2(F,F,BUF,BUF(IOFFBF),NBAS,NBAST,ISPIN)
         CALL FILTER(F,NBAS*NBAS,1.D-10)
         IF(PRINT)THEN
          WRITE(6,5000)LABEL(IXYZ)(8:8)
5000      FORMAT(T3,A,' MO basis dipole integrals ')
          WRITE(6,'((2I5,F10.5,20X,2I5,F10.5))')((I,J,F(I,J),
     &          I=1,NBAS),J=1,NBAS)
         ENDIF  
C
C FORM SYMMETRY PACKED OCCUPIED-OCCUPIED PART.
C
         ITHRU=0
         DO 10 IRREPR=1,NIRREP
          IRREPL=DIRPRD(IRREPR,IRREPX)
          DO 20 J=1,POP(IRREPR,ISPIN)
           DO 21 I=1,POP(IRREPL,ISPIN)
            ITHRU=ITHRU+1
            BUF(ITHRU)=F(I+IOFFO(IRREPL),J+IOFFO(IRREPR))
21         CONTINUE
20        CONTINUE
10       CONTINUE
         CALL PUTLST(BUF,IXYZ,1,1,IRREPX,475+ISPIN)
C
C FORM SYMMETRY PACKED VIRTUAL-VIRTUAL PART.
C
         ITHRU=0
         DO 110 IRREPR=1,NIRREP
          IRREPL=DIRPRD(IRREPR,IRREPX)
          DO 120 J=1,VRT(IRREPR,ISPIN)
           DO 121 I=1,VRT(IRREPL,ISPIN)
            ITHRU=ITHRU+1
            BUF(ITHRU)=F(I+IOFFV(IRREPL),J+IOFFV(IRREPR))
c            BUF(ITHRU)=0.0D0
121        CONTINUE
120       CONTINUE
110      CONTINUE
         CALL PUTLST(BUF,IXYZ,1,1,IRREPX,477+ISPIN)
C
C FORM SYMMETRY PACKED VIRTUAL-OCCUPIED PART.
C
         ITHRU=0
         DO 210 IRREPR=1,NIRREP
          IRREPL=DIRPRD(IRREPR,IRREPX)
          DO 220 I=1,POP(IRREPR,ISPIN)
           DO 221 A=1,VRT(IRREPL,ISPIN)
            ITHRU=ITHRU+1
            BUF(ITHRU)=F(A+IOFFV(IRREPL),I+IOFFO(IRREPR))
c            BUF(ITHRU)=0.0D0
221        CONTINUE
220       CONTINUE
210      CONTINUE
         CALL PUTLST(BUF,IXYZ,1,1,IRREPX,479+ISPIN)
C
1002    CONTINUE
1001   CONTINUE
1000  CONTINUE
C
      RETURN
      END 
