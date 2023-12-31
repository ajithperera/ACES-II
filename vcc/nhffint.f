      SUBROUTINE NHFFINT(ICORE,MAXCOR,IUHF,SINGLE)
C
C THIS SUBROUTINE COMPUTES THE PART OF THE F INTERMEDIATES WHICH
C  IS DUE TO OFF-DIAGONAL TERMS IN THE FOCK MATRIX AND THEN AUGMENTS
C  THE F INTERMEDIATES WITH THESE VALUES.
C
C
C       f(E,A)  - (1/2) SUM F(E,M)*T1(A,M) (FOR CASE 'FAE')
C                        M
C
C       f(M,I)  + (1/2) SUM F(E,I)*T1(E,M) (FOR CASE 'FMI')
C                        M
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION HALF,MHALF,ZILCH,ONE
      LOGICAL SINGLE,ROHF4,ITRFLG,ROHFMB
      DIMENSION ICORE(MAXCOR),IOFFT(8,2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /ROHF/ ROHF4,ITRFLG
      COMMON /FLAGS/ IFLAGS(100)
      DATA HALF /0.5/
      DATA MHALF /-0.5/
      DATA ZILCH /0.0/
      DATA ONE /1.0/
C
      ROHFMB=.FALSE.
      IF((IFLAGS(38).NE.0 .OR. IFLAGS(11).EQ.2) .AND.
     &                         IFLAGS( 2).LE.4) ROHFMB=.TRUE.
C
      IF(SINGLE)CALL GETT1(ICORE,MAXCOR,MXCOR,IUHF,IOFFT)
C
C LOOP OVER SPIN CASES
C
      DO 5 ISPIN=1,1+IUHF
C
C DETERMINE SIZE OF THE SYMMETRY-PACKED FOCK MATRIX [f(m,e)].
C
       SIZFOK=NT(ISPIN)
C
C DO F(EA) INTERMEDIATES.
C
       SIZTAR=NFEA(ISPIN)
       I000=1
       I010=I000+SIZFOK*IINTFP
       I020=I010+SIZTAR*IINTFP
       IOFFFME=I000
       IOFFTAR=I010
       IOFFT1 =IOFFT(1,ISPIN)
       IF(ROHFMB.AND..NOT.ROHF4)THEN
        CALL ZERO(ICORE(I000),SIZFOK)
       ELSE
        CALL GETLST(ICORE(I000),1,1,1,2+ISPIN,93)
       ENDIF
       CALL ZERO(ICORE(I010),SIZTAR)
       IF(SINGLE)THEN
        DO 10 IRREP=1,NIRREP
         NVRT=VRT(IRREP,ISPIN)
         NOCC=POP(IRREP,ISPIN)
         CALL XGEMM('N','T',NVRT,NVRT,NOCC,MHALF,ICORE(IOFFFME),
     &              NVRT,ICORE(IOFFT1),NVRT,ZILCH,ICORE(IOFFTAR),
     &              NVRT)
         IOFFT1=IOFFT1+NOCC*NVRT*IINTFP
         IOFFFME=IOFFFME+NOCC*NVRT*IINTFP
         IOFFTAR=IOFFTAR+NVRT*NVRT*IINTFP
10      CONTINUE
       ENDIF
       I030=I020+SIZTAR*IINTFP
       IF(ROHFMB.AND.ROHF4) THEN
         CALL GETLST(ICORE(I020),1,1,1,ISPIN,92)
         CALL SAXPY(SIZTAR,ONE,ICORE(I010),1,ICORE(I020),1)
       ELSE
         CALL GETLST(ICORE(I020),1,1,1,2+ISPIN,92)
         CALL SAXPY(SIZTAR,ONE,ICORE(I010),1,ICORE(I020),1)
         CALL GETLST(ICORE(I010),1,1,1,ISPIN,92)
         CALL SAXPY(SIZTAR,ONE,ICORE(I010),1,ICORE(I020),1)
       ENDIF
       CALL PUTLST(ICORE(I020),1,1,1,ISPIN,92)
C
C DO F(MI) INTERMEDIATES
C
       SIZTAR=NFMI(ISPIN)
       I000=1
       I010=I000+SIZFOK*IINTFP
       I020=I010+SIZTAR*IINTFP
       IOFFFME=I000
       IOFFTAR=I010
       IOFFT1=IOFFT(1,ISPIN)
       IF(ROHFMB.AND..NOT.ROHF4)THEN
        CALL ZERO(ICORE(I000),SIZFOK)
       ELSE
        CALL GETLST(ICORE(I000),1,1,1,2+ISPIN,93)
       ENDIF
       CALL ZERO(ICORE(I010),SIZTAR)
       IF(SINGLE)THEN
        DO 20 IRREP=1,NIRREP
         NVRT=VRT(IRREP,ISPIN)
         NOCC=POP(IRREP,ISPIN)
         CALL XGEMM('T','N',NOCC,NOCC,NVRT,HALF,ICORE(IOFFFME),
     &              NVRT,ICORE(IOFFT1),NVRT,ZILCH,ICORE(IOFFTAR),
     &              NOCC)
         IOFFT1=IOFFT1+NOCC*NVRT*IINTFP
         IOFFFME=IOFFFME+NOCC*NVRT*IINTFP
         IOFFTAR=IOFFTAR+NOCC*NOCC*IINTFP
20      CONTINUE
       ENDIF
       I030=I020+SIZTAR*IINTFP
       IF(ROHFMB.AND.ROHF4) THEN
         CALL SCOPY(SIZTAR,ICORE(I010),1,ICORE(I020),1)
       ELSE
         CALL GETLST(ICORE(I020),1,1,1,2+ISPIN,91)
         CALL SAXPY(SIZTAR,ONE,ICORE(I010),1,ICORE(I020),1)
         CALL GETLST(ICORE(I010),1,1,1,ISPIN,91)
         CALL SAXPY(SIZTAR,ONE,ICORE(I010),1,ICORE(I020),1)
       ENDIF
       CALL PUTLST(ICORE(I020),1,1,1,ISPIN,91)
       IF(.NOT.ROHF4)THEN
C
C DO F(ME) INTERMEDIATES
C
        SIZTAR=SIZFOK
        I000=1
        I010=I000+SIZFOK*IINTFP
        I020=I010+SIZTAR*IINTFP
        CALL GETLST(ICORE(I000),1,1,1,2+ISPIN,93)
        IF(.NOT.SINGLE)THEN
c         CALL UPDMOI(1,SIZTAR,ISPIN,93,0,0)
        ELSE
         CALL GETLST(ICORE(I010),1,1,1,ISPIN,93)
         CALL SAXPY(SIZTAR,ONE,ICORE(I010),1,ICORE(I000),1)
        ENDIF
        CALL PUTLST(ICORE(I000),1,1,1,ISPIN,93)
       ENDIF
5     CONTINUE
C
C GO BACK HOME
C
      RETURN
      END
