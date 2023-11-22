      SUBROUTINE OOINOO(IOO,DOO,ICORE,MXCOR,IUHF)
C
C THIS ROUTINE COMPUTES THE CONTRIBUTIONS:
C
C      IOO(I,J) = - SUM D(M,N) * <MI||NJ> + D(m,n) * <Im|Jn> (ISPIN=1)
C
C      IOO(i,j) = - SUM D(m,n) * <mi||nj> + D(M,N) * <iM|jN> (ISPIN=2)
C
      IMPLICIT INTEGER (A-Z)
      CHARACTER*4 TYPE
      DOUBLE PRECISION IOO,DOO
      DIMENSION DOO(1),IOO(1),ICORE(MXCOR)
      COMMON /SYM2/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP2/ IRPDPD(8,22)
      COMMON /SYMPOP/ IRP_DM(8,22),ISYTYP(2,500),ID(18)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      common /dropgeo/ ndrgeo
      COMMON /SHIFT/ ISHIFT 
C
      IF(IUHF.EQ.0)THEN
C
C SPIN ADAPTED RHF CODE
C
       LISTW=13 + ISHIFT 
       DO 120 IRREP=1,NIRREP
        DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW)) 
        NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW)) 
        I000=1
        I010=I000+IINTFP*NUMDIS*DISSIZ
        I020=I010+IINTFP*DISSIZ
        I030=I020+IINTFP*DISSIZ
        CALL GETLST(ICORE(I000),1,NUMDIS,1,IRREP,LISTW)
        CALL SPINAD1(IRREP,POP(1,1),DISSIZ,ICORE(I000),
     &               ICORE(I010),ICORE(I020))
        CALL VMINUS(ICORE(I000),NUMDIS*DISSIZ)
        CALL DOT24(IRREP,IOO,DOO,ICORE(I000),ICORE(I010),DISSIZ,
     &             POP(1,1),POP(1,1),POP(1,1),POP(1,1),POP(1,1),
     &             POP(1,1),'STST')
120    CONTINUE
      ELSEIF(IUHF.NE.0)THEN
       DO 10 ISPIN=1,2
        LISTW1=10+ISPIN + ISHIFT 
        LISTW2=13 + ISHIFT 
        DO 20 IRREP=1,NIRREP
         NDSZ1=IRPDPD(IRREP,ISYTYP(1,LISTW1))
         NDSZF=IRPDPD(IRREP,20+ISPIN)
         NDIS1=IRPDPD(IRREP,ISYTYP(2,LISTW1))
         NDISF=IRPDPD(IRREP,20+ISPIN)
         NDSZ2=IRPDPD(IRREP,ISYTYP(1,LISTW2))
         NDIS2=IRPDPD(IRREP,ISYTYP(2,LISTW2))
         IOFFI=1+(ISPIN-1)*NFMI(1)
C
C DO FIRST PART - D(m,n) * I(mi,nj)
C
         I000=1
         I010=I000+IINTFP*NDSZF*NDISF
         I020=I010+NFMI(ISPIN)
         IOFFD=1+(ISPIN-1)*NFMI(1)
         CALL GETLST(ICORE(I000),1,NDIS1,2,IRREP,LISTW1)
         CALL SYMEXP(IRREP,POP(1,ISPIN),NDSZ1,ICORE(I000))
         CALL SYMEXP2(IRREP,POP(1,ISPIN),NDSZF,NDSZ1,NDISF,ICORE(I000),
     &                ICORE(I000))
         CALL VMINUS(ICORE(I000),NDSZF*NDISF)
         CALL DOT24(IRREP,IOO(IOFFI),DOO(IOFFD),ICORE(I000),
     &              ICORE(I010),NDSZF,POP(1,ISPIN),POP(1,ISPIN),
     &              POP(1,ISPIN),POP(1,ISPIN),POP(1,ISPIN),POP(1,ISPIN),
     &             'STST')
        
C NOW DO SECOND PART 
C
C         D(m,n) * <Im|Jn> [ISPIN=1]
C 
C         D(M,N) * <Mi|Nj> [ISPIN=2]
C
         I000=1
         I010=I000+IINTFP*NDSZ2*NDIS2
         I020=I010+NFMI(ISPIN)
         IOFFD=1+(2-ISPIN)*NFMI(1)
         CALL GETLST(ICORE(I000),1,NDIS2,2,IRREP,LISTW2)
         CALL VMINUS(ICORE(I000),NDSZ2*NDIS2)
         IF(ISPIN.EQ.1)THEN
          TYPE='TSTS'
         ELSE
          TYPE='STST'
         ENDIF
         CALL DOT24(IRREP,IOO(IOFFI),DOO(IOFFD),ICORE(I000),
     &              ICORE(I010),NDSZ2,POP(1,ISPIN),POP(1,ISPIN),
     &              POP(1,1),POP(1,2),POP(1,1),POP(1,2),TYPE)
20      CONTINUE
10     CONTINUE
      ENDIF
      RETURN
      END
