         
      SUBROUTINE UAIINVVC(IRREPX,IOFFX,NX,FABA,UAIA,ICORE,
     &                    MXCOR,IUHF)
C
C THIS ROUTINE COMPUTES THE CONTRIBUTIONS:
C
C      FD(A,B) =  SUM U(E,M) (<MA||EB> -  <EA||MB>) (ISPIN=1)
C
C      FD(a,b) =  SUM U(e,m) (<ma||eb> - <ea||mb>) (ISPIN=2)
C
C NOTE THAT THIS ROUTINE DEALS WITH COMPLEX PERTURBATIONS AND
C THAT IN THIS CASE THE COMPLEX CONJUGATE IS THE NEGATIVE OF
C THE ORIGINAL (U(E,M)* = - U(E,M). THERE IS NO SPIN COUPLING.
C
CEND
C

      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,UAIA,FABA,ONEM,AZERO
      DIMENSION FABA(1000),UAIA(1000),ICORE(MXCOR)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA ONE/1.0D0/,AZERO,ONEM/0.D0,-1.D0/
C
      IF(IUHF.EQ.0)THEN
C
C SPIN ADAPTED RHF CODE
C
C FOR RHF THIS REDUCES TO 
C
C - P_(A,B) SUM U(E,M) <Ma||Be>
C
C   - P_(A,B) SUM U(E,M) <Em||Ab>
C
       LISTW=30
       DISSIZ=IRPDPD(IRREPX,ISYTYP(1,LISTW)) 
       NUMDIS=IRPDPD(IRREPX,ISYTYP(2,LISTW)) 
       I000=1
       I010=I000+IINTFP*NUMDIS*DISSIZ
       I020=I010+IINTFP*DISSIZ
       IF(I020.LT.MXCOR) THEN 
        CALL GETLST(ICORE(I000),1,NUMDIS,1,IRREPX,LISTW)
        DO 121 IPERT=IOFFX+1,NX
         CALL GETLST(UAIA,IPERT,1,1,IRREPX,182)
         CALL ZERO(FABA,IRPDPD(IRREPX,19))
         CALL XGEMM('N','N',DISSIZ,1,NUMDIS,ONEM,ICORE(I000),
     &              DISSIZ,UAIA,NUMDIS,AZERO,FABA,DISSIZ)
         CALL ASYMMET(IRREPX,VRT(1,1),FABA,ICORE(I010),
     &                IRPDPD(IRREPX,19))
         CALL GETLST(ICORE(I010),IPERT,1,1,IRREPX,178)
         CALL SAXPY(IRPDPD(IRREPX,19),ONE,FABA,1,ICORE(I010),1)
         CALL PUTLST(ICORE(I010),IPERT,1,1,IRREPX,178)
121     CONTINUE  
       ELSE
        IOFFSET=1
        IOFFU=1
        DO 122 IRREP=1,NIRREP
         NOCC=POP(IRREP,1)
         NVRT=VRT(DIRPRD(IRREP,IRREPX),1)
         I010=I000+IINTFP*DISSIZ*NVRT
         I020=I010+IINTFP*DISSIZ
         IF(I020.GE.MXCOR) CALL INSMEM('UAIINVVC',I010,MXCOR)
         DO 122 IOCC=1,NOCC
          CALL GETLST(ICORE(I000),IOFFSET,NVRT,1,IRREPX,LISTW)
          IOFFSET=IOFFSET+NVRT
          DO 123 IPERT=IOFFX+1,NX
           CALL GETLST(UAIA,IPERT,1,1,IRREPX,182)
           CALL ZERO(FABA,IRPDPD(IRREPX,19))
           CALl XGEMM('N','N',DISSIZ,1,NVRT,ONEM,ICORE(I000),
     &                DISSIZ,UAIA(IOFFU),NVRT,AZERO,FABA,DISSIZ)
           CALL ASYMMET(IRREPX,VRT(1,1),FABA,ICORE(I010),
     &                  IRPDPD(IRREPX,19))
           CALL GETLST(ICORE(I010),IPERT,1,1,IRREPX,178)
           CALL SAXPY(IRPDPD(IRREPX,19),ONE,FABA,1,ICORE(I010),1) 
           CALL PUTLST(ICORE(I010),IPERT,1,1,IRREPX,178)
123       CONTINUE
          IOFFU=IOFFU+NVRT
122      CONTINUE
        ENDIF
C
      ELSEIF(IUHF.NE.0)THEN
       call errex
       DO 10 ISPIN=1,2
        LISTW1=26+ISPIN
        LISTW2=31-ISPIN
c        DO 20 IRREP=1,NIRREP
         NDSZ1=IRPDPD(IRREP,ISYTYP(1,LISTW1))
         NDSZF=IRPDPD(IRREP,18+ISPIN)
         NDIS1=IRPDPD(IRREP,ISYTYP(2,LISTW1))
         NDSZ2=IRPDPD(IRREP,ISYTYP(1,LISTW2))
         NDIS2=IRPDPD(IRREP,ISYTYP(2,LISTW2))
C
C DO FIRST PART - SD(m,n) * I(mi,nj)
C
         I000=1
         I010=I000+IINTFP*NDSZF*NDIS1
         CALL GETLST(ICORE(I000),1,NDIS1,2,IRREP,LISTW1)
         CALL SYMEXP2(IRREP,VRT(1,ISPIN),NDSZF,NDSZ1,NDIS1,
     &                ICORE(I000),ICORE(I000))
         DO 21 IPERT=IOFFX+1,NX
         CALL GETLST(UAIA,IPERT,1,1,IRREPX,181+ISPIN)
         CALL ZERO(FABA,IRPDPD(IRREPX,18+ISPIN))
         CALL DDOT24(IRREPX,IRREP,FABA,UAIA,ICORE(I000),
     &              ICORE(I010),NDSZF,VRT(1,ISPIN),VRT(1,ISPIN),
     &              VRT(1,ISPIN),VRT(1,ISPIN),VRT(1,ISPIN),
     &              POP(1,ISPIN),'TSTS')
         CALL SYMMET6(IRREPX,VRT(1,ISPIN),FABA,ICORE(I010),
     &                IRPDPD(IRREPX,18+ISPIN))
         CALL GETLST(ICORE(I010),IPERT,1,1,IRREPX,177+ISPIN)
         CALL SAXPY(IRPDPD(IRREPX,18+ISPIN),ONE,FABA,1,ICORE(I010),1)
         CALL PUTLST(ICORE(I010),IPERT,1,1,IRREPX,177+ISPIN)
21       CONTINUE
C
C NOW DO SECOND PART 
C
C        SD(m,n) * <Im|Jn> [ISPIN=1]
C 
C        SD(M,N) * <Mi|Nj> [ISPIN=2]
C
         I000=1
         I010=I000+IINTFP*NDSZ2*NDIS2
         CALL GETLST(ICORE(I000),1,NDIS2,2,IRREP,LISTW2)
         DO 22 IPERT=IOFFX+1,NX
         CALL GETLST(UAIA,IPERT,1,1,IRREPX,184-ISPIN)
         CALL ZERO(FABA,IRPDPD(IRREPX,18+ISPIN))
         IF(ISPIN.EQ.1) THEN
         CALL DDOT24(IRREPX,IRREP,FABA,UAIA,ICORE(I000),
     &              ICORE(I010),NDSZ2,VRT(1,ISPIN),VRT(1,ISPIN),
     &              VRT(1,ISPIN),VRT(1,3-ISPIN),VRT(1,ISPIN),
     &              POP(1,3-ISPIN),'TSTS')
         ELSE
         CALL DDOT24(IRREPX,IRREP,FABA,UAIA,ICORE(I000),
     &              ICORE(I010),NDSZ2,VRT(1,ISPIN),VRT(1,ISPIN),
     &              VRT(1,3-ISPIN),VRT(1,ISPIN),POP(1,3-ISPIN),
     &              VRT(1,ISPIN),'STST')
         ENDIF
         CALL SYMMET6(IRREPX,VRT(1,ISPIN),FABA,ICORE(I010),
     &                IRPDPD(IRREPX,18+ISPIN))
         CALL GETLST(ICORE(I010),IPERT,1,1,IRREPX,177+ISPIN)
         CALL SAXPY(IRPDPD(IRREPX,18+ISPIN),ONE,FABA,1,ICORE(I010),
     &              1)
         CALL PUTLST(ICORE(I010),IPERT,1,1,IRREPX,177+ISPIN)
22       CONTINUE
20      CONTINUE
10     CONTINUE
      ENDIF
      RETURN
      END
