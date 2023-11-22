       
      SUBROUTINE SINDEN(X, FAB, FIJ, IRREP0,ISPIN)
C
C FORM SINGLE EXCITATION DENOMINATOR ARRAY FOR VARIOUS DPD IRREP
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER POP,VRT,DIRPRD
      DIMENSION X(1),IOFFO(8),IOFFV(8),SCR(1), FAB(1,1), FIJ(1,1)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
C
      integer findex
      findex(i,j)=(i*(i-1))/2+j

      NBAS  = NOCCO(1)+NVRTO(1)
      NOCCA = NOCCO(1)
      NOCCB = NOCCO(2)
      NVRTA = NVRTO(1)
      NVRTB = NVRTO(2)
C
      IF(ISPIN.EQ.1)THEN
C   
        CALL GETREC(20, "JOBARC", "FOCKIJA ", NOCCA*NOCCA*IINTFP, FIJ)
        call output(FIJ, 1, NOCCA, 1, NOCCA, NOCCA, NOCCA, 1)
        CALL GETREC(20, "JOBARC", "FOCKABA ", NVRTA*NVRTA*IINTFP, FAB)
        
#ifdef _DEBUG_LVL0
        Write(6,*) "The vrt-vrt AFock matrix"
        call output(FAB, 1, NVRTA, 1, NVRTA, NVRTA, NVRTA, 1)
        Write(6,*) "The occ-occ AFock matrix"
        call output(FIJ, 1, NOCCA, 1, NOCCA, NOCCA, NOCCA, 1)
#endif

        CALL GETREC(20,'JOBARC','SCFEVALA',NBAS*IINTFP,SCR)
      ELSE
        CALL GETREC(20,'JOBARC','SCFEVALB',NBAS*IINTFP,SCR)

        CALL GETREC(20, "JOBARC", "FOCKIJB ", NOCCB*NOCCB*IINTFP, FIJ)
        CALL GETREC(20, "JOBARC", "FOCKABB ", NVRTB*NVRTB*IINTFP, FAB)
C
#ifdef _DEBUG_LVL0
        Write(6,*) "The vrt-vrt BFock matrix"
        call output(FAB, 1, NVRTB, 1, NVRTB, NVRTB, NVRTB, 1)
        Write(6,*) "The occ-occ BFock matrix"
        call output(FIJ, 1, NOCCB, 1, NOCCB, NOCCB, NOCCB, 1)
#endif
      ENDIF
cjp add corrections for non-hf character of higher references
      if(isbwcc) then
      do i=1,nbas
        scr(i)=scr(i)+fockcontr(findex(i,i),ispin)
      enddo
      endif
cjp end
   
C
C COMPUTE OCCUPIED AND VIRTUAL OFFSETS
C
      IOFFO(1)=1
      IOFFV(1)=NOCCO(ISPIN)+1
      DO 10 IRREP=2,NIRREP
       IOFFO(IRREP)=IOFFO(IRREP-1)+POP(IRREP-1,ISPIN)
       IOFFV(IRREP)=IOFFV(IRREP-1)+VRT(IRREP-1,ISPIN)
10    CONTINUE
C
      IOFF=1
      DO 20 IRRR=1,NIRREP
       IRRL=DIRPRD(IRRR,IRREP0)
       IOFFO0=IOFFO(IRRR)
       IOFFV00=IOFFV(IRRL)
       IOFFV0=IOFFV00
cjp avoid name clash ->ioccx
       DO 30 IOCCx=1,POP(IRRR,ISPIN)
        DO 40 IVRT=1,VRT(IRRL,ISPIN)
         X(IOFF)=SCR(IOFFV0)-SCR(IOFFO0)
         IOFF=IOFF+1
         IOFFV0=IOFFV0+1
40      CONTINUE
        IOFFO0=IOFFO0+1
        IOFFV0=IOFFV00
30     CONTINUE
20    CONTINUE
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The number of elements", ioff, NIRREP
      Write(6, "(5(1x,F10.5))") (X(I), I=1, IOff)
      write(6,*)
#endif
C
      Do IRREP0 = 1, NIRREP
      CALL BUILDX(X, FAB, FIJ, IRREP0, NOCCO(ISPIN), NVRTO(ISPIN), 
     &            ISPIN)
      Enddo

      STOP
      RETURN
      END
