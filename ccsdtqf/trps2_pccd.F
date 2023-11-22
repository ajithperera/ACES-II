      SUBROUTINE TRPS2_PCCD(T1A,T1B,FOVA,FOVB,D1T1A,D1T1B,S1A,S1B,
     1                 EVAL,SCR1,SCR2,SCR3,
     1                 GOOOV,BUF,
     1                 T2ABA,T2ABBI,T2ABBJ,T2AAAI,T2AAAJ,  
     1                 IADBLK,LENBLK,IADT2,LENT2,
     1                 IADV,LENV,OOOVAD,LENINT,
     1                 ICORE,
     1                 LNOOOV,NOCA,NOCB,NVRTA,NVRTB,NFREE,
     1                 INT1,INT2,NONHF,IUHF,LCCSD_T,FINAL,
     1                 EXTRNL_CCSD,PCCD)
      IMPLICIT INTEGER (A-Z)
      LOGICAL SKIP
      LOGICAL INT1,INT2,NONHF, FINAL
      LOGICAL IJEQL,NONEQL
      LOGICAL IJKEQL2,JKEQL2,IJEQL2,NONEQL2
      LOGICAL CCSDT4,CCSDT,LCCSD_T
      LOGICAL TRIPNI,TRIPNI1,TRIPIT,T3STOR
      LOGICAL EXTRNL_CCSD
      LOGICAL PCCD
      DOUBLE PRECISION T1A(1),T1B(1),FOVA(1),FOVB(1),
     1                 D1T1A(1),D1T1B(1),S1A(1),S1B(1),
     1                 SCR1(1),SCR2(1),SCR3(1),EVAL(NOCA+NVRTA,2)
      DOUBLE PRECISION FACT1,FACT2
C     GOOOV ijka gammas.
      DOUBLE PRECISION GOOOV(1),BUF(1)
      DOUBLE PRECISION T2ABA(1),T2ABBJ(1),T2ABBI(1)
      DOUBLE PRECISION T2AAAI(1),T2AAAJ(1)
      DOUBLE PRECISION ICORE(1)
C
      DOUBLE PRECISION E4TAAA,E4TBBB,E4TBAA,E4TABB,E5STA,E5STB,
     1                 E5ST2A,E5ST2B,E4STA,E4STB,E4ST2A,E4ST2B
      DOUBLE PRECISION E5TEST,SDOT, E4, DDOT,tmp,E6TEST
      DOUBLE PRECISION DOOA,DOOB,DVVA,DVVB
      DOUBLE PRECISION DIJK, suma, sumb
C
#include "maxbasfn.par"
C
C The three most relevent papers, for all the standard triple excitation
C method. Watss et al., JCP, 98, 8718, 1993 (JCP2)
C for L-CCSD(T), Crawford and Stanton, IJQC,70,601,Eqns. A16 (IJQC)
C
      DIMENSION OOOVAD(8,4),LENINT(8,4)
      DIMENSION LEN(8,8),IADBLK(8),LENBLK(8)
      DIMENSION                        IADT2(8) ,LENT2(8),
     1                                 IADV(8)  ,LENV(8)
      DIMENSION IADW(8),LENW(8)
      DIMENSION IADW2(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /INFO/   NOCCO(2),NVRTO(2)
      COMMON /FLAGS/  IFLAGS(100)
      EQUIVALENCE(ICLLVL,IFLAGS( 2))
      EQUIVALENCE(IDRLVL,IFLAGS( 3))
      EQUIVALENCE(IREFNC,IFLAGS(11))
      EQUIVALENCE(IQRHFP,IFLAGS(32))
      EQUIVALENCE(IQRHFM,IFLAGS(33))
      EQUIVALENCE(INONHF,IFLAGS(38))
      EQUIVALENCE(IORBTP,IFLAGS(39))
      EQUIVALENCE(INONHF,IFLAGS(38))
      EQUIVALENCE(IORBTP,IFLAGS(39))
C
C     T3 ARRAYS
C
      COMMON /TRIPLES/ TRIPNI,TRIPNI1,TRIPIT,T3STOR
      COMMON /ACTORB/ ABSVRT(MAXBASFN,8,2),ABSOCC(MAXBASFN,8,2)
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
      COMMON /T3ENRG/ E4TAAA,E4TBBB,E4TBAA,E4TABB,E5STA,E5STB,
     1                E4STA ,E4STB,E4TBAA2, E4TBAA3
      COMMON /DENST3/ DOOA(MAXBASFN),DOOB(MAXBASFN),
     1                DVVA(MAXBASFN),DVVB(MAXBASFN)
      COMMON /LISWI/  LWIC11,LWIC12,LWIC13,LWIC14,
     1                LWIC15,LWIC16,LWIC17,LWIC18,
     1                LWIC21,LWIC22,LWIC23,LWIC24,
     1                LWIC25,LWIC26,LWIC27,LWIC28,
     1                LWIC31,LWIC32,LWIC33,
     1                LWIC34,LWIC35,LWIC36,
     1                LWIC37,LWIC38,LWIC39,LWIC40,LWIC41,LWIC42
      COMMON /GAMLIS/ LISGO1,LISGO2,LISGO3,LISGO4,LISGV1,LISGV2,
     1                LISGV3,LISGV4
      COMMON /TESTEN/ E5TEST,E6TEST
      COMMON /AUXIO / DISTSZ(8,100),NDISTS(8,100),INIWRD(8,100),LNPHYR,
     1                NRECS,LUAUX
      COMMON /T3IOOF/ IJKPOS(8,8,8,2),IJKLEN(36,8,4),IJKOFF(36,8,4),
     1                NCOMB(4)
      COMMON /T3METH/ CCSDT4,CCSDT
C
      INDEX(I) = I*(I-1)/2
C
      WRITE(6,1015)
 1015 FORMAT(' @TRPS2-I, Spin case AAB ')
C
C     Initialize D1T1 and "S1" increments.
C
      CALL ZERO(D1T1A,NTAA)
      CALL ZERO(D1T1B,NTBB)
      CALL ZERO(S1A  ,NTAA)
      CALL ZERO(S1B  ,NTBB)

#ifdef _DEBUG_LVL0
      iprint=1
#endif
C
C     Get eigenvalues (more precisely, diagonal Fock matrix elements)
C
      CALL GETREC(20,'JOBARC','SCFEVALA',(NOCA+NVRTA)*IINTFP,
     1            EVAL(1,1))
      CALL GETREC(20,'JOBARC','SCFEVALB',(NOCB+NVRTB)*IINTFP,
     1            EVAL(1,2))
C
C     Determine the lengths and dimensions of the NIRREP**2 T3 AAB arrays.
C
C      WRITE(6,1070)
C 1070 FORMAT(' @TRPS2-I, Triples symmetry block information ')
      DO  100 IRPIJK=1,NIRREP
      DO   90 IRPC  =1,NIRREP
      IRPAB = DIRPRD(IRPIJK,IRPC)
      LEN(IRPC,IRPIJK) = IRPDPD(IRPAB,1) * VRT(IRPC,2)
C      WRITE(6,1080) IRPIJK,IRPC,IRPAB,LEN(IRPC,IRPIJK)
 1080 FORMAT('   FOR IJK BLOCK ',I4,/,
     1       '       C   BLOCK ',I4,' AB BLOCK ',I4,' LENGTH IS ',I6)
   90 CONTINUE
  100 CONTINUE
C
C     Compute addresses and lengths of OOOV quantities.
C
      DO  410 ISPIN=1,4
      DO  400 IRREP=1,NIRREP
      IF(ISPIN.EQ.1)THEN
      LENINT(IRREP,ISPIN) = IRPDPD(IRREP, 3) * IRPDPD(IRREP,16)
      ENDIF
      IF(ISPIN.EQ.2)THEN
      LENINT(IRREP,ISPIN) = IRPDPD(IRREP, 4) * IRPDPD(IRREP,17)
      ENDIF
      IF(ISPIN.EQ.3)THEN
      LENINT(IRREP,ISPIN) = IRPDPD(IRREP,14) * IRPDPD(IRREP,11)
      ENDIF
      IF(ISPIN.EQ.4)THEN
      LENINT(IRREP,ISPIN) = IRPDPD(IRREP,14) * IRPDPD(IRREP,18)
      ENDIF
  400 CONTINUE
  410 CONTINUE
C
      OOOVAD(1,1) = 1
      OOOVAD(1,2) = 1
      OOOVAD(1,3) = 1
      OOOVAD(1,4) = 1
      DO  430 ISPIN=1,4
      IF(NIRREP.GT.1)THEN
      DO  420 IRREP = 2,NIRREP
      OOOVAD(IRREP,ISPIN) = OOOVAD(IRREP-1,ISPIN) +
     1                      LENINT(IRREP-1,ISPIN)
  420 CONTINUE
      ENDIF
  430 CONTINUE
c      WRITE(6,1200)
 1200 FORMAT(' @TRPS2-I, OOOV ARRAY LENGTHS AND ADDRESSES ')
      DO  450 ISPIN=1,4
      DO  440 IRREP=1,NIRREP
c      WRITE(6,1210) ISPIN,IRREP,LENINT(IRREP,ISPIN),OOOVAD(IRREP,ISPIN)
 1210 FORMAT(' SPIN, IRREP, LENGTH, ADDRESS ',2I4,I6,I6)
  440 CONTINUE
  450 CONTINUE
C
      IF(IDRLVL.GT.0)THEN
C
      IF(IUHF.EQ.0)THEN
      IGPOS4 = 1
      ELSE
      IGPOS1 = 1
      IGPOS3 = 3
      IGPOS4 = 4
      ENDIF
      DO  9461 IRREP=1,NIRREP
      CALL GETTRN(GOOOV((IGPOS4-1)*LNOOOV + OOOVAD(IRREP,4)),BUF,
     1            IRPDPD(IRREP,14),IRPDPD(IRREP,18),
     1            2,IRREP,LISGO4)
      IF(IUHF.GT.0)THEN
      CALL GETTRN(GOOOV((IGPOS1-1)*LNOOOV + OOOVAD(IRREP,1)),BUF,
     1            IRPDPD(IRREP, 3),IRPDPD(IRREP,16),
     1            2,IRREP,LISGO1)
      CALL GETTRN(GOOOV((IGPOS3-1)*LNOOOV + OOOVAD(IRREP,3)),BUF,
     1            IRPDPD(IRREP,14),IRPDPD(IRREP,11),
     1            2,IRREP,LISGO3)
      CALL SYMTR3(IRREP,VRT(1,1),POP(1,2),IRPDPD(IRREP,11),
     1            IRPDPD(IRREP,14),
     1            GOOOV((IGPOS3-1)*LNOOOV + OOOVAD(IRREP,3)),
     1            SCR1,SCR2,SCR3)
      ENDIF
 9461 CONTINUE
      ENDIF

      IF (EXTRNL_CCSD) THEN

         IF(IUHF.EQ.0)THEN
            IGPOS4 = 1
         ELSE
            IGPOS1 = 1
            IGPOS3 = 3
            IGPOS4 = 4
         ENDIF

         DO IRREP=1,NIRREP
            CALL GETTRN(GOOOV((IGPOS4-1)*LNOOOV + OOOVAD(IRREP,4)),
     &                  BUF,IRPDPD(IRREP,14),IRPDPD(IRREP,18),
     &                  2,IRREP,110)

           call checksum("AAAA_IJKA",
     &           GOOOV((IGPOS4-1)*LNOOOV + OOOVAD(IRREP,4)),
     &           IRPDPD(IRREP,14)*IRPDPD(IRREP,18))

            IF (IUHF.GT.0)THEN
               CALL GETTRN(
     &              GOOOV((IGPOS1-1)*LNOOOV + OOOVAD(IRREP,1)),BUF,
     &              IRPDPD(IRREP, 3),IRPDPD(IRREP,16),2,IRREP,107)
              CALL GETTRN(
     &             GOOOV((IGPOS3-1)*LNOOOV + OOOVAD(IRREP,3)),BUF,
     &             IRPDPD(IRREP,14),IRPDPD(IRREP,11),2,IRREP,109)
              CALL SYMTR3(IRREP,VRT(1,1),POP(1,2),IRPDPD(IRREP,11),
     &             IRPDPD(IRREP,14),
     &             GOOOV((IGPOS3-1)*LNOOOV + OOOVAD(IRREP,3)),
     &             SCR1,SCR2,SCR3)
            ENDIF
         ENDDO
      ENDIF 

C     COMPUTE NUMBER OF FREE DOUBLE WORDS
      MCORE = NFREE / IINTFP
cYAU      WRITE(6,1230) NFREE,MCORE,IINTFP
 1230 FORMAT(' @TRPS2-I, Number of free integer words ',I10,/,
     1       '           Number of free double  words ',I10,/,
     1       '           Ratio                        ',I10)
C
      E4TBAA = 0.0D+00
      suma   = 0.0D+00
      sumb   = 0.0D+00
C
      ISPIN1 = 1
      ISPIN2 = 2
C
      DO 1000 IRPIJK=1,NIRREP
C
      IJKVAL = 0
C
      DO  990 IRPK  =1,NIRREP
      IF(POP(IRPK,ISPIN2).EQ.0) GOTO 990
C
      KLOW  = 1
      KHIGH = POP(IRPK,ISPIN2)
C
      DO  980 IRPJ=1,NIRREP
      IF(POP(IRPJ,ISPIN1).EQ.0) GOTO 980
      IRPJK = DIRPRD(IRPJ,IRPK)
      IRPI  = DIRPRD(IRPJK,IRPIJK)
      IF(IRPI.GT.IRPJ) GOTO 980
C
      IF(POP(IRPI,ISPIN1).EQ.0) GOTO 980
C
      IJEQL  = .FALSE.
      NONEQL = .FALSE.
      IF(IRPI.EQ.IRPJ)THEN
      IJEQL  = .TRUE.
      JLOW  = 2
      JHIGH = POP(IRPJ,ISPIN1)
      ILOW  = 1
      NIJ = (POP(IRPJ,ISPIN1) * (POP(IRPJ,ISPIN1)-1))/2
      ELSE
      NONEQL = .TRUE.
      JLOW  = 1
      JHIGH = POP(IRPJ,ISPIN1)
      ILOW  = 1
      IHIGH = POP(IRPI,ISPIN1)
      NIJ =  POP(IRPI,ISPIN1) * POP(IRPJ,ISPIN1)
      ENDIF
C
      NIK = POP(IRPI,ISPIN1) * POP(IRPK,ISPIN2)
      NJK = POP(IRPJ,ISPIN1) * POP(IRPK,ISPIN2)
C
      IF(IJEQL.AND.POP(IRPJ,ISPIN1).LT.2) GOTO 980
C
      IRPIJ =  DIRPRD(IRPI,IRPJ)
      IRPIK =  DIRPRD(IRPI,IRPK)
C
C     COMPUTE ADDRESSES AND LENGTHS OF ABC SYMMETRY BLOCKS FOR
C     THIS IJK. THIS LOOP IS EFFECTIVELY OVER IRREPS OF C (IE
C     COLUMNS OF THE EVENTUAL TARGET).
C
C     MAIN CORE VECTOR ADDRESSING : (DOUBLE PRECISION WORDS)
C
C           I000   SYMMETRY PACKED T3 FOR THIS IJK BLOCK
C           I010                   W
C           I020   UNEXPANDED T2(  ,IJ) AAAA FOR THIS IJ BLOCK
C           I030   EXPANDED   T2(  ,IJ) AAAA
C           I040   T2(  ,JK) ABAB FOR THIS JK BLOCK
C           I050   T2(  ,IK) ABAB FOR THIS IK BLOCK
C           ISTART CURRENT STARTING ADDRESS OF REUSABLE CORE
C
      I000 = 1
      DO  600 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
      IADBLK(IRREP) = 1
      ELSE
      IADBLK(IRREP) = IADBLK(IRREP-1) + LEN(IRREP-1,IRPIJK)
      ENDIF
      LENBLK(IRREP) = LEN(IRREP,IRPIJK)
  600 CONTINUE
C
      LENABC=0
      DO IRREP=1,NIRREP
         LENABC = LENABC + LENBLK(IRREP)
      END DO
C
      I010 = I000 + LENABC
      DO IRREP=1,NIRREP
         IRPA = DIRPRD(IRPIJK,IRREP)
         LENW(IRREP) = VRT(IRPA,1) * IRPDPD(IRREP,13)
      END DO
C
      IADW(1) = 1
      IF (NIRREP.GT.1) THEN
         DO IRREP=2,NIRREP
            IADW(IRREP) = IADW(IRREP-1) + LENW(IRREP-1)
         END DO
      END IF
C
      IOFFT2 = I010
      DO  641 IRREP=1,NIRREP
      IOFFT2 = IOFFT2 + LENW(IRREP)
  641 CONTINUE
C
      I020 = IOFFT2
C     GET VARIOUS "SMALL" BLOCKS OF T2 AND <IJ//KA> IN CORE, ORGANIZED
C     APPROPRIATELY.
C
C     BLOCK OF AAAA T(  ,IJ) WHERE I AND J ARE DIFFERENT SYMMETRY
C     BLOCKS. LENGTH IS IRPDPD(IRPIJ,1) * POP(IRPI,1) * POP(IRPJ,1).
C
      CALL GETLST(ICORE(I020),IOFFOO(IRPJ,IRPIJ,1)+1,NIJ,1,IRPIJ,44)
C
C     EXPAND INTO A DIFFERENT REGION OF CORE SO WE KEEP UNTOUCHED T2
C     AS WELL AS EXPANDED ONE. LENGTH IS DSZEXP * POP(IRPI,1) * POP(IRPJ,1)
C
      DSZ   = IRPDPD(IRPIJ,   ISPIN1)
      DSZEXP= IRPDPD(IRPIJ,18+ISPIN1)
      I030 = I020 + DSZ    * NIJ
      I040 = I030 + DSZEXP * NIJ
      CALL SYMEXP2(IRPIJ,VRT(1,1),DSZEXP,DSZ,NIJ,ICORE(I030),
     1             ICORE(I020))
C
C     BLOCK OF ABAB T(  ,JK). LENGTH IS IRPDPD(IRPJK,13) * NJK
C
      I050 = I040 + IRPDPD(IRPJK,13) * NJK
      CALL GETLST(ICORE(I040),IOFFOO(IRPK,IRPJK,5)+1,NJK,1,IRPJK,46)
C
C     BLOCK OF ABAB T(  ,IK). LENGTH IS IRPDPD(IRPIK,13) * NIK
C
      CALL GETLST(ICORE(I050),IOFFOO(IRPK,IRPIK,5)+1,NIK,1,IRPIK,46)
C
      I060 = I050 + IRPDPD(IRPIK,13) * NIK
C
      LNVOIJ  = IRPDPD(IRPIJ, 9)
      LNVOIK  = IRPDPD(IRPIK,12)
      LNVOJK  = IRPDPD(IRPJK,12)
      LNVOIK2 = IRPDPD(IRPIK,11)
      LNVOJK2 = IRPDPD(IRPJK,11)
C
C     SET A STARTING ADDRESS FOR REST OF CORE, THAT PART WHICH CAN BE
C     REUSED AFTER D1T2.
C
      I070 = I060 + LENABC
C 
C     If this is a CCSD(T) or QCISD(T) gradient calculation, allocate
C     space for "disconnected" triple amplitudes.
C
CSSS      IF ((IDRLVL.GT.0).AND.((ICLLVL.EQ.21).OR.(ICLLVL.EQ.22))) THEN
         I080 = I070 + LENABC
         I090 = I080 + (I020 - I010)
         I100 = I090 + LENABC
         I110 = I100 + LENABC
         I120 = I110 + LENABC
         ISTART = I120
CSSS      ELSE
CSSS         I090 = I070
CSSS         I100 = I090 + LENABC
CSSS         ISTART = I100
CSSS      END IF
      IF (LCCSD_T) THEN
C
C IL00 is L3 (symmetry packed)
C
          IL00 = ISTART
          IL10 = IL00 + LENABC
          IL90 = IL10 + (I020-I010)
          ISTART = IL90 + LENABC
      ENDIF
      NLEFT = MCORE - ISTART
C
C The following memory test is not correct. We use ICORE(ISTART) 
C to store quantities of length (a<b,c) for gradients and iterative
C triples. This might be the cause for mysterious crashes of large
C triples gradient calculations for low symmetry molecules.
C Ajith Perera 11/2005.
CSSS      IF (ISTART.GE.MCORE) THEN
CSSS         WRITE(*,*) '@TRPS2: Insufficient memory to continue.'
CSSS         CALL INSMEM('TRPS2',ISTART,MCORE)
CSSS         STOP 1
CSSS      END IF
      IF(ISTART+LENABC.GE.MCORE)THEN                                  
         WRITE(6,9010)                                                
         CALL INSMEM('TRPS2',ISTART+LENABC,MCORE)                    
      ENDIF                             
C
C     COMPUTE DENOMINATOR ARRAY
C
      CALL MKD32(ICORE(I060),EVAL,IADBLK,IRPIJK,NOCA,NOCB,NVRTA,NVRTB)
C
      DO  649    K=KLOW,KHIGH
C
      CALL GETABA(T2ABA,IRPK,K,.FALSE.)
C
      DO  648    J=JLOW,JHIGH
C
      CALL GETABB(T2ABBJ,IRPJ,J,.FALSE.)
      CALL GETAAA(T2AAAJ,IRPJ,J,1,.FALSE.)
C
      IF(IJEQL) IHIGH = J-1
C
      DO  647    I=ILOW,IHIGH
C
C     See if we can skip this triple.
C
c      CALL TSTLAB(IRPI,IRPJ,IRPK,I,J,K,SKIP)
c      IF(SKIP) GOTO 547
C
      CALL GETABB(T2ABBI,IRPI,I,.FALSE.)
      CALL GETAAA(T2AAAI,IRPI,I,1,.FALSE.)
C
      CALL ZERO(ICORE(I000),LENABC)
C
      LENGTHW = 0
      DO 9647 IRREP=1,NIRREP
      LENGTHW=LENGTHW + LENW(IRREP)
 9647 CONTINUE
      CALL ZERO(ICORE(I010),LENGTHW)
C
      IF (.NOT. FINAL) THEN

      CALL D1T2NEW( ICORE(I000),ICORE(I010),
     1           ICORE(I030),ICORE(I040),ICORE(I050),
     1 ICORE(ISTART),IUHF,
     1 IADBLK,LENBLK,IADW,LENW,IADT2,LENT2,IADV,LENV,
     1 I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPJK,IRPIK,IRPIJK,
     1 DSZEXP,IRPDPD(IRPJK,13),IRPDPD(IRPIK,13),
     1 SCR1,SCR2,SCR3,NLEFT,T2ABA,T2ABBJ,T2ABBI,
     1 T2AAAJ,T2AAAI)
C
c      IJK = IJKPOS(IRPI,IRPJ,IRPK,2)
c      IOFF = IJKOFF(IJK,IRPIJK,ISPIN1+1)
c      IF(IJEQL ) IJKVAL = IOFF + (K-1)*NIJ + INDEX(J-1) + I
c      IF(NONEQL) IJKVAL = IOFF + (K-1)*NIJ + (J-1)*POP(IRPI,ISPIN1) + I
      IJKVAL = IJKVAL + 1
C
C     CCSDT or non-HF if semicanonical orbitals are not used.
C     (i) Virtual orbital sum.
C
      IF(ICLLVL.EQ.18.OR.(INONHF.EQ.1.AND.IORBTP.EQ.0))THEN
      CALL T3FT323V(ICORE(I000),ICORE(I010),ICORE(ISTART),NFREE,IUHF,
     1              1,2,LENABC,LENGTHW,IRPIJK,IJKVAL,IADBLK,IADW)
      ENDIF
      IF(CCSDT4.OR.CCSDT)THEN
C     CCSDT abef intermediate contribution.
c     write(6,*) ' before t3wt323p '
c      call sumblk(icore(i000),lenabc)
      CALL T3WT323P(ICORE(I000),ICORE(I010),ICORE(ISTART),NFREE,IUHF,
     1              1,2,LENABC,LENGTHW,IRPIJK,IJKVAL,IADBLK,IADW)
c     write(6,*) ' after  t3wt323p '
c      call sumblk(icore(i000),lenabc)
CC     CCSDT maei intermediate contribution.
      CALL T3WT323R(ICORE(I000),ICORE(I010),ICORE(ISTART),NFREE,1,2,
     1              IADBLK,IADW,IRPI,IRPJ,IRPK,I,J,K,IRPIJK,IUHF)
c     write(6,*) ' after  t3wt323r '
c      call sumblk(icore(i000),lenabc)
      ENDIF
C
C     CCSDT or non-HF if semicanonical orbitals are not used.
C     (ii) Occupied orbital sum.
C
      IF(ICLLVL.EQ.18.OR.(INONHF.EQ.1.AND.IORBTP.EQ.0))THEN
      CALL T3FT323O(ICORE(I000),ICORE(ISTART),NFREE,IUHF,1,2,LENABC,
     1              IRPI,IRPJ,IRPK,IRPIJK,I,J,K)
      ENDIF
      IF(CCSDT4.OR.CCSDT)THEN
      CALL T3WT323HH(ICORE(I000),ICORE(ISTART),NFREE,IUHF,1,2,LENABC,
     1              IRPI,IRPJ,IRPK,IRPIJK,I,J,K)
c     write(6,*) ' after  t3wt323hh '
c      call sumblk(icore(i000),lenabc)
      ENDIF
C
      CALL EXPSC2(ICORE(I000),ICORE(I010),IADBLK,IADW,IRPIJK)
C
C     Remove denominators (D3T3 ---> T3)
C
      DIJK = EVAL(ABSOCC(I,IRPI,1),1) + EVAL(ABSOCC(J,IRPJ,1),1)
     1     + EVAL(ABSOCC(K,IRPK,2),2)
      CALL RMD314(ICORE(I000),ICORE(I060),LENABC,DIJK)
CSSS      Call checksum("T3(D)-(I)", ICORE(I000), Lenabc)
      ELSE
C
C If this iterative NONHF-CCSD(T), converged amplitudes are already
C on the disk.
C
      IJKVAL = IJKVAL + 1
#ifdef _NOSKIP
      IF(ICLLVL.EQ.18.OR.(INONHF.EQ.1.AND.IORBTP.EQ.0))THEN
      CALL ZERO(ICORE(I000), LENABC)
      CALL T3FT323V(ICORE(I000),ICORE(I010),ICORE(ISTART),NFREE,IUHF,
     1              1,2,LENABC,LENGTHW,IRPIJK,IJKVAL,IADBLK,IADW)
      CALL T3FT323O(ICORE(I000),ICORE(ISTART),NFREE,IUHF,1,2,LENABC,
     1              IRPI,IRPJ,IRPK,IRPIJK,I,J,K)
      ENDIF

      IF (ISTART + LENABC .GT. NLEFT) THEN
         CALL INSMEM("@-TRPS2", ISTART+LENABC, NLEFT)
      ENDIF 
      CALL GETLIST(ICORE(ISTART),IJKVAL,1,1,IRPIJK,1 + 1 + 4)
      E4TBAA = E4TBAA + Ddot(LENABC, ICORE(ISTART), 1, ICORE(I000), 1) 
#endif 
       CALL ZERO(ICORE(I000),LENABC)
       CALL ZERO(ICORE(I010),LENGTHW)
       CALL D1T2NEW( ICORE(I000),ICORE(I010),
     1           ICORE(I030),ICORE(I040),ICORE(I050),
     1 ICORE(ISTART),IUHF,
     1 IADBLK,LENBLK,IADW,LENW,IADT2,LENT2,IADV,LENV,
     1 I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPJK,IRPIK,IRPIJK,
     1 DSZEXP,IRPDPD(IRPJK,13),IRPDPD(IRPIK,13),
     1 SCR1,SCR2,SCR3,NLEFT,T2ABA,T2ABBJ,T2ABBI,
     1 T2AAAJ,T2AAAI)

       CALL EXPSC2(ICORE(I000),ICORE(I010),IADBLK,IADW,IRPIJK)
       CALL GETLIST(ICORE(ISTART),IJKVAL,1,1,IRPIJK,1 + 1 + 4)

       E4TBAA = E4TBAA + Ddot(LENABC, ICORE(ISTART), 1, ICORE(I000), 1) 
       E6TEST = E6TEST + Ddot(LENABC, ICORE(ISTART), 1, ICORE(I000), 1)
CSSS       WRITE(6,*) ' @TRPS2-I, E6ST from disconnected triples ',E6TEST
CSSS
CSSS      Write(*,*)
CSSS      Call checksum("T3(C)+T3D-(D)", ICORE(ISTART), Lenabc)

CSSS       CALL E4TIJK(ICORE(ISTART),ICORE(I060),LENABC,DIJK,E4TBAA,
CSSS     &             __FILE__)
      DIJK = EVAL(ABSOCC(I,IRPI,1),1) + EVAL(ABSOCC(J,IRPJ,1),1)
     1     + EVAL(ABSOCC(K,IRPK,2),2)
      CALL RMD314(ICORE(I000),ICORE(I060),LENABC,DIJK)
CSSS      Call checksum("T3(C)-(I)", ICORE(I000), Lenabc)
CSSS       CALL DAXPY(LENABC, -1.0D0, ICORE(I000), 1, ICORE(ISTART), 1)
      CALL DCOPY(LENABC, ICORE(ISTART), 1, ICORE(I070), 1)
      ENDIF 

C For L based CCSD(T), we also need to compute L3(ijk,abc)= L2(ij,{le}a)<bc||{le}k>
C and have both T3 and L3 contributions together. The current implementation
C requires two arrays of o^3v^3 vectors in memory.
C There is the possibility for an ilternative implementation
C T3(ijk,abc)o<ek||bc> = I(ij,ae) and E4 = I(ij,ae)oL(ij,ae). The memory
C requirments for the second schemem is slightly less but it is not as stright
C forward. The alternative implementaion can be considered later if the
C memory requirements become excessive. (Eqns. 9 in IJQC) Ajith Perera, 10/05.
C
       IF (LCCSD_T) THEN
C
C Get various "small" blocks of L2 and <IJ||KA> in core, organized
c appropriately.
C Block of AAAA L(  ,IJ) where I and J are different symmetry
C blocks. LENGTH is IRPDPD(IRPIJ,1) * POP(IRPI,1) * POP(IRPJ,1).
C
         CALL ZERO(ICORE(IL00), LENABC)
         CALL ZERO(ICORE(IL10), LENGTHW)
C
         CALL GETLST(ICORE(I020),IOFFOO(IRPJ,IRPIJ,1)+1,NIJ,1,IRPIJ,
     &               144)
C
C Expand into a different region of core so we keep untouched L2
C as well as expanded one. LENGTH is DSZEXP * POP(IRPI,1) * POP(IRPJ,1)
C
         CALL SYMEXP2(IRPIJ,VRT(1,1),DSZEXP,DSZ,NIJ,ICORE(I030),
     1                ICORE(I020))
C
C Block of abab T(  ,JK). LENGTH is irpdPD(IRPJK,13) * NJK
C
         CALL GETLST(ICORE(I040),IOFFOO(IRPK,IRPJK,5)+1,NJK,1,IRPJK,
     &               146)
C
C Block of ABAB T(  ,IK). LENGTH is IRPDPD(IRPIK,13) * NIK
C
         CALL GETLST(ICORE(I050),IOFFOO(IRPK,IRPIK,5)+1,NIK,1,IRPIK,
     &               146)
C
         CALL GETABA(T2ABA,IRPK,K,.TRUE.)
         CALL GETABB(T2ABBJ,IRPJ,J,.TRUE.)
         CALL GETAAA(T2AAAJ,IRPJ,J,1,.TRUE.)
         CALL GETABB(T2ABBI,IRPI,I,.TRUE.)
         CALL GETAAA(T2AAAI,IRPI,I,1,.TRUE.)

         CALL D1T2NEW(ICORE(IL00),ICORE(IL10),
     1                ICORE(I030),ICORE(I040),ICORE(I050),
     1                ICORE(ISTART),IUHF,
     1                IADBLK,LENBLK,IADW,LENW,IADT2,LENT2,IADV,LENV,
     1                I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPJK,IRPIK,IRPIJK,
     1                DSZEXP,IRPDPD(IRPJK,13),IRPDPD(IRPIK,13),
     1                SCR1,SCR2,SCR3,NLEFT,T2ABA,T2ABBJ,T2ABBI,
     1                T2AAAJ,T2AAAI)
C
         CALL EXPSC2(ICORE(IL00),ICORE(IL10),IADBLK,IADW,IRPIJK)
         
         CALL RMD314(ICORE(IL00),ICORE(I060),LENABC,DIJK)

         CALL E4TIJKL(ICORE(I000),ICORE(IL00),ICORE(I060),LENABC,DIJK,
     &                E4TBAA, __FILE__)
         
         CALL GETABA(T2ABA,IRPK,K,.FALSE.)
         CALL GETABB(T2ABBJ,IRPJ,J,.FALSE.)
         CALL GETAAA(T2AAAJ,IRPJ,J,1,.FALSE.)
         CALL GETABB(T2ABBI,IRPI,I,.FALSE.)
         CALL GETAAA(T2AAAI,IRPI,I,1,.FALSE.)
         
C
C Reload the T2 vectors back in to the appropriate locations, so any methods 
C other than CCSD(T) can  smoothly proceed.
C
         CALL GETLST(ICORE(I020),IOFFOO(IRPJ,IRPIJ,1)+1,NIJ,1,IRPIJ,
     &               44)
         CALL SYMEXP2(IRPIJ,VRT(1,1),DSZEXP,DSZ,NIJ,ICORE(I030),
     1                ICORE(I020))
         CALL GETLST(ICORE(I040),IOFFOO(IRPK,IRPJK,5)+1,NJK,1,IRPJK,
     &               46)
         CALL GETLST(ICORE(I050),IOFFOO(IRPK,IRPIK,5)+1,NIK,1,IRPIK,
     &               46)
      ELSE
C
C Compute T(CCSD) energy
C
         IF (.NOT. FINAL)
     &   CALL E4TIJK(ICORE(I000),ICORE(I060),LENABC,DIJK,E4TBAA,
     &                __FILE__)

         IF (EXTRNL_CCSD) THEN
#ifdef _DEBUG_LVL0
         if (iprint.eq.1) write(6,"(2a)") " In trps2.F,Entered",
     &                                    " EXTRNL_CC block"
         iprint=-1
#endif 
         ENDIF 
C
C This Endif is for LCCSD_T
C
      ENDIF
C
C     If this is a CCSD(T) or QCISD(T) gradient calculation, compute the
C     disconnected triple amplitudes.  For LCCSD(T) gradients, need the
C     disconnected lambda triple amplitudes, but NOT the disconnected
C     t-triple amplitudes.
C
CAP      IF(IDRLVL.GT.0.AND.(ICLLVL.EQ.21.OR.ICLLVL.EQ.22))THEN

      IF((IDRLVL.GT.0.AND.
     1  (ICLLVL.EQ.21.OR.ICLLVL.EQ.22).OR.EXTRNL_CCSD))THEN

      IF (.NOT. FINAL) THEN
      CALL ZERO(ICORE(I070),LENABC)
      CALL ZERO(ICORE(I010),IOFFT2 - I010)
      CALL S1S223(ICORE(ISTART),ICORE(I070),ICORE(I010),IADBLK,IADW,
     1            1,2,I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,
     1            IUHF,1,LCCSD_T)
      IF(NONHF)THEN
      CALL ZERO(ICORE(I110),LENABC)
      CALL ZERO(ICORE(I010),IOFFT2 - I010)
      CALL S1S223(ICORE(ISTART),ICORE(I110),ICORE(I010),IADBLK,IADW,
     1            1,2,I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,
     1            IUHF,2,LCCSD_T)

      CALL VADD(ICORE(I070),ICORE(I070),ICORE(I110),LENABC, 1.0D+00)
      ENDIF

#ifdef _DEBUG_NONHFTRP
      CALL RMD314(ICORE(I070),ICORE(I060),LENABC,DIJK)
      CALL DAXPY(LENABC, 1.0D0, ICORE(I000), 1, ICORE(I070),1)
      CALL E4TCTDIJK(ICORE(I000), ICORE(I070), ICORE(I060),
     &               LENABC,DIJK, E5TEST,__FILE__)
CSSS      Write(*,*)
CSSS      Call checksum("T3(C)+T3D-(I)", ICORE(I070), Lenabc)
CSSS      WRITE(6,*) ' @TRPS2-I, E5ST from disconnected triples ',E5TEST
#endif 
      ENDIF 
      IF (.NOT. FINAL)
     & CALL RMD314(ICORE(I070),ICORE(I060),LENABC,DIJK)
C
C     If (RHF), form all disconnected alpha amplitudes at I100.
      IF (IUHF.EQ.0) THEN
c YAU : old
c        CALL ICOPY(LENABC*IINTFP,ICORE(I070),1,ICORE(I100),1)
c YAU : new
         CALL DCOPY(LENABC,ICORE(I070),1,ICORE(I100),1)
c YAU : end
         CALL SC1SC2(ICORE(I100),ICORE(I070),LENABC/IINTFP,
     &               IADBLK,IRPIJK)
      END IF
c Endif for the (IDRLVL.GT.0.AND.(ICLLVL.EQ.21.OR.ICLLVL.EQ.22)       

CAP      END IF
      END IF
C
C     ICORE(I000)     A<B,C T3(C) 
C     ICORE(I070)     A<B,C T3(D)/L3(D)
C     ICORE(IL00)     A<B,C L3(C)
COLD     ICORE(IL70)     A<B,C L3(D)
C
C     If this is a QCISD(T) calculation, multiply T3(C) by 2.
C
      IF (ICLLVL.EQ.21) CALL XSCAL(LENABC,2.0D+00,ICORE(I000),1)
C
C     Expand alpha*t3(c)
C
CAP      IF(INT1.OR.INT2.OR.IDRLVL.GT.0)THEN

      IF(INT1.OR.INT2.OR.IDRLVL.GT.0.OR.EXTRNL_CCSD)THEN

      CALL SYMTRW2(ICORE(I000),ICORE(I010),ICORE(ISTART),
     1             IADBLK,IADW2,IRPIJK,1,2)

      IF (LCCSD_T) 
C     Expand alpha*l3(c)
     &   CALL SYMTRW2(ICORE(IL00),ICORE(IL10),ICORE(ISTART),
     &                IADBLK,IADW2,IRPIJK,1,2)
 
CAP      ENDIF
      ENDIF 
C
C     If this is rhf, compute all alpha amplitudes and put at ISTART,
C     compute E4TAAA, which must later be divided by 9. Then form effective
C     amplitude at I000.
C
      IF (IUHF.EQ.0) THEN
c YAU : old
c        CALL ICOPY(LENABC*IINTFP,ICORE(I000),1,ICORE(I090),1)
c YAU : new
         IF (LCCSD_T) THEN
            CALL DCOPY(LENABC,ICORE(IL00),1,ICORE(IL90),1)
            CALL SC1SC2(ICORE(IL90),ICORE(IL00),LENABC,IADBLK,IRPIJK)
         ENDIF 
         CALL DCOPY(LENABC,ICORE(I000),1,ICORE(I090),1)
         CALL SC1SC2(ICORE(I090),ICORE(I000),LENABC,IADBLK,IRPIJK)
c YAU : end
CSSS Is this a duplicate
         IF (ICLLVL.EQ.21) CALL XSCAL(LENABC,0.5D+00,ICORE(I090),1)
         IF (LCCSD_T) THEN
            CALL E4TIJKL(ICORE(I090), ICORE(IL90), ICORE(I060),
     &                   LENABC,DIJK,E4TAAA, __FILE__)
            CALL VADD(ICORE(IL00),ICORE(IL00),ICORE(IL90),LENABC,1.0D0)
         ELSE
            CALL E4TIJK(ICORE(I090),ICORE(I060),LENABC,DIJK,
     &                  E4TAAA, __FILE__)
         ENDIF
C
C This should be OK for iterated T3 since it is still at I000.
C
         IF (ICLLVL.EQ.21) CALL XSCAL(LENABC,2.0D+00,ICORE(I090),1)
         CALL VADD(ICORE(I000),ICORE(I000),ICORE(I090),LENABC,1.0D+00)
      END IF

      IF (PCCD) CALL DZERO(ICORE(I000),LENABC)
C
C     Include T3 in T1. In the case of CCSD(T) gradients this
C     WT3 (S(i,a)) can be added to Lambda (ie. these terms are
C     remains constant) but if we are doing L-CCSD(T) gradients, these
C     terms are the right hand side of a linear equations 
C     that detemine a new set of amplitudes (designated as epsilon
C     Crawford and Stanton, IJQC,70,601,Eqns. A16 and Eqn. 36 in JCP2).
C     See below for the formation S(ij,ab).

      IF(INT1)THEN
      CALL  T1T32N(D1T1A,D1T1B,S1A,S1B,
     1             ICORE(I000),ICORE(I010),
     1             ICORE(I020),ICORE(I040),ICORE(I050),
     1             SCR1,SCR2,SCR3,
     1             IADBLK,IADW2,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,
     1             IRPIJK,I,J,K,1,2,
     1             IRPDPD(IRPIJ,1),IRPDPD(IRPJK,13),IRPDPD(IRPIK,13),
     1             NONHF,LCCSD_T)
#ifdef _DEBUG_LVLM
      Call checksum("D1T1A", D1T1A, NTAA)
      Call checksum("D1T1B", D1T1A, NTBB)
#endif 
      ENDIF
C
      IF(IUHF.EQ.0)THEN
       CALL VADD(ICORE(I000),ICORE(I000),ICORE(I090),LENABC,-1.0D+00)
       IF (LCCSD_T)
     &    CALL VADD(ICORE(IL00),ICORE(IL00),ICORE(IL90),LENABC,-1.0D+00)
        IF(ICLLVL.EQ.21)THEN
         CALL SSCAL(LENABC,0.5D+00,ICORE(I090),1)
        ENDIF
      ENDIF
C
C     Form T3(C) contribution to GAMMA(IJ,AB). 
C     If this is a QCISD(T) or CCSD(T) gradient calculation, calculate
C     the contribution of T3(C) to GAMMA(IJ,AB). Gamma(IJ,AB)= <T1||T3>
C     in the case of CCSD(T) and Gamma(IJ,AB)= <L1||T3> for L-CCSD(T). 
C     (Eqn. 40 in JCP2 and part of the first terms in Eqn. A22 in IJQC)
C     For external CCSD this contribution is included in FORM_HT3C3 
C     (see below).

      IF(IDRLVL.GT.0.AND.(ICLLVL.EQ.21.OR.ICLLVL.EQ.22))THEN
      CALL T2FT32(ICORE(I000),ICORE(I010),ICORE(ISTART),
     1            IADBLK,IADW2,NLEFT,
     1            IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,
     1            ICLLVL,IUHF,1,LCCSD_T)
      IF(NONHF)THEN
      CALL T2FT32(ICORE(I000),ICORE(I010),ICORE(ISTART),
     1            IADBLK,IADW2,NLEFT,
     1            IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,
     1            ICLLVL,IUHF,2,LCCSD_T)
      ENDIF
      ENDIF
C
C     If this is a QCISD(T) calculation, divide 2*T3(C) by 2.
C
      IF(ICLLVL.EQ.21)THEN
      CALL SSCAL(LENABC     ,0.5D+00,ICORE(I000),1)
      CALL SSCAL(IOFFT2-I010,0.5D+00,ICORE(I010),1)
      ENDIF
C
C     Form alpha*t3(d) and expand it.
C
CAP      IF(IDRLVL.GT.0.AND.(ICLLVL.EQ.21.OR.ICLLVL.EQ.22))THEN

      IF((IDRLVL.GT.0.AND.
     1  (ICLLVL.EQ.21.OR.ICLLVL.EQ.22).OR.EXTRNL_CCSD))THEN

      IF(ICLLVL.EQ.21)THEN
      CALL SSCAL(LENABC,2.0D+00,ICORE(I070),1)
      CALL SSCAL(LENABC,2.0D+00,ICORE(I100),1)
      ENDIF
      CALL SYMTRW2(ICORE(I070),ICORE(I080),ICORE(ISTART),
     1             IADBLK,IADW2,IRPIJK,1,2)
      ENDIF

CAP      ENDIF 
C
C     We now have :
C                     I000 T3(C)
C                     IL00 L3(C)
C                     I010 expanded T3(C)
C                     IL10 expanded L3(C)
C                     I070 alpha*T3(D)/L3(D)
COLD                     IL70 alpha*L3(D)
C                     I080 expanded alpha*T3(D)/L3(D)
COLD                     IL80 expanded alpha*L3(D)
C                     I090 T3(C) all alpha if rhf
C                     IL90 L3(C) all alpha if rhf
C                     I100 alpha*T3(D)/L3(D)
C
CAP      IF(IDRLVL.GT.0)THEN

      IF (IDRLVL.GT.0 .OR. EXTRNL_CCSD) THEN
C
C     Put t3(c)+alpha*t3(d) or l3(c)+alpha*l3(d) at I070 and I080.
C
      IF(ICLLVL.EQ.21.OR.ICLLVL.EQ.22.OR.EXTRNL_CCSD)THEN

         IF (LCCSD_T) THEN
            CALL VADD(ICORE(I070),ICORE(IL00),ICORE(I070),LENABC
     &         ,1.0D+00)
            CALL VADD(ICORE(I080),ICORE(IL10),ICORE(I080),IOFFT2-I010
     &         ,1.0D+00)
            IF(IUHF.EQ.0)THEN
               CALL VADD(ICORE(I100),ICORE(IL90),ICORE(I100),LENABC
     &            ,1.0D+00)
            ENDIF

         ELSE

            IF (.NOT. FINAL) THEN

            CALL VADD(ICORE(I070),ICORE(I000),ICORE(I070),LENABC
     &         ,1.0D+00)
            CALL VADD(ICORE(I080),ICORE(I010),ICORE(I080),IOFFT2-I010
     &         ,1.0D+00)
            IF(IUHF.EQ.0)THEN
               CALL VADD(ICORE(I100),ICORE(I090),ICORE(I100),LENABC
     &            ,1.0D+00)
            ENDIF

            ENDIF 

         ENDIF
      ENDIF
C
      IF(ICLLVL.EQ.4.OR.ICLLVL.EQ.9.OR.ICLLVL.EQ.11)THEN
      JFLAG = 0
      CALL DOT32(DOOA,DOOB,ICORE(I000),ICORE(I000),
     1           LENABC,IRPI,IRPJ,IRPK,
     1           I,J,K,NOCCO(1),NOCCO(2),JFLAG)
C
      CALL DVT32(DVVA,DVVB,ICORE(I000),ICORE(I010),
     1                     ICORE(I000),ICORE(I010),IADBLK,IADW2,
     1           IRPIJK,NVRTO(1),NVRTO(2),JFLAG)
      IF(IUHF.EQ.0)THEN
      JFLAG = 1
      FACT2 = 1.0D+00 / DSQRT(3.0D+00)
      CALL SSCAL(LENABC,FACT2,ICORE(I090),1)
      CALL DOT32(DOOA,DOOB,ICORE(I090),ICORE(I090),
     1           LENABC,IRPI,IRPJ,IRPK,
     1           I,J,K,NOCCO(1),NOCCO(2),JFLAG)
C
      CALL DVT32(DVVA,DVVB,ICORE(I090),ICORE(I010),
     1                     ICORE(I090),ICORE(I010),IADBLK,IADW2,
     1           IRPIJK,NVRTO(1),NVRTO(2),JFLAG)
      ENDIF
      ENDIF
C
      IF(ICLLVL.EQ.21.OR.ICLLVL.EQ.22)THEN
      JFLAG = 0
      CALL DOT32(DOOA,DOOB,ICORE(I000),ICORE(I070),
     1           LENABC,IRPI,IRPJ,IRPK,
     1           I,J,K,NOCCO(1),NOCCO(2),JFLAG)
      
#ifdef _DEBUG_LVLM
      call checksum("DOOA", DOOA, NOCCO(1)*NOCCO(1))
#endif

      CALL DVT32(DVVA,DVVB,ICORE(I000),ICORE(I010),
     1                     ICORE(I070),ICORE(I080),IADBLK,IADW2,
     1           IRPIJK,NVRTO(1),NVRTO(2),JFLAG)

#ifdef _DEBUG_LVLM
      call checksum("DVVA", DVVA, NVRTO(1)*NVRTO(1))
#endif
C
      IF(IUHF.EQ.0)THEN
      JFLAG = 1
      FACT2 = 1.0D+00 / DSQRT(3.0D+00)
      CALL SSCAL(LENABC,FACT2,ICORE(I090),1)
      CALL SSCAL(LENABC,FACT2,ICORE(I100),1)
      CALL DOT32(DOOA,DOOB,ICORE(I090),ICORE(I100),
     1           LENABC,IRPI,IRPJ,IRPK,
     1           I,J,K,NOCCO(1),NOCCO(2),JFLAG)
C
      CALL DVT32(DVVA,DVVB,ICORE(I090),ICORE(I010),
     1                     ICORE(I100),ICORE(I080),IADBLK,IADW2,
     1           IRPIJK,NVRTO(1),NVRTO(2),JFLAG)
      ENDIF
      ENDIF
C
C     GRADIENTS : GAMMA EVALUATION.
C
C     Put t3(c)+(alpha/2)*t3(d) at I000 and I010
C
      IF(ICLLVL.EQ.21.OR.ICLLVL.EQ.22.OR.EXTRNL_CCSD)THEN
         IF (LCCSD_T) THEN
C     In the LCCSD(T) case, put l3(c)+alpha*l3(d) at IL00 and IL10
             call dcopy(lenabc,icore(i070),1,icore(il00),1)
             call dcopy(IOFFT2-I010,icore(i080),1,icore(il10),1)
         ELSE
      IF (.NOT. FINAL) THEN

            CALL VADD(ICORE(I000),ICORE(I000),ICORE(I070),LENABC
     &         ,1.0D+00)
            CALL VADD(ICORE(I010),ICORE(I010),ICORE(I080),IOFFT2-I010
     &         ,1.0D+00)
      CALL SSCAL(LENABC     ,0.5D+00,ICORE(I000),1)
      CALL SSCAL(IOFFT2-I010,0.5D+00,ICORE(I010),1)

      ENDIF 

      ENDIF
      ENDIF
C
      IF (LCCSD_T) THEN
C Scale Lambda3 and T3(c) because there is an implicit factor of two
C for gamma
         call dscal(lenabc,0.5d0,icore(il00),1)
         call dscal(IOFFT2-I010,0.5d0,icore(il10),1)
         call dscal(lenabc,0.5d0,icore(i000),1)
         call dscal(IOFFT2-I010,0.5d0,icore(i010),1)
C L3 - T2
         CALL GVT32(ICORE(IL00),ICORE(IL10),ICORE(ISTART),IADBLK,IADW2,
     &      IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,SCR1,SCR2,
     &      SCR3,IUHF,ICORE(I030),ICORE(I040),ICORE(I050),DSZEXP,
     &      IRPDPD(IRPJK,13),IRPDPD(IRPIK,13))
         CALL GOT32(GOOOV,LNOOOV,OOOVAD,ICORE(IL00),ICORE(IL10),
     &      ICORE(ISTART),IADBLK,IADW2,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,
     &      IRPJK,I,J,K,IUHF,T2ABA,T2ABBJ,T2ABBI,T2AAAJ,T2AAAI,IGPOS4)

C     Load Lambda_2 vectors into memory
         CALL GETLST(ICORE(I020),IOFFOO(IRPJ,IRPIJ,1)+1,NIJ,1,IRPIJ,
     &               144)
         CALL SYMEXP2(IRPIJ,VRT(1,1),DSZEXP,DSZ,NIJ,ICORE(I030),
     1                ICORE(I020))
         CALL GETLST(ICORE(I040),IOFFOO(IRPK,IRPJK,5)+1,NJK,1,IRPJK,
     &               146)
         CALL GETLST(ICORE(I050),IOFFOO(IRPK,IRPIK,5)+1,NIK,1,IRPIK,
     &               146)
         CALL GETABA(T2ABA,IRPK,K,.TRUE.)
         CALL GETABB(T2ABBJ,IRPJ,J,.TRUE.)
         CALL GETAAA(T2AAAJ,IRPJ,J,1,.TRUE.)
         CALL GETABB(T2ABBI,IRPI,I,.TRUE.)
         CALL GETAAA(T2AAAI,IRPI,I,1,.TRUE.)
C     Add these contributions with t3(c) to the gamma matrices
         CALL GVT32(ICORE(I000),ICORE(I010),ICORE(ISTART),IADBLK,IADW2,
     &      IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,SCR1,SCR2,
     &      SCR3,IUHF,ICORE(I030),ICORE(I040),ICORE(I050),DSZEXP,
     &      IRPDPD(IRPJK,13),IRPDPD(IRPIK,13))
         CALL GOT32(GOOOV,LNOOOV,OOOVAD,ICORE(I000),ICORE(I010),
     &      ICORE(ISTART),IADBLK,IADW2,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,
     &      IRPJK,I,J,K,IUHF,T2ABA,T2ABBJ,T2ABBI,T2AAAJ,T2AAAI,IGPOS4)
C     Reload t2 vectors
         CALL GETABA(T2ABA,IRPK,K,.FALSE.)
         CALL GETABB(T2ABBJ,IRPJ,J,.FALSE.)
         CALL GETAAA(T2AAAJ,IRPJ,J,1,.FALSE.)
         CALL GETABB(T2ABBI,IRPI,I,.FALSE.)
         CALL GETAAA(T2AAAI,IRPI,I,1,.FALSE.)
         CALL GETLST(ICORE(I020),IOFFOO(IRPJ,IRPIJ,1)+1,NIJ,1,IRPIJ,
     &               44)
         CALL SYMEXP2(IRPIJ,VRT(1,1),DSZEXP,DSZ,NIJ,ICORE(I030),
     1                ICORE(I020))
         CALL GETLST(ICORE(I040),IOFFOO(IRPK,IRPJK,5)+1,NJK,1,IRPJK,
     &               46)
         CALL GETLST(ICORE(I050),IOFFOO(IRPK,IRPIK,5)+1,NIK,1,IRPIK,
     &               46)
C Rescale Lambda3 and T3
         call dscal(lenabc,2.0d0,icore(il00),1)
         call dscal(IOFFT2-I010,2.0d0,icore(il10),1)
         call dscal(lenabc,2.0d0,icore(i000),1)
         call dscal(IOFFT2-I010,2.0d0,icore(i010),1)

      ELSE
        IF (.NOT. EXTRNL_CCSD) THEN
#ifdef _DEBUG_LVLM
          call checksum("T32(C)---",icore(i000),lenabc)
          call checksum("T32(d)---",icore(i010),iofft2-i010)
#endif 
         CALL GVT32(ICORE(I000),ICORE(I010),ICORE(ISTART),IADBLK,IADW2,
     &      IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,SCR1,SCR2,
     &      SCR3,IUHF,ICORE(I030),ICORE(I040),ICORE(I050),DSZEXP,
     &      IRPDPD(IRPJK,13),IRPDPD(IRPIK,13))
         CALL GOT32(GOOOV,LNOOOV,OOOVAD,ICORE(I000),ICORE(I010),
     &      ICORE(ISTART),IADBLK,IADW2,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,
     &      IRPJK,I,J,K,IUHF,T2ABA,T2ABBJ,T2ABBI,T2AAAJ,T2AAAI,IGPOS4)
        ENDIF 

      ENDIF
C
CAP      ENDIF
      ENDIF 

      IF (EXTRNL_CCSD) THEN
         CALL GETLST(ICORE(I020),IOFFOO(IRPJ,IRPIJ,1)+1,NIJ,1,
     &               IRPIJ,14)
         CALL SYMEXP2(IRPIJ,VRT(1,1),DSZEXP,DSZ,NIJ,ICORE(I030),
     &                ICORE(I020))
         CALL GETLST(ICORE(I040),IOFFOO(IRPK,IRPJK,5)+1,NJK,1,
     &               IRPJK,16)
         CALL GETLST(ICORE(I050),IOFFOO(IRPK,IRPIK,5)+1,NIK,1,  
     &               IRPIK,16)
         CALL SYMTRW2(ICORE(I000),ICORE(I010),ICORE(ISTART),
     &                IADBLK,IADW2,IRPIJK,1,2)
#ifdef _DEBUG_LVLM
         call checksum("T32(C)---",icore(i000),lenabc)
         call checksum("T32(d)---",icore(i010),iofft2-i010)
#endif 
         CALL T32ABCI2(ICORE(I000),ICORE(I010),ICORE(ISTART),
     &                 IADBLK,IADW2,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,
     &                 IRPJK,IRPIJK,I,J,K,SCR1,SCR2,SCR3,IUHF,
     &                 ICORE(I030),ICORE(I040),ICORE(I050),DSZEXP,
     &                 IRPDPD(IRPJK,13),IRPDPD(IRPIK,13))
         CALL GETLST(ICORE(I020),IOFFOO(IRPJ,IRPIJ,1)+1,NIJ,1,
     &               IRPIJ,44)
         CALL SYMEXP2(IRPIJ,VRT(1,1),DSZEXP,DSZ,NIJ,ICORE(I030),
     &                ICORE(I020))
         CALL GETLST(ICORE(I040),IOFFOO(IRPK,IRPJK,5)+1,NJK,1,
     &               IRPJK,46)
         CALL GETLST(ICORE(I050),IOFFOO(IRPK,IRPIK,5)+1,NIK,1,  
     &               IRPIK,46)

         CALL GETWABA(T2ABA,IRPK,K,.FALSE.)
         CALL GETWABB(T2ABBJ,IRPJ,J,.FALSE.)
         CALL GETWAAA(T2AAAJ,IRPJ,J,1,.FALSE.)
         CALL GETWABB(T2ABBI,IRPI,I,.FALSE.)
         CALL GETWAAA(T2AAAI,IRPI,I,1,.FALSE.)

         CALL T32IJKA2(GOOOV,LNOOOV,OOOVAD,ICORE(I000),ICORE(I010),
     &                 ICORE(ISTART),IADBLK,IADW2,IRPI,IRPJ,IRPK,
     &                 IRPIJ,IRPIK,IRPJK,I,J,K,IUHF,T2ABA,T2ABBJ,
     &                 T2ABBI,T2AAAJ,T2AAAI,IGPOS4)
         CALL GETABA(T2ABA,IRPK,K,.FALSE.)
         CALL GETABB(T2ABBJ,IRPJ,J,.FALSE.)
         CALL GETAAA(T2AAAJ,IRPJ,J,1,.FALSE.)
         CALL GETABB(T2ABBI,IRPI,I,.FALSE.)
         CALL GETAAA(T2AAAI,IRPI,I,1,.FALSE.)
      ENDIF 
C
      IF(INT2)THEN
C
C     --- D2T2 = F T3 for non-Hartree-Fock cases ---
C
      IF((ICLLVL.EQ.13.AND.NONHF).OR.
     1   (ICLLVL.GE.14.AND.ICLLVL.LE.18).OR.ICLLVL.EQ.33.OR.
     &                                      ICLLVL.EQ.34)THEN
      CALL T2FT32(ICORE(I000),ICORE(I010),ICORE(ISTART),
     1            IADBLK,IADW2,NLEFT,
     1            IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,
     1            ICLLVL,IUHF,0,.FALSE.)
      ENDIF
C
      IF (EXTRNL_CCSD .AND. NONHF) THEN
         CALL T2FT32(ICORE(I000),ICORE(I010),ICORE(ISTART),
     1               IADBLK,IADW2,NLEFT,
     1               IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,
     1               ICLLVL,IUHF,0,.FALSE.)
      ENDIF 
C
C     --- D2T2 = W T3 ---
C
      CALL T2T32(ICORE(I000),ICORE(I010),ICORE(ISTART),
     1           IADBLK,IADW2,IRPI,IRPJ,IRPK,
     1           IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,IUHF,
     1           SCR1,SCR2,SCR3,.FALSE.)

      CALL T2T32O(ICORE(I000),ICORE(I010),ICORE(ISTART),IADBLK,IADW2,
     1            IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,
     1            IUHF,SCR1,SCR2,SCR3,.FALSE.)
      IF (LCCSD_T) THEN
C     --- D2L2 = L3 W ---
         CALL T2T32(ICORE(IL00),ICORE(IL10),ICORE(ISTART),
     1      IADBLK,IADW2,IRPI,IRPJ,IRPK,
     1      IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,IUHF,
     1      SCR1,SCR2,SCR3,.TRUE.)
C
         CALL T2T32O(ICORE(IL00),ICORE(IL10),ICORE(ISTART),IADBLK,IADW2,
     1      IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,
     1      IUHF,SCR1,SCR2,SCR3,.TRUE.)
      ENDIF
      ENDIF

C This is no longer needed since T3 -> T2 are computed by default since
C INT1 and INT2 is true for extrnl_ccsd
CSSS      IF (EXTRNL_CCSD) THEN
CSSS         CALL FORM_HT3C2(ICORE(I000),ICORE(I010),SCR1,SCR2,
CSSS     &                   SCR3,ICORE(ISTART),
CSSS     &                   IADBLK,IADW,IADW2,NLEFT,IRPI,IRPJ,
CSSS     &                   IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,
CSSS     &                   ICLLVL,IUHF,0,.FALSE.,NONHF)
CSSS      ENDIF 
C
c      IJK = IJKPOS(IRPI,IRPJ,IRPK,2)
c      IOFF = IJKOFF(IJK,IRPIJK,ISPIN1+1)
c      IF(IJEQL ) IJKVAL = IOFF + (K-1)*NIJ + INDEX(J-1) + I
c      IF(NONEQL) IJKVAL = IOFF + (K-1)*NIJ + (J-1)*POP(IRPI,ISPIN1) + I
C     write(6,*) ' @TRPS2-I, sym, record ',IRPIJK,IJKVAL
      IF(T3STOR)THEN
      CALL PUTLIST_PCCD(ICORE(I000),IJKVAL,1,1,IRPIJK,ISPIN1 + 1,PCCD)
      ENDIF
#ifdef _DEBUG_LVLM
      Write(*,*)
      call checksum("TRPS2-1", ICORE(I000), Lenabc)
#endif
C
C     In RHF, form AAA (A<B<C) from AAB and write to disk.
C
      IF(IUHF.EQ.0.AND.T3STOR)THEN
c      write(6,*) ' @trps2-i, writing triples check sums '
c      call sumblk(icore(i000),distsz(irpijk,2))
c      IF(IUHF.EQ.0)THEN
C
C     Check if i<j<k. Need only check k>j.
C
c      write(6,7199) i,j,k,irpi,irpj,irpk
c 7199 FORMAT(6I6)
      IF(IRPK.GT.IRPJ.OR.(IRPK.EQ.IRPJ.AND.K.GT.J))THEN
c      write(6,*) i,j,k
C
      IJKEQL2 = .FALSE.
      IJEQL2  = .FALSE.
      JKEQL2  = .FALSE.
      NONEQL2 = .FALSE.
      IF(IRPI.EQ.IRPJ.AND.IRPI.EQ.IRPK) IJKEQL2 = .TRUE.
      IF(IRPJ.EQ.IRPK.AND.IRPI.NE.IRPK)  JKEQL2 = .TRUE.
      IF(IRPI.EQ.IRPJ.AND.IRPI.NE.IRPK)  IJEQL2 = .TRUE.
      IF(IRPI.NE.IRPJ.AND.IRPI.NE.IRPK.AND.IRPJ.NE.IRPK)
     1                                  NONEQL2 = .TRUE.
C
      IJK  = IJKPOS(IRPI,IRPJ,IRPK,1)
      IOFF = IJKOFF(IJK,IRPIJK,1)
      IF(IJKEQL2) IJKAAA = IOFF + ((K-1)*(K-2)*(K-3))/6 + INDEX(J-1) + I
      IF(IJEQL2 ) IJKAAA = IOFF + (K-1)*NIJ             + INDEX(J-1) + I
      IF(JKEQL2 ) IJKAAA = IOFF + (INDEX(K-1)+ J-1)*POP(IRPI,1)  + I
      IF(NONEQL2) IJKAAA = IOFF + (K-1)*NIJ + (J-1)*POP(IRPI,1)  + I
c      write(6,*) ' @trps2-i, calling symcontw '
      call zero(icore(istart),distsz(irpijk,1))
      CALL SYMCONTW(ICORE(istart),ICORE(I000),IADBLK,1,IRPIJK)
c      call sumblk(icore(istart),distsz(irpijk,1))
c      call sumblk(icore(i000),distsz(irpijk,2))
c      write(6,*) i000,i090,distsz(irpijk,1),distsz(irpijk,2)
      CALL PUTLIST_PCCD(ICORE(istart),IJKAAA,1,1,IRPIJK,1,PCCD)
#ifdef _DEBUG_LVLM
      Write(*,*) 
      call checksum("TRPS2-2", ICORE(istart), distsz(irpijk,1))
#endif
      ENDIF
      ENDIF
C
  547 CONTINUE
C
  647 CONTINUE
  648 CONTINUE
  649 CONTINUE
C
  980 CONTINUE
  990 CONTINUE
C
 1000 CONTINUE

      IF (PCCD) Then
         Write(6,"(20x,a)") " ----------Warning------------"
         Write(6,"(2a)")  " A pCCD or pCCD like calculation is being",
     +                    " performed and T3 is set to zero"
         Write(6,"(20x,a)") " -----------------------------"
      Endif

      IF (.NOT. FINAL) THEN
      WRITE(6,3000) E4TBAA
 3000 FORMAT(' @TRPS2-I, E4TBAA ',F20.12)
      IF(IUHF.EQ.0)THEN
      E4TAAA = E4TAAA / 9.0D+00
      WRITE(6,3002) E4TAAA
      ENDIF

      ELSE

      WRITE(6,3001) E4TBAA
 3001 FORMAT(' @TRPS2-I, E4TBAA+E5ST+E4ST ',F20.12)
      IF(IUHF.EQ.0)THEN
      E4TAAA = E4TAAA / 9.0D+00
      WRITE(6,3002) E4TAAA
 3002 FORMAT(' @TRPS2-I, E4TAAA ',F20.12)
      ENDIF
     
      ENDIF 

      IF (EXTRNL_CCSD) THEN
         CALL GETLST(SCR1,1,1,1,3,90)
#ifdef _DEBUG_LVLM
         call checksum("D1T1A:",SCR1,NTAA)
         call checksum("D1T1A:",D1T1A,NTAA)
#endif 
         CALL VADD(SCR1,SCR1,D1T1A,NTAA,1.0D+00)
         CALL PUTLST(SCR1,1,1,1,3,90)
#ifdef _DEBUG_LVLm
         call checksum("D1T1A:",SCR1,NTAA)
#endif   
         IF (IUHF .GT. 0) THEN
            CALL GETLST(SCR1,1,1,1,3+IUHF,90)
#ifdef _DEBUG_LVLM
         call checksum("D1T1B:",SCR1,NTBB)
         call checksum("D1T1B:",D1T1B,NTBB)
#endif 
            CALL VADD(SCR1,SCR1,D1T1B,NTBB,1.0D+00)
            CALL PUTLST(SCR1,1,1,1,3+IUHF,90)
#ifdef _DEBUG_LVLM
         call checksum("D1T1B:",SCR1,NTBB)
#endif 
         ELSE
            CALL VADD(SCR1,SCR1,D1T1B,NTBB,1.0D+00)
            CALL PUTLST(SCR1,1,1,1,3,90)
         ENDIF 
         DO IRREP=1,NIRREP
            CALL PUTTRN(GOOOV((IGPOS4-1)*LNOOOV + OOOVAD(IRREP,4)),
     &           BUF,IRPDPD(IRREP,14),IRPDPD(IRREP,18),
     &           2,IRREP,110)

           call checksum("AAAA_IJKA",
     &           GOOOV((IGPOS4-1)*LNOOOV + OOOVAD(IRREP,4)),
     &           IRPDPD(IRREP,14)*IRPDPD(IRREP,18))

            IF (IUHF.GT.0)THEN
               CALL PUTTRN(GOOOV((IGPOS1-1)*LNOOOV + OOOVAD(IRREP,1)),
     &              BUF,IRPDPD(IRREP, 3),IRPDPD(IRREP,16),
     &              2,IRREP,107)
              CALL SYMTR3(IRREP,POP(1,2),VRT(1,1),IRPDPD(IRREP,11),
     &             IRPDPD(IRREP,14),
     &             GOOOV((IGPOS3-1)*LNOOOV + OOOVAD(IRREP,3)),
     &             SCR1,SCR2,SCR3)
              CALL PUTTRN(GOOOV((IGPOS3-1)*LNOOOV + OOOVAD(IRREP,3)),
     &             BUF,IRPDPD(IRREP,14),IRPDPD(IRREP,11),
     &             2,IRREP,109)
            ENDIF
         ENDDo
#ifdef _DEBUG_LVL0
      Write(6,"(a)") " ---After trps2,Lists:61-63,90(3,4)---"
CSSS      call check_ints(icore(istart),nleft,iuhf,.false.)
      call check_t32w(icore(istart),nleft,iuhf)
      call check_t2(icore(istart),nleft,iuhf)
#endif
      ENDIF 
C
      IF(IDRLVL.GT.0)THEN
C
      DO  9463 IRREP=1,NIRREP
      CALL PUTTRN(GOOOV((IGPOS4-1)*LNOOOV + OOOVAD(IRREP,4)),BUF,
     1            IRPDPD(IRREP,14),IRPDPD(IRREP,18),
     1            2,IRREP,LISGO4)
      IF(IUHF.GT.0)THEN
      CALL PUTTRN(GOOOV((IGPOS1-1)*LNOOOV + OOOVAD(IRREP,1)),BUF,
     1            IRPDPD(IRREP, 3),IRPDPD(IRREP,16),
     1            2,IRREP,LISGO1)
      CALL SYMTR3(IRREP,POP(1,2),VRT(1,1),IRPDPD(IRREP,11),
     1            IRPDPD(IRREP,14),
     1            GOOOV((IGPOS3-1)*LNOOOV + OOOVAD(IRREP,3)),
     1            SCR1,SCR2,SCR3)
      CALL PUTTRN(GOOOV((IGPOS3-1)*LNOOOV + OOOVAD(IRREP,3)),BUF,
     1            IRPDPD(IRREP,14),IRPDPD(IRREP,11),
     1            2,IRREP,LISGO3)
      ENDIF
 9463 CONTINUE
      ENDIF
C
      IF((ICLLVL.GE.13.AND.ICLLVL.LE.18).OR.ICLLVL.EQ.33.OR.
     &                                      ICLLVL.EQ.34)THEN
      CALL GETLST(SCR1,1,1,1,3,90)
      CALL   VADD(SCR1,SCR1,D1T1A,NTAA,1.0D+00)
      CALL PUTLST(SCR1,1,1,1,3,90)
C      IF(IUHF.EQ.1)THEN
      CALL GETLST(SCR1,1,1,1,3+IUHF,90)
      CALL   VADD(SCR1,SCR1,D1T1B,NTBB,1.0D+00)
      CALL PUTLST(SCR1,1,1,1,3+IUHF,90)
C      ENDIF
      ENDIF
C
C     If this is CCSD(T) or QCISD(T) gradient, write D1T1 triples contribution
C     to disk.
C
      IF(IDRLVL.GT.0.AND.(ICLLVL.EQ.21.OR.ICLLVL.EQ.22))THEN
C The list 104 is used eleswhere.
CSSS        CALL GETLST(SCR1,1,1,1,1,104)
        CALL GETLST(SCR1,1,1,1,1,4)
        CALL   VADD(SCR1,SCR1,D1T1A,NTAA,1.0D+00)
        CALL PUTLST(SCR1,1,1,1,1,4)
        IF(NONHF)THEN
          CALL GETLST(SCR1,1,1,1,5,160)
          CALL   VADD(SCR1,SCR1,S1A,NTAA,0.5D+00)
          CALL PUTLST(SCR1,1,1,1,5,160)
        ENDIF
C
        IF(IUHF.EQ.0)THEN
          CALL GETLST(SCR1,1,1,1,1,4)
          CALL   VADD(SCR1,SCR1,D1T1B,NTBB,1.0D+00)
          CALL PUTLST(SCR1,1,1,1,1,4)
        ELSE
          CALL GETLST(SCR1,1,1,1,2,4)
          CALL   VADD(SCR1,SCR1,D1T1B,NTBB,1.0D+00)
          CALL PUTLST(SCR1,1,1,1,2,4)
          IF(NONHF)THEN
            CALL GETLST(SCR1,1,1,1,6,160)
            CALL   VADD(SCR1,SCR1,S1B,NTBB,0.5D+00)
            CALL PUTLST(SCR1,1,1,1,6,160)
          ENDIF
        ENDIF
      ENDIF
C
      IF(ICLLVL.EQ.9.OR.ICLLVL.EQ.13.OR.ICLLVL.EQ.14.OR.ICLLVL.EQ.15.OR.
     1   ICLLVL.EQ.16.OR.ICLLVL.EQ.17.OR.ICLLVL.EQ.18) RETURN
      IF(ICLLVL.EQ.33.OR.ICLLVL.EQ.34) RETURN
      IF(ICLLVL.EQ.4.AND.IREFNC.LE.1.AND..NOT.NONHF) RETURN
C
C COMPUTE CCSD(T) OR QCISD(T) FIFTH-ORDER PIECE OF ENERGY
C
C For L-based CCSD instead of T1, we use L1 to compute the
C 5th order triples contribution. Ajith Perera, 10/05.
C
      IF (.NOT. FINAL) THEN
C
C For iterative NON-HF triples, the total contribution is
C already computed)

      IF (LCCSD_T) THEN 
         CALL GETLST(T1A,1,1,1,1,190)
         IF (IUHF.EQ.0) THEN
             CALL GETLST(T1B,1,1,1,1,190)
         ELSE
             CALL GETLST(T1B,1,1,1,2,190)
         END IF
      ELSE
         CALL GETLST(T1A,1,1,1,1,90) 
         IF (IUHF.EQ.0) THEN
             CALL GETLST(T1B,1,1,1,1,90)
         ELSE
             CALL GETLST(T1B,1,1,1,2,90)
         END IF
      ENDIF 
C
      E5ST2A = 0.0D0
      DO AI = 1, NTAA
         E5ST2A = E5ST2A + ( T1A(AI) * D1T1A(AI) )
      END DO
      E5STA = E5STA + E5ST2A
C
      E5ST2B = 0.0D0
      DO AI = 1, NTBB
         E5ST2B = E5ST2B + ( T1B(AI) * D1T1B(AI) )
      END DO
      E5STB = E5STB + E5ST2B
C
      WRITE(*,'(a,f20.15)') ' @TRPS2-I, E5ST A ',E5STA
      WRITE(*,'(a,f20.15)') ' @TRPS2-I, E5ST B ',E5STB
C
      IF (NONHF) THEN
C
C        COMPUTE NON-HF CCSD(T) OR QCISD(T) FOURTH-ORDER PIECE.
C
         CALL GETLST(FOVA,1,1,1,3,93)
         IF (IUHF.EQ.0) THEN
            CALL GETLST(FOVB,1,1,1,3,93)
         ELSE
            CALL GETLST(FOVB,1,1,1,4,93)
         END IF
C
         E4ST2A = 0.0D0
         DO AI = 1, NTAA
            E4ST2A = E4ST2A + ( FOVA(AI) * S1A(AI) )
         END DO
         E4STA = E4STA + E4ST2A
C
         E4ST2B = 0.0D0
         DO AI = 1, NTBB
            E4ST2B = E4ST2B + ( FOVB(AI) * S1B(AI) )
         END DO
         E4STB = E4STB + E4ST2B
C
         WRITE(*,'(a,f20.15)') ' @TRPS2-I, E4ST A ',E4STA
         WRITE(*,'(a,f20.15)') ' @TRPS2-I, E4ST B ',E4STB
C
C     END IF (NONHF)
      END IF

      END IF 
 9010 FORMAT(' @TRPS2-I, Insufficient memory to continue. ')
#ifdef _DEBUG_LVL0
      if (.not. extrnl_ccsd) then
      Write(6,"(a)") " ---After trps2,Lists:61-63,90(3,4)---"
CSSS      call check_ints(icore(istart),nleft,iuhf,.false.)
CSSS      If (Idrlvl .ge.0) call check_t32w_grad(icore(istart),nleft,Iuhf)
CSSS      call check_t32w(icore(istart),nleft,Iuhf)
      call check_t2(icore(istart),nleft,iuhf)
      call check_t2grd(icore(istart),nleft,iuhf)
C      irrepx=1
C      Lenab=IDSYMSZ(IRREPX,ISYTYP(1,63),ISYTYP(2,63))
C      iend = istart+lenab
C      nleft = nleft - iend 
C      call Form_t2aaaa(icore(istart),icore(iend),nleft,lenab)
      endif 
#endif

      RETURN
      END
