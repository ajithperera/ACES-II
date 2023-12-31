      SUBROUTINE HEFF01(ICORE,MAXCOR,IUHF)
C
C FORMATION OF EFFECTIVE HAMILTONIAN FOR (0,1) SECTOR.  THE
C SPIN ORBITAL EQUATIONS ARE:
C
C   HEFF(ij) = -[F(ij) - F(mj) T(im) + F(me) T(jm,ie) 
C        AA        AA      IA    AI      NN    AN AN 
C
C              +  T(mn,ie) W(mn,je)]
C                   NN AN    NN AN
C  
CEND
      IMPLICIT INTEGER (A-H,O-Z)
      LOGICAL RHF
      DOUBLE PRECISION ONE,ONEM,ZILCH
      DIMENSION ICORE(MAXCOR)

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /FSSYMPOP/ FSDPDAN(8,22),FSDPDNA(8,22),FSDPDAA(8,22),
     &                  FSDPDIN(8,22),FSDPDNI(8,22),FSDPDII(8,22),
     &                  FSDPDAI(8,22),FSDPDIA(8,22)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /FSSYM/ POPA(8,2),POPI(8,2),VRTA(8,2),VRTI(8,2),
     &               NTAN(2),NTNA(2),NTAA(2),
     &               NTIN(2),NTNI(2),NTII(2),
     &               NTAI(2),NTIA(2),
     &               NF1AN(2),NF1AA(2),NF1IN(2),NF1II(2),NF1AI(2),
     &               NF2AN(2),NF2AA(2),NF2IN(2),NF2II(2),NF2AI(2)
C
      DATA ZILCH /0.0/
      DATA ONE    /1.0/
      DATA ONEM /-1.0/
C
      RHF=IUHF.EQ.0
C
C LOOP OVER SPIN CASES OF HEFF
C
      I000=1
      I010=I000+(NF1AA(1)+MIN(IUHF,1)*NF1AA(2))*IINTFP
       call izero(icore(i000),i010-i000)
      IOFFZ0=I000
      DO 10 ISPIN=1,IUHF+1
       I020=I010+MAX(NT(1),NT(2),NF1AI(1),NF1AI(2))*IINTFP
       I030=I020+NF1AI(ISPIN)*IINTFP
C
C INITIALIZE HEFF(0,1) WITH F(ij) INTERMEDIATES
C                             AA
C
       CALL FSGETT1(ICORE(IOFFZ0),ISPIN,91,'AA',20+ISPIN)
C
C ADD F(mj) * T(im) CONTRIBUTION.  ONLY ONE SPIN CASE.
C       IA      AI
C
       CALL FSGETT1(ICORE(I010),ISPIN,91,'IA',20+ISPIN)
       CALL FSGETT1(ICORE(I020),ISPIN,94,'AI',20+ISPIN)
       IOFFZ=IOFFZ0
       IOFFF=I010
       IOFFT=I020
       DO 100 IRREP=1,NIRREP
        NUMINACT=POPI(IRREP,ISPIN)
        NUMACT  =POPA(IRREP,ISPIN)
        CALL XGEMM('N','N',NUMACT,NUMACT,NUMINACT,ONEM,ICORE(IOFFT),
     &             NUMACT,ICORE(IOFFF),NUMINACT,ONE,ICORE(IOFFZ),
     &             NUMACT)
        IOFFZ=IOFFZ+NUMACT*NUMACT*IINTFP
        IOFFF=IOFFF+NUMINACT*NUMACT*IINTFP
        IOFFT=IOFFT+NUMINACT*NUMACT*IINTFP
100    CONTINUE 
C
C NOW DO F(me) * T(jm,ie).  
C          NN      AN AN
C
C        T(IJ,EM) * F(EM) + T(IJ,em) * F(em)  [ALPHA HEFF]
C          AA NN      NN      AA NN      NN
C
C        T(ij,em) * F(em) + T(ij,EM) * F(EM)  [ALPHA HEFF]
C          AA NN      NN      AA NN      NN
C
C        [2 * T(IJ,EM) - T(JI,EM)] * F(EM) [SPIN ADAPTED RHF]
C               AA NN      AA NN       NN
C
C USE RESORTED T AMPLITUDES.  NOTE THAT ONLY IRREP=1 CONTRACTION NEEDED.
C TWO SPIN CASES FOR UHF REFERENCES, SPIN ADAPTED CODE FOR RHF.
C
C RHF HBAR:
C
       IF(RHF)THEN
        LISTT =4
        DISSIZ=FSDPDAA(1,ISYTYP(1,LISTT))
        NUMDIS=IRPDPD (1,ISYTYP(2,LISTT))
        I030=I020+NUMDIS*DISSIZ*IINTFP
        I040=I030+MAX(NUMDIS,DISSIZ)*IINTFP
        I050=I040+MAX(NUMDIS,DISSIZ)*IINTFP
        I060=I050+MAX(NUMDIS,DISSIZ)*IINTFP
        CALL GETLST(ICORE(I010),1,1,1,1,93)
C
C READ T(Jm,Ie) AMPLITUDES INTO T(JI,me) AND REARRANGE TO T(IJ,em).
C        AN AN                    AA NN                     AA NN 
C
        CALL FSGET(ICORE(I020),1,NUMDIS,1,1,LISTT,'AANN')
        CALL SYMTR1 (1,POP(1,1),VRT(1,1),DISSIZ,ICORE(I020),
     &               ICORE(I030),ICORE(I040),ICORE(I050))
        CALL SYMTR3 (1,POPA(1,1),POPA(1,1),DISSIZ,NUMDIS,ICORE(I020),
     &               ICORE(I030),ICORE(I040),ICORE(I050))
        NROW=NF1AA(1)
        NCOL=1
        NSUM=NT(1)
        CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(I020),NROW,
     &             ICORE(I010),NSUM,ONE,ICORE(IOFFZ0),NROW)
C
C OPEN-SHELL HBAR:
C
       ELSE
C
C DO       T(JM,IE) * F(EM)
C            AN AN      NN 
C
        LISTT =ISPIN
        DISSIZ=FSDPDAA(1,ISYTYP(1,LISTT))
        NUMDIS=IRPDPD (1,ISYTYP(2,LISTT))
        I030=I020+NUMDIS*DISSIZ*IINTFP
        I040=I030+MAX(NUMDIS,DISSIZ)*IINTFP
        I050=I040+MAX(NUMDIS,DISSIZ)*IINTFP
        I060=I050+MAX(NUMDIS,DISSIZ)*IINTFP
C
C READ T(JM,IE) INTO T(JI,ME) AND REORDER TO T(IJ,EM)
C        AN AN         AA NN                   AA NN
C 
        CALL GETLST(ICORE(I010),1,1,1,ISPIN,93)
        CALL FSGET (ICORE(I020),1,NUMDIS,1,1,LISTT,'AANN') 
        CALL SYMTR1(1,POP(1,ISPIN),VRT(1,ISPIN),DISSIZ,ICORE(I020),
     &              ICORE(I030),ICORE(I040),ICORE(I050))
        CALL SYMTR3(1,POPA(1,ISPIN),POPA(1,ISPIN),DISSIZ,NUMDIS,
     &              ICORE(I020),ICORE(I030),ICORE(I040),ICORE(I050))
        NROW=NF1AA(ISPIN)
        NCOL=1
        NSUM=NT(ISPIN)
        CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(I020),NROW,
     &             ICORE(I010),NSUM,ONE,ICORE(IOFFZ0),NROW)
        IF(ISPIN.EQ.1)THEN
C
C DO       T(Jm,Ie) * F(em)
C            AN AN      NN 
C
         LISTT =4
         DISSIZ=FSDPDAA(1,ISYTYP(1,LISTT))
         NUMDIS=IRPDPD (1,ISYTYP(2,LISTT))
         I030=I020+NUMDIS*DISSIZ*IINTFP
         I040=I030+MAX(NUMDIS,DISSIZ)*IINTFP
         I050=I040+MAX(NUMDIS,DISSIZ)*IINTFP
         I060=I050+MAX(NUMDIS,DISSIZ)*IINTFP
         CALL GETLST(ICORE(I010),1,1,1,2,93) 
C
C READ T(Jm,Ie) INTO T(JI,me) AND REORDER TO T(IJ,em)
C        AN AN         AA NN                   AA NN
C
         CALL FSGET (ICORE(I020),1,NUMDIS,1,1,LISTT,'AANN')
         CALL SYMTR1(1,POP(1,2),VRT(1,2),DISSIZ,ICORE(I020),
     &               ICORE(I030),ICORE(I040),ICORE(I050))
         CALL SYMTR3(1,POPA(1,1),POPA(1,1),DISSIZ,NUMDIS,ICORE(I020),
     &               ICORE(I030),ICORE(I040),ICORE(I050))
         NROW=NF1AA(1)
         NCOL=1
         NSUM=NT(2)
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(I020),NROW,
     &              ICORE(I010),NSUM,ONE,ICORE(IOFFZ0),NROW)
        ELSE
C
C DO     F(EM) * T(Mj,Ei)
C          NN      NA NA
C
         LISTT =3
         DISSIZ=IRPDPD(1,ISYTYP(1,LISTT))
         NUMDIS=FSDPDAA(1,ISYTYP(2,LISTT))
         I030=I020+NUMDIS*DISSIZ*IINTFP
         I040=I030+MAX(NUMDIS,DISSIZ)*IINTFP
         I050=I040+MAX(NUMDIS,DISSIZ)*IINTFP
         I060=I050+MAX(NUMDIS,DISSIZ)*IINTFP
         CALL GETLST(ICORE(I010),1,1,1,1,93) 
C
C READ T(Mj,Ei) INTO T(ME,ji) AND REORDER TO T(EM,ij)
C        NA NA         NN AA                   NN AA
C
         CALL FSGET (ICORE(I020),1,NUMDIS,1,1,LISTT,'NNAA')
         CALL SYMTR1(1,POPA(1,2),POPA(1,2),DISSIZ,ICORE(I020),
     &               ICORE(I030),ICORE(I040),ICORE(I050))
         CALL SYMTR3(1,POP(1,1),VRT(1,1),DISSIZ,NUMDIS,ICORE(I020),
     &               ICORE(I030),ICORE(I040),ICORE(I050))
         NROW=1
         NCOL=NF1AA(2)
         NSUM=NT(1)
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(I010),NROW,
     &              ICORE(I020),NSUM,ONE,ICORE(IOFFZ0),NROW)
C
C END OF ISPIN=1/ISPIN=2 IF..THEN..ELSE
C
        ENDIF
C
C END OF RHF/UHF IF..THEN..ELSE
C
       ENDIF
C
       IOFFZ0=IOFFZ0+NF1AA(1)*IINTFP
C        
10    CONTINUE
C
C NOW AUGMENT HBAR WITH THE T2->T1 CONTRIBUTION
C
      CALL T2INH01(ICORE,MAXCOR,IUHF)
      IOFF=1 
      DO 50 ISPIN=1,1+IUHF
       CALL VMINUS(ICORE(IOFF),NF1AA(ISPIN))
       CALL FSPUTT1(ICORE(IOFF),ISPIN+2,91,'AA',20+ISPIN) 
       IOFF=IOFF+NF1AA(ISPIN)*IINTFP
50    CONTINUE 
      RETURN
      END
