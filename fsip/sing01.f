      SUBROUTINE SING01(HEFF,ICORE,MAXCOR,IUHF)
C
C SINGLES AMPLITUDE EQUATION FOR (0,1) FOCK SPACE SECTOR.  THE SPIN
C  ORBITAL EQUATIONS ARE:
C
C   T1(ij) =   F(ij) - F(mj) T(im) + F(me) T(jm,ie) 
C      AI        AI      II    AI      NN    IN AN 
C
C              +  T(mn,ie) W(mn,je) - F(mj) HEFF(im)
C                   NN AN    NN IN      AI         AA
C
C NOTE THAT THIS IS VERY SIMILAR TO THE EQUATION FOR THE EFFECTIVE
C HAMILTONIAN AND THE CODE IS CONSEQUENTLY A PRETTY STRAIGHTFORWARD 
C (AND EXPLICIT) HACKUP.
C  
CEND
      IMPLICIT INTEGER (A-H,O-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      LOGICAL RHF
      DIMENSION ICORE(MAXCOR),HEFF(1)

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
      DATA ONE/1.0/, ONEM/-1.0/, ZILCH /0.0/
C
      RHF=IUHF.EQ.0
C
C LOOP OVER SPIN CASES OF AMPLITUDES
C
      I000=1
      I010=I000+(NF1AI(1)+MIN(IUHF,1)*NF1AI(2))*IINTFP
      IOFFZ0=I000
      IOFFH0=1
      DO 10 ISPIN=1,IUHF+1
       I020=I010+MAX(NT(1),NT(2),NF1II(1),NF1II(2))*IINTFP
       I030=I020+NF1AI(ISPIN)*IINTFP
C
C INITIALIZE T1 WITH F(ij) INTERMEDIATES
C                      AI
C
       CALL FSGETT1(ICORE(IOFFZ0),ISPIN,91,'AI',20+ISPIN)
C
C ADD F(mj) * T(im) AND H(im) * T(mj) CONTRIBUTIONS.  ONLY ONE SPIN CASE.
C       II      AI        AA      AI
C
       CALL FSGETT1(ICORE(I010),ISPIN,91,'II',20+ISPIN)
       CALL FSGETT1(HEFF ,2+ISPIN,91,'AA',20+ISPIN)
       CALL FSGETT1(ICORE(I020),ISPIN,94,'AI',20+ISPIN)
       IOFFZ=IOFFZ0
       IOFFF=I010
       IOFFT=I020
       IOFFH=IOFFH0
       DO 100 IRREP=1,NIRREP
        NUMINACT=POPI(IRREP,ISPIN)
        NUMACT  =POPA(IRREP,ISPIN)
        CALL XGEMM('N','N',NUMACT,NUMINACT,NUMINACT,ONEM,ICORE(IOFFT),
     &             NUMACT,ICORE(IOFFF),NUMINACT,ONE,ICORE(IOFFZ),
     &             NUMACT)
        CALL XGEMM('N','N',NUMACT,NUMINACT,NUMACT,ONEM,HEFF(IOFFH),
     &             NUMACT,ICORE(IOFFT),NUMACT,ONE,ICORE(IOFFZ),
     &             NUMACT)
        IOFFT=IOFFT+NUMINACT*NUMACT*IINTFP
        IOFFF=IOFFF+NUMINACT*NUMINACT*IINTFP
        IOFFH=IOFFH+NUMACT*NUMACT*IINTFP
        IOFFZ=IOFFZ+NUMACT*NUMINACT*IINTFP
100    CONTINUE 
C
C NOW DO F(me) * T(jm,ie).  
C          NN      IN AN
C
C        T(JM,IE) * F(EM) + T(Jm,Ie) * F(em)  [ALPHA HEFF]
C          IN AN      NN      IN AN      NN
C
C        T(jm,ie) * F(em) + T(jM,iE) * F(EM)  [BETA HEFF ]
C          IN AN      NN      IN AN      NN
C
C        [2 * T(Jm,Ie) - T(Mj,Ie)] * F(EM) [SPIN ADAPTED RHF]
C               IN AN      NI AN       NN
C
C USE RESORTED T AMPLITUDES.  NOTE THAT ONLY IRREP=1 CONTRACTION NEEDED.
C TWO SPIN CASES FOR UHF REFERENCES, SPIN ADAPTED CODE FOR RHF.
C
C RHF HBAR:
C
       IF(RHF)THEN
        LISTT =4
        DISSIZ=FSDPDIA(1,ISYTYP(1,LISTT))
        NUMDIS=IRPDPD (1,ISYTYP(2,LISTT))
        I030=I020+NUMDIS*DISSIZ*IINTFP
        I040=I030+MAX(NUMDIS,DISSIZ)*IINTFP
        I050=I040+MAX(NUMDIS,DISSIZ)*IINTFP
        I060=I050+MAX(NUMDIS,DISSIZ)*IINTFP
        CALL GETLST(ICORE(I010),1,1,1,1,93)
C
C READ T(Jm,Ie) AMPLITUDES INTO T(JI,me) AND REARRANGE TO T(IJ,em).
C        IN AN                    IA NN                     AI NN
C
        CALL FSGET(ICORE(I020),1,NUMDIS,1,1,LISTT,'IANN')
        CALL SYMTR1 (1,POP(1,1),VRT(1,1),DISSIZ,ICORE(I020),
     &               ICORE(I030),ICORE(I040),ICORE(I050))
        CALL SYMTR3 (1,POPI(1,1),POPA(1,1),DISSIZ,NUMDIS,ICORE(I020),
     &               ICORE(I030),ICORE(I040),ICORE(I050))
        NROW=NF1AI(1)
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
C            IN AN      NN 
C
        LISTT =ISPIN
        DISSIZ=FSDPDIA(1,ISYTYP(1,LISTT))
        NUMDIS=IRPDPD (1,ISYTYP(2,LISTT))
        I030=I020+NUMDIS*DISSIZ*IINTFP
        I040=I030+MAX(NUMDIS,DISSIZ)*IINTFP
        I050=I040+MAX(NUMDIS,DISSIZ)*IINTFP
        I060=I050+MAX(NUMDIS,DISSIZ)*IINTFP
C
C READ T(JM,IE) INTO T(JI,ME) AND REORDER TO T(IJ,EM)
C        IN AN         IA NN                   AI NN
C 
        CALL GETLST(ICORE(I010),1,1,1,ISPIN,93)
        CALL FSGET (ICORE(I020),1,NUMDIS,1,1,LISTT,'IANN') 
        CALL SYMTR1(1,POP(1,ISPIN),VRT(1,ISPIN),DISSIZ,ICORE(I020),
     &              ICORE(I030),ICORE(I040),ICORE(I050))
        CALL SYMTR3(1,POPI(1,ISPIN),POPA(1,ISPIN),DISSIZ,NUMDIS,
     &              ICORE(I020),
     &              ICORE(I030),ICORE(I040),ICORE(I050))
        NROW=NF1AI(ISPIN)
        NCOL=1
        NSUM=NT(ISPIN)
        CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(I020),NROW,
     &             ICORE(I010),NSUM,ONE,ICORE(IOFFZ0),NROW)
        IF(ISPIN.EQ.1)THEN
C
C DO       T(Jm,Ie) * F(em)
C            IN AN      NN 
C
         LISTT =4
         DISSIZ=FSDPDIA(1,ISYTYP(1,LISTT))
         NUMDIS=IRPDPD (1,ISYTYP(2,LISTT))
         I030=I020+NUMDIS*DISSIZ*IINTFP
         I040=I030+MAX(NUMDIS,DISSIZ)*IINTFP
         I050=I040+MAX(NUMDIS,DISSIZ)*IINTFP
         I060=I050+MAX(NUMDIS,DISSIZ)*IINTFP
         CALL GETLST(ICORE(I010),1,1,1,2,93) 
C
C READ T(Jm,Ie) INTO T(JI,me) AND REORDER TO T(IJ,em)
C        IN AN         IA NN                   AI NN
C
         CALL FSGET (ICORE(I020),1,NUMDIS,1,1,LISTT,'IANN')
         CALL SYMTR1(1,POP(1,2),VRT(1,2),DISSIZ,ICORE(I020),
     &               ICORE(I030),ICORE(I040),ICORE(I050))
         CALL SYMTR3(1,POPI(1,1),POPA(1,1),DISSIZ,NUMDIS,ICORE(I020),
     &               ICORE(I030),ICORE(I040),ICORE(I050))
         NROW=NF1AI(1)
         NCOL=1
         NSUM=NT(2)
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(I020),NROW,
     &              ICORE(I010),NSUM,ONE,ICORE(IOFFZ0),NROW)
        ELSE
C
C DO     F(EM) * T(Mj,Ei)
C          NN      NI NA
C
         LISTT =3
         DISSIZ=IRPDPD(1,ISYTYP(1,LISTT))
         NUMDIS=FSDPDIA(1,ISYTYP(2,LISTT))
         I030=I020+NUMDIS*DISSIZ*IINTFP
         I040=I030+MAX(NUMDIS,DISSIZ)*IINTFP
         I050=I040+MAX(NUMDIS,DISSIZ)*IINTFP
         I060=I050+MAX(NUMDIS,DISSIZ)*IINTFP
         CALL GETLST(ICORE(I010),1,1,1,1,93) 
C
C READ T(Mj,Ei) INTO T(ME,ji) AND REORDER TO T(EM,ij)
C        NI NA         NN IA                   NN AI
C
         CALL FSGET (ICORE(I020),1,NUMDIS,1,1,LISTT,'NNIA')
         CALL SYMTR1(1,POPI(1,2),POPA(1,2),DISSIZ,ICORE(I020),
     &               ICORE(I030),ICORE(I040),ICORE(I050))
         CALL SYMTR3(1,POP(1,1),VRT(1,1),DISSIZ,NUMDIS,ICORE(I020),
     &               ICORE(I030),ICORE(I040),ICORE(I050))
         NROW=1
         NCOL=NF1AI(2)
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
       IOFFZ0=IOFFZ0+NF1AI(1)*IINTFP
C        
10    CONTINUE
C
C NOW AUGMENT HBAR WITH THE T2->T1 CONTRIBUTION
C
      CALL T2INT101(ICORE,MAXCOR,IUHF)
C
C UPON RETURN, THE T1 VECTOR IS STORED AT THE BOTTOM OF CORE
C
      DO 1000 ISPIN=1,1+IUHF
       IOFFT=I000+(ISPIN-1)*NF1AI(1)*IINTFP
       CALL FSPUTT1(ICORE(IOFFT),ISPIN+2,94,'AI',20+ISPIN)
1000  CONTINUE
      RETURN
      END