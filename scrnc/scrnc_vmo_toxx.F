
      SUBROUTINE SCRNC_VMO_TOXX(ICORE,MAXCOR,IUHF,LIST_SRC,LIST_TAR,
     +                         IRREPX,ICASE,TRANS)
C
C This routine drives the transformation of (PQ,RS) Coulomb integral
C to (XX,RS)  (X is AO index)
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION X,SNRM2,SDOT,ABLAD,ZILCH,ONE
      double precision sum,ddot 
      LOGICAL TAU, TRANS
      CHARACTER*2 SPLAB(3)
      CHARACTER*5 SIDE 
      DIMENSION ICORE(MAXCOR)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)

      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),ioffv(8,2),ioffo(8,2),
     &             IRPDPDAO(8),IRPDPDAOMO_OCCBK(8,2),
     &             IRPDPDAOMO_VRTBK(8,2),IRPDPDAOMO_OCCKB(8,2),
     &             IRPDPDAOMO_VRTKB(8,2),
     &             IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)

      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/INFO/NOCCO(2),NVRTO(2)
C
      DATA ZILCH,ONE /0.0D0,1.0D0/
      DATA SPLAB /'AA','BB','AB'/
C
      NNM1O2(I)=(I*(I-1))/2

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " Entered - scrnc_vmo_toao"
      sum=0.0d0
#endif

      NMO=NOCCO(1)+NVRTO(1)
      IONE=1
      CALL GETREC(20,'JOBARC','NBASTOT ',IONE,NAO)
C
C LOOP OVER SYMMETRIES OF BJ of (AI,BJ)
C
      DO 3 ISPIN=1, (IUHF+1)
       DO 4 IRREPBJ=1,NIRREP  
        IRREPAI=DIRPRD(IRREPBJ,IRREPX)
        IF (TRANS) THEN
           NUMBJ =IRPDPD(IRREPBJ,8+ISPIN)
           NUMAI =IRPDPDAO(IRREPAI)
           NUMAO1=NUMAI 
           NUMAO2=IRPDPDAO(IRREPBJ)

           I000=1
           I010=I000+NUMAI*NUMBJ*IINTFP
           I020=I010+NUMAO1*NUMAO2*IINTFP
           SIDE="LEFT "
        ELSE
           NUMAI=IRPDPD(IRREPAI,8+ISPIN)
           NUMBJ=IRPDPD(IRREPBJ,8+ISPIN)
           NUMAO=IRPDPDAO(IRREPAI)

           I000=1
           I010=I000+NUMAI*NUMBJ*IINTFP
           I020=I010+NUMAO*NUMBJ*IINTFP
           SIDE="LEFT "
        ENDIF 

        I030=I020+MAX(NAO*NMO,NT(1))*IINTFP
        I040=I030+MAX(NAO*NMO,NT(2))*IINTFP
        IEND=I040
        IF (TRANS) THEN
            I050 = I040+NUMAI*NUMBJ*IINTFP
            IEND = I050
        ENDIF 
        
        IF(IEND.GE.MAXCOR) CALL INSMEM('SCRNC_VMO_TOAO',IEND,MAXCOR)

        IF (TRANS) THEN
           CALL GETLST(ICORE(I040),1,NUMBJ,1,IRREPBJ,LIST_SRC+ISPIN)
           CALL TRANSP(ICORE(I040),ICORE(I000),NUMBJ,NUMAI)
        ELSE
           CALL GETLST(ICORE(I000),1,NUMBJ,1,IRREPBJ,LIST_SRC+ISPIN)
        ENDIF 

#ifdef _DEBUG_LVL0
        If (ispin .eq.1) then
        Call Checksum("SCRNC_VMO_TOAO(AI|BJ):",ICORE(I000),
     +                 NUMBJ*NUMAI)
        if (trans) then
        write(6,"(a,a,i2)") "scrnc_vmo_toxx(AI,XX);","irrepbj:",
     +  irrepbj
        call output(icore(i000),1,numbj,1,numai,numbj,numai,1)
        else
        write(6,"(a,a,i2)") "scrnc_vmo_toxx(AI,BJ);","irrepbj:",
     +  irrepbj
        call output(icore(i000),1,numai,1,numbj,numai,numbj,1)
        endif
        endif
#endif 
        CALL SCRNC_I1TRAN(ICASE,ICORE(I000),ICORE(I010),ICORE(I020),
     &                    ICORE(I030),ICORE(I040),MAXCOR,NAO,NMO,
     &                    ISPIN,IUHF,IRREPBJ,IRREPAI,.TRUE.,SIDE)
#ifdef _DEBUG_LVL0
        if (ispin .eq. 1) then
C        Call Checksum("SCRNC_VMO_TOAO(XI|BJ):",ICORE(I010),
C     +                 NUMAO*NUMBJ)
C        sum = sum + ddot(NUMAO*NUMBJ,ICORE(I010),1,ICORE(I010),1)
C        write(6,*) "sum =", sum
C         call output(icore(i010),1,NUMAO,1,NUMBJ,NUMAO,NUMBJ,1)
        If (Trans) then 
        write(6,"(a,a,i2)") "@-scrnc_vmo_toao (XX|XX);","irrepbj:",
     +  irrepbj
        call output(icore(i010),1,NUMAO1,1,NUMAO2,NUMAO1,NUMAO2,1)
        Else
        write(6,"(a,a,i2)") "@-scrnc_vmo_toao (XX|BJ);","irrepbj:",
     +  irrepbj
        call output(icore(i010),1,NUMAO,1,NUMBJ,NUMAO,NUMBJ,1)
        endif
        endif 
#endif 
        IF (TRANS) THEN
        CALL PUTLST(ICORE(I010),1,NUMAO1,1,IRREPAI,LIST_TAR+ISPIN)
        ELSE
        CALL PUTLST(ICORE(I010),1,NUMBJ,1,IRREPBJ,LIST_TAR+ISPIN)
        ENDIF 
4      CONTINUE
3     CONTINUE
C
      RETURN
      END
