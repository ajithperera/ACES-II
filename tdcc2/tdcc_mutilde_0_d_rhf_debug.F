C
      SUBROUTINE Tdcc_mutilde_0_d_rhf_debug(MIJBAR,MABBAR,MIABAR,
     &                                      ICORE,MAXCOR,IUHF,IRREPX,
     &                                      ISIDE,MBAR00)
C
C SOLVES FOR THE TWO-PARTICLE PART OF THE CANONICALLY TRANSFORMED DIPOLE
C MOMENT OPERATOR
C
C           _
C           M = exp(-T) M exp (T)
C
C   _                     _                  _
C   M(ab,ij) = T(ae,ij) * M(be) - T(ab,im) * M(mj)
C
C
C THE CONTRACTIONS ARE CARRIED OUT USING RING-ORDERED LISTS :
C                             _                  _
C   RHS(ai,bj) = - T(ai,bm) * M(mj) + T(ai,ej) * M(be) [Iside = 1]
C
C                             _                  _
C   LHS(ai,bj) = - L(ia,mb) * M(jm) + L(ia,je) * M(eb) [Iside = 2]
C
C This code is RHF reference only. Also notice that we can straightly
C use the same code for antisymmetric perturbations. 
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,SNRM2,HALF,sdot
      DOUBLE PRECISION MABBAR(*),MIJBAR(*),MIABAR(*)
      DOUBLE PRECISION MBAR00
      LOGICAL RHF
C
      DIMENSION ICORE(MAXCOR),IDID(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMLOC/ ISYMOFF(8,8,25)
      COMMON /FLAGS/IFLAGS(100)
      COMMON /FLAGS2/IFLAGS2(500)
C
      DATA ONE  /1.0/
      DATA ONEM /-1.0/
      DATA ZILCH/0.0/
      DATA HALF/0.5/

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a)") "---Entered tdcc_mutilde_0_d_rhf_debug---"
        Write(6,*)
#endif
C
C 08/94 Ajith 
C
      LSTTMP=302
      LISTT0=100*(ISIDE-1)+34
C
C LOOP OVER IRREPS 
C
      DO 100 IRREPZR=1,NIRREP
       IRREPZL=DIRPRD(IRREPZR,IRREPX)
       IRREPTR=IRREPZL
       IRREPTL=IRREPZL
       LISTT=LISTT0-1+4
       DISSYZ=IRPDPD(IRREPZL,ISYTYP(1,LISTT))
       NUMDSZ=IRPDPD(IRREPZR,ISYTYP(2,LISTT))
       DISSYT=IRPDPD(IRREPTL,ISYTYP(1,LISTT))
       NUMDST=IRPDPD(IRREPTR,ISYTYP(2,LISTT))
       MAXT=MAX(DISSYT,NUMDST,DISSYZ,NUMDSZ)
       Ndim_aibj= DISSYZ*NUMDSZ
C
       I000=1
       I010=I000+IINTFP*NUMDSZ*DISSYZ
       I020=I010+IINTFP*NUMDST*DISSYT
       I030=I020+IINTFP*IRPDPD(IRREPX,22)
       I040=I030+IINTFP*IRPDPD(IRREPX,20)
       I050=I040+IINTFP*MAXT
       I060=I050+IINTFP*MAXT
       I070=I060+IINTFP*MAXT
       Call Dzero(Icore(I000),Ndim_aibj)
C 
C 07/95 Memory checks were added by Ajith
C   
       IF (I070 .GE. MAXCOR) CALL INSMEM('OPRBAR2R', I070, MAXCOR)
C
C Disconnected term for left hand state [L(a,i)*Mbar(b,j)].
C But this is a connected when the totall contribution is considered
C and has to be retained (in both Hedio's, John's and quadratic).
C
        IF (IRREPZR .EQ. IRREPX .AND. ISIDE .EQ. 2) THEN
           CALL GETLST(ICORE(I010),1,1,1,1,190)
           CALL XGEMM ('N','N',DISSYZ,NUMDSZ,1,ONE,ICORE(I010),DISSYT,
     &                  MIABAR,1,ONE,ICORE(I000),DISSYT)
        ENDIF
#ifdef _DEBUG_LVLM
      Write(6,*)
      call checksum("Tdcc_mutilde_0_d_rhf_debug,mu^t_ab(1) :",
     +               Icore(I000),Ndim_aibj,S)
#endif
       CALL GETLST(ICORE(I010),1,NUMDST,1,IRREPTR,LISTT)
C
C Unlinked term for the left-hand state - MBAR00*L(ab,ij). See the
C comment above. IEOMPROP = 0 John's original implementation.
C IEOMPROP = 1 Hideo's modification.
C
        IF(IRREPX .EQ. 1 .AND. ISIDE .EQ.2) THEN
           CALL SAXPY(DISSYT*NUMDST,HALF*MBAR00,ICORE(I010),1,
     &                ICORE(I000),1)
        ENDIF
#ifdef _DEBUG_LVLM
      Write(6,*)
      call checksum("Tdcc_mutilde_0_d_rhf_debug,mu^t_ab(2) :",
     +               Icore(I000),Ndim_aibj,S)
#endif
C
C CARRY OUT  -T(AI,bm) * Q(mj) CONTRACTION
C
        IOFFQ0=I020
        IOFFT0=I010
        IOFFZ0=I000
        DO 110 IRREPM=1,NIRREP
         IRREPJ=DIRPRD(IRREPM,IRREPX)
         IRREPB=DIRPRD(IRREPM,IRREPTR)
         NUMM=POP(IRREPM,2)
         NUMJ=POP(IRREPJ,2)
         NUMB=VRT(IRREPB,2)
         NROW=DISSYT*NUMB
         NCOL=NUMJ
         NSUM=NUMM
         IOFFT=IOFFT0+IINTFP*DISSYT*(ISYMOFF(IRREPM,IRREPTR,10)-1)
         IOFFZ=IOFFZ0+IINTFP*DISSYZ*(ISYMOFF(IRREPJ,IRREPZR,10)-1)
         IF(ISIDE.EQ.1)THEN
          IOFFM=1+(ISYMOFF(IRREPJ,IRREPX,22)-1)
          CALL XGEMM('N','N',NROW,NCOL,NSUM,ONEM,ICORE(IOFFT),NROW,
     &               MIJBAR(IOFFM),NSUM,ZILCH,ICORE(IOFFZ),NROW)
         ELSEIF(ISIDE.EQ.2)THEN
          IOFFM=1+(ISYMOFF(IRREPM,IRREPX,22)-1)
          CALL XGEMM('N','T',NROW,NCOL,NSUM,ONEM,ICORE(IOFFT),NROW,
     &               MIJBAR(IOFFM),NCOL,ONE,ICORE(IOFFZ),NROW)
         ENDIF
110     CONTINUE
#ifdef _DEBUG_LVLM
      Write(6,*)
      call checksum("Tdcc_mutilde_0_d_rhf_debug,mu^t_ab(3) :",
     +               Icore(I000),Ndim_aibj,S)
#endif
C
C NOW TRANSPOSE KET INDICES
C
        CALL SYMTR1(IRREPZR,VRT(1,2),POP(1,2),DISSYT,ICORE(I000),
     &              ICORE(I040),ICORE(I050),ICORE(I060))
        CALL SYMTR1(IRREPTR,VRT(1,2),POP(1,2),DISSYT,ICORE(I010),
     &              ICORE(I040),ICORE(I050),ICORE(I060))
C
C NOW CARRY OUT T(AI,je) * F(eb) CONTRACTION
C
        IOFFQ0=I030
        IOFFT0=I010
        IOFFZ0=I000
        DO 120 IRREPE=1,NIRREP
         IRREPB=DIRPRD(IRREPE,IRREPX)
         IRREPJ=DIRPRD(IRREPE,IRREPTR)
         NUMJ=POP(IRREPJ,2)
         NUME=VRT(IRREPE,2)
         NUMB=VRT(IRREPB,2)
         NROW=DISSYT*NUMJ
         NCOL=NUMB
         NSUM=NUME
         IOFFT=IOFFT0+IINTFP*DISSYT*(ISYMOFF(IRREPE,IRREPTR,17)-1)
         IOFFZ=IOFFZ0+IINTFP*DISSYZ*(ISYMOFF(IRREPB,IRREPZR,17)-1)
         IF(ISIDE.EQ.1)THEN
          IOFFM=1+(ISYMOFF(IRREPB,IRREPX,19)-1)
          CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFT),NROW,
     &               MABBAR(IOFFM),NSUM,ONE,ICORE(IOFFZ),NROW)
         ELSEIF(ISIDE.EQ.2)THEN
          IOFFM=1+(ISYMOFF(IRREPE,IRREPX,19)-1)
          CALL XGEMM('N','T',NROW,NCOL,NSUM,ONE,ICORE(IOFFT),NROW,
     &               MABBAR(IOFFM),NCOL,ONE,ICORE(IOFFZ),NROW)
         ENDIF
120     CONTINUE
#ifdef _DEBUG_LVLM
      Write(6,*)
      call checksum("Tdcc_mutilde_0_d_rhf_debug,mu^t_ab(4) :",
     +               Icore(I000),Ndim_aibj,S)
#endif
C
C  TRANSPOSE KET INDICES OF TARGET AND WRITE TO RING LIST
C
        CALL SYMTR1(IRREPZR,POP(1,2),VRT(1,2),DISSYT,ICORE(I000),
     &              ICORE(I040),ICORE(I050),ICORE(I060))
        CALL SYMTR1(IRREPTR,POP(1,2),VRT(1,2),DISSYT,ICORE(I010),
     &              ICORE(I040),ICORE(I050),ICORE(I060))
C         
        CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPZR,LSTTMP)
100    CONTINUE
C
C SYMMETRIZE INCREMENTS 
C
       CALL IZERO(IDID,NIRREP)
       DO 200 IRREPR=1,NIRREP
         IF(IDID(IRREPR).EQ.0)THEN
         IRREPL=DIRPRD(IRREPR,IRREPX)
         NUMDIS=IRPDPD(IRREPR,ISYTYP(2,LSTTMP))
         DISSIZ=IRPDPD(IRREPL,ISYTYP(1,LSTTMP))
         I000=1
         I010=I000+IINTFP*NUMDIS*DISSIZ
C
C 07/95 Memory checks were added by Ajith       
C  
         IF (I010 .GE. MAXCOR) CALL INSMEM('OPRBAR2R', I020, MAXCOR)
C
         IF(IRREPX.EQ.1)THEN
          CALL GETLST(ICORE(I000),1,DISSIZ,1,IRREPL,LSTTMP)
          CALL MPMT  (ICORE(I000),DISSIZ)
          CALL PUTLST(ICORE(I000),1,DISSIZ,1,IRREPL,LSTTMP)
         ELSE
          I020=I010+IINTFP*NUMDIS*DISSIZ
C
C 07/95 Memory checks were added by Ajith       
C  
          IF (I020 .GE. MAXCOR) CALL INSMEM('OPRBAR2R', I020, MAXCOR)
C
          CALL GETLST(ICORE(I000),1,DISSIZ,1,IRREPL,LSTTMP)
          CALL TRANSP(ICORE(I000),ICORE(I010),DISSIZ,NUMDIS)
          CALL GETLST(ICORE(I000),1,NUMDIS,1,IRREPR,LSTTMP)
          CALL SAXPY (NUMDIS*DISSIZ,ONE,ICORE(I010),1,ICORE(I000),1)
          CALL PUTLST(ICORE(I000),1,NUMDIS,1,IRREPR,LSTTMP)
          CALL TRANSP(ICORE(I000),ICORE(I010),NUMDIS,DISSIZ)
          CALL PUTLST(ICORE(I010),1,DISSIZ,1,IRREPL,LSTTMP)
         ENDIF
         IDID(IRREPR)=1
         IDID(IRREPL)=1
        ENDIF
200    CONTINUE
C
C NOW SSTRNG THESE AND WRITE THEM TO THE FINAL DIPOLE LISTS
C
       ISIZE=IDSYMSZ(IRREPX,ISYTYP(1,LSTTMP),ISYTYP(2,LSTTMP))
       I000=1
       I010=I000+IINTFP*ISIZE
       I020=I010+IINTFP*ISIZE
C
C 07/95 Memory checks were added by Ajith       
C  
       IF (I020 .GE. MAXCOR) CALL INSMEM('OPRBAR2R', I020, MAXCOR)
C
       CALL GETALL(ICORE(I000),ISIZE,IRREPX,LSTTMP)
       CALL SSTGEN(ICORE(I000),ICORE(I010),ISIZE,VRT(1,1),POP(1,1),
     &             VRT(1,1),POP(1,1),ICORE(I020),IRREPX,'1324')
C
C  PUT RESULT DIRECTLY ON A PROPER LIST
C
#ifdef _DEBUG_LVL0
      Isize_pphh=idsymsz(irrepx,13,14)
      call checksum("Tdcc_mutilde_0_d_rhf_debug,mu^t_ab(f) :",
     +               Icore(I010),Isize_pphh,S)
#endif
CSSS       CALL PUTALL(ICORE(I010),ISIZE,IRREPX,326)
C
C       CALL SCOPY(ISIZE,ICORE(I010),1,MABIJBAR,1)
C
90    CONTINUE
      RETURN
      END
