      SUBROUTINE PCCD_E4S(ICORE,MAXCOR,IUHF,SE4)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP1,POP2,VRT1,VRT2
      DIMENSION ICORE(MAXCOR)
      LOGICAL LAMBDA 
      COMMON/SYM/ POP1(8),POP2(8),VRT1(8),VRT2(8),NTAA,NTBB,NF1AA,
     &            NF1BB,NF2AA,NF2BB
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /NHFREF/ NONHF
C
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &              DIRPRD(8,8)
      DATA ONE /1.0/

      LAMBDA = .FALSE.
      NOCCA=NOCCO(1)
      NOCCB=NOCCO(2)
      NVRTA=NVRTO(1)
      NVRTB=NVRTO(2)
      NWAA=NTAA
      NWBB=NTBB
C
      ITAA=MAXCOR+1-NWAA*IINTFP
      MXCOR=MAXCOR-NWAA*IINTFP
      IF(IUHF.EQ.0) THEN
       ITBB=ITAA
      ELSE
       ITBB=ITAA-NWBB*IINTFP
       MXCOR=MXCOR-NWBB*IINTFP
      ENDIF
C
      CALL GETLST(ICORE(ITAA),1,1,1,3,90)
      IF(IUHF.EQ.1) CALL GETLST(ICORE(ITBB),1,1,1,4,90)

#ifdef _DEBUG_LVLM
       Call checksum("T1-AA   :",Icore(itaa),nwaa)
       If (iuhf.ne.0) Call checksum("T1-BB   :",Icore(itbb),nwbb)
#endif
      CALL T2T1AB2(ICORE(ITAA),ICORE,MXCOR,POP1,POP2,VRT1,VRT2,1,
     &             IUHF,LAMBDA,0)
C
      IF(IUHF.EQ.1) THEN
      CALL T2T1AB2(ICORE(ITBB),ICORE,MXCOR,POP2,POP1,VRT2,VRT1,2,
     &             IUHF,LAMBDA,0)
      ENDIF
C
C  WITHIN THE SPIN ADAPTED RHF CODE NO CALL TO T2T1AA1 IS NECESSARY
C
       CALL T2T1AB1(ICORE(ITAA),ICORE,MXCOR,POP1,POP2,VRT1,VRT2,1,
     &              IUHF,LAMBDA,0)
C
       IF(IUHF.EQ.1) THEN
       CALL T2T1AB1(ICORE(ITBB),ICORE,MXCOR,POP2,POP1,VRT2,VRT1,2,
     &              IUHF,LAMBDA,0)
       ENDIF

#ifdef _DEBUG_LVL0
       Call checksum("T1-AA    :",Icore(itaa),nwaa)
       If (iuhf.ne.0) Call checksum("T1-BB    :",Icore(itbb),nwbb)
#endif

      IF(NONHF)THEN
       CALL GETLST(ICORE,1,1,1,3,93)
       CALL SAXPY(NWAA,ONE,ICORE,1,ICORE(ITAA),1)
       IF(IUHF.NE.0)THEN
        CALL GETLST(ICORE,1,1,1,4,93)
        CALL SAXPY(NWBB,ONE,ICORE,1,ICORE(ITBB),1)
       ENDIF
      ENDIF

#ifdef _DEBUG_LVL0
       Call checksum("T1-AA   :",Icore(itaa),nwaa)
       If (iuhf.ne.0) Call checksum("T1-BB    :",Icore(itbb),nwbb)
#endif
C
C CALCULATE NEW AMPLITUDES AND MBPT(4) SINGLE ENERGY
C
      SE4=0.0D0
      IDEN=1
      IE=IDEN+IINTFP*NWAA
      ID=IE+IINTFP*(NOCCA+NVRTA)
      IEND=ID+IINTFP*NWAA

      IF(IEND.LT.MXCOR) THEN
#ifdef _DEBUG_LVLM
        call checksum("T1AA-i :", ICORE(Itaa), NWAA)
#endif
       CALL PCCD_E4SENG(ICORE(ITAA),ICORE(ID),ICORE(IDEN),ICORE(IE),
     &                  NWAA,POP1,VRT1,NOCCA,NVRTA,SE4,1)
       CALL PUTLST(ICORE(ID),1,1,1,3,90)

#ifdef _DEBUG_LVL0
        call checksum("T1AA-in:", ICORE(ID), NWAA)
#endif
      ENDIF

      IF(IUHF.EQ.0) THEN
       SE4=2.0D0*SE4
      ELSE
       IDEN=1
       IE=IDEN+IINTFP*NWBB
       ID=IE+IINTFP*(NOCCB+NVRTB)
       IEND=ID+IINTFP*NWBB
       IF(IEND.LT.MXCOR) THEN
#ifdef _DEBUG_LVL0
        call checksum("T1BB-i:", ICORE(ITBB), NWBB)
#endif 
        CALL PCCD_E4SENG(ICORE(ITBB),ICORE(ID),ICORE(IDEN),ICORE(IE),
     &                   NWBB,POP2,VRT2,NOCCB,NVRTB,SE4,2)
         CALL PUTLST(ICORE(ID),1,1,1,4,90)
       ENDIF

      ENDIF
      Write(6,"(a,F15.10)") "The energy from singels: ", Se4

      RETURN
      END
