

      SUBROUTINE MKT1DBL(W, WT, Z, T1ALPA, T1BETA, DISSYW, DISSYZ,
     &                   DISTMP,NUMSYW, NUMSYZ, NUMTMP, NOCCSQ, POP1,
     &                   POP2, VRT1, VRT2, LISTW, LISTZ, FACT, ITRANS, 
     &                   INCREM, IRREPWL, IRREPWR, IRREPQL, IRREPQR,
     &                   IRREPX)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER DISSYW,DISSYZ,DISTMP,DIRPRD,POP1,POP2,VRT1,VRT2
      DIMENSION W(DISSYW,NUMSYW),WT(NUMSYW,DISSYW),Z(DISTMP,NUMTMP),
     &          T1ALPA(1),T1BETA(1),POP1(8),POP2(8),VRT1(8),VRT2(8),
     &          IWOFF(8),IZOFF(8),ITOFFA(8),ITOFFB(8)
C     
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FILES/ LUOUT,MOINTS
C
      DATA AZERO,ONE /0.D0,1.D0/
C
      IF (MIN(NUMSYZ, DISSYZ) .EQ. 0) RETURN
C
C Get the transpose of the integral list and transpose last tow indices
C Note that in CC code they used the symmetry of the two electron 
C integrals, which is not true for Hbar elements.
C
      CALL GETLST(Z, 1, NUMSYW, 2, IRREPWR, LISTW)
      CALL TRANSP(Z, W, NUMSYW, DISSYW)
      CALL SCOPY(NUMSYW*DISSYW, W, 1, Z, 1)
C
      CALL SYMTRA(IRREPWL, VRT2, POP1, NUMSYW, Z, WT)
C
C Zero target array
C
      CALL ZERO(Z, DISTMP*NOCCSQ)
C
C Now integrals are ordered (AJ,ME). This case NUMSYW = DISTMP
C                                      I
C Perform first multiplication, SUM E T(ALPHA)*Hbar(AJ,ME)
C                                      E
      IOFFZ = 1
      IWOFF(1) = 1
      IOFF  = 1
C 
      DO 5000 IRREPE = 2, NIRREP   
C
         IRREP = IRREPE - 1
         IRREPM = DIRPRD(IRREP, IRREPWL)
C
         IWOFF(IRREPE) = IWOFF(IRREPE - 1) + POP1(IRREPM)*VRT2(IRREP)
C
 5000 CONTINUE
C 
      DO 5100 IRREPI = 1, NIRREP   
C
         IRREPE = DIRPRD(IRREPI, IRREPX)
C
         ITOFFA(IRREPE) = IOFF
         IOFF = IOFF + VRT2(IRREPE)*POP2(IRREPI)
C
 5100 CONTINUE
C
      DO 100 IRREPI = 1, NIRREP
C     
         IRREPTI = IRREPI
         IRREPTE = DIRPRD(IRREPTI, IRREPX)
         IRREPWE = IRREPTE
         IRREPWM = DIRPRD(IRREPWE, IRREPWL)
C
         NVRTWE = VRT2(IRREPWE)
         NOCCTJ = POP2(IRREPTI) 
         NOCCWM = POP1(IRREPWM)
C
         IOFFW = IWOFF(IRREPWE)
         IOFFT = ITOFFA(IRREPTE)
C
         IF(NOCCTJ .NE. 0 .AND. NOCCWM .NE. 0 .AND. NVRTWE .NE. 0) THEN
C
            CALL XGEMM('N', 'N', NUMSYW*NOCCWM, NOCCTJ, NVRTWE, ONE, 
     &                  WT(1,IOFFW), NUMSYW*NOCCWM, T1ALPA(IOFFT),
     &                  NVRTWE, AZERO, Z(1,IOFFZ), NUMSYW*NOCCWM)
       ENDIF
C
       IOFFZ = IOFFZ + NOCCWM*NOCCTJ
C
 100  CONTINUE
C
C Transpose again last two indicies. Now we need to know the 
C symmetry of the right hand side (MI). We know the overall
C symmetry of the intermediate product which is IRREPX and
C the left hand symmetry. Standard CC code these difficulties
C do not occur.
C
      IRREPJM = DIRPRD(IRREPWR, IRREPX)
      CALL SYMTRA(IRREPJM, POP1, POP2, NUMSYW, Z, WT)
C
C Zero traget array
C
      CALL ZERO(Z, DISTMP*NUMTMP)
C
C Perform second multiplication 
C
      IOFFW = 1
      IZOFF(1) = 1
      IOFF = 1
C
C Set up a offset array for the taget
C     
      DO 5200 IRREPBO = 2, NIRREP
C
         IRREP = IRREPBO - 1
         IRREPI = DIRPRD(IRREP, IRREPQR)
C
         IZOFF(IRREPBO) = IZOFF(IRREPBO - 1) + POP2(IRREPI)*VRT1(IRREP)
C
 5200 CONTINUE
C
      DO 5300 IRREPM = 1, NIRREP
C
         IRREPBO = DIRPRD(IRREPM, IRREPX)
C
         ITOFFB(IRREPBO) = IOFF
         IOFF =  IOFF + POP1(IRREPM)*VRT1(IRREPBO)         
C
 5300 CONTINUE
C
      DO 200 IRREPM = 1, NIRREP
C
         IRREPZM = IRREPM
         IRREPZJ = DIRPRD(IRREPZM, IRREPJM)
         IRREPTM = IRREPZM 
         IRREPTB = DIRPRD(IRREPTM, IRREPX)
C
         NVRTA = VRT1(IRREPTB)
         NOCCM = POP1(IRREPTM)
         NOCCJ = POP2(IRREPZJ)
C
         IOFFZ = IZOFF(IRREPTB)
         IOFFT = ITOFFB(IRREPTB)
C
         IF(NOCCJ .NE. 0 .AND. NOCCM .NE. 0 .AND. NVRTA .NE. 0) THEN
C
            CALL XGEMM('N', 'T', NUMSYW*NOCCJ, NVRTA, NOCCM, FACT, 
     &                  WT(1,IOFFW), NUMSYW*NOCCJ, T1BETA(IOFFT),
     &                  NVRTA, AZERO, Z(1,IOFFZ), NUMSYW*NOCCJ)
         ENDIF
C
         IOFFW = IOFFW + NOCCJ*NOCCM
C
 200  CONTINUE
C
C Transpose again the last two indicies and save results on list
C
      IF (IFLAGS(1) .GE. 20) THEN

         NSIZE = DISSYZ*NUMSYZ

         CALL HEADER('Checksum @-T1T1IND1 per sym. block', 0, LUOUT)
     
         WRITE(LUOUT, *) ' SPIN CASE =  ', SDOT(NSIZE, Z, 1, Z, 1)
      ENDIF
C
      CALL SYMTRA(IRREPQR, POP2, VRT1, NUMSYW, Z, W)
C
C First case, no transpose of results are required 
C
      IF (ITRANS .EQ. 0) THEN
C
         IF (INCREM .EQ. 1) THEN
C Update  
            CALL GETLST(Z, 1, NUMSYZ, 1, IRREPQR, LISTZ)
            CALL VADD(Z, Z, W, NUMSYZ*DISSYZ, ONE)
            CALL PUTLST(Z, 1, NUMSYZ, 1, IRREPQR, LISTZ)
C
         ELSE
C     
C Initialize the list
C
            CALL PUTLST(W, 1, NUMSYZ, 1, IRREPQR, LISTZ)
C
         ENDIF
C
      ELSE
C
C Transpose the results first 
C
         CALL TRANSP(W, Z, DISSYZ, NUMSYZ)
C
         IF (INCREM .EQ. 1) THEN
C Update 
            CALL GETLST(W, 1, NUMSYZ, 1, IRREPQR, LISTZ)
            CALL VADD(W, W, Z, NUMSYZ*DISSYZ, ONE)
            CALL PUTLST(W, 1, NUMSYZ, 1, IRREPQR, LISTZ)
C
         ELSE
C
C Initilize the list
C
            CALL PUTLST(Z, 1, NUMSYZ, 1, IRREPQR, LISTZ)
C
         ENDIF
C
      ENDIF
C
      RETURN
      END
