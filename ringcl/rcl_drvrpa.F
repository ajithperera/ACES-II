      SUBROUTINE RCL_DRVRPA(W,MAXCOR,IUHF)
C
C THIS ROUTINE DRIVES THE drCCD Hbar diagonalizations 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSIZ
      LOGICAL CIS,EOM,RPA,EOM_SDRCCD,EOM_SRCCD,EOM_SFDRCCD,
     +        EOM_SFRCCD,EOM_SDXRCCD,EOM_SFDXRCCD
      DOUBLE PRECISION W(MAXCOR)

      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREP0(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/EXCITE/EOM,CIS,RPA,EOM_SRCCD,EOM_SDRCCD,EOM_SFRCCD,
     &              EOM_SFDRCCD,EOM_SDXRCCD,EOM_SFDXRCCD

      DATA TWO /2.0D0/
C
      IF (CIS) THEN
         LISTDIAG = 42
      ELSE
         LISTDIAG = 56
      ENDIF 
         
      Write(6,*)
      IF (CIS) THEN
         Write(6,"(a)") " CIS H(ai,bj) is fully diagonalized" 
      ELSEIF (EOM_SFRCCD) THEN
         Write(6,"(a,a)") " EOM(SF)-RCCD (RPA) Heff(ai,bj) is fully",
     +                  " diagonalized" 
      ELSEIF (EOM_SFDRCCD) THEN
         Write(6,"(a,a)") " EOM(SF)-DRCCD (DRPA) Heff(ai,bj) is fully",
     +                  " diagonalized" 
      ELSEIF (EOM_SDRCCD) THEN
         Write(6,"(a,a)") " EOM(S)-DRCCD Heff(ai,bj) is fully",
     +                    " diagonalized" 
      ELSEIF (EOM_SRCCD) THEN
         Write(6,"(a,a)") " EOM(S)-RCCD Heff(ai,bj) is fully",
     +                    " diagonalized" 
      ENDIF 

C LOOP OVER SYMMETRY BLOCKS, READ THE MATRIX IN AND DIAGONALIZE IT
C
      DO IRREP=1,NIRREP

C RHF block

         IF (IUHF.EQ.0)THEN
            NUMDIS=IRPDPD(IRREP,ISYTYP(2,56))
            MATDIM = NUMDIS
            I000=1
            I010=I000+NUMDIS*NUMDIS
            IEND=I010
            IF (IEND .GT. MAXCOR) CALL INSMEM("RCL_DRVRPA",IEND,
     +                                         MAXCOR)
            CALL GETLST(W(I000),1,NUMDIS,1,IRREP,LISTDIAG)
         ELSE

C UHF AA and BB blocks 

            NUMAA=IRPDPD(IRREP,ISYTYP(1,19))
            NUMBB=IRPDPD(IRREP,ISYTYP(1,20))
            MATDIM=NUMAA+NUMBB

            I000=1
            I010=I000+MATDIM*MATDIM
            IEND=I010
            IF (IEND .GT. MAXCOR) CALL INSMEM("RCL_DRVRPA",IEND,
     +                                         MAXCOR)
            IOFFAAAA=1
            IOFFBBBB=IOFFAAAA+MATDIM*NUMAA+NUMAA
            IOFFBBAA=IOFFAAAA+NUMAA
            IOFFAABB=IOFFAAAA+MATDIM*NUMAA

            IF (CIS) THEN
               LISTW=23
            ELSE
               LISTW=54
            ENDIF 
         
            NUMDIS=IRPDPD(IRREP,ISYTYP(1,LISTW))        
            IOFF=IOFFAAAA
            DO 11 IDIS=1,NUMDIS
               CALL GETLST(W(IOFF),IDIS,1,1,IRREP,LISTW)
               IOFF=IOFF+MATDIM
11          CONTINUE

            IF (CIS) THEN
               LISTW=24
            ELSE
               LISTW=55
            ENDIF 

            NUMDIS=IRPDPD(IRREP,ISYTYP(1,LISTW))        
            IOFF=IOFFBBBB
            DO 12 IDIS=1,NUMDIS
               CALL GETLST(W(IOFF),IDIS,1,1,IRREP,LISTW)
               IOFF=IOFF+MATDIM
12          CONTINUE

            IF (CIS) THEN
               LISTW1=17
               LISTW2=17
            ELSE 
               LISTW1=56
               LISTW2=57
            ENDIF 

            DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))        
            NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))        
            IOFF1=IOFFBBAA
            IOFF2=IOFFAABB

            DO 13 IDIS=1,NUMDIS
               CALL GETLST(W(IOFF1),IDIS,1,1,IRREP,LISTW1)
               CALL GETLST(W(IOFF2),IDIS,1,1,IRREP,LISTW2)
               IOFF1=IOFF1+MATDIM
               IOFF2=IOFF2+MATDIM 
13          CONTINUE

            NUMDIS=MATDIM
         ENDIF
C       
#ifdef _DEBUG_LVLM
        Write(6,"(a)") " The Hbar matrix for the requested model"
        call output(W,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
#endif 
C Do a full non-symmetric diagonalization

         I000 = 1
         I010 = I000 + MATDIM * MATDIM
         I020 = I010 + MATDIM 
         I030 = I020 + MATDIM
         I040 = I030 + MATDIM * MATDIM
         I050 = I040 + MATDIM * MATDIM
         IEND = I050
         IF (IEND .GT. MAXCOR) CALL INSMEM("RCL_DRVRPA",IEND,
     +                                       MAXCOR)

          CALL DGEEV("V","V",MATDIM,W(I000),MATDIM,W(I010),W(I020),
     +                W(I030),MATDIM,W(I040),MATDIM,W(I050),
     +                4*MATDIM,IERROR)
#ifdef _DEBUG_LVL0  
       Write(6,*) 
       Write(6,"(a)") " The real eigenvalues"
       Bohrtoev=27.2113961D0
CSSS       Bohrtoev=1.0D0

       Write(6,"(5(1x,F15.8))") (W(I010+I-1)*Bohrtoev,i=1,MATDIM)
       Write(6,"(a)") "The imaginary eigenvalues"
       Write(6,"(5(1x,F15.8))") (W(I020+I-1)*Bohrtoev,i=1,MATDIM)
#endif 
      ENDDO 

      RETURN 
      END
