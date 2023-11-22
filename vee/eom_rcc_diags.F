#include "flags.h"
      SUBROUTINE EOM_RCC_DIAGS(W,MAXCOR,IUHF,IMULT)
C
C THIS ROUTINE DRIVES THE drCCD Hbar diagonalizations 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSIZ
      LOGICAL MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
      LOGICAL METH, CIS,RPA,EOMCC,CISD,FULDIAG,INCORE,READGUES
      LOGICAL EOM_S_RCCD,EOM_S_DRCCD,EOM_SF_RCCD,EOM_SF_DRCCD
      LOGICAL EOM_S_DXRCCD,EOM_SF_DXRCCD
      LOGICAL ADD_DBLS_CORR
      DOUBLE PRECISION W(MAXCOR)

      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREP0(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/REFTYPE/MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
      COMMON/METH/CIS,RPA,EOMCC,CISD,FULDIAG,INCORE,READGUES
#include "flags.com"

      DATA TWO /2.0D0/
C
      IF (CIS) THEN
         LISTDIAG = 42
      ELSE
         LISTDIAG = 56
      ENDIF 
      EOM_SF_RCCD   = (IFLAGS(h_IFLAGS_excite) .EQ.17) 
      EOM_SF_DRCCD  = (IFLAGS(h_IFLAGS_excite) .EQ.18)
      EOM_S_RCCD    = (IFLAGS(h_IFLAGS_excite) .EQ.15) 
      EOM_S_DRCCD   = (IFLAGS(h_IFLAGS_excite) .EQ.13)
      EOM_S_DXRCCD  = (IFLAGS(h_IFLAGS_excite) .EQ.19)
      EOM_SF_DXRCCD = (IFLAGS(h_IFLAGS_excite) .EQ.20)
         
      Write(6,*)
      IF (CIS) THEN
         Write(6,"(a)") " CIS H(ai,bj) is fully diagonalized" 
      ELSEIF (RPA) THEN
         IF (EOM_SF_RCCD) THEN
            Write(6,"(a,a)") " EOM(SF)-RCCD (RPA) Heff(ai,bj) is",
     +                       " fully diagonalized"    
         ELSEIF(EOM_SF_DRCCD) THEN
            Write(6,"(a,a)") " EOM(SF)-RCCD (DRPA) Heff(ai,bj) is",
     +                       " fully diagonalized"    
         ELSEIF (EOM_S_RCCD) THEN
            Write(6,"(a,a)") " EOM(S)-RCCD Heff(ai,bj) is fully",
     +                       " diagonalized" 
         ELSEIF (EOM_S_DRCCD) THEN
            Write(6,"(a,a)") " EOM(S)-DRCCD Heff(ai,bj) is fully ",
     +                    " diagonalized" 
         ELSEIF (EOM_S_DXRCCD) THEN
            Write(6,"(a,a)") " EOM(S)-DXRCCD Heff(ai,bj) is fully ",
     +                       " diagonalized"
         ELSEIF (EOM_SF_DXRCCD) THEN
            Write(6,"(a,a)") " EOM(SF)-DXRCCD Heff(ai,bj) is fully ",
     +                       " diagonalized"
         ENDIF 
      ENDIF 

C LOOP OVER SYMMETRY BLOCKS, READ THE MATRIX IN AND DIAGONALIZE IT
C
      JROOT = 0
      DO IRREP=1,NIRREP

C RHF block

         IF (IUHF.EQ.0)THEN
            NUMDIS=IRPDPD(IRREP,ISYTYP(2,56))
            MATDIM = NUMDIS
            I000=1
            I010=I000+NUMDIS*NUMDIS
            CALL GETLST(W(I000),1,NUMDIS,1,IRREP,LISTDIAG)
         ELSE

C UHF AA and BB blocks 

            NUMAA=IRPDPD(IRREP,ISYTYP(1,19))
            NUMBB=IRPDPD(IRREP,ISYTYP(1,20))
            MATDIM=NUMAA+NUMBB

            I000=1
            I010=I000+MATDIM*MATDIM

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
        Write(6,"(a)") " The drCCD Hbar matrix"
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
          IF (IEND .GT. MAXCOR) CALL INSMEM("@--eom_rcc_diags",IEND,
     +                                       MAXCOR)

          CALL DGEEV("V","V",MATDIM,W(I000),MATDIM,W(I010),W(I020),
     +                W(I030),MATDIM,W(I040),MATDIM,W(I050),
     +                4*MATDIM,IERROR)
 
          CALL SORT_RCC_EIGS(W(I010),W(I040),W(I030),MATDIM)

       Write(6,*) 
       Write(6,"(a)") " The sorted real eigenvalues (eV)"
       Bohrtoev=27.2113961D0
       Write(6,"(5(1x,F15.8))") (W(I010+I-1)*Bohrtoev,i=1,MATDIM)
       Write(6,*) 
       Write(6,"(a)") "The sorted imaginary eigenvalues (eV)"
       Write(6,"(5(1x,F15.8))") (W(I020+I-1)*Bohrtoev,i=1,MATDIM)

#ifdef _DEBUG_LVLM
       Write(6,"(a)") "The sorted left eigenvectors"
       call output(W(I030),1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
       Write(6,"(a)") "The sorted right eigenvectors"
       call output(W(I040),1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
#endif 

       I100 = 1
       I110 = I100 + MATDIM
       I120 = I110 + MATDIM*MATDIM
       I130 = I120 + MATDIM*MATDIM
       I140 = I130 + MATDIM*8
       IEND = I140
       IF (IEND .GT. MAXCOR) CALL INSMEM("@-eom_rcc_diags",IEND,MAXCOR)
       MEMLEFT = MAXCOR - IEND 

       CALL DCOPY(MATDIM,W(I010),1,W(I100),1)
       CALL DCOPY(MATDIM*MATDIM,W(I030),1,W(I110),1)
       CALL DCOPY(MATDIM*MATDIM,W(I040),1,W(I120),1)

       CALL DZERO(W(I130),MATDIM*8)
       ADD_DBLS_CORR = .TRUE.

       IF (ADD_DBLS_CORR) THEN
       CALL NEWLST(IRREP,W(IEND),LENGTH,IUHF)
       CALL FORM_DBLS_CORRECTNS_2RCC(W(I100),W(I110),W(I120),
     +                               W(I130),W(IEND),
     +                               MEMLEFT,IUHF,MATDIM,IRREP,
     +                               IMULT,EOM_SF_RCCD,EOM_SF_DRCCD)
       ENDIF 
       
       IF (IUHF .EQ. 0) CALL PROCESS_RHF_CCD_ROOTS(W(I120),W(I110),
     +                       W(I100),W(I130),MATDIM,IRREP,1,RPA,
     +                       EOM_SF_RCCD,EOM_SF_DRCCD,
     +                       EOM_S_RCCD,EOM_S_DRCCD,IMULT,JROOT)

       IF (IUHF .EQ. 1) CALL PROCESS_UHF_CCD_ROOTS(W(I120),W(I110),
     +                       W(I100),W(I130),MATDIM,NUMAA,NUMBB,IRREP,
     +                       IUHF,RPA,EOM_SF_RCCD,EOM_SF_DRCCD,
     +                       EOM_S_RCCD,EOM_S_DRCCD,JROOT)

      ENDDO 

      RETURN 
      END
