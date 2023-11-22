#include "flags.h"
      SUBROUTINE SCRNC_FORM_RCC_HBAR(ICORE,MAXCOR,IUHF,IMULT)
C
C THIS ROUTINE FORMS THE SINGLE-SINGLE BLOCK OF EITHER THE REGULAR
C OR THE CCSD EFFECTIVE HAMILTONIAN
C
C       H(ai,bj) = F(ab) d(ij) - F(ij) d(ab) - W(ai,bj)
C
C WHERE F AND W ARE THE ONE AND TWO PARTICLE PARTS OF EXP(-T) H EXP(T).
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL CIS,RPA,EOMCC,CISD,FULDIAG,INCORE,READGUES
      LOGICAL MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
      LOGICAL NONHF 
      LOGICAL HEFF,TDA,EOM_SFRCCD,EOM_SFDRCC
      DOUBLE PRECISION ONE,ONEM,TWO,HALF
      DIMENSION ICORE(MAXCOR),IOFFO(2),IOFFV(2)
C
      COMMON/FLAGS/IFLAGS(100)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREP0(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/INFO/NOCCO(2),NVRTO(2)
C 
      DATA ONE,ONEM /1.0D0,-1.0D0/
      DATA TWO,HALF /2.0D0,0.50D0/
#ifdef _DEBUG_LVL0
      Write(6,"(1x,a)") "-----Entering scrnc_form_rcc_hbar------"
#endif

      HEFF  = .FALSE.
      TDA   = .FALSE.
      CIS   = .FALSE.
      RPA   = .TRUE.
      NONHF =  (IFLAGS(h_IFLAGS_nonhf)  .GT. 0)

C Occ-Occ block

      I000=1
      I010=I000+NFMI(1)*IINTFP
      IF(IUHF.NE.0)THEN
       I020=I010+NFMI(2)*IINTFP
      ELSE
       I020=I010
      ENDIF

C Vrt-Vrt  block

      I030=I020+NFEA(1)*IINTFP
      IF(IUHF.NE.0)THEN
       I040=I030+NFEA(2)*IINTFP
      ELSE
       I040=I030
      ENDIF

      IF (CIS) THEN
C This simply do CIS and is usefull for debugging.
          HEFF = .FALSE.
          TDA  = .TRUE. 

      ELSEIF(RPA) THEN

C This is drRPA or RPA, but using Hbar from drCCD or rCCD with HF diagonals.
C These methods are designated as EOM(Sf)-RCCD or  EOM(Sf)-DRCCD

          HEFF = .FALSE.
          TDA  = .FALSE. 
      ENDIF 

      IF (HEFF) THEN

C If Hbar(i,j) and Hbar(a,b) must be used; load them here. 
         
          CALL GETLST(ICORE(I000),1, 1, 1, 1, 91)
          CALL GETLST(ICORE(I020),1, 1, 1, 1, 92)
          IF (IUHF .EQ. 1) THEN
            CALL GETLST(ICORE(I010),1, 1, 1, 2, 91)
            CALL GETLST(ICORE(I030),1, 1, 1, 2, 92)
          ENDIF
      ELSE 

C For NON-HF we need to F(a,b) and F(i,j) elements. For EOM(S)-RCC
C or DRCC methods these have already been added when we cnstruct
C the Hbar(a,b) and Hbar(i,j). But for EOM(SF)-RCC or -DRCC we need to
C add them here.

          IF (NONHF) THEN
             IF (RPA .OR. EOM_SFRCCD .OR. EOM_SFDRCCD) THEN
                CALL GETLST(ICORE(I000),1, 1, 1, 3, 91)
                CALL GETLST(ICORE(I020),1, 1, 1, 3, 92)
                IF (IUHF .EQ. 1) THEN
                  CALL GETLST(ICORE(I010),1, 1, 1, 4, 91)
                  CALL GETLST(ICORE(I030),1, 1, 1, 4, 92)
                ENDIF
             ELSE
                CALL IZERO(ICORE(I000),I040-1)
             ENDIF 
          ELSE
             CALL IZERO(ICORE(I000),I040-1)
          ENDIF 
       ENDIF

C The Diagonal elements. First alpha occupied and followed by virtuals 

       NBAS=NOCCO(1)+NVRTO(1)
       CALL GETREC(20,'JOBARC','SCFEVALA',IINTFP*NBAS,ICORE(I040))

       IOFFE =I040
       IOFFTO=I000
       IOFFTV=I020

       DO 10 IRREP=1,NIRREP
          NOCC=POP(IRREP,1)
          CALL SAXPY(NOCC,ONE,ICORE(IOFFE),1,ICORE(IOFFTO),NOCC+1)
          IOFFE =IOFFE +NOCC*IINTFP
          IOFFTO=IOFFTO+NOCC*NOCC*IINTFP
10     CONTINUE

       DO 11 IRREP=1,NIRREP
          NVRT=VRT(IRREP,1)
          CALL SAXPY(NVRT,ONE,ICORE(IOFFE),1,ICORE(IOFFTV),NVRT+1)
          IOFFE =IOFFE +NVRT*IINTFP
          IOFFTV=IOFFTV+NVRT*NVRT*IINTFP
11     CONTINUE
C
       IF(IUHF.EQ.1)THEN
         CALL GETREC(20,'JOBARC','SCFEVALB',IINTFP*NBAS,ICORE(I040))
         IOFFE =I040
         IOFFTO=I010
         IOFFTV=I030

        DO 12 IRREP=1,NIRREP
           NOCC=POP(IRREP,2)
           CALL SAXPY(NOCC,ONE,ICORE(IOFFE),1,ICORE(IOFFTO),NOCC+1)
           IOFFE =IOFFE+NOCC*IINTFP
           IOFFTO=IOFFTO+NOCC*NOCC*IINTFP
12      CONTINUE

        DO 13 IRREP=1,NIRREP
           NVRT=VRT(IRREP,2)
           CALL SAXPY(NVRT,ONE,ICORE(IOFFE),1,ICORE(IOFFTV),NVRT+1)
           IOFFE =IOFFE+NVRT*IINTFP
           IOFFTV=IOFFTV+NVRT*NVRT*IINTFP
13        CONTINUE
       ENDIF        
C
      IOFFO(1)=I000
      IOFFO(2)=I010
      IOFFV(1)=I020
      IOFFV(2)=I030

C Read in Hbar(ib,aj) elements; ordered as (bi,aj)

      DO ISPIN=1,1+IUHF
         IF(IUHF.EQ.0)THEN
           ISIZE=ISYMSZ(ISYTYP(1,56),ISYTYP(1,56))
           I050=I040+ISIZE*IINTFP
           I060=I050+ISIZE*IINTFP
           IEND=I060
           IF (IEND .GT. MAXCOR) CALL INSMEM("RCL_HEFFXC",
     &                                       IEND,MAXCOR)
           IF (TDA) THEN
              CALL GETALL(ICORE(I040),ISIZE,1,23)
              CALL VMINUS(ICORE(I040),ISIZE)
              CALL GETALL(ICORE(I050),ISIZE,1,18)
              CALL SAXPY (ISIZE,ONE,ICORE(I050),1,ICORE(I040),1)
           ELSE 
              IF (IMULT .EQ. 1) THEN
                 CALL GETALL(ICORE(I040),ISIZE,1,56)
              ELSE IF (IMULT .EQ. 2) THEN
                 CALL GETALL(ICORE(I040),ISIZE,1,54)
              ENDIF 
           ENDIF 
         ELSE
           IF (TDA) THEN
               LISTAA=22+ISPIN
           ELSE 
               LISTAA=53+ISPIN
           ENDIF 
           ISIZE=ISYMSZ(ISYTYP(1,LISTAA),ISYTYP(1,LISTAA))
           I050=I040+ISIZE*IINTFP
           I060=I050+ISIZE*IINTFP
           IEND=I060 
           IF (IEND .GT. MAXCOR) CALL INSMEM("RCL_HEFFXC",
     &                                        IEND,MAXCOR)
           CALL GETALL(ICORE(I040),ISIZE,1,LISTAA)
           IF (TDA) CALL VMINUS(ICORE(I040),ISIZE)
          enDIF

C Reorder; (ai,bj) -> (ab,ij)

          CALL SSTGEN(ICORE(I040),ICORE(I050),ISIZE,
     &                VRT(1,ISPIN),POP(1,ISPIN),VRT(1,ISPIN),
     &                POP(1,ISPIN),ICORE(I060),1,'1324')

C Add the one particle contributions 

          DISSZ1=IRPDPD(1,18+ISPIN)
          NUMDS1=IRPDPD(1,20+ISPIN)
          CALL SCRNC_ADDONEH(ICORE(I050),ICORE(I040),ISIZE,
     &                       ICORE(IOFFV(ISPIN)),
     &                       ICORE(IOFFO(ISPIN)),
     &                       DISSZ1,NUMDS1,ISPIN)

C Reorder (ab,ij) -> (ai,bj)

          CALL SSTGEN(ICORE(I050),ICORE(I040),ISIZE,
     &                VRT(1,ISPIN),VRT(1,ISPIN),POP(1,ISPIN),
     &                POP(1,ISPIN),ICORE(I060),1,'1324')
     
          IF (TDA) THEN
            LISTDUMP = 42
          ELSE
            LISTDUMP = 56
          ENDIF 
          IF (IUHF .NE. 0) LISTDUMP = LISTAA 
          CALL PUTALL(ICORE(I040),ISIZE,1,LISTDUMP) 

      ENDDO 

      RETURN
      END
