C
      SUBROUTINE UNWIND(ICORE, MAXCOR, IRREPX, IUHF, ILIST, NUMPERT, 
     &                   ICOUNT, BETAFLAG)
C
C This routines read a list (ILIST) which has both T1 and T2 amplitudes
C ordered in the following manner and split it in to T1 and T2 two separate 
C lists. Also, the lists are rewritten to be consistent with the routines
C used in quadratic term.
C
C  T1(AA) - T1(BB) - T2(AB) - T2(BB) - T2(AA) [UHF] 
C  T1[AA] - T2[AB]                            [RHF]
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION ICORE(MAXCOR)
      INTEGER DISSYT2, DISSYAB, DISSYAA, POP, VRT, DIRPRD
      LOGICAL BETAFLAG
C
      COMMON/MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/FILES/LUOUT, MOINTS
      COMMON/FLAGS/IFLAGS(100)
C   
C Common blocks used in the quadratic term
C 
      COMMON /QADINTI/ INTIMI, INTIAE, INTIME
      COMMON /QADINTG/ INGMNAA, INGMNBB, INGMNAB, INGABAA,
     &                 INGABBB, INGABAB, INGMCAA, INGMCBB,
     &                 INGMCABAB, INGMCBABA, INGMCBAAB,
     &                 INGMCABBA
      COMMON /APERTT1/ IAPRT1AA             
      COMMON /BPERTT1/ IBPRT1AA
      COMMON /APERTT2/ IAPRT2AA1, IAPRT2BB1, IAPRT2AB1, IAPRT2AB2, 
     &                 IAPRT2AB3, IAPRT2AB4, IAPRT2AA2, IAPRT2BB2,
     &                 IAPRT2AB5 
C
      LEN  = 0
      I000 = 1
      IPRTLVL = IFLAGS(1)
C
C Compute the length of the full vector.
C
      DO 10 ISPIN = 1, IUHF + 1
         LEN = LEN + IRPDPD(IRREPX, 8 + ISPIN)
 10   CONTINUE  
C
      DO 20 ISPIN = 3 , 3 - 2*IUHF, -1
         LSTTYP = 43 + ISPIN
         LEN = LEN + IDSYMSZ(IRREPX, ISYTYP(1, LSTTYP), 
     &                       ISYTYP(2, LSTTYP))
 20   CONTINUE
C
C Allocate space for the full vector.
C
      I010 = I000 + LEN*IINTFP
C      
C Load the full list into the core memory.
C
      CALL GETLST(ICORE(I000), ICOUNT, 1, 1, IRREPX, ILIST)
C
C Print out the perturb T amplitudes if this is a debugging process.
C
      IF (IPRTLVL .GE. 40) THEN
         CALL PRTEXTLST(ICORE(I000), IRREPX, LEN, IUHF, .TRUE., LUOUT,
     &                  DIRPRD, IRPDPD, NIRREP)
      ENDIF
C
C Compute the length of individual T1 and T2 vectors and
C also allocate the space for individual T1 and T2 vectors.
C
      IF (IUHF .NE. 0) THEN
C
         LENT1AA = IRPDPD(IRREPX, 9)
         LENT1BB = IRPDPD(IRREPX, 10)
         LENT2AA = IDSYMSZ(IRREPX, ISYTYP(1, 44), ISYTYP(2, 44))
         LENT2BB = IDSYMSZ(IRREPX, ISYTYP(1, 45), ISYTYP(2, 45))
         LENT2AB = IDSYMSZ(IRREPX, ISYTYP(1, 46), ISYTYP(2, 46))
C
      ELSE
C
         LENT1AA = IRPDPD(IRREPX, 9)
         LENT2AB = IDSYMSZ(IRREPX, ISYTYP(1, 46), ISYTYP(2, 46))
C
      ENDIF         
C     
C Create the new lists by writing out the appropriate portions
C of the SCR array. First do the T1(ALPHA) and T1(BETA) vectors.
C When BETAFLAG = .TRUE. we write the T1(BETA) vectors in to the
C appropriate list, unless we only write T1(ALPHA) vectors in to the
C disk.
C
      DO 30 ISPIN = 1, IUHF + 1
C         
         IF (BETAFLAG) THEN
C
            IOFF = I000 + (ISPIN - 1)*LENT1AA*IINTFP   
            CALL PUTLST (ICORE(IOFF), 1, 1, 1, ISPIN, IBPRT1AA)
C     
         ELSE
C
            IOFF = I000 + (ISPIN - 1)*LENT1AA*IINTFP
            CALL PUTLST (ICORE(IOFF), 1, 1, 1, ISPIN, IAPRT1AA)
         ENDIF
C
 30   CONTINUE
C
C T2(ALPHA) and T2(BETA) are ordered as (A<B,I<J), (a<b,i<j) and (Ab,Ij).
C List IAPRT2AA2, IAPRT2BB2 and IAPRT2AB5 expects perturb T2 amplitudes 
C as ordered above. Before switch ordering write these vectors to the 
C appropriate lists. For T2 lists we overwrite T2(ALPHA) by T2(BETA)
C so we do not need to use the BETAFLAG.
C
C First do the ABAB spin case (Both UHF and RHF case)
C
      IF (IUHF .NE.0 ) THEN
         LENT1 = (LENT1AA + LENT1BB)*IINTFP 
      ELSE
         LENT1 = LENT1AA*IINTFP
      ENDIF
C
      LIST2 = IAPRT2AB5
      IOFFT2 =  LENT1 + 1
C     
      DO 60 IRREP = 1, NIRREP
C
         IRREPTIJ = IRREP
         IRREPTAB = DIRPRD(IRREPTIJ, IRREPX)
C     
         DISSYT2 = IRPDPD(IRREPTAB, ISYTYP(1, LIST2))
         NUMSYT2 = IRPDPD(IRREPTIJ, ISYTYP(2, LIST2))
C
         CALL PUTLST(ICORE(IOFFT2), 1, NUMSYT2, 1, IRREPTIJ, LIST2)
C     
         IOFFT2 = IOFFT2 + DISSYT2*NUMSYT2*IINTFP
C     
 60   CONTINUE
C
      LENT2 =  LENT1 + LENT2AB*IINTFP
C     
C Now do the AAAA and BBBB spin cases. (Only for UHF)
C
      IF (IUHF .NE. 0) THEN
C
         DO 40 ISPIN = 2, IUHF, -1
C         
            LIST2  = (IAPRT2AA2 - 1) + ISPIN 
            IOFFT2 =  LENT2 - (ISPIN - 2)*LENT2BB*IINTFP + 1
C     
            DO 50 IRREP = 1, NIRREP
C     
               IRREPTIJ = IRREP
               IRREPTAB = DIRPRD(IRREPTIJ, IRREPX)
C     
               DISSYT2 = IRPDPD(IRREPTAB, ISYTYP(1, LIST2))
               NUMSYT2 = IRPDPD(IRREPTIJ, ISYTYP(2, LIST2))
C     
               CALL PUTLST(ICORE(IOFFT2), 1, NUMSYT2, 1, IRREPTIJ,
     &                     LIST2)
C     
               IOFFT2 = IOFFT2 + DISSYT2*NUMSYT2*IINTFP
C     
 50         CONTINUE
C   
 40      CONTINUE
C
      ENDIF
C      
C Now we have T2(ALPHA), T2(BETA) written in to the list IAPRT2AA2 
C - IAPRT2AB5. These lists are ordered (A<B,I<J), (a<b,i<j) and (Ab,Ij). 
C To create T2 lists IAPRT2AA1 - IAPRT2AB4 we need reorder IAPRT2AA2 
C - IAPRT2AB5 list. Now we can overwrite on the original T1 and T2
C vectors, in the following calls we can use the whole available core
C memory.
C
C For RHF runs we need to generate T2 vectors for AAAA spin case
C from the T2 vectors of ABAB spin case
C
      IF (IUHF .EQ. 0) THEN
C
         DO 70 IRREP = 1, NIRREP
C     
            IRREPTIJ = IRREP
            IRREPTAB = DIRPRD(IRREPTIJ, IRREPX)
C
            DISSYAB = IRPDPD(IRREPTAB, ISYTYP(1, LIST2))
            NUMSYAB = IRPDPD(IRREPTIJ, ISYTYP(2, LIST2))
            DISSYAA = IRPDPD(IRREPTAB, ISYTYP(1, LIST2 - 2))
            NUMSYAA = IRPDPD(IRREPTIJ, ISYTYP(2, LIST2 - 2))
C
            I000 = 1
            I010 = I000 + MAX(DISSYAB, DISSYAA)*MAX(NUMSYAA, NUMSYAB)
     &             *IINTFP
            I020 = I010 + MAX(DISSYAB, DISSYAA)*MAX(NUMSYAA, NUMSYAB)
     &             *IINTFP
            I030 = I020 + MAX(NUMSYAB, NUMSYAA)*IINTFP
            I040 = I030 + MAX(NUMSYAB, NUMSYAA)*IINTFP
C
            CALL GETLST(ICORE(I000), 1, NUMSYAB, 1, IRREPTIJ, LIST2)
C
            CALL ASSYM2A(IRREPTAB, VRT(1, 1), DISSYAB, NUMSYAB,
     &                   ICORE(I000), ICORE(I020), ICORE(I030))
C
            CALL SQSYM(IRREPTAB, VRT(1, 1), DISSYAA, DISSYAB, NUMSYAB,
     &                 ICORE(I010), ICORE(I000))
C
            CALL TRANSP(ICORE(I010), ICORE(I000), NUMSYAB, DISSYAA)
            CALL SQSYM(IRREPTIJ, POP(1, 1), NUMSYAA, NUMSYAB, DISSYAA,
     &                 ICORE(I010), ICORE(I000))
            CALL TRANSP(ICORE(I010), ICORE(I000), DISSYAA, NUMSYAA)
C     
            CALL PUTLST(ICORE(I000), 1, NUMSYAA, 1, IRREPTIJ, LIST2-2)
C
 70      CONTINUE
C
      ENDIF
C
C If it is a debugging process print out the contents in lists
C just created
C    
      IF (IPRTLVL .GE. 40) THEN
         DO 5 ISPIN = 1, IUHF + 1
            LIST2 = (IAPRT2AA2 - 1) + ISPIN
            IF (BETAFLAG) THEN
               CALL PRTLIST(ICORE(I000), MAXCOR, ISPIN, IBPRT1AA, 
     &                      IRREPX, LUOUT, .TRUE., DIRPRD, IRPDPD,
     &                      ISYTYP, NIRREP)
            ELSE
               CALL PRTLIST(ICORE(I000), MAXCOR, ISPIN, IAPRT1AA, 
     &                      IRREPX, LUOUT, .TRUE., DIRPRD, IRPDPD,
     &                      ISYTYP, NIRREP)
            ENDIF
C     
            CALL PRTLIST(ICORE(I000), MAXCOR, ISPIN, LIST2, IRREPX,
     &                  LUOUT, .FALSE., DIRPRD, IRPDPD, ISYTYP, NIRREP)
 5       CONTINUE
C
         CALL PRTLIST(ICORE(I000), MAXCOR, 3, IAPRT2AB5, IRREPX, 
     &                LUOUT, .FALSE., DIRPRD, IRPDPD, ISYTYP, NIRREP)
      ENDIF
            
C
C Generate T2 lists for the list numbers IAPRT2AA1 - IAPRT2AB4.
C     
      CALL DBLRESORT(ICORE(I000), MAXCOR, IRREPX, IUHF, 'T')
C
C This will finish the resorting perturb T2 vectors. Now we have
C all the appropriate lists to proceed.
C
      RETURN
      END
