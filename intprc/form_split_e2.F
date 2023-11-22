#include "flags.h"
      SUBROUTINE FORM_SPLIT_E2(W,MAXCOR,IUHF)
C
C THIS ROUTINE FORMS THE FIRST ORDER T2 VECTOR IN SYMMETRY 
C  PACKED FORM AND WRITES IT OUT.  ALSO COMPUTE THE CONTRIBUTION
C  TO THE CORRELATION ENERGY.
C
C      ARRAYS : W - USED TO HOLD INTEGRAL ARRAYS.  MUST
C                   BE DIMENSIONED TO LARGEST IRREP OF THIS
C                   ARRAY.
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION W(MAXCOR),DENOM(1),EAA,EAB,ENERGY,
     &       SDOT,ETAA,ESCF,EAAC,EAAXC,EAAX
      DOUBLE PRECISION ESING
      LOGICAL NONHF
      LOGICAL POLSCF
      CHARACTER*2 SPCASE(3)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFEA(2),NFMI(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),D(18)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /SHIFT/  ISHIFT,NDRGEO

      DATA SPCASE /'AA','BB','AB'/

      ETAA =0.0

      NONHF=(IFLAGS(38)+IFLAGS(77).NE.0)

      IF(NONHF)THEN
         WRITE(6,1011)
1011   FORMAT(T3,'@FORMT2-I, Singles contribution will be calculated.')
       IF(IFLAGS(39).EQ.0)THEN
        WRITE(6,1012)
1012    FORMAT(T3,'@FORMT2-W, MBPT(2) energies are correct only ',
     &           'for semicanonical orbitals.')
       ENDIF
      ENDIF 
      WRITE(6,*) 
      Write(6,"(3x,a,a)")"MBPT(2) energies separated into Coulomb and ",
     &                "Exchange contributions."

      CALL GETREC(20,'JOBARC','SCFENEG ',IINTFP,ESCF)
      WRITE(6,1005)
      WRITE(6,1004)ESCF

      DO ISPIN=1,1+IUHF

         EAA   = 0.0D0
         EAAC  = 0.0D0 
         EAAXC = 0.0D0

         LISTWC=213+ISPIN 
         LISTW =13+ISPIN 
         LEN_WC=IDSYMSZ(1,ISYTYP(1,LISTWC),ISYTYP(2,LISTWC))
         LEN_W =IDSYMSZ(1,ISYTYP(1,LISTW),ISYTYP(2,LISTW))

         I000   = 1
         I020   = I000 + LEN_W
         I030   = I020 + LEN_WC
         I040   = I030 + LEN_WC
         I050   = I040 + LEN_WC
         I060   = I050 + LEN_WC
         I070   = I060 + LEN_W
         IF(I070.GT.MAXCOR)CALL INSMEM('@-FORM_SPLIT_E2',I070,MAXCOR)
C
         CALL GETALL(W(I000), LEN_W, 1, LISTW)
         CALL GETALL(W(I020), LEN_WC, 1,LISTWC)

#ifdef _DEBUG_LVLM
      call checksum("W(A<B,I<J)",W(I000),LEN_W)
      call checksum("WC(AB,IJ)",W(I020),LEN_WC)
      Write(6,*) "The length len_w and len_wc", LEN_W, LEN_WC
#endif 
         INCRM1 = 0

         DO IRREPR = 1, NIRREP
            IRREPL = DIRPRD(IRREPR,1)

            IF (ISPIN .EQ. 1) THEN
                LEN1 = IRPDPD(IRREPL,19)
                LEN2 = IRPDPD(IRREPR,21)
            ELSE
                LEN1 = IRPDPD(IRREPL,20)
                LEN2 = IRPDPD(IRREPR,22)
            ENDIF
            CALL TRANSP(W(I020+INCRM1),W(I030+INCRM1),LEN2,LEN1)
            INCRM1 = INCRM1 + LEN1 * LEN2

         ENDDO 

CSSS         Write(6,*) "After transpose"
CSSS         call checksum("WC(AB,IJ)",W(I030),LEN_WC)

         IOFF1  = I030
         IOFF2  = I040

         DO IRREPR = 1, NIRREP

            IRREPL = DIRPRD(IRREPR,1)
            NUMDISC=IRPDPD(IRREPR,18+ISPIN)

            NDSSIZC=IRPDPD(IRREPL,20+ISPIN)
            NDSSIZW=IRPDPD(IRREPL,2+ISPIN)

CSSS      call output(W(IOff1),1,NDSSIZC,1,NUMDISC,NDSSIZC,NUMDISC,1)
            CALL SQSYM(IRREPL,POP(1,ISPIN),NDSSIZW,NDSSIZC,NUMDISC,
     &                 W(IOFF2),W(IOFF1))
CSSS      write(6,*) 
CSSS      call output(W(ioff2),1,NDSSIZW,1,NUMDISC,NDSSIZW,NUMDISC,1)
            IOFF1 = IOFF1 + NUMDISC * NDSSIZC  
            IOFF2 = IOFF2 + NUMDISC * NDSSIZW 
         ENDDO

         INCRM1 = 0
         DO IRREPR = 1, NIRREP
            IRREPL = DIRPRD(IRREPR,1)

            IF (ISPIN .EQ. 1) THEN
                LEN1 = IRPDPD(IRREPL,19)
                LEN2 = IRPDPD(IRREPR,3)
            ELSE
                LEN1 = IRPDPD(IRREPL,20)
                LEN2 = IRPDPD(IRREPR,4)
            ENDIF

CSSS           call output(W(I040+INCRM1),1,len2,1,len1,len2,len1,1)
            CALL TRANSP(W(I040+INCRM1),W(I050+INCRM1),LEN1,LEN2)
CSSS            call checksum("Intermediate",W(I050+INCRM1),LEN1*LEN2)
CSSS            call output(W(I050+INCRM1),1,len1,1,len2,len1,len2,1)
            INCRM1 = INCRM1 + LEN1 * LEN2
         ENDDO 

         IOFF1  = I050
         IOFF2  = I060

         DO IRREPR = 1, NIRREP
             IRREPL = DIRPRD(IRREPR,1)

             NUMDISC=IRPDPD(IRREPR,2+ISPIN)
             NDSSIZC=IRPDPD(IRREPL,18+ISPIN)
             NDSSIZW=IRPDPD(IRREPL,ISPIN)

CSSS             call output(W(IOFF1), 1, NDSSIZC, 1, NUMDISC, NDSSIZC, 
CSSS     &                   NUMDISC, 1)

             CALL SQSYM(IRREPL,VRT(1,ISPIN),NDSSIZW,NDSSIZC,NUMDISC,
     &                  W(IOFF2),W(IOFF1))

CSSS            call checksum("Intermediate22",W(ioff2),NUMDISC*NDSSIZW)
CSSS             call output(W(IOFF2), 1, NDSSIZW, 1, NUMDISC, NDSSIZW, 
CSSS     &                   NUMDISC, 1)

             IOFF1 = IOFF1 +  NUMDISC *  NDSSIZC
             IOFF2 = IOFF2 +  NUMDISC *  NDSSIZW 
         ENDDO

         CALL PUTALL(W(I060), LEN_W, 1, LISTW)

#ifdef _DEBUG_LVLM
      call checksum("W(A<B,I<J)",W(I050),LEN_W)
#endif 
         DO IRREP = 1, NIRREP
            LISTW = 13 + ISPIN
            LISTT = 43 + ISPIN
            NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
            DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))

            I030 = I020 + NUMDIS * DISSIZ 
            I040 = I030 + NUMDIS * DISSIZ 
            CALL GETLST(W(I020),1,NUMDIS,2,IRREP,LISTW)
            CALL GETLST(W(I030),1,NUMDIS,2,IRREP,LISTT)
            EAAC=EAAC+SDOT(NUMDIS*DISSIZ,W(I020),1,W(I030),1)
CSSS            write(6,"(a,1X,f12.7)") "The Coulomb Energy", EAAC
         ENDDO 

#ifdef _DEBUG_LVLM
      write(6,"(a,1X,f12.7)") "The Coulomb Energy", EAAC
#endif 
        IF(IUHF.EQ.0)EAAC=2.0*EAAC
        WRITE(LUOUT,999)SPCASE(ISPIN),EAAC
        ETAA=ETAA+EAAC

        CALL PUTALL(W(I000), LEN_W, 1, LISTW)

#ifdef _DEBUG_LVLM
      write(6,*) "Restored the list 14,15"
      call checksum("W(A<B,I<J)",W(I000),LEN_W)
#endif 
         DO IRREP = 1, NIRREP
            LISTW = 13 + ISPIN
            LISTT = 43 + ISPIN
            NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
            DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))

            I000 = 1
            I010 = I000 + NUMDIS * DISSIZ
            I020 = I010 + NUMDIS * DISSIZ
            CALL GETLST(W(I000),1,NUMDIS,2,IRREP,LISTW)
            CALL GETLST(W(I010),1,NUMDIS,2,IRREP,LISTT)
            EAAXC=EAAXC+SDOT(NUMDIS*DISSIZ,W(I000),1,W(I010),1)
         ENDDO  

#ifdef _DEBUG_LVLM
      write(6,"(a,1X,f12.7)") "The Exchng. Energy", EAAXC
#endif 
         IF(IUHF.EQ.0) THEN
           EAAX=(2.0D0*EAAXC-EAAC)
         ELSE
           EAAX=(EAAXC-EAAC)
         ENDIF
             
         WRITE(LUOUT,1000)SPCASE(ISPIN),EAAX
         ETAA=ETAA+EAAX

      ENDDO 

      LISTW=16 
      LISTD=50
      LISTT=46
      EAB  =0.0

      DO IRREP=1,NIRREP
         NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
         DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))

         I000 = 1
         I010 = I000 + NUMDIS * DISSIZ 
         I020 = I010 + NUMDIS * DISSIZ 
         CALL GETLST(W(I000),1,NUMDIS,2,IRREP,LISTW)
         CALL GETLST(W(I010),1,NUMDIS,1,IRREP,LISTT)
         EAB=EAB+SDOT(NUMDIS*DISSIZ,W(I000),1,W(I010),1)

      ENDDO 

      WRITE(LUOUT,1006)SPCASE(3),EAB
      ENERGY=EAB+ETAA
C
C CALCULATE SINGLES CONTRIBUTION TO SECOND-ORDER ENERGY
C  FOR SEMICANONICAL ORBITALS
C
      ESING=0.0

      IF(NONHF)THEN

        DO ISPIN=1,1+IUHF
         I000=1
         I010=I000+NT(ISPIN)
         I020=I010+NT(ISPIN)

         CALL GETLST(W(I000),1,1,1,ISPIN,90)
#ifdef _DEBUG_LVLM1
         Write(6,*) "T1 amps"
         call output(W(I000), 1, NT(ispin), 1, 1, NT(ispin), 1)
#endif
         CALL GETLST(W(I010),1,1,1,ISPIN+2,93)

#ifdef _DEBUG_LVLM1
         Write(6,*) "FAI amps"
         call output(W(I010), 1, NT(ispin), 1, 1, NT(ispin), 1)
#endif
         ESING=ESING+SDOT(NT(ISPIN),W(I000),1,W(I010),1)

         If ((Iflags(h_IFLAGS_calc) .EQ. 45) .OR. 
     &       (Iflags(h_IFLAGS_calc) .EQ. 46)) 
     &       ESING = 2.0D0*ESING 
   
        ENDDO 
         WRITE(LUOUT,1003)ESING
      ENDIF
C
      WRITE(LUOUT,1002)ENERGY
      WRITE(LUOUT,1001)ENERGY+ESCF+ESING
C
      WRITE(6,1005)
 999  FORMAT(T15,'E2_CO(',A2,')',T37,'=',F18.12,' a.u.')
1000  FORMAT(T15,'E2_EX(',A2,')',T37,'=',F18.12,' a.u.')
1001  FORMAT(T15,'Total MBPT(2) energy',T37,'=',F18.12,' a.u.')
1002  FORMAT(T15,'E2(TOT)',T37,'=',F18.12,' a.u.')
1003  FORMAT(T15,'E2(SINGLE) ',T37,'=',F18.12,' a.u.')
1004  FORMAT(T15,'E(SCF)',T37,'=',F18.12,' a.u.')
1005  FORMAT(T15,46('-'))
1006  FORMAT(T15,'E2(',A2,')',T37,'=',F18.12,' a.u.')

       RETURN
       END
