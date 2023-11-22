C
      SUBROUTINE WRTDBLPRJ(SCR, DMAXCOR, IRREPX, ICASE, IOFFSET)
C     
      IMPLICIT INTEGER (A - Z)
      DOUBLE PRECISION SCR(DMAXCOR)
      DIMENSION LS2OUT(2, 2), IOFFSET(8, 3)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      LS2OUT(1,1) = 444
      LS2OUT(1,2) = 446
      LS2OUT(2,1) = 446
      LS2OUT(2,2) = 445
C
C Let's calcualte the length of the projection vector
C
      LENGTH = IOFFSET(NIRREP, ICASE)
C
      IBEGIN = LENGTH
C
      IF (ICASE .EQ. 1) THEN
         ASPIN = 1
         BSPIN = 2
      ELSEIF (ICASE .EQ. 2) THEN
         ASPIN = 1
         BSPIN = 1
      ELSEIF (ICASE .EQ. 3) THEN
         ASPIN = 2
         BSPIN = 2
      ENDIF 
C
      LISTS2EX = LS2OUT(ASPIN, BSPIN)
C
      DO 60 RIRREP = 1, NIRREP
C
         LIRREP = DIRPRD(RIRREP, IRREPX)
C
         IF (ICASE .EQ. 1) THEN
            DISSYS = IRPDPD(LIRREP, ISYTYP(1, LISTS2EX))
            NUMDSS = IRPDPD(RIRREP, ISYTYP(2, LISTS2EX))
            DISSYA = DISSYS
            NUMDSA = NUMDSS
         ELSE
            DISSYA = IRPDPD(LIRREP, ISYTYP(1, LISTS2EX))
            NUMDSA = IRPDPD(RIRREP, ISYTYP(2, LISTS2EX))
            DISSYS = IRPDPD(LIRREP, 18 + ASPIN)
            NUMDSS = IRPDPD(RIRREP, 20 + ASPIN)
         ENDIF
C
         IF (DISSYA*NUMDSA .GT. 0) THEN
C
            I000 = 1 
            I010 = IBEGIN + DISSYS*NUMDSS
            INEED = DISSYS*NUMDSS
            IAVIL = DMAXCOR - IBEGIN
            IF (I010 .GT. DMAXCOR) CALL INSMEM("@-WRTDBLPART", INEED, 
     &                                          IAVIL)
            IOFF = IOFFSET(RIRREP, ICASE) 
C
            IF (ICASE .EQ. 1) THEN
               CALL PUTLST(SCR(I000 + IOFF), 1, NUMDSS, 1, RIRREP,
     &                     LISTS2EX)
            ELSE
C
C SQUEEZE ARRAY IN PROPER FORM A<B, I<J
C
               CALL SQSYM(LIRREP, VRT(1,ASPIN), DISSYA, DISSYS,
     &                    NUMDSS,SCR(I010), SCR(I000 + IOFF))
C
               CALL TRANSP(SCR(I010), SCR(I000 + IOFF),NUMDSS,DISSYA)
               CALL SQSYM(RIRREP,POP(1,ASPIN), NUMDSA, NUMDSS,
     &                    DISSYA, SCR(I010), SCR(I000 + IOFF))
               CALL TRANSP(SCR(I010), SCR(I000 + IOFF), DISSYA,
     &                     NUMDSA)
               CALL PUTLST(SCR(I000 + IOFF), 1, NUMDSA, 1, RIRREP,
     &                     LISTS2EX)    
            ENDIF   
C
         ENDIF
C
 60   CONTINUE
C
      RETURN
      END
