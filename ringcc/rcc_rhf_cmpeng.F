      SUBROUTINE RCC_RHF_CMPENEG(ICORE,MAXCOR,IUHF,LIST_IN,E_S,E_T)

      IMPLICIT INTEGER(A-Z)
    
      DOUBLE PRECISION E_S, E_T,DDOT
      DIMENSION ICORE(MAXCOR)

      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
C
      E_S = 0.0D0
      E_T = 0.0D0

      DO IMULT= 1, 2
         IF (IMULT .EQ. 1) THEN
             LISTT = LIST_IN + 3
             LISTW = 197
          ELSE IF (IMULT .EQ. 2) THEN
             LISTT = LIST_IN + 1
             LISTW = 198
          ENDIF 

#ifdef _DEBUG_LVL0
      write(6,"(a,i4)") "The input T2 list: ", listt 
#endif 
          DO IRREP=1,NIRREP
             DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
             NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
             NSIZE = DISSYT * NUMSYT 

             IF (MIN(NUMSYT,DISSYT).NE.0) THEN
                I000 = 1
                I010 = I000 + NSIZE * IINTFP
                I020 = I010 + NSIZE * IINTFP
                IEND = I020
                IF (IEND .GE. MAXCOR) CALL INSMEM("rcc_rhf_compeneg",
     +                                             IEND,MAXCOR)

                CALL GETLST(ICORE(I000),1,NUMSYT,1,IRREP,LISTT)
#ifdef _DEBUG_LVL0
                call checksum("T2s :",icore(i000),Nsize)
#endif 
                CALL GETLST(ICORE(I010),1,NUMSYT,1,IRREP,LISTW)
               
                IF (IMULT .EQ. 1) E_S = E_S + DDOT(NSIZE,ICORE(I010),
     +                                            1,ICORE(I000),1)
                IF (IMULT .EQ. 2) E_T = E_T + DDOT(NSIZE,ICORE(I010),
     +                                             1,ICORE(I000),1)
             ENDIF 
          ENDDO 

      ENDDO

#ifdef _DEBUG_LVL0
      WRITE(6,"(a,F15.10)") "E_S = ", E_S
      WRITE(6,"(a,F15.10)") "E_T = ", E_T
#endif 

      RETURN
      END 
