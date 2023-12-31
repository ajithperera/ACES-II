      SUBROUTINE RCL_EXPANDED_ABIJ(ICORE,MAXCOR,IUHF)

      IMPLICIT INTEGER (A-Z)

      DIMENSION ICORE(MAXCOR)

      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),IJUNK(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)

      DO ISPIN=1,1+IUHF
         DO IRREP=1,NIRREP
            NCOL_ILJ=IRPDPD(IRREP,ISYTYP(2,13+ISPIN))
            NROW_ALB=IRPDPD(IRREP,ISYTYP(1,13+ISPIN))
            NCOL_IAJ=IRPDPD(IRREP,ISYTYP(2,60+ISPIN))
            NROW_AAB=IRPDPD(IRREP,ISYTYP(1,60+ISPIN))

            I000 = 1
            I010 = I000 + NCOL_IAJ * NROW_AAB * IINTFP
            I020 = I010 + NCOL_IAJ * NROW_AAB * IINTFP
            IEND = I020
            IF (IEND .GT. MAXCOR) CALL INSMEM("rcl_expanded_abij",
     +                                         IEND,MAXCOR)

            CALL GETLST(ICORE(I000),1,NCOL_ILJ,1,IRREP,13+ISPIN)
C W(A<B;I<J) -> W(A<B;IJ)
            CALL SYMEXP(IRREP,POP(1,ISPIN),NROW_ALB,ICORE(I000))
C W(A<B;IJ) -> W(AB,IJ)
            CALL SYMEXP2(IRREP,VRT(1,ISPIN),NROW_AAB,NROW_ALB,
     +                   NCOL_IAJ,ICORE(I010),ICORE(I000))
#ifdef _DEBUG_LVLM
      write(6,"(a)") " <A<B||I<J> ->  (AB,IJ)"
      call output(icore(i010),1,NROW_AAB,1,NCOL_IAJ,NROW_AAB,
     +            NCOL_IAJ,1)
#endif
            CALL PUTLST(ICORE(I010),1,NCOL_IAJ,1,IRREP,2+ISPIN)
        ENDDO
      ENDDO

      RETURN 
      END
