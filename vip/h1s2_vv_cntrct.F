      SUBROUTINE H1S2_VV_CNTRCT(S2IN,HBAR,S2EX,NUMDIS,DISSIZ,L1,
     &                          XIRREP,ISPIN,MIXSPIN,FACT)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      INTEGER DISSIZ,XIRREP,DIRPRD,POP,VRT
      DIMENSION S2IN(DISSIZ*NUMDIS),HBAR(L1),S2EX(DISSIZ*NUMDIS)

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)

      DATA IONE,ONE /1,1.0D0/

C H(A,E)*S2IN(JIP,E) -> S2(JIP,A)

         IOFF_AE  = IONE
         IOFF_JIP = IONE
         IOFF_JIA = IONE

         DO IRREPE = 1, NIRREP
            IRREPA = DIRPRD(IRREPE,1)
            
            NROW = DISSIZ
            NCOL = POP(IRREPA,MIXSPIN)
            NSUM = POP(IRREPE,MIXSPIN)
             
            IF (NROW*NCOL*NSUM .NE. 0) THEN
                CALL XGEMM('N','T',NROW,NCOL,NSUM,FACT,
     &                      HBAR(IOFF_AE),NROW,
     &                      S2IN(IOFF_JIP),NCOL,ONE,
     &                      S2EX(IOFF_JIA),NROW)

             ENDIF
             IOFF_AE   = IOFF_AE  + NSUM * NCOL
             IOFF_JIP  = IOFF_JIP + NROW * NSUM
             IOFF_JIA  = IOFF_JIA + NROW * NCOL
          ENDDO
C      ENDDO

      RETURN
      END 

