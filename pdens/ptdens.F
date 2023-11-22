      SUBROUTINE PTDENS(IRREPXL,IRREPXR,IRREPX,DOO,DVV,DVO,DOV,
     &                  SCR,MXCOR,IUHF,FACT,R0,L0,ZNORM,
     &                  LSTGRL,LSTGTL,LSTGRLOF,LSTGTLOF,LSTTMP,
     &                  LSTTMPOF,LSTR1,LSTL1,LSTR1OFF,LSTL1OFF,
     &                  LISTR2,LISTL2,LSTR2RS,LSTL2RS,LSTT1,
     &                  LSTT1OFF)
C
C DRIVER FOR FORMATION OF GENERAL ONE-PARTICLE REDUCED DENSITY MATRIX
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION L0
      DIMENSION SCR(MXCOR),LENVV(2),LENOO(2),LENVO(2)
      DIMENSION DOO(*),DVV(*),DVO(*),DOV(*)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      DATA ONEM /-1.0D0/
C
      DO ISPIN=1,1+IUHF
         LENVV(ISPIN)=IRPDPD(IRREPX,18+ISPIN)
         LENOO(ISPIN)=IRPDPD(IRREPX,20+ISPIN)
         LENVO(ISPIN)=IRPDPD(IRREPX, 8+ISPIN)
      ENDDO
C
C INITIALIZE THE DENSITY MATRIX WITH GROUND STATE
C
C
      R0 =0.0D0
      CALL ZERO(DOO,LENOO(1)+IUHF*LENOO(2))
      CALL ZERO(DVV,LENVV(1)+IUHF*LENVV(2))
      CALL ZERO(DVO,LENVO(1)+IUHF*LENVO(2))
      CALL ZERO(DOV,LENVO(1)+IUHF*LENVO(2))
C
C This is ground state density <0|(1+L){p+q}|0>. The case that I am 
C interested in L0 = 1 and R0 = 1.
C
      IF (IRREPXR.EQ.1) THEN
         CALL GGSDEN(IRREPXL,DOO,DVV,DVO,DOV, LSTL1OFF,LSTL1,0,90,
     &               LISTL2,44,LSTGTL,LSTGTLOF,L0,SCR,MXCOR,IUHF)
       CALL SSCAL(LENOO(1)+IUHF*LENOO(2),R0,DOO,1)
       CALL SSCAL(LENVV(1)+IUHF*LENVV(2),R0,DVV,1)
       CALL SSCAL(LENVO(1)+IUHF*LENVO(2),R0,DVO,1)
       CALL SSCAL(LENVO(1)+IUHF*LENVO(2),R0,DOV,1)
      ENDIF
C                      +
C NOW CALCULATE <0|L {p q exp(T)} R_n|0> CONTRIBUTIONS TO THE
C DENSITY, WHERE n=1 AND 2.
C
C EVALUATE L2*R2 G-TYPE INTERMEDIATES AND PLACE ON LISTS 91 AND 92
C
      IF (LSTGRL.NE.-1) THEN
          CALL GFORMG_MOD_4PDENS(IRREPXL,IRREPXR,LISTL2,LISTR2,
     &                          LSTGRL,SCR,
     &                          MXCOR,LSTGRLOF,ONEM,IUHF)
      ENDIF
C
C EVALUATE L2*T2 G-TYPE INTERMEDIATES AND PLACE ON LISTS 491 AND 492
C
      CALL GFORMG(IRREPXL,1,LISTL2,44,LSTGTL,SCR,MXCOR,LSTGTLOF,
     &            ONEM,IUHF)
C
C OCCUPIED-VIRTUAL TERMS
C
CSSS      Print*, "The occ-vrt dens"
CSSS      call output(dov, 1, lenvo, 1, 1, lenvo, 1, 0)
C

      CALL GDENSOV(IRREPXL,IRREPXR,IRREPX,DOV,R0,SCR,MXCOR,IUHF,
     &             LSTGRL,LSTGTL,LSTGRLOF,LSTGTLOF,LSTTMP,LSTTMPOF,
     &             LSTR1,LSTL1,LSTR1OFF,LSTL1OFF,LISTR2,LISTL2,
     &             LSTR2RS,LSTL2RS)
#ifdef _DEBUG_LVLM
          Print*, "----Afetr GDENSOV----"
      Print*, "The occ-vrt density A and B (R1*L2)"
      call output(dov(1), 1, lenvo(1), 1, 1, lenvo(1), 1, 0)
      call output(dov(1+lenvo(1)), 1, lenvo(2), 1, 1, lenvo(2), 1, 0)
#endif
C
C VIRTUAL-VIRTUAL AND OCCUPIED-OCCUPIED TERMS
C
      CALL GDENS1(IRREPXL,IRREPXR,IRREPX,DOO,DVV,DOV,R0,SCR,MXCOR,
     &            IUHF,LSTGRL,LSTGTL,LSTGRLOF,LSTGTLOF,LSTTMP,
     &            LSTTMPOF,LSTR1,LSTL1,LSTR1OFF,LSTL1OFF,LISTR2,
     &            LISTL2,LSTR2RS,LSTL2RS,LSTT1,LSTT1OFF)
C
#ifdef _DEBUG_LVLM
      Print*, "----Afetr GDENS1----"
      Print*, "The occ-occ density A and B"
      call output(doo(1), 1, lenoo(1), 1, 1, lenoo(1), 1, 0)
      call output(doo(1+lenoo(1)), 1, lenoo(2), 1, 1, lenoo(2), 1, 0)
      Write(6,*)
      Print*, "The vrt-vrt density A and B"
      call output(dvv(1), 1, lenvv(1), 1, 1, lenvv(1), 1, 0)
      call output(dvv(1+lenvv(1)), 1, lenvv(2), 1, 1, lenvv(2), 1, 0)
      Write(6,*)
      Print*, "The occ-vrt density A and B"
      call output(dov(1), 1, lenvo(1), 1, 1, lenvo(1), 1, 0)
      call output(dov(1+lenvo(1)), 1, lenvo(2), 1, 1, lenvo(2), 1, 0)
      Write(6,*)
      Print*, "The vrt-ort density A and B"
      call output(dvo(1), 1, lenvo(1), 1, 1, lenvo(1), 1, 0)
      call output(dvo(1+lenvo(1)), 1, lenvo(2), 1, 1, lenvo(2), 1, 0)
#endif
C VIRTUAL-OCCUPIED TERMS

      CALL GDENSVO(IRREPXL,IRREPXR,IRREPX,DVO,DOV,R0,L0,ZNORM,SCR,
     &             MXCOR,IUHF,LSTGRL,LSTGTL,LSTGRLOF,LSTGTLOF,LSTTMP,
     &             LSTTMPOF,LSTR1,LSTL1,LSTR1OFF,LSTL1OFF,LISTR2,
     &             LISTL2,LSTR2RS,LSTL2RS)
C
#ifdef _DEBUG_LVLM
      Print*, "----Afetr GDENSVO----"
      Print*, "The occ-occ density A and B"
      call output(doo(1), 1, lenoo(1), 1, 1, lenoo(1), 1, 0)
      call output(doo(1+lenoo(1)), 1, lenoo(2), 1, 1, lenoo(2), 1, 0)
      Write(6,*)
      Print*, "The vrt-vrt density A and B"
      call output(dvv(1), 1, lenvv(1), 1, 1, lenvv(1), 1, 0)
      call output(dvv(1+lenvv(1)), 1, lenvv(2), 1, 1, lenvv(2), 1, 0)
      Write(6,*)
      Print*, "The occ-vrt density A and B"
      call output(dov(1), 1, lenvo(1), 1, 1, lenvo(1), 1, 0)
      call output(dov(1+lenvo(1)), 1, lenvo(2), 1, 1, lenvo(2), 1, 0)
      Write(6,*)
      Print*, "The vrt-occ density A and B"
      call output(dvo(1), 1, lenvo(1), 1, 1, lenvo(1), 1, 0)
      call output(dvo(1+lenvo(1)), 1, lenvo(2), 1, 1, lenvo(2), 1, 0)
#endif
C
      RETURN
      END
