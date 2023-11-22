      SUBROUTINE FAOPOP(ALASKA)
      IMPLICIT INTEGER (A-Z)
      LOGICAL SPHHRM,MOEQAO, ALASKA
      COMMON /BASTYP/ SPHHRM,MOEQAO
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /AOPOPS/ AOPOP(8),MOPOP(8),NAO,NAOSPH,NMO
      COMMON /SZAOGM/ AOGMSZ(2,4,100)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /AOOFST/ INDOCC(8,2),INDVRT(8,2)
      COMMON /INFO/   NOCCO(2),NVRTO(2)
      IONE=1
      NMO=NOCCO(1)+NVRTO(1)
      CALL GETREC(20,'JOBARC','NAOBASFN',IONE,NAO) 
      CALL GETREC(20,'JOBARC','NBASTOT ',IONE,NAOSPH) 
      IF (ALASKA) NAO = NAOSPH
c----------------------------------------------------------------
c---   To correct for the drop--mo cases  ----  Mar. 94,  KB ----- 
c----------------------------------------------------------------
      naosph0 = naosph 
      CALL GETREC(20,'JOBARC','NUMDROPA',IONE,NDROP) 
      IF (NDROP.NE.0) NAOSPH=NAOSPH + NDROP
c----------------------------------------------------------------
      CALL GETREC(20,'JOBARC','NUMBASIR',NIRREP,MOPOP)
      IF(MOEQAO .OR. ALASKA)THEN
       CALL GETREC(20,'JOBARC','NUMBASIR',NIRREP,AOPOP)
      ELSE
       CALL GETREC(20,'JOBARC','FAOBASIR',NIRREP,AOPOP)
      ENDIF
      RETURN
      END
