      SUBROUTINE PUTS_HHP(S, NIP, ISPIN, IUHF, LISTS1, LISTS2)
C
C  THIS ROUTINE PUTS THE COMPLETE S-VECTOR CONTAINED IN S, ONTO THE 
C  PROPER LISTS. THE ELEMENTS IN S ARE ORDERED AS FOLLOWS: 
C
C    SINGLES(ISPIN), DOUBLES (ISPIN, ISPIN), DOUBLES(ISPIN, 3-ISPIN)  [ UHF ]
C
C    SINGLES(ISPIN), DOUBLES (ISPIN, 3-ISPIN)                         [ RHF ]
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION S
      DIMENSION S(NIP), LISTS2(2,2)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SINFO/ NS(8), SIRREP
C
      CALL PUTLST(S(1), 1, 1, 1, ISPIN, LISTS1)
#ifdef _DEBUG_LVLM
      write(6,"(a)") "@-puts_hhp-S"
      write(6,"(6(1x,f12.6))")(s(i),i=1,Nip) 
      Call checksum("@putshhs:",S,Pop(Sirrep,Ispin))
      write(6,"(6(1x,F10.5))") (s(i),i=1,Pop(Sirrep,Ispin))
#endif
      IOFF = 1 + POP(SIRREP, ISPIN) 
      IF (IUHF .NE. 0) THEN
         CALL GETLEN2(LENAA, IRPDPD(1, ISYTYP(1, LISTS2(ISPIN, ISPIN))),
     $      VRT(1,ISPIN),NS)
         CALL PUTALLS2_HHP(S(IOFF),LENAA, VRT(1,ISPIN), NS, 1,
     $      LISTS2(ISPIN, ISPIN))
         IOFF = IOFF + LENAA
      ENDIF
#ifdef _DEBUG_LVLM
      call checksum("s2aaaa  :",s(1+pop(sirrep,ispin)),lenaa)
#endif 
      MSPIN = 3 - ISPIN
      CALL GETLEN2(LENAB, IRPDPD(1, ISYTYP(1, LISTS2(ISPIN, MSPIN))),
     $      VRT(1,MSPIN+IUHF-1),NS)
      CALL PUTALLS2_HHP(S(IOFF),LENAB, VRT(1,MSPIN+IUHF-1), NS, 1,
     $      LISTS2(ISPIN,MSPIN))
#ifdef _DEBUG_LVLM
      call checksum("s2abab  :",s(1+pop(sirrep,ispin)+lenaa),lenab)
      Write(6,"(a,i6)"), "Length =", Lenab+Pop(Sirrep,Ispin)
      Length=Lenab+Pop(Sirrep,Ispin)
      Call checksum("@putshhp:",S,length)
#endif
C
      RETURN
      END
