      SUBROUTINE CHECK_IPVECS(IUHF, ISPIN, ICORE, MAXCOR, LS1OUT,
     $                        LS2OUT)
C
C THE LISTS LS1OUT AND LS2OUT ARE PUT TO ZERO
C
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR), LS2OUT(2,2)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SINFO/ NS(8), SIRREP
C
      DATA IONE /1/

      IF ((IUHF .EQ. 0) .AND. (ISPIN .NE. 1)) THEN
         RETURN
      ENDIF
C
      NS1 = NS(SIRREP)*POP(SIRREP,ISPIN) 
      IOFF = IONE
      JOFF = IONE
      CALL GETLST(ICORE(IOFF),1,1,1,ISPIN,LS1OUT)
      call checksum("chk_ips1:",icore,NS(SIRREP)*POP(SIRREP,ISPIN))
C
      DO 50 KSPIN = 1, 1+IUHF
      IF (IUHF. NE.0) THEN
         MSPIN = KSPIN
      ELSE
         MSPIN = 3-ISPIN
      ENDIF
         LISTS2EX = LS2OUT(ISPIN, MSPIN)
         IOFF = IOFF + NS1
         JOFF = JOFF + NS1 
         LEN  = 0
         DO 100 IRREPR = 1 , NIRREP
            IRREPL = DIRPRD(IRREPR,SIRREP)
            DISSYS = IRPDPD(IRREPL,ISYTYP(1, LISTS2EX))
            NUMDSS = VRT(IRREPR, MSPIN) * NS(SIRREP)
            IF ((NUMDSS * DISSYS) .GT. 0) THEN
               CALL GETLST(ICORE(IOFF), 1, NUMDSS, 1, IRREPL, LISTS2EX)
               IOFF = IOFF + NUMDSS * DISSYS
            ENDIF
            LEN=LEN+NUMDSS * DISSYS
 100     CONTINUE

#ifdef _DEBUG_LVL0
      call checksum("chk_ips2:",icore(joff),len)
      Joff = Joff+Len 
#endif 

 50   CONTINUE
#ifdef _DEBUG_LVL0
      call checksum("chck_ipf:",icore,Ioff-1)
      write(6,*)
#endif 
      RETURN
      END
