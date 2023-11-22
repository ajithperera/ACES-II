      SUBROUTINE S2PPS2(ICORE,MAXCOR,IUHF, ISIDE, ISPIN)
C
C   THIS SUBROUTINE CALCULATES THE CONTRIBUTION OF S2 TO
C   S2 INVOLVING THE PP BLOCK OF THE EFFECTIVE ONE-PARTICLE 
C   INTEGRALS.
C
C CONTRACTION :
C
C  Z(ji,ap) = SUM [Hbar(a,e) * S2(ji,ep)]  S2 is spin adapted RHF
C             e 
C  Ispin = spin-label of e and Mixspin = spin-lable of j,e
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      LOGICAL SKIP
C
      DIMENSION ICORE(MAXCOR), NUMSZS(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SINFO/ NS(8), SIRREP
      COMMON/SLISTS/LS1IN, LS1OUT, LS2IN(2,2), LS2OUT(2,2)
C
      DATA ONE  /1.0/
      DATA ZILCH/0.0/
      DATA ONEM/-1.0/

      Write(6,"(a)") "-----------------S2PPS2---------------------"

      DO MIXSPIN = 1, 1+IUHF

         LISTS2EX = LS2OUT(ISPIN, MIXSPIN + 1 - IUHF)
         LISTS2IN = LS2IN(ISPIN,  MIXSPIN + 1 - IUHF)

         DO XIRREP=1,NIRREP
            IRREPA0 = DIRPRD(XIRREP,SIRREP)

            DISSYS  =IRPDPD(XIRREP,ISYTYP(1,LISTS2IN))
            DISSYSEX=IRPDPD(XIRREP,ISYTYP(1,LISTS2EX))
            NUMDSS  =VRT(IRREPA0,MIXSPIN)*NS(SIRREP)

            IF (DISSYSEX .NE. DISSYS) THEN
               WRITE(6,"(2A)") " Can not proceeed dut to inconsitency",
     &                         " in internal array dimenssions"
               CALL ERREX
            ENDIF 
            I000=1
            I010=I000+IINTFP*DISSYS*NUMDSS
            I020=I010+IINTFP*DISSYSEX*NUMDSS
            IEND=I020+NFEA(MIXSPIN)
            IF(IEND .GT.MAXCOR) CALL INSMEM("S2PPS2",IEND,MAXCOR)

            CALL GETLST(ICORE(I000),1,NUMDSS,1,XIRREP,LISTS2IN)
            CALL GETLST(ICORE(I020), 1,1,1,MIXSPIN,92)

#ifdef _DEBUG_LVL0
      call checksum("-> s2   :",icore(i000),dissys*numdss)
      call checksum("->hs2   :",icore(i010),dissyex*numdss)
      call checksum("->fae   :",icore(i020),nfea(mixspin))
#endif
C  
C PERFORM MATRIX MULTIPLICATION, ISIDE=1,S2IN(IJ,PE)*H(E,A)->S2EX(IJ,PA)=S2EX(JI,PA)
C                                ISIDE=2,S2IN(IJ,PE)*H(A,E)->S2EX(IJ,PA)=S2EX(JI,PA)

             CALL DZERO(ICORE(I010),DISSYSEX*NUMDSS)
             IOFF_S2IN = I000
             IOFF_S2EX = I010
             IOFF_HBAR = I020

             DO IRREPA = 1,NIRREP
                IRREPE = DIRPRD(IRREPA,1)
                IRREPP = DIRPRD(IRREPE,XIRREP) 

                NROW = DISSYS*NS(IRREPP)
                NCOL = VRT(IRREPA,MIXSPIN)
                NSUM = VRT(IRREPE,MIXSPIN)
CSSS                call output(icore(ioff_s2in),1,nrow,1,nsum,nrow,nsum,1)
CSSS                call output(icore(ioff_hbar),1,nsum,1,ncol,nsum,ncom,1)

                IF (NROW*NCOL*NSUM .NE. 0) THEN

                IF (ISIDE.EQ.1) THEN
                    CALL XGEMM("N","N",NROW,NCOL,NSUM,ONE,
     &                          ICORE(IOFF_S2IN),NROW,ICORE(IOFF_HBAR),
     &                          NSUM,ONE,ICORE(IOFF_S2EX),NROW)

                 ELSE
                    CALL XGEMM('N','T',NROW,NCOL,NSUM,ONE,
     &                          ICORE(IOFF_S2IN),NROW,ICORE(IOFF_HBAR),
     &                          NSUM,ONE,ICORE(IOFF_S2EX),NROW)
                ENDIF
                ENDIF
CSSS                call output(icore(ioff_s2ex),1,nrow,1,ncol,nrow,ncol,1)
                IOFF_S2IN = IOFF_S2IN + NROW*NSUM*IINTFP
                IOFF_S2EX = IOFF_S2EX + NROW*NCOL*IINTFP
                IOFF_HBAR = IOFF_HBAR + NCOL*NSUM*IINTFP
             ENDDO 

#ifdef _DEBUG_LVL0
      call checksum("s2->s2  :",icore(i010),numdss*dissysex)
#endif
             CALL GETLST(ICORE(I000),1,NUMDSS,1,XIRREP,LISTS2EX) 
             CALL DAXPY(DISSYSEX*NUMDSS,ONE,ICORE(I010),1,ICORE(I000),1)
             CALL PUTLST(ICORE(I000),1,NUMDSS,1,XIRREP,LISTS2EX)
#ifdef _DEBUG_LVL0
      call checksum("s2->s2f :",icore(i000),numdss*dissys)
#endif
         ENDDO     
      ENDDO 
C  
      RETURN
      END

