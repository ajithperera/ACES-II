
#include "flags.h"
#include "disttype.h"

      SUBROUTINE MKT1AT2(WORK, MAXCOR)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DOUBLE PRECISION WORK(MAXCOR)
 
C
      INTEGER POP, VRT, NT, NFMI, NFEA, IRREP
      COMMON /SYM/ POP(8,2), VRT(8,2), NT(2), NFMI(2), NFEA(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22), ISYTYP(2,500), ID(18)
C
      INTEGER IRPDPD, ISYTYP, ID
      INTEGER ISYMSZ, ITOTALAA, ITOTALBB, ITOTALAB, LENABAASQ,
     &        LENABBBSQ, LENIJAASQ, LENIJBBSQ, DISSZAAI, DISSZBBI,
     &        DISSZAAT, DISSZBBT, DISSZWAB
C
      INTEGER I000, I010, I020, I030, I040, I050, I060, IEND, I
C
      ITOTALAA = ISYMSZ(9, 9)
      ITOTALBB = ISYMSZ(10, 10)
      ITOTALAB = ISYMSZ(10, 9)
      ISQTOTAA = ISYMSZ(1, 3)
      ISQTOTBB = ISYMSZ(2, 4)
C
   
      I000 = 1
      I010 = I000 + ITOTALAA
      I020 = I010 + ITOTALBB
      I030 = I020 + ITOTALAB
C
      CALL GETALL(WORK(I000), ITOTALAA, 1, 200)
      CALL GETALL(WORK(I010), ITOTALAA, 1, 201)
      CALL GETALL(WORK(I020), ITOTALAA, 1, 202)

      I040 = I030 + ITOTALAA
      I050 = I040 + ITOTALBB
      I060 = I050 + ITOTALAB
      I070 = I060 + ISQTOTAA
      I080 = I070 + ISQTOTBB
      I090 = I080 + ISQTOTAA
      I100 = I090 + ISQTOTBB
      I110 = I100 + ITOTALAB
      IEND = I110 + MAX(ITOTALAA, ITOTALBB, ITOTALAB)
C
      IF(IEND.GT.MAXCOR)CALL INSMEM('MKT1AT2', IEND, MAXCOR)
C
#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,*) "T2AA from disk"
      write(6,*)
      Write(6, "(5(F11.9,3x))") (WORK(I000+I-1), I=1, ITOTALAA)
      Write(6,*)
      Write(6,*) "T2BB from disk"
      write(6,*)
      Write(6, "(5(F11.9,3x))") (WORK(I010+I-1), I=1, ITOTALBB)
      Write(6,*)
      Write(6,*) "T2AB from disk"
      write(6,*)
      Write(6, "(5(F11.9,3x))") (WORK(I020+I-1), I=1, ITOTALAB)
#endif
   
      CALL SSTGEN(WORK(I000), WORK(I030), ITOTALAA, VRT(1, 1), 
     &            POP(1, 1), VRT(1, 1), POP(1, 1), WORK(I110), 1,
     &            "2413")
      CALL SSTGEN(WORK(I010), WORK(I040), ITOTALBB, VRT(1, 2), 
     &            POP(1, 2), VRT(1, 2), POP(1, 2), WORK(I110), 1, 
     &            "2413")
      CALL SSTGEN(WORK(I020), WORK(I050), ITOTALAB, VRT(1, 1), 
     &            POP(1, 1), VRT(1, 2), POP(1, 2), WORK(I110), 1,
     &            "1324")

#ifdef _DEBUG_LVLM
        Write(6,*)
        Write(6,*) "Reordered T2AA"
        write(6,*)
        Write(6, "(5(F11.9,3x))") (WORK(I030+I-1), I=1, ITOTALAA)
        Write(6,*)
        Write(6,*) "Reordered T2BB"
        write(6,*)
        Write(6, "(5(F11.9,3x))") (WORK(I040+I-1), I=1, ITOTALBB)
        Write(6,*)
        Write(6,*) "Reordered  T2AB"
        write(6,*)
        Write(6, "(5(F11.9,3x))") (WORK(I050+I-1), I=1, ITOTALAB)
#endif
C
C Dump the reorderd RPA T2 amplitudes to a file (OEP work).
C
      ITOTLWAB = ISYMSZ(13, 14)
C
      I120 =  IEND 
      I130 =  I120 + ITOTLWAB
      IEND =  I130 + ITOTLWAB
      IF(IEND.GT.MAXCOR)CALL INSMEM('MKT1AT2', IEND, MAXCOR)
      
      CALL GETALL(WORK(I120), ITOTALAB, 1, 16)
      IOFFB = I120
      IOFFE = I130
      DO IRREP = 1, NIRREP
         DISSZWAB = IRPDPD(IRREP, 13)
         NUMDZWIJ = IRPDPD(IRREP, 14)
         CALL SYMTRA(IRREP, POP(1, 1), POP(1,2), DISSZWAB, WORK(IOFFB),
     &               WORK(IOFFE))
         CALL DSCAL(DISSZWAB*NUMDZWIJ, 2.0D0, WORK(IOFFB), 1)
         CALL DAXPY(DISSZWAB*NUMDZWIJ, -1.0D0, WORK(IOFFE), 1, 
     &              WORK(IOFFB), 1)
         IOFFB = IOFFB + DISSZWAB*NUMDZWIJ
         IOFFE = IOFFE + DISSZWAB*NUMDZWIJ
      ENDDO

      Ecor = Ddot(ITOTALAB, WORK(I120), 1, WORK(I050), 1)
      Write(6, "(a,F14.10)")"The T(Ij,Ab)W(Ij,Ab) energy  = ", Ecor
      Call dmp_rpat2vecs(Work(I030), I060-I030, "RPAT2AMP")

      ISTARTAA = I030
      ISTARTBB = I040
      JSTARTAA = I030
      JSTARTBB = I040
      KSTARTAA = I060
      KSTARTBB = I070
      LSTARTAA = I080
      LSTARTBB = I090
      JOFFAA   = 0
      JOFFBB   = 0
      KOFFAA   = 0
      KOFFBB   = 0

      DO IRREPAB = 1, NIRREP
         DISSZAAI = IRPDPD(IRREPAB, 21)
         NUMDZAAI = IRPDPD(IRREPAB, 19) 
         DISSZBBI = IRPDPD(IRREPAB, 22)
         NUMDZBBI = IRPDPD(IRREPAB, 20)

         DISSZAAT = IRPDPD(IRREPAB, 3)
         NUMDZAAT = IRPDPD(IRREPAB, 1) 
         DISSZBBT = IRPDPD(IRREPAB, 4)
         NUMDZBBT = IRPDPD(IRREPAB, 2)
         DISSZABT = IRPDPD(IRREPAB, 14)
         NUMDZABT = IRPDPD(IRREPAB, 13)
C
#ifdef _DEBUG_LVLM1
	 Write(6,*) NUMDZAAT, DISSZAAT, DISSZAAI
	 Write(6,*) NUMDZBBT, DISSZBBT, DISSZBBI
         Write(6,*) ISTARTAA, ISTARTBB, JSTARTAA, JSTARTBB
#endif

         CALL ASSYM2(IRREPAB, VRT(1, 1), DISSZAAI, WORK(ISTARTAA))
         CALL ASSYM2(IRREPAB, VRT(1, 2), DISSZBBI, WORK(ISTARTBB))
C
#ifdef _DEBUG_LVLM1
        itmpaa =  Irpdpd(irrepab, 1)*Irpdpd(irrepab, 21)
        itmpbb =  Irpdpd(irrepab, 2)*Irpdpd(irrepab, 22)
        Write(6,*)
        write(6,*)  "ASYMM T2AA"
        Write(6, "(5(F11.9,3x))") (WORK(istartaa+I-1), I=1, itmpaa)
        Write(6,*)
        Write(6,*) "ASYMM T2BB"
        write(6,*)
        Write(6, "(5(F11.9,3x))") (WORK(istartbb+I-1), I=1, itmpbb)
#endif
         ISTARTAA = ISTARTAA + DISSZAAI*NUMDZAAI
         ISTARTBB = ISTARTBB + DISSZBBI*NUMDZBBI

#ifdef _DEBUG_LVLM1
	 Write(6,*) NUMDZAAT, DISSZAAT, DISSZAAI
	 Write(6,*) NUMDZBBT, DISSZBBT, DISSZBBI
         Write(6,*) ISTARTAA, ISTARTBB, JSTARTAA, JSTARTBB
#endif

         DO IDIS = 1, NUMDZAAT 
            CALL SQSYM(IRREPAB, POP(1, 1), DISSZAAT, DISSZAAI, 1,
     &              WORK(KSTARTAA+KOFFAA), WORK(JSTARTAA+JOFFAA))  
           JOFFAA = JOFFAA + DISSZAAI
           KOFFAA = KOFFAA + DISSZAAT
           
         ENDDO
         DO IDIS = 1,  NUMDZBBT
            CALL SQSYM(IRREPAB, POP(1, 2), DISSZBBT, DISSZBBI, 1,
     &              WORK(KSTARTBB+KOFFBB), WORK(JSTARTBB+JOFFBB))  
           JOFFBB = JOFFBB + DISSZBBI
           KOFFBB = KOFFBB + DISSZBBT
         ENDDO
C
#ifdef _DEBUG_LVLM1
        itmpaa =  Irpdpd(irrepab, 1)*Irpdpd(irrepab, 3)
        itmpbb =  Irpdpd(irrepab, 2)*Irpdpd(irrepab, 4)
        Write(6,*)
        write(6,*)  "SQSYM T2AA"
        Write(6, "(5(F11.9,3x))") (WORK(KSTARTAA+I-1), I=1, itmpaa)
        Write(6,*)
        Write(6,*) "SQSYM T2BB"
        write(6,*)
        Write(6, "(5(F11.9,3x))") (WORK(KSTARTBB+I-1), I=1, itmpbb)
#endif
C         
         CALL TRANSP(WORK(KSTARTAA), WORK(LSTARTAA), NUMDZAAT,
     &               DISSZAAT)
         CALL TRANSP(WORK(KSTARTBB), WORK(LSTARTBB), NUMDZBBT,
     &               DISSZBBT)
C
         JSTARTAA = JSTARTAA + NUMDZAAI*DISSZAAI
         JSTARTBB = JSTARTBB + NUMDZBBI*DISSZBBI 
C
         KSTARTAA = KSTARTAA + DISSZAAT*NUMDZAAT
         KSTARTBB = KSTARTBB + DISSZBBT*NUMDZBBT
C
         LSTARTAA = LSTARTAA + DISSZAAT*NUMDZAAT
         LSTARTBB = LSTARTBB + DISSZBBT*NUMDZBBT
         
         JOFFAA = 0
         JOFFBB = 0
         KOFFAA = 0
         KOFFBB = 0

      ENDDO
#ifdef _DEBUG_LVLM1
        itmpaa =  ISYMSZ(1, 3)
        itmpbb =  ISYMSZ(2, 4)
        itmpcc =  ISYMSZ(13, 14)
       
        Write(6,*) "Final T2AA"
        write(6,*)
        Write(6, "(5(F11.9,3x))") (WORK(I080+I-1), I=1, itmpaa)
        Write(6,*)
        Write(6,*) "Final T2BB"
        write(6,*)
        Write(6, "(5(F11.9,3x))") (WORK(I090+I-1), I=1, itmpbb)
        Write(6,*)
#endif
C
      CALL PUTALL(WORK(I080), ISQTOTAA, 1, 44)
      CALL PUTALL(WORK(I090), ISQTOTBB, 1, 45)
      CALL PUTALL(WORK(I050), ITOTALAB, 1, 46)
C
#ifdef _DEBUG_LVLM1
        itmpaa =  ISYMSZ(1, 3)
        itmpbb =  ISYMSZ(2, 4)
        itmpcc =  ISYMSZ(13, 14)
       
        Write(6,*) "Final T2AA"
        write(6,*)
        Write(6, "(5(F11.9,3x))") (WORK(I080+I-1), I=1, itmpaa)
        Write(6,*)
        Write(6,*) "Final T2BB"
        write(6,*)
        Write(6, "(5(F11.9,3x))") (WORK(I090+I-1), I=1, itmpbb)
        Write(6,*)
        Write(6,*) "final  T2AB"
        write(6,*)
        Write(6, "(5(F11.9,3x))") (WORK(I050+I-1), I=1, itmpcc)
#endif
C
#ifdef _RCC
      Write(6,*)
      Write(6, "(a)") "The antisymmetrized r-CCD correlation energy"
      Write(6,*)
#else
      Write(6,*)
      Write(6, "(a)") "The antisymmetrized dr-CCD correlation energy"
      Write(6,*)
#endif
      CALL CMPENG(WORK, MAXCOR, 43, 0, ECORR, ETOT, ETOTT2, 1,
     &            1)

      CALL DDMPTGSS(WORK, MAXCOR, 1, 0, 'TGUESS  ')

      RETURN
      END
