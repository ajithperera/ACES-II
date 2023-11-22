      SUBROUTINE DG6TOAO(ICORE, MAXCOR, IUHF, ISPIN,
     &                   LISTMO, NAO, IRREPX)
C
C TRANSFORMS THE G(AB,CI) AMPLITUDES TO THE AO-BASIS
C
C    G(AB,CI)   --> G(MU NU,SIGMA I) 
C
CEND
C Once again this is a modified version of a routine written
C by Jurgen Gauss. Modifications in order to adapt to ACES II 
C (Gnv) version was done by A. Perera, 05/2005. 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP,VRT
C
      DIMENSION ICORE(MAXCOR)
      LOGICAL EXIT, DONE
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),IOFFV(8,2),IOFFO(8,2),
     &             IRPDPDAO(8),IRPDPDAOS(8),ISTART(8,8),ISTARTMO(8,3)
      COMMON/AOMOMIX/IRPDPDAM(8,2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/SYMLOC/ISYMOFF(8,8,25)
      COMMON/G6OOC/IPOPS(8,8),IPOPE(8,8),NDONE(8,8),NSTARTM,
     &             ISTARTI,NSIZG,IOFFSET,DONE
C
      IONE=1
C
C NMO:    NUMBER OF MOS
C NAO:    NUMBER OF AOS USED IN THE CONTRACTION
C NAOMAX: MAXIMUM NUMBER OF AOS IN ONE IRREP
C
      NMO   = NOCCO(1) + NVRTO(1)
      NAOMAX= 0
C
      DO IRREP = 1, NIRREP
         NAOMAX = MAX(NAOMAX, IAOPOP(IRREP))
      ENDDO
C
C Loop over symmetries and form G(munu,lambda,i) and the IOFF is
C the offset for the said quantitiy. 
C
      NSIZG  = 0
      NREADT = 0
      IDONE  = 0
      IOFF   = 1
      EXIT   = .FALSE.
C
      LISTG = LISTMO + ISPIN
C
C loop over irreps of right hand side, X refer to the AO index.
C
      DO 100 IRREPXI = ISTARTI, NIRREP
C
         IRREPAB = DIRPRD(IRREPXI, IRREPX)
         IRREPXX = IRREPAB
         IRREPAI = IRREPXI
C
         ISTARTD = 1
C
         NUMAB = IRPDPD(IRREPAB, ISYTYP(1,LISTG))
         NUMAI = IRPDPD(IRREPAI, ISYTYP(2,LISTG))
C
         NUMXX = IRPDPDAO(IRREPXX)
C
         IF (ISPIN .LE. 2) THEN
             NUMXI = IRPDPDAM(IRREPXI, ISPIN)
             NUMABX= IRPDPD(IRREPAB, 18+ISPIN) 
         ELSE
             NUMXI = IRPDPDAM(IRREPXI,ISPIN-2)
             NUMABX= NUMAB
         ENDIF
C
C NAO, NMO, NUMABX, NUMXI, NUMAB and NUMAI (NUMXI=NUMAI) are 
C total number of AO, total number of MOs, number of MO pairs 
C in irrep AB (expanded), number of MO pairs in irrep AI (expanded)
C number of MO pairs in irrepa (A<B when applies) and number of 
C MO pairs in irrepi. NUMXX is the number of AO pairs in irrepxx. 
 
         NUMSET = 2
         INEED  = IINTFP*(NUMSET*NAO*NMO + NAO*NAO +
     &            NUMABX*NAOMAX)
         MEMLEFT= (MAXCOR-IOFF-INEED)/IINTFP
C
         MEM1 = 0
         MEM2 = 0
         NUMD = 0
C
C Loop over the slowest index I of G(AB|CI)
C
         DO 10 IRREPI =1, NIRREP
C
            IF (ISPIN .LE. 2) THEN
              ISPINI = ISPIN
              ISPINA = ISPIN
            ELSE
              ISPINI = ISPIN - 2
              ISPINA = 3 - ISPINI
            ENDIF
C
C Skip this irrep after updating ISTARTD if we are done 
C with this irrep. 
C
            IF (NDONE(IRREPI,IRREPAI) .EQ. 1) THEN
               ISTARTD = ISTARTD + POP(IRREPI,ISPINI)*
     &                           VRT(DIRPRD(IRREPAI,IRREPI),
     &                           ISPINA)
               GO TO 10 
            ENDIF
C
            NUMI   = POP(IRREPI,ISPINI)
            IRREPA = DIRPRD(IRREPAI,IRREPI)
            NUMA   = VRT(IRREPA,ISPINA)
            NUMX   = IAOPOP(IRREPA)
C
C NUMXXX: mu*nu*lambda and NUMAAA=A*B*C (expanded).
C
            NUMXXX = NUMXX*NUMX
            NUMAAA = NUMABX*NUMA
C
C First distribution to be treated (offset for the CI s.). The
C NSTARTM refer to I and NUMA is all Cs for given I.  
C
           ISTARTD = ISTARTD + (NSTARTM - 1)*NUMA
C
C How many Is can we handle with the available memory. NREAD
C is the number of Is can be handled in the irrep I. The 
C (NUMXXX + NUMAAA) can only be zero when both NUMXXX and 
C NUMAAA are individualy zero. That means there are no mu*nu*
C lambda and A*B*Cs in this irrep. 
C
           IF ((NUMXXX + NUMAAA) .NE. 0) THEN
C
C I read the minimum number of I s possible.
C
                NREAD = MIN(NUMI+1-NSTARTM,MEMLEFT/
     &                     (NUMXXX+NUMAAA))
           ELSE
C
C I read whatever left to do
C
                NREAD = NUMI + 1 - NSTARTM
           ENDIF
C
C Increment counter NREADT
C
           NREADT = NREADT + NREAD
C
C NREAD is zero either beacuse all CIs has been read or there is 
C not enough memory. Either case we need to get out of here. If
C we get in here because we can not keep everything along with
C all ABC for given occupied orbital, we simply do not have enough
C memory, we transfer the control to 20, and then get out with
C error message saying there is not enough memory (see the call 
C to INSMEM below).  
C
           IF (NREAD.EQ.0 .AND. NUMI.NE.0) THEN
               EXIT = .TRUE.
               GO TO 20 
C
C If we get out from here to 20, NUMD=0, but we do not
C go that far in the code to test for NUMD=0 test (see below),
C because we simply do not have memory to do anuthing. 
C
           ENDIF
C
C IPOPS is the offset for I (occupied) and IPOPE is the offset for 
C C (virtual).
C
           IPOPS(IRREPI, IRREPXI) = NSTARTM
           IPOPE(IRREPI, IRREPXI) = NSTARTM + NREAD -  1
C
C NUMD is the number of CI distributions read (NREAD correspond 
C to I and NUMA correspond to C). 
C
           NUMD    = NUMD + NUMA*NREAD
           MEMLEFT = MEMLEFT - NREAD*(NUMXXX+NUMAAA)
C
C MEM1 and MEM2: length of A,B,C,I and mu,nu,lam,I
C
           MEM1 = MEM1 + NUMAAA*NREAD
           MEM2 = MEM2 + NUMXXX*NREAD
C
C We have more CI distributions to read. Get out of this loop
C and read G(ab,ci)s and process them and then return to complete
C this block. Note that we return to the caller before we come
C back. When the condition is true we have more to read in this
C spin type, IRREPXI and IRREPI loops. Update the NSTARTM and 
C return the call back to spin-type loop after reading a block 
C of A*B*NUMD of G(AB,CI). 
C
           IF ((NSTARTM + NREAD - 1) .NE. NUMI) THEN
                NSTARTM = NSTARTM + NREAD
                EXIT    = .TRUE.
                GO TO 20
           ENDIF
C
C If we come here we are done for this spin-type and IRREPI. So,
C start the next IRREPI, after tagging the NDONE array and updating
C NSTARTM go back to the next value of IRREPI. 
C
           NDONE(IRREPI, IRREPAI)= 1
           NSTARTM               = 1
C
 10      CONTINUE
 20      CONTINUE
C
         NSIZG = NSIZG + MEM2
C
         I000 = IOFF
         I005 = I000 + MEM2*IINTFP
         I010 = I005 + MEM1*IINTFP
         I020 = I010 + NAO*NMO*IINTFP
         I030 = I020 + NAO*NMO*IINTFP
         I040 = I030 + NUMABX*NAOMAX*IINTFP
         IEND = I040 + NAO*NAO*IINTFP 
C
         IF (ISPIN .EQ.3 .OR. IUHF .EQ.0) THEN
            ITMP1 = I010
            ITMP2 = ITMP1 + IINTFP*MAX(NUMAB,NUMD)
            ITMP3 = ITMP2 + IINTFP*MAX(NUMAB,NUMD)
            IEND2 = ITMP3 + IINTFP*MAX(NUMAB,NUMD)
            IEND  = MAX(IEND, IEND2)
         ENDIF

         IF (IEND .GE. MAXCOR) THEN
             IF (NREADT .EQ. 0) CALL INSMEM("@-DG6TOAO",IEND,MAXCOR)
         ENDIF 
C
C If there are no more CI distributions left to do in this spin 
C type and more spin types left, then go back to the caller (that has
C the loop over spin types) otherwise continue the loop over IRREPXI
C for the next IRREPXI starting from (ISTARTI + 1). 
C NUMD = 0 EXIT=TRUE: go back the loop over spin-type block
C NUMD = 0 EXIT=FALSE: Go back to the loop over IRREPXI
C NUMD != 0: READ!!!!!
C
         IF (NUMD .EQ. 0) THEN
            CALL ZERO(ICORE(I000), MEM2)
            IOFF = IOFF + IINTFP*MEM2
            IF( .NOT. EXIT) ISTARTI = ISTARTI + 1
            IF (EXIT) RETURN
            GO TO 100
         ENDIF 
C
C Read in G(ab,ci) amplitudes 
C
         IF (ISPIN .NE. 3) THEN
C
             CALL GETLST(ICORE(I005),ISTARTD,NUMD,1,IRREPAI,
     &                   LISTG)
C
         ELSE
C
C Treat the ISPIN=3 seperately. The reason is the its left hand 
C side is stored as CI instead if IC as in the other 3 cases. 
C 
             IOFFF = I005
             NDIS  = 0
C
             DO IRREPI = 1, NIRREP
                NUMI  = POP(IRREPI, ISPINI)
                IRREPA= DIRPRD(IRREPI, IRREPAI)
                NUMA  = VRT(IRREPA, ISPINA)
C
                DO I = IPOPS(IRREPI,IRREPXI),IPOPE(IRREPI,IRREPXI)
                   DO IA = 1, NUMA
                      IDIS = ISYMOFF(IRREPA,IRREPAI,18) + 
     &                       (IA-1)*NUMI+I-1
                      CALL GETLST(ICORE(IOFFF),IDIS,1,1, IRREPAI, 
     &                            LISTG)
                      IOFFF = IOFFF + NUMABX*IINTFP
                      NDIS  = NDIS + 1
                   ENDDO
                ENDDO
             ENDDO
C
             IF (NUMD .NE. NDIS) THEN
                 WRITE(*,*) ' Something is very wrong '
                 CALL ERREX
             ENDIF 
C
C Interchange first two and last two indices for ispin is 3. 
C
             CALL SYMTR3(IRREPAB,VRT(1,1),VRT(1,2),NUMABX,NUMD,
     &                   ICORE(I005),ICORE(ITMP1),ICORE(ITMP2),
     &                   ICORE(ITMP3))
C
         ENDIF
C
C Expand left side for ispin 1 or 2.
C
        IF (ISPIN .LE. 2) THEN
C
          CALL SYMEXP2(IRREPAB,VRT(1,ISPIN),NUMABX,NUMAB,
     &                 NUMD,ICORE(I005),ICORE(I005))
        ENDIF
C
C Spin adaptation for RHF
C
        IF (IUHF .EQ. 0) THEN
C
           CALL SPINAD3(IRREPAB,VRT(1,1),NUMAB,NUMD,ICORE(I005),
     &                  ICORE(ITMP1),ICORE(ITMP2))
C
        ENDIF
C
          CALL DG6TRAN(ICORE(I005),ICORE(I000),ICORE(I010),
     &                 ICORE(I020),ICORE(I030),ICORE(I040),
     &                 NAO,NMO,NUMABX,NUMD,NUMXX,NUMXI,
     &                 ISPIN,IRREPAI,IRREPAB,IUHF)
C
         IOFF = IOFF + IINTFP*MEM2
C
         IF( .NOT. EXIT) ISTARTI = ISTARTI + 1
C
         IF (EXIT) RETURN
C
 100  CONTINUE 
C
      IF (.NOT. EXIT) DONE = .TRUE.
C
      RETURN
      END       
