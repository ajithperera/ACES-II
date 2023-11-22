      Subroutine Setup_4nos(Work,Maxcor,Irrepx,Ioffset,Iside,IUhf)

      Implicit Double Precision(A-H,O-Z)
      Dimension Work(Maxcor)
      Dimension Length(8)
      Dimension Length_sd(8,2,3)

#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "sympop.com"
  
      Data Ione /1/
  
      DO IRREP = 1, NIRREP
        LENGTH(IRREP) = 0
        DO ISPIN = 1, 1 + IUHF
                        LENGTH(IRREP) = LENGTH(IRREP) +
     &                                  IRPDPD(IRREP,8+ISPIN)
           LENGTH_SD(IRREP,1,ISPIN) = IRPDPD(IRREP,8+ISPIN)
        END DO

        DO ISPIN = 3, 3-2*IUHF, -1
                               LSTTYP = 43 + ISPIN
                        LENGTH(IRREP) = LENGTH(IRREP) +
     &                                  IDSYMSZ(IRREP,
     &                                  ISYTYP(1,LSTTYP),
     &                                  ISYTYP(2,LSTTYP))
           LENGTH_SD(IRREP,2,ISPIN) = IDSYMSZ(IRREP,
     &                                ISYTYP(1,LSTTYP),
     &                                ISYTYP(2,LSTTYP))
        END DO
      END DO

      LISTPT2   = 372+Iside
      LISTR1    = 490
      LISTR2_AA = 444
      LISTR2_BB = 445
      LISTR2_AB = 446

      IOFF_T1_A  = IONE
      IOFF_T2_AB = IONE + (LENGTH_SD(IRREPX,1,1))*IINTFP
      IOFF_T2_AA = IONE + (LENGTH_SD(IRREPX,1,1) +
     &                     LENGTH_SD(IRREPX,2,1))*IINTFP

      IF (IUHF .NE. 0) THEN
          IOFF_T1_B  = IONE + (LENGTH_SD(IRREPX,1,1))*IINTFP
          IOFF_T2_AB = IONE + (LENGTH_SD(IRREPX,1,1) + 
     &                         LENGTH_SD(IRREPX,1,2))*IINTFP
          IOFF_T2_AA = IONE + (LENGTH_SD(IRREPX,1,1) +
     &                         LENGTH_SD(IRREPX,1,2) + 
     &                         LENGTH_SD(IRREPX,2,3))*IINTFP
          IOFF_T2_BB = IONE + (LENGTH_SD(IRREPX,1,1) +
     &                         LENGTH_SD(IRREPX,1,2) + 
     &                         LENGTH_SD(IRREPX,2,3) +
     &                         LENGTH_SD(IRREPX,2,2))*IINTFP
      ENDIF
C    
      ILOC_PERT = IONE
      IEND_PERT = ILOC_PERT + LENGTH(IRREPX)
      MEMLEFT   = (MAXCOR - IEND_PERT)
C
      IF (IEND_PERT .GE. MAXCOR) CALL INSMEM("setup",ILOC_PERT,
     &                                        MAXCOR)
      CALL GETLST(WORK(ILOC_PERT),IOFFSET,1,1,IRREPX,LISTPT2)

      RNORM = DDOT(LENGTH(IRREPX),WORK(ILOC_PERT),1,
     &             WORK(ILOC_PERT),1)

      Write(6,"(a,F15.9)") " The Rnorm: ", Rnorm 
      ILEN_T1_A  = LENGTH_SD(IRREPX,1,1)
      ILEN_T2_AA = LENGTH_SD(IRREPX,2,1)
      ILEN_T2_AB = LENGTH_SD(IRREPX,2,3)

      IF (IUHF .NE. 0) THEN
         ILEN_T1_B  = LENGTH_SD(IRREPX,1,2)
         ILEN_T2_BB = LENGTH_SD(IRREPX,2,2)
      ENDIF

      CALL PUTLST(WORK(IOFF_T1_A),1,1,1,1,LISTR1)
      IF (IUHF .NE. 0) CALL PUTLST(WORK(IOFF_T1_B),1,1,1,2,LISTR1)

      DO IRREP = 1, NIRREP
         IRREPR = IRREP
         IRREPL = DIRPRD(IRREPR, IRREPX)

         LENGTH_OOAA = IRPDPD(IRREPR, 3)
         LENGTH_VVAA = IRPDPD(IRREPL, 1)
         LENGTH_AAAA = LENGTH_OOAA*LENGTH_VVAA
         IF (IUHF .NE . 0) THEN
             LENGTH_OOBB = IRPDPD(IRREPR, 4)
             LENGTH_VVBB = IRPDPD(IRREPL, 2)
             LENGTH_BBBB = LENGTH_OOBB*LENGTH_VVBB
         ENDIF

         LENGTH_OOAB = IRPDPD(IRREPR, 14)
         LENGTH_VVAB = IRPDPD(IRREPL, 13)
         LENGTH_ABAB = LENGTH_OOAB*LENGTH_VVAB

         CALL PUTLST(WORK(IOFF_T2_AA),1,LENGTH_OOAA,1,IRREPR,
     &               LISTR2_AA)
         CALL PUTLST(WORK(IOFF_T2_AB),1,LENGTH_OOAB, 1,IRREPR,
     &               LISTR2_AB)
         IF (IUHF .NE. 0) CALL PUTLST(WORK(IOFF_T2_BB),1,LENGTH_OOBB,
     &                                1,IRREPR,LISTR2_BB)

         IOFF_T2_AA = IOFF_T2_AA +  LENGTH_AAAA
         IOFF_T2_AB = IOFF_T2_AB +  LENGTH_ABAB
         IF (IUHF .NE. 0) IOFF_T2_BB = IOFF_T2_BB + LENGTH_BBBB

      END DO

      RETURN
      END
   
      
