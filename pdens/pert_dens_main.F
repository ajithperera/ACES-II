      SUBROUTINE PERT_DENS_MAIN(ICORE, MAXCOR, IUHF)
C
#include "flags.h"
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
      INTEGER POP, VRT, DIRPRD
      LOGICAL SPN_SPN, POLAR, JFC, JPSO, JSD
C
#include "machsp.com"
#include "mxatms.par"
#include "jodaflags.com"
C
      DIMENSION ICORE(MAXCOR), LENGTH(8), LENGTH_SD(8, 2, 3), 
     &          POL_TENSOR(3, 3), SPNC_TENSOR(6*MXATMS*6*MXATMS),
     &          NPERT(8)
C
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /SPN_SPN_TYPE/ JFC, JPSO, JSD
C
      DATA IONE /1/      
C
#ifdef _MKLST
      CALL MAKE_LISTS4_PERTRB_DENS(ICORE, MAXCOR, IUHF)
#endif 
      DO IRREP = 1, NIRREP    
C
        LENGTH(IRREP) = 0
C
        DO ISPIN = 1, 1 + IUHF
                        LENGTH(IRREP) = LENGTH(IRREP) + 
     &                                  IRPDPD(IRREP, 8 + ISPIN)
           LENGTH_SD(IRREP, 1, ISPIN) = IRPDPD(IRREP, 8 + ISPIN)
        END DO
C
        DO ISPIN = 3, 3-2*IUHF, -1
                               LSTTYP = 43 + ISPIN
                        LENGTH(IRREP) = LENGTH(IRREP) + 
     &                                  IDSYMSZ(IRREP,
     &                                  ISYTYP(1, LSTTYP), 
     &                                  ISYTYP(2, LSTTYP))
           LENGTH_SD(IRREP, 2, ISPIN) = IDSYMSZ(IRREP,
     &                                  ISYTYP(1, LSTTYP), 
     &                                  ISYTYP(2, LSTTYP))
        END DO    
      END DO
C 
#ifdef _DEBUG_LVLM2
      print*, "Length of the perturbed T2 vector"  
      print*, (length(i), i=1, nirrep)
      print*, "Length of the perturbed singles"
      print*, ((length_sd(irrep, 1, i), irrep =1, nirrep), 
     &          i=1, iuhf+1)
      print*, "Length of perturrbed doubles"
      print*, ((length_sd(irrep, 2, i), irrep =1, nirrep), 
     &          i=3, 3-2*IUHF, -1 )
#endif
#ifdef _TDENS
      CALL GETREC(0, 'JOBARC', 'NTOTPERT', IRECLEN, NTPERT)  
      IF (IRECLEN .GT. 0) THEN
          CALL GETREC(20, 'JOBARC', 'NTOTPERT', 1, NTPERT)
          LISTPT2 = 373 
          CALL ZERO(SPNC_TENSOR, 9*MXATMS*MXATMS)
      ELSE
          NTPERT  = 3
          LISTPT2 = 373
          CALL ZERO(POL_TENSOR, 9)
      ENDIF 
      CALL GETREC(-1, 'JOBARC', 'NPERTIRR', 8, NPERT)
          LISTR1 = 490
       LISTR2_AA = 444 
       LISTR2_BB = 445
       LISTR2_AB = 446
#endif
#ifdef _DEBUG_LVLM
       print*, "The total number of perturbations and symmetries", 
     &          NTPERT
       Print*,  (NPERT(I), I=1, NIRREP)
       Print*, "The perturbed amplitudes are read from list:", 
     &          LISTPT2
#endif
#ifdef _TEST1
         LISTPT2 = 374
          LISTR1 = 90      
       LISTR2_AA = 44
       LISTR2_BB = 45
       LISTR2_AB = 46
#endif
C
      DO IRREPX = 1, NIRREP
C
         DO IPERT = 1,  NPERT(IRREPX)

            IOFF_T1_A  = IONE
            IOFF_T2_AB = IONE + (LENGTH_SD(IRREPX, 1, 1))*IINTFP
            IOFF_T2_AA = IONE + (LENGTH_SD(IRREPX, 1, 1) +
     &                           LENGTH_SD(IRREPX, 1, 2))*IINTFP

            IF (IUHF .NE. 0) THEN
              IOFF_T1_B  = IONE + (LENGTH_SD(IRREPX,1,1))*IINTFP
              IOFF_T2_AB = IONE + (LENGTH_SD(IRREPX,1,1) +
     &                             LENGTH_SD(IRREPX,1,2))*IINTFP
              IOFF_T2_AA = IONE + (LENGTH_SD(IRREPX,1,1) +
     &                             LENGTH_SD(IRREPX,1,2) +
     &                             LENGTH_SD(IRREPX,2,3))*IINTFP
              IOFF_T2_BB = IONE + (LENGTH_SD(IRREPX,1,1) +
     &                             LENGTH_SD(IRREPX,1,2) +
     &                             LENGTH_SD(IRREPX,2,3) +
     &                             LENGTH_SD(IRREPX,2,2))*IINTFP
            ENDIF
C          
            IOFFSET = IPERT
C
            ILOC_PERT = IONE
            IEND_PERT = ILOC_PERT + LENGTH(IRREPX)*IINTFP
              MEMLEFT = (MAXCOR - IEND_PERT)/IINTFP 
C
            IF (IEND_PERT .GE. MAXCOR) CALL INSMEM("@-WAVEFN_ANALYSE",
     &                                             ILOC_PERT, MAXCOR)
#ifdef _GETPTS
#ifdef _DEBUG_LVLM0
        print*, "Reading the perturbed amplitudes"
#endif
            cALL GETLST(ICORE(ILOC_PERT), IOFFSET, 1, 1, IRREPX,
     &                  LISTPT2)

            CALL GET_RNORM(IRREPX, LENGTH(IRREPX), ICORE(ILOC_PERT),
     &                     ICORE(IEND_PERT), MEMLEFT, IUHF, RNORM)
#endif
C
#ifdef _DEBUG_LVLM2
            nsizec = length(irrepx)
            print*, "The length of the vector: ", nsizec, ILOC_PERT,
     &      ioffset
            call checkSum("Comp_pert_dens", icore(iloc_pert), nsizec)
            print*, "The pert. T amps"
            call output(icore(iloc_pert),1,nsizec, 1, 1,nsizec,1,1)
#endif
C 
            ILEN_T1_A = LENGTH_SD(IRREPX, 1, 1)
            ILEN_T2_AA = LENGTH_SD(IRREPX, 2, 1)
            ILEN_T2_AB = LENGTH_SD(IRREPX, 2, 3)

            IF (IUHF .NE. 0) THEN
               ILEN_T1_B = LENGTH_SD(IRREPX, 1, 2)
              ILEN_T2_BB = LENGTH_SD(IRREPX, 2, 2)
             ENDIF
C
#ifdef _DEBUG_LVL0
            print*, "The singles length AI and ai: ",
     &      ilen_T1_A, ilen_T1_B
            write(6,*)
            print*, "The offsets for T 1s: ", Ioff_t1_a,
     &      ioff_t1_b
            print*, "The perturbed T1s AI"
            call output(icore(ioff_t1_a), 1,nt(1), 1,1,nt(1),1,1)
            call checksumP("T1AA", icore(ioff_t1_a),nt(1))
            print*, "The perturbed T1s ai"
            call output(icore(ioff_t1_b), 1,nt(2), 1,1,nt(2),1,1)
            call checksumP("T1BB", icore(ioff_t1_b),nt(2))
#endif
#ifdef _WPTS
            CALL PUTLST(ICORE(IOFF_T1_A), 1, 1, 1, 1, LISTR1)
            IF (IUHF .NE. 0) CALL PUTLST(ICORE(IOFF_T1_B), 1, 1, 1,
     &                                   2, LISTR1) 
#endif
     
            IOFF_T1_A = IOFF_T1_A + ILEN_T1_A*IINTFP    
            IF (IUHF .NE . 0) IOFF_T1_B = IOFF_T1_B + ILEN_T1_B*
     &                                    IINTFP
            DO IRREP = 1, NIRREP
               IRREPR = IRREP
               IRREPL = DIRPRD(IRREPR, IRREPX) 
#ifdef _DEBUG_LVLM2
               Print*, IRREPR, IRREPL,IRREPX
#endif

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
C
#ifdef _DEBUG_LVLM
              print*, "The doubles length, AA, BB and AB: ", 
     &        length_aaaa,  length_bbbb,  length_abab,length_ooaa,
     &        length_vvaa, length_oobb, length_vvbb,length_ooab,
     &        length_vvab
              write(6,*)
              Write(6,*) "The offsets for T2 s: ", ioff_t2_aa,
     &        ioff_t2_bb, ioff_t2_ab
              write(6,*)
              print*, "The perturbed T2 AA s: "
              call output(icore(ioff_t2_aa),1,length_ooaa, 1,
     &                    length_vvaa,length_ooaa,length_aa,1)
              write(6,*)
              call checksumP("T2AAAA", icore(ioff_t2_aa),LENGTH_AAAA)
              if (iuhf .eq. 1) then
              print*, "The perturbed T2 BB s: "
              call output(icore(ioff_t2_bb),1,length_oobb, 1,
     &                    length_vvbb,length_oobb,length_vvbb,1)
              write(6,*)
              call checksumP("T2BBBB", icore(ioff_t2_bb),LENGTH_BBBB)
              endif 
              print*, "The perturbed T2 AB s: "
              call output(icore(ioff_t2_ab),1,length_ooab, 1,
     &                    length_vvab,length_ooab,length_vvab,1)
              Write(6,*)
              call checksumP("T2ABAB", icore(ioff_t2_ab),LENGTH_ABAB)
#endif
C
#ifdef _WPTS
              CALL PUTLST(ICORE(IOFF_T2_AA), 1, LENGTH_OOAA, 1,
     &                    1, LISTR2_AA)
              CALL PUTLST(ICORE(IOFF_T2_AB), 1, LENGTH_OOAB, 1,
     &                    1, LISTR2_AB)
              IF (IUHF .NE. 0) CALL PUTLST(ICORE(IOFF_T2_BB), 1,
     &                                     LENGTH_OOBB, 1, 1,
     &                                     LISTR2_BB)
#endif
              IOFF_T2_AB = IOFF_T2_AB +  LENGTH_ABAB*IINTFP
              IOFF_T2_AA = IOFF_T2_AA +  LENGTH_AAAA*IINTFP
              IF (IUHF .NE. 0) IOFF_T2_BB = IOFF_T2_BB + 
     &                                      LENGTH_BBBB*IINTFP

            END DO
C
#ifdef _TEST
            CALL GENRSPNSDEN(ICORE, MAXCOR, IUHF) 
#endif
C
#ifdef _TDENS
#ifdef _DEBUG_LVLM2
            print*, "entering  GEN_PERTRB_RSPNS_DENS"
#endif
            CALL GEN_PERTRB_RSPNS_DENS(ICORE, MAXCOR, IRREPX, IUHF)
#endif
            CALL PRP_DEN4PROPS(ICORE, MAXCOR, NAO, NATOMS, NCENTER, 
     &                         IPERT, IRREPX, IUHF)
C
               I000 = IONE
              ILEFT = I000 + NAO*NAO*IINTFP
            MEMLEFT = MAXCOR - ILEFT
                                     
            CALL COMP_2NDORD_PROPS(ICORE(I000), ICORE(ILEFT), NTPERT,
     &                            MEMLEFT, NAO, NATOMS, IPERT, IUHF,
     &                            SPN_SPN, POLAR, SPNC_TENSOR, 
     &                            POL_TENSOR)
         ENDDO
      ENDDO
C
      IF (SPN_SPN) THEN
#ifdef _DEBUG_LVLM2
         CALL OUTPUT(SPNC_TENSOR, 1, NTPERT, 1, NTPERT, 
     &               NTPERT, NTPERT, 1)
#endif
         CALL SYMMET2(SPNC_TENSOR, NTPERT)
         CALL DSCAL(NTPERT*NTPERT, 2.0D0, SPNC_TENSOR, 1)
#ifdef _DEBUG_LVL0
         CALL OUTPUT(SPNC_TENSOR, 1, NTPERT, 1, NTPERT, 
     &               NTPERT, NTPERT, 1)
#endif
         CALL PRNT_SPNSPN_TNSR(SPNC_TENSOR, NTPERT, JFC, JSD, JPSO,
     &                         .FALSE., 1)
      ELSE IF (POLAR) THEN
CSSS         CALL SYMMET2(POL_TENSOR, 3) 
         CALL DSCAL(9, 2.0D0, POL_TENSOR, 1)
         CALL NEATPRINT(6, POL_TENSOR, 3, 3, 3, 3)
      ENDIF

      RETURN 
      END
















