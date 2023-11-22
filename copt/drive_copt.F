#include "flags.h"

      Subroutine Drive_copt(Icore, Maxcor, Iuhf)
      
      Implicit Double Precision (A-H, O-Z)

      Dimension Icore(Maxcor)
      Logical COPT,BRKCNV,BRUECK,SING1,AOLOG,ORDER2
      Logical MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,
     &                QCISD,LAMBDA,NONHF,ROHF4,ITRFLG,UCC
      Integer T1_CYCLE
C
      COMMON /AOLOG/ AOLOG
      COMMON/SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &              QCISD,UCC
#include "flags.com"
#include "flags2.com"
#include "machsp.com"
#include "sym.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
#include "istart.com"

#ifdef _DEBUG_LVL0
      CALL AMPSUM(ICORE(1),MAXCOR,IUHF,0,.TRUE.,'T')
#endif
C
C Standard EOM-MBPT(2) ends here. Start the PEOM_X. What is it. The
C first-step of it to solve the Brueckner MBPT(2). That is equivalent
C to solving Q2H(bar,T2(1)) = 0.0 Q1Hbar(T2(1),T1)=0.
C 
      CALL SETMET
      CALL SETLOG
C
      write(6,"(3x,a,a)") "Perform a canonical (not parametrized) CC ",
     +                 "like calculations."
      write(6,*)
      Fae_scale    = 0.0D0
      Fmi_scale    = 0.0D0
      Wmnij_scale  = 1.0D0
      Wmbej_scale  = 1.0D0
      Gae_scale    = 1.0D0
      Gmi_scale    = 1.0D0

      AOLOG  = (Iflags(h_IFLAGS_abcdtype)  .EQ. 2) 
      BRUECK = (Iflags(h_IFLAGS_brueckner) .EQ. 1)
      COPT   = (Iflags2(h_IFLAGS2_copt)    .EQ. 1 .AND.
     &          Iflags(h_IFLAGS_calclevel) .EQ. 1) 

      IF (COPT) THEN
C
C Initialize the singles lists (list 90; 1,2,3,4; 3 and 4 are
C T1 increments). The SETLST form all Wmbej intermediate lists.
C SING1 does not need to be .TRUE. unless T1 guess is read 
C from somwhere else. The INMBPT at this point is not necessary
C (perhaps) but sure does do no harm. 
C
         SING1 = .FALSE.
C
C First Initialize some list that do not get initialized by
C intilization calls (these are normally done INSING)
C      
         DO ISPIN = 1, IUHF+1
             CALL UPDMOI(1,NFMI(ISPIN),ISPIN,  91,0,0)
             CALL UPDMOI(1,NFEA(ISPIN),ISPIN,  92,0,0)
             CALL UPDMOI(1,NT(ISPIN),  ISPIN,  93,0,0)
         END DO
C
C For non-HF methods, first-order T1 (= f(a,i)/ea-ei) is stored in
C list 90 (1,2) location. The firs-order T2 is computed for all the
C references. We simply built the 2nd order T1 
C 
	 CALL UPDMOI(1,NT(1),3,90,0,0)
         IF (IUHF.EQ.1) THEN
            CALL UPDMOI(1,NT(2),4,90,0,0)
         ENDIF
C
         CALL INMBPT(ICORE,MAXCOR,IUHF)
         CALL RNABIJ(ICORE,MAXCOR,IUHF,'T')
         CALL SETLST(ICORE,MAXCOR,IUHF)
C           
         DIFF     = 0.1D0
         T1_CYCLE = 0
C
C Iterate the T1 equation for a fixed T2 untill the desired 
C convergence is achived.
C
         TOL = 10.0**(-IFLAGS(h_IFLAGS_cc_conv))

         ORDER2    =  (Iflags(h_IFLAGS_calclevel) .EQ. 1)

         IF (ORDER2) THEN 

#ifdef _DEBUG_LVL0
      write(6,"(2x,a,a)") " The COPT-MBPT(2) T1 contributions are", 
     +                 " being computed"
#endif
            call aces_list_memset(3,90,0)
            call aces_list_memset(1,91,0)
            call aces_list_memset(1,92,0)
            call aces_list_memset(1,93,0)
C            
            call zersym(icore(inext),51)
            call zersym(icore(inext),52)
            call zersym(icore(inext),54)
            call zersym(icore(inext),56)
            call zersym(icore(inext),58)

            if (iuhf .ne. 0) then
               call aces_list_memset(4,90,0)
               call aces_list_memset(2,91,0)
               call aces_list_memset(2,92,0)
               call aces_list_memset(2,93,0)
               call zersym(icore(inext),53)
               call zersym(icore(inext),55)
               call zersym(icore(inext),57)
               call zersym(icore(inext),59)
            endif

            CALL COMPUTE_T1_4T2(ICORE,MAXCOR,IUHF,DIFF)
C

#ifdef _DEBUG_LVL0
      write(6,"(1x,a)") " The COPT-MBPT(2) T1 contributions computed"
#endif 
         ELSE

         DO WHILE (DIFF .GT. TOL .AND.
     &             T1_CYCLE .LE. IFLAGS(h_IFLAGS_cc_maxcyc))
          
            call aces_list_memset(3,90,0)
            call aces_list_memset(1,91,0)
            call aces_list_memset(1,92,0)
            call aces_list_memset(1,93,0)
            call aces_list_memset(3,93,0)
      
            call zersym(icore,51)
            call zersym(icore,52)
            call zersym(icore,54)
            call zersym(icore,56)
            call zersym(icore,58)

            if (iuhf .ne. 0) then
                call aces_list_memset(4,90,0)
                call aces_list_memset(2,91,0)
                call aces_list_memset(2,92,0)
                call aces_list_memset(2,93,0)
                call aces_list_memset(4,93,0)
                call zersym(icore,53)
                call zersym(icore,55)
                call zersym(icore,57)
                call zersym(icore,59)
            endif

            CALL COMPUTE_T1_4T2(ICORE,MAXCOR,IUHF,DIFF)
            T1_CYCLE = T1_CYCLE + 1
     
         ENDDO 

         IF  (T1_CYCLE .GT. IFLAGS(h_IFLAGS_cc_maxcyc) .OR. 
     &        DIFF .GT. IFLAGS(h_IFLAGS_cc_conv)) THEN
              Write(6,"(A,A)") "T1 iterations did not converged",
     &                         " in alloted number of cyclces."
              call ERREX
         ELSE

              Write(6,"(A,1x,i3,1x,A)")"T1 iterations converged in",
     &                           T1_CYCLE, "cyclces."
         ENDIF

         ENDIF
C
C It is important to have the Brueckner key-word set. So, the test
C Brueckner is not absolutely necessary. 
C
         IF (BRUECK) THEN
            CALL BRUECKIT(ICORE,MAXCOR,IUHF,BRKCNV)
         END IF
            
         IF (BRKCNV) THEN
            CALL INMBPT(ICORE,MAXCOR,IUHF)
            CALL INSING(10,IUHF)
            CALL SETLST(ICORE,MAXCOR,IUHF)
            CALL INITIN(ICORE,MAXCOR,IUHF)
            CALL INITSN(ICORE,MAXCOR,IUHF)
            CALL INCOR(I0,ICRSIZ,IUHF)
            MAXCOR = ICRSIZ
            CALL RNABIJ(ICORE,MAXCOR,IUHF,'T')
            CALL GENINT(ICORE,MAXCOR,IUHF,2,DO_HBAR_4LCCSD)
            CALL INIT2(IUHF)

C The FEACONT_MODF, FMICONT_MODF and FMECONT_MODF are simply
C modifications to the orgianal routines (as the name indicates).
C They compute only FT1 + FT2 to T1 equation. The original 
C also compute 

            CALL FEACONT_MODF(ICORE,MAXCOR,IUHF)
            CALL FMICONT_MODF(ICORE,MAXCOR,IUHF)
            CALL FMECONT(ICORE,MAXCOR,IUHF,1)
C 
C This will compute the WT2 to the T1 equation. 


C Built Hbar(p,q) and Hbar(mb,ej) and Hbar(mn,ij) with the 
C MBPT(2) T2. The Hbar(ij,ka), Hbar(ab,ci) and Hbar(ab,cd) is
C handlled in lambda code.
C Hbar(p,q)  :  f(p,q) + T1*W + Tau2*W2 terms.
C Hbar(pq,rs):  W2 + T1*W + Tau*W2 terms (note that during LCCD
C iteration only W2 terms are invloved).
C
C This may not be necessary; but safer than being sorry. All
C the (mb,ej) intermediates are set to zero inside genint. So,
C it is not done here.

            IF (IUHF .NE. 0) THEN
               CALL ZERSYM(ICORE,51)
               CALL ZERSYM(ICORE,52)
            ELSE
               CALL ZERSYM(ICORE,53)
            ENDIF
C
            CALL GENINT(ICORE,MAXCOR,IUHF,2,DO_HBAR_4LCCSD)
            IF (IUHF .EQ. 0) CALL RESET(ICORE,MAXCOR,IUHF)
C
C I do not know whether these are necessry. There is no L1 and L2
C for LCCSD. I believe our EOM codes looks for these lists for
C the left hand vectors. So, this is simply a trick to get the
C left side (may be gradients to work).

             DO ISPIN = 3, 3-2*IUHF, -1
                CALL MVSYM(ICORE,MAXCOR,IUHF,1,43+ISPIN,60+ISPIN)
                CALL MVSYM(ICORE,MAXCOR,IUHF,1,43+ISPIN,143+ISPIN)
             ENDDO
             DO ISPIN=1,IUHF+1
                CALL UPDMOI(1,NT(ISPIN),ISPIN,190,0,0)
                CALL UPDMOI(1,NT(ISPIN),ISPIN+2,190,0,0)
             END DO

             CALL GETLST(ICORE, 1, 1, 1, 1, 90)
             CALL PUTLST(ICORE, 1, 1, 1, 1, 190)

             IF (IUHF .NE. 0) THEN
                 CALL GETLST(ICORE, 1, 1, 1, 2, 90)
                 CALL PUTLST(ICORE, 1, 1, 1, 2, 190)
             ENDIF
         ENDIF 
C
      ENDIF 

      RETURN
      END

