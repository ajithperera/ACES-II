#include "flags.h"
       SUBROUTINE DRIVE_CORE_EE_STATE_SP(SCR,MAXCOR,IRREPX,IUHF,IROOT)
C
C The core EEs are obtained by setting the R(ij,ab) = 0 when i and
C j do  not belong to the core region. The present version works
C state specifc fashion.
C
C                           all; each
C LIST 444:    C(IJ,AB )       A<B ; I<J     AA AA
C      445:    C(ij,ab )       a<b ; i<j     BB BB
C      446:    C(Ij,Ab )       A,b ; I,j     AB A
C 
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION SCR(MAXCOR)
C
#include "maxbasfn.par"
C
      LOGICAL PROJECT_SINGLES 
      CHARACTER*12 STRINGI(2), STRINGJ(2), STRINGA(2),STRINGB(2)
      DIMENSION LS2OUT(2,2), IMAP_A(MAXBASFN),IMAP_B(MAXBASFN) 
      DIMENSION OCA_OFF(8),OCB_OFF(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON/ROOTS/EIGVAL(100,8),EIGVAL_T(100,8),OSCSTR(100,8),
     &             BGN(100,8),BGN_IRP(100,8),END(100,8),
     &             END_IRP(100,8),NATURE
      COMMON/EXCITATION/IORB_A(20),JORB_B(20),IORB_AA(20),JORB_AA(20),
     &                  IORB_BB(20),JORB_BB(20),IORB_AB(20),JORB_AB(20)
#include "info.com"
#include "flags.com"
#include "flags2.com"
C
      PROJECT_SINGLES = (Iflags2(h_IFLAGS2_prjct_single) .NE. 0)

C#ifdef _DEBUG_LVL0
      Write(6,"(a,a)") " Entering the core projection for specific", 
     &                 " core specified via %EXCITE*."
      Write(6,"(a,a)") " At entry, root and core indices: "
      Write(6,"(a)")   " RootID,IA,JB,IAA,JAA,IBB,JBB,IAB,JAB"
      Write(6,"(2x,9I4)")
     &         IROOT,IORB_A(IROOT),JORB_B(IROOT),
     &         IORB_AA(IROOT),JORB_AA(IROOT),IORB_BB(IROOT),
     &         JORB_BB(IROOT),IORB_AB(IROOT),JORB_AB(IROOT)
C#endif 
C
      ICHECK = MAX(IORB_AA(IROOT),JORB_AA(IROOT),IORB_BB(IROOT),
     +             JORB_BB(IROOT),IORB_AB(IROOT),JORB_AB(IROOT))

      IF (ICHECK .EQ. 0) THEN
         WRITE(6,"(2A)")" %EXCITE is used but no double amplitudes",
     +                  " are specified. Singles specification"
         WRITE(6,"(2A)")" is used to built a core projection scheme",
     +                  " for doubles."
         WRITE(6,*) 
         IF (IUHF .NE. 0) THEN
            IORB_AA(IROOT) = IORB_A(IROOT)
            JORB_AA(IROOT) = IORB_A(IROOT)
            IORB_BB(IROOT) = JORB_B(IROOT)
            JORB_BB(IROOT) = JORB_B(IROOT)
            IORB_AB(IROOT) = IORB_A(IROOT)
            JORB_AB(IROOT) = JORB_B(IROOT)
         ELSE
            IORB_AB(IROOT) = IORB_A(IROOT)
            JORB_AB(IROOT) = IORB_A(IROOT)
         ENDIF 
         
      ENDIF 

      LS1OUT      = 490
      LS2OUT(1,1) = 444
      LS2OUT(1,2) = 446
      LS2OUT(2,1) = 446
      LS2OUT(2,2) = 445
C
      STRINGI(2)='i  [i_SYM]  '
      STRINGJ(2)='j  [j_SYM]  '
      STRINGA(2)='a  [a_SYM]  '
      STRINGB(2)='b  [b_SYM]  '
      STRINGI(1)='I  [I_SYM]  '
      STRINGJ(1)='J  [J_SYM]  '
      STRINGA(1)='A  [A_SYM]  '
      STRINGB(1)='B  [B_SYM]  '
     
      NBAS   = NOCCO(1) + NVRTO(1)
      CALL GETREC(20,"JOBARC","ORBMAP_A",NBAS,IMAP_A)
      CALL GETREC(20,"JOBARC","ORBMAP_A",NBAS,IMAP_B)
      IF (IUHF.NE.0) CALL GETREC(20,"JOBARC","ORBMAP_B",NBAS,
     &                           IMAP_B)
      
      OCA_OFF(1) = 0
      OCB_OFF(1) = 0

      DO I = 1, (NIRREP-1)
         OCA_OFF(I+1) = OCA_OFF(I) + POP(I,1)
         OCB_OFF(I+1) = OCB_OFF(I) + POP(I,2)
      ENDDO

      CALL ZERO(SCR,NSIZEC) 
C  FIRST CONSIDER SINGLE EXCITATION COEFFICIENTS

      dO 500 SSPIN = 1, 1+IUHF
        CALL GETLST(SCR, 1, 1, 1, SSPIN, LS1OUT)
        ICOUNT = 1
        DO IIRREP = 1, NIRREP
          AIRREP=DIRPRD(IIRREP,IRREPX)
          DO 1 I = 1, POP(IIRREP, SSPIN)
            DO 2 A = 1, VRT(AIRREP,SSPIN)

                 IF (PROJECT_SINGLES) THEN
                 IF (SSPIN .EQ. 1) THEN
                    II = IMAP_A(I+OCA_OFF(IIRREP))
                 ELSE
                    II = IMAP_B(I+OCB_OFF(AIRREP))
                 ENDIF

#ifdef _DEBUG_LVLM
      Write(6,*) "Printing from 1"
      Write(6,"(a,1x,4I4)")"iirrep,airrep,i,a:",iirrep,airrep,i,a
      Write(6,"(1x,3i3,F15.8)") IORB_A(IROOT),II,A,SCR(ICOUNT) 
#endif 
                 IF (IORB_A(IROOT) .EQ. II) THEN
                    SCR(ICOUNT) = 1.0D0
                 ELSE 
                    SCR(ICOUNT) = 0.0D0
                 ENDIF 
                 ELSE
                   SCR(ICOUNT) = 1.0D0
                 ENDIF 
#ifdef _DEBUG_LVLM
      Write(6,*) "Printing from 2"
      Write(6,"(1x,3i3,F15.8)") IORB_A(IROOT),II,A,SCR(ICOUNT) 
#endif 
              ICOUNT = ICOUNT + 1
    2       CONTINUE
    1     CONTINUE
        ENDDO

        CALL PUTLST(SCR, 1, 1, 1, SSPIN, LS1OUT)
  500 CONTINUE
C
      DO ICASE = 1, 1 + 2 * IUHF
        IF (ICASE .EQ. 1) THEN
          ASPIN = 1
          BSPIN = 2
        ELSEIF (ICASE .EQ. 2) THEN
          ASPIN = 1
          BSPIN = 1
        ELSEIF (ICASE .EQ. 3) THEN
          ASPIN = 2
          BSPIN = 2
        ENDIF

      LISTS2EX = LS2OUT(ASPIN, BSPIN)

      DO 10 RIRREP = 1, NIRREP
        LIRREP = DIRPRD(RIRREP,IRREPX)
        IF (ICASE .EQ. 1) THEN
          DISSYS = IRPDPD(LIRREP, ISYTYP(1,LISTS2EX))
          NUMDSS = IRPDPD(RIRREP, ISYTYP(2,LISTS2EX))
          DISSYA = DISSYS
          NUMDSA = NUMDSS
        ELSE
          DISSYA = IRPDPD(LIRREP, ISYTYP(1,LISTS2EX))
          NUMDSA = IRPDPD(RIRREP, ISYTYP(2,LISTS2EX))
          DISSYS = IRPDPD(LIRREP, 18 + ASPIN)
          NUMDSS = IRPDPD(RIRREP, 20 + ASPIN)
        ENDIF

       IF (DISSYA * NUMDSA .GT. 0) THEN
          I000 = 1
          I010 = I000 + DISSYS*NUMDSS *IINTFP
          IF (ICASE .EQ. 1) THEN
            CALL GETLST(SCR(I000), 1, NUMDSS, 1, RIRREP,
     $       LISTS2EX)
            CALL ZERO(SCR(I000),DISSYS*NUMDSS)
          ELSE
            CALL GETLST(SCR(I000), 1, NUMDSA, 1, RIRREP,
     $       LISTS2EX)
            CALL SYMEXP(RIRREP, POP(1,ASPIN),DISSYA,SCR(I000))
            CALL SYMEXP2(LIRREP,VRT(1,ASPIN),DISSYS, DISSYA,
     &         NUMDSS, SCR(I000), SCR(I000))
            CALL ZERO(SCR(I000),DISSYS*NUMDSS)
          ENDIF

          ICOUNT = I000
          DO 20 JIRREP = 1, NIRREP
            IIRREP = DIRPRD(JIRREP, RIRREP)
            DO 30 J= 1, POP(JIRREP,BSPIN)
              DO 40 I = 1, POP(IIRREP, ASPIN)
                DO 50 BIRREP = 1, NIRREP
                  AIRREP = DIRPRD(LIRREP, BIRREP)
                  DO 60 B = 1, VRT(BIRREP, BSPIN)
                    DO 70 A = 1, VRT(AIRREP, ASPIN)

                       IF (ICASE .EQ. 1) THEN

                          II = IMAP_A(I+OCA_OFF(IIRREP))
                          JJ = IMAP_A(J+OCB_OFF(JIRREP))
#ifdef _DEBUG_LVLM
      Write(6,*) "Printing from here"
      Write(6,"(1x,4i3,F15.8)") IORB_AB(IROOT),JORB_AB(IROOT),II,JJ,
     &                          SCR(ICOUNT),ICOUNT 
#endif  
                         IF ((IORB_AB(IROOT) .EQ. II  .OR. 
     &                        IORB_AB(IROOT) .EQ. JJ) .OR. 
     &                       (IORB_AB(IROOT) .EQ. II  .AND. 
     &                        IORB_AB(IROOT) .EQ. JJ)) 
     &                        SCR(ICOUNT) = 1.0D0
                         IF ((JORB_AB(IROOT) .EQ. II  .OR. 
     &                        JORB_AB(IROOT) .EQ. JJ) .OR. 
     &                       (JORB_AB(IROOT) .EQ. II  .AND. 
     &                        JORB_AB(IROOT) .EQ. JJ))
     &                        SCR(ICOUNT) = 1.0D0
#ifdef _DEBUG_LVLM
      Write(6,"(1x,4i3,F15.8)") IORB_AB(IROOT),JORB_AB(IROOT),II,JJ,
     &                           SCR(ICOUNT)
#endif 
                       ELSE IF (ICASE .EQ. 2) THEN 

                          II = IMAP_A(I+OCA_OFF(IIRREP))
                          JJ = IMAP_A(J+OCA_OFF(JIRREP))
#ifdef _DEBUG_LVLM
      Write(6,"(1x,4i3,F15.8)") IORB_AB,JORB_AB,II,JJ,
     &                           SCR(ICOUNT)
#endif 
                         IF ((IORB_AA(IROOT) .EQ. II  .OR.
     &                        IORB_AA(IROOT) .EQ. JJ) .OR.
     &                       (IORB_AA(IROOT) .EQ. II  .AND.
     &                        IORB_AA(IROOT) .EQ. JJ))
     &                        SCR(ICOUNT) = 1.0D0
                         IF ((JORB_AA(IROOT) .EQ. II  .OR.
     &                        JORB_AA(IROOT) .EQ. JJ) .OR.
     &                       (JORB_AA(IROOT) .EQ. II  .AND.
     &                        JORB_AA(IROOT) .EQ. JJ))
     &                        SCR(ICOUNT) = 1.0D0

                       ELSE IF (ICASE .EQ. 3) THEN 

                          II = IMAP_B(I+OCB_OFF(IIRREP))
                          JJ = IMAP_B(J+OCB_OFF(JIRREP))
#ifdef _DEBUG_LVLM
      Write(6,"(1x,4i3,F15.8)") IORB_AB,JORB_AB,I,J,
     &                           SCR(ICOUNT)
#endif 
                         IF ((IORB_BB(IROOT) .EQ. II  .OR.
     &                        IORB_BB(IROOT) .EQ. JJ) .OR.
     &                       (IORB_BB(IROOT) .EQ. II  .AND.
     &                        IORB_BB(IROOT) .EQ. JJ))
     &                        SCR(ICOUNT) = 1.0D0
                         IF ((JORB_BB(IROOT) .EQ. II  .OR.
     &                        JORB_BB(IROOT) .EQ. JJ) .OR.
     &                       (JORB_BB(IROOT) .EQ. II  .AND.
     &                        JORB_BB(IROOT) .EQ. JJ))
     &                        SCR(ICOUNT) = 1.0D0
                       ENDIF 

                    ICOUNT = ICOUNT + 1
   70             CONTINUE
   60           CONTINUE
   50         CONTINUE
   40       CONTINUE
   30     CONTINUE
   20   CONTINUE

      ENDIF

      IF (ICASE .EQ. 1) THEN
         CALL PUTLST(SCR(I000), 1, NUMDSS, 1, RIRREP,
     $               LISTS2EX)
      ELSE
C
C SQUEEZE ARRAY IN PROPER FORM A<B, I<J
C
         CALL SQSYM(LIRREP, VRT(1,ASPIN),DISSYA,DISSYS,
     &              NUMDSS,SCR(I010), SCR(I000))
         CALL TRANSP(SCR(I010), SCR(I000),NUMDSS,DISSYA)
         CALL SQSYM(RIRREP,POP(1,ASPIN), NUMDSA, NUMDSS,
     &              DISSYA, SCR(I010), SCR(I000))
         CALL TRANSP(SCR(I010), SCR(I000),DISSYA,NUMDSA)
         CALL PUTLST(SCR(I000), 1, NUMDSA, 1, RIRREP,
     $               LISTS2EX)
      ENDIF

   10 CONTINUE
      ENDDO
C
      CALL LOADVEC1(IRREPX,SCR,MAXCOR,IUHF,490,0,443,NSIZEC,
     &              .FALSE.)
C
C  PUT NORMALIZED PROJECTION VECTORS BACK ON LIST
C     
#ifdef _DEBUG_LVLM
      Write(6,"(a)")" The checksum of the state spe. core projection"
      Call checksum(" @-DRIVE_CORE_EE",SCR,NSIZEC,S1) 
#endif 
      CALL UPDATES(IRREPX,SCR,444,0,490,IUHF)

C
      RETURN
      END
                






