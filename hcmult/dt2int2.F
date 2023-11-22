      SUBROUTINE DT2INT2(ICORE,MAXCOR,IUHF,IRREPX,ISIDE,
     &                   LISTIJKL,LISTABCD,
     &                   LISTAIBJ,LISTT2,LISTT2IN,LISTT2RS,
     &                   LSTT2RNG)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM
      DIMENSION ICORE(MAXCOR)
      LOGICAL DISCO,INCREM,WSPIN,THREEBOD,HBAR_4LCCSD
      LOGICAL  MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2,ADC2
      COMMON/FLAGS/IFLAGS(100)
      COMMON/FLAGS2/IFLAGS2(500)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /REFTYPE/ MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
C
#include "parcc.com"
C
      DATA ONE,ONEM /1.0D0,-1.0D0/

#ifdef _DEBUG_VEE
      write(*,*) '@VEE: Entered ', __FILE__
#endif

      HBAR_4LCCSD = .FALSE.
      HBAR_4LCCSD = (IFLAGS(2) .EQ. 6 .OR. IFLAGS2(117) .EQ. 7)
      ADC2 = .FALSE.
      ADC2 =  (IFLAGS(2) .EQ. 1 .AND. IFLAGS(87) .EQ. 3 .AND.
     &         IFLAGS2(117) .EQ. 10)
C
C SET UP LADDER CONTRACTIONS
C
      IF (.NOT. CC2) THEN

          IF (ADC2) THEN
             CALL DLADDER(ICORE,MAXCOR,IUHF,IRREPX,1,LISTT2, 
     &                    11,LISTT2IN,ISIDE)
          ELSE
CSSS#ifdef _DCC_FLAG

             If (Ispar) Then
C Here I need the W(mn,ij) without the T2 term for pCC. It is
C constructed as a post_vcc and stored in 250-253 lists for pCC. 
c 
                 CALL DLADDER(ICORE,MAXCOR,IUHF,IRREPX,1,LISTT2, 
     &                        251,LISTT2IN,ISIDE)
              Else
CSSS#else
                 CALL DLADDER(ICORE,MAXCOR,IUHF,IRREPX,1,LISTT2, 
     &                       LISTIJKL,LISTT2IN,ISIDE)
              Endif 
          Endif 
CSSS#endif 
C#ifdef _NOSKIP
      IF(ISIDE.EQ.1 .AND. .NOT. (CC2 .OR. ADC2)) 
     &  CALL MODIFYT2(ICORE,MAXCOR,IUHF,IRREPX,ONE)

      IF(IFLAGS(93).EQ.2)THEN
       CALL DRAOLAD(ICORE,MAXCOR,IUHF,.FALSE.,.FALSE.,IRREPX,0,
     &              443,460,213,463)
c     &              443,460,413,463)
      ELSEIF(ISYTYP(1,233).EQ.5)THEN
       CALL DLADDER2(ICORE,MAXCOR,IUHF,IRREPX,6,LISTT2,
     &               LISTABCD,LISTT2IN,ISIDE)
      ELSE
       CALL DLADDER(ICORE,MAXCOR,IUHF,IRREPX,6,LISTT2,
     &              LISTABCD,LISTT2IN,ISIDE)
      ENDIF
      IF(ISIDE.EQ.1 .AND. .NOT. (CC2 .OR. ADC2)) 
     &   CALL MODIFYT2(ICORE,MAXCOR,IUHF,IRREPX,ONEM)

CSSS#ifdef _DCC_FLAG
      If (Ispar) Then
C This is the extra HH and PP ladder terms. In lambda code this is done 
C using  V(mn,ij). It is directly computed here. Note that HH and PP both 
C are identical to the final R(ij,ab). 
C 
         IF (.NOT. ADC2) CALL PDCC_W1LAD(ICORE,MAXCOR,IUHF,ISIDE,IRREPX)
      Else 
CSSS#else
         IF (.NOT. ADC2) CALL W1LAD  (ICORE,MAXCOR,IUHF,ISIDE,IRREPX)
      Endif 
CSSS#endif

      IF (.NOT. ADC2) CALL MODAIBC(ICORE,MAXCOR,IUHF,ONEM)

      IF(ISIDE.EQ.2)THEN
        CALL L1W1(ICORE,MAXCOR,IUHF,IRREPX,443,90,26,460)
      ELSEIF(ISIDE.EQ.1)THEN
       IF(.NOT. ADC2) CALL T1W1(ICORE,MAXCOR,IUHF,IRREPX,443,90,26,460)
      ENDIF 

      IF (.NOT. ADC2) CALL MODAIBC(ICORE,MAXCOR,IUHF,ONE)

C#endif
      ENDIF 
C
C SET UP RING CONTRACTIONS.  THIS IS SLIGHTLY COMPLICATED.
C
CSSS#ifdef _DCC_FLAG
      If (Ispar) Then
C
C The pCC needs scalled F intermediates for the DD terms. These terms
C are stored in 191, 192 (10,11) instead of 1 and 2 positions. 
C 
        CALL PDCC_DFT2INT2(ICORE,MAXCOR,IUHF,IRREPX,LSTT2RNG,LISTT2RS,
     &                    ISIDE)
      Else
CSSS#else
        CALL DFT2INT2(ICORE,MAXCOR,IUHF,IRREPX,LSTT2RNG,LISTT2RS,ISIDE)
      Endif
CSSS#endif

#ifdef _DEBUG_LVL0
        Write(6,*) "F contributions"
        call Check_vee_lists(icore,Maxcor,Irrepx,Iuhf,Iside)
#endif

C
C THREE BODY TERM
C
      THREEBOD= .NOT. (CC2 .OR. ADC2)
      THREEBOD = .TRUE.

      IF (HBAR_4LCCSD) THEN
         CALL MODIAJK(ICORE,MAXCOR,IUHF,ONEM)
         CALL MODAIBC(ICORE,MAXCOR,IUHF,ONEM)
      ENDIF

      IF(THREEBOD)THEN
C
       IF(ISIDE.EQ.1)THEN
C
C RIGHT-HAND EIGENVECTOR THREE-BODY CONTRIBUTIONS
C
        CALL FORMQ1(IRREPX,ICORE,MAXCOR,IUHF)

CSSS#ifdef _DCC_FLAG
        If (Ispar) Then

C The pCC need scalled G(AE) and G(MI) for the DD terms.

          CALL GFORMG2(IRREPX,1,444,14,400,ICORE,MAXCOR,0,
     &                 Gae_scale,Gmi_scale,IUHF)
        Else
CSSS#else
          CALL GFORMG2(IRREPX,1,444,14,400,ICORE,MAXCOR,0,ONE,ONE,
     &                 IUHF)
        Endif 
CSSS#endif
        CALL CNT3BOD2(IRREPX,ICORE,MAXCOR,IUHF,
     &               LISTT2RS)

       ELSEIF(ISIDE.EQ.2)THEN
C
        CALL GFORMG (IRREPX,1,444,44,400,ICORE,MAXCOR,0,
     &               ONE,ONE,IUHF)
        CALL GINC1L(IRREPX,ICORE,MAXCOR,IUHF)

CSSS#ifdef _DCC_FLAG
        If (Ispar) Then
C The pCC need scalled G(AE) and G(MI) for the DD terms. Reform
C the with the scalling factor.
 
           CALL GFORMG (IRREPX,1,444,44,400,ICORE,MAXCOR,0,
     &                   Gae_scale,Gmi_scale,IUHF)
        Endif 
CSSS#endif 
        CALL GINC2L(IRREPX,ICORE,MAXCOR,IUHF,LISTT2RS)
       ENDIF
      ENDIF

      IF (HBAR_4LCCSD) THEN
         CALL MODIAJK(ICORE,MAXCOR,IUHF,ONE)
         CALL MODAIBC(ICORE,MAXCOR,IUHF,ONE)
      ENDIF
C
      IF(IUHF.EQ.0)CALL ZEROLIST(ICORE,MAXCOR,441)
      DISCO =ISIDE.EQ.2
      INCREM=.TRUE.
      WSPIN =IUHF.EQ.0
      IF(IUHF.EQ.0)THEN

       IF (.NOT. ADC2) THEN
          CALL DRHFRNG(ICORE,MAXCOR,IUHF,IRREPX,DISCO,INCREM,WSPIN,
     &                 LSTT2RNG,LISTAIBJ,LISTT2IN,LISTT2RS,ISIDE)
       ELSE 
          CALL DRHFRNG_MODF(ICORE,MAXCOR,IUHF,IRREPX,DISCO,INCREM,WSPIN,
     &                      LSTT2RNG,LISTT2IN,LISTT2RS,ISIDE)
       ENDIF 

      ELSE

       IF (.NOT. ADC2) THEN
          CALL DUHFRNG(ICORE,MAXCOR,IUHF,IRREPX,DISCO,INCREM,
     &                 LSTT2RNG,LISTAIBJ,LISTT2IN,LISTT2RS,ISIDE)
       ELSE
          CALL DUHFRNG_MODF(ICORE,MAXCOR,IUHF,IRREPX,DISCO,INCREM,
     &                      LSTT2RNG,LISTT2IN,LISTT2RS,ISIDE)
       ENDIF

      ENDIF
C
      RETURN
      END
