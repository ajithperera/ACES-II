      SUBROUTINE TDEE_DR2INR2(ICORE,MAXCOR,IUHF,IRREPX,ISIDE,
     &                        LISTIJKL,LISTABCD,LISTAIBJ,LISTT2,
     &                        LISTT2IN,LISTT2RS,LSTT2RNG)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM
      DIMENSION ICORE(MAXCOR)
      DIMENSION SOURCE(2),TARGET(2),TARRNG(2)
      LOGICAL DISCO,INCREM,WSPIN,THREEBOD,HBAR_4LCCSD
      LOGICAL MBPT2,CC,CCD,RPA,DRPA,LCCD,LCCSD,CC2,ADC2
      COMMON/FLAGS/IFLAGS(100)
      COMMON/FLAGS2/IFLAGS2(500)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /REFTYPE/ MBPT2,CC,CCD,RPA,DRPA,LCCD,LCCSD,CC2
      DATA TARGET /334,344/
      DATA SOURCE /314,324/
      DATA TARRNG /350,360/
      DATA SORRNG /304/
      DATA GINLST /300/

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
             CALL TDEE_DLADDER(ICORE,MAXCOR,IUHF,IRREPX,1,
     &                         SOURCE(ISIDE),11,TARGET(ISIDE),
     &                         ISIDE)
          ELSE

             If (Ispar) Then
C Here I need the W(mn,ij) without the T2 term for pCC. It is
C constructed as a post_vcc and stored in 250-253 lists for pCC. 
c 
                 CALL TDEE_DLADDER(ICORE,MAXCOR,IUHF,IRREPX,1,
     &                            SOURCE(ISIDE),251,TARGET(ISIDE),
     &                            ISIDE)
              Else
                 CALL TDEE_DLADDER(ICORE,MAXCOR,IUHF,IRREPX,1,
     &                             SOURCE(ISIDE),LISTIJKL,
     &                             TARGET(ISIDE),ISIDE)
              Endif 
          Endif 
      
      IF(ISIDE.EQ.1 .AND. .NOT. (CC2 .OR. ADC2)) 
     &  CALL TDEE_MODIFYT2(ICORE,MAXCOR,IUHF,IRREPX,ONE)

      IF(IFLAGS(93).EQ.2)THEN

       CALL TDEE_DRAOLAD(ICORE,MAXCOR,IUHF,.FALSE.,.FALSE.,IRREPX,0,
     &                   (SOURCE(ISIDE)-1),(TARGET(ISIDE)-1),374,384,
     &                    ISIDE)

      ELSEIF(ISYTYP(1,233).EQ.5)THEN
       CALL TDEE_DLADDER2(ICORE,MAXCOR,IUHF,IRREPX,6,SOURCE(ISIDE),
     &                    LISTABCD,TARGET(ISIDE),ISIDE)
      ELSE
       CALL TDEE_DLADDER(ICORE,MAXCOR,IUHF,IRREPX,6,SOURCE(ISIDE),
     &                   LISTABCD,TARGET(ISIDE),ISIDE)
      ENDIF

#ifdef _DEBUG_LVL0
      Call check_target_lists(icore,maxcor,Irrepx,Iuhf,Iside)
#endif
C
      IF(ISIDE.EQ.1 .AND. .NOT. (CC2 .OR. ADC2)) 
     &   CALL TDEE_MODIFYT2(ICORE,MAXCOR,IUHF,IRREPX,ONEM)

      If (Ispar) Then
C This is the extra HH and PP ladder terms. In lambda code this is done 
C using  V(mn,ij). It is directly computed here. Note that HH and PP both 
C are identical to the final R(ij,ab). 
C 
         IF (.NOT. ADC2) CALL TDEE_PDCC_W1LAD(ICORE,MAXCOR,IUHF,ISIDE,
     &                                        IRREPX)
      Else 
         IF (.NOT. ADC2) CALL TDEE_W1LAD(ICORE,MAXCOR,IUHF,ISIDE,IRREPX)
      Endif 
#ifdef _DEBUG_LVL0
      Call check_target_lists(icore,maxcor,Irrepx,Iuhf,Iside)
#endif

      IF (.NOT. ADC2) CALL MODAIBC(ICORE,MAXCOR,IUHF,ONEM)

      IF(ISIDE.EQ.2)THEN

        CALL TDEE_L1W1(ICORE,MAXCOR,IUHF,IRREPX,(SOURCE(ISIDE)-1),
     &                 90,26,(TARGET(ISIDE)-1))

      ELSEIF(ISIDE.EQ.1)THEN
       IF(.NOT. ADC2) CALL TDEE_T1W1(ICORE,MAXCOR,IUHF,IRREPX,
     &                               (SOURCE(ISIDE)-1),90,26,
     &                               (TARGET(ISIDE)-1))
      ENDIF 

      IF (.NOT. ADC2) CALL MODAIBC(ICORE,MAXCOR,IUHF,ONE)

      ENDIF 
#ifdef _DEBUG_LVL0
      Call check_target_lists(icore,maxcor,Irrepx,Iuhf,Iside)
#endif
C
C SET UP RING CONTRACTIONS.  THIS IS SLIGHTLY COMPLICATED.
C
      If (Ispar) Then
C
C The pCC needs scalled F intermediates for the DD terms. These terms
C are stored in 191, 192 (10,11) instead of 1 and 2 positions. 

        CALL TDEE_PDCC_DFT2INT2(ICORE,MAXCOR,IUHF,IRREPX,SORRNG,
     &                          TARRNG(ISIDE),ISIDE)
      Else
        CALL TDEE_DFT2INT2(ICORE,MAXCOR,IUHF,IRREPX,SORRNG,
     &                     TARRNG(ISIDE),ISIDE)
      Endif
#ifdef _DEBUG_LVL0
      Call check_target_lists(icore,maxcor,Irrepx,Iuhf,Iside)
#endif
C
C THREE BODY TERM
C
      THREEBOD= .NOT. CC2 .OR. .NOT. ADC2
   
      IF (HBAR_4LCCSD) THEN
         CALL MODIAJK(ICORE,MAXCOR,IUHF,ONEM)
         CALL MODAIBC(ICORE,MAXCOR,IUHF,ONEM)
      ENDIF

      IF(THREEBOD)THEN
C
       IF(ISIDE.EQ.1)THEN
C
C RIGHT-HAND EIGENVECTOR THREE-BODY CONTRIBUTIONS
C CREARE TWO LISTS of VV and OO 394,396 (IRREPX)
C
        CALL TDEE_FORMQ1(IRREPX,ICORE,MAXCOR,IUHF,ISIDE)

        If (Ispar) Then

C The pCC need scalled G(AE) and G(MI) for the DD terms.

          CALL TDEE_GFORMG2(IRREPX,1,SOURCE(ISIDE),14,GINLST,
     &                      ICORE,MAXCOR,0,Gae_scale,Gmi_scale,IUHF)
        Else
          CALL TDEE_GFORMG2(IRREPX,1,SOURCE(ISIDE),14,GINLST,
     &                      ICORE,MAXCOR,0,ONE,ONE,IUHF)
        Endif 

        CALL TDEE_CNT3BOD2(IRREPX,ICORE,MAXCOR,IUHF,TARRNG(ISIDE))

       ELSEIF(ISIDE.EQ.2)THEN
C
        CALL TDEE_GFORMG(IRREPX,1,SOURCE(ISIDE),44,GINLST,
     &                   IcORE,MAXCOR,0,ONE,ONE,IUHF)
        CALL TDEE_GINC1L(IRREPX,ICORE,MAXCOR,IUHF)

        If (Ispar) Then
C The pCC need scalled G(AE) and G(MI) for the DD terms. Reform
C the with the scalling factor.
C 
           CALL TDEE_GFORMG(IRREPX,1,SOURCE(ISIDE),44,GINLST,
     &                      ICORE,MAXCOR,0,Gae_scale,Gmi_scale,IUHF)
        Endif 
        CALL TDEE_GINC2L(IRREPX,ICORE,MAXCOR,IUHF,TARRNG(ISIDE))
       ENDIF
      ENDIF

      IF (HBAR_4LCCSD) THEN
         CALL MODIAJK(ICORE,MAXCOR,IUHF,ONE)
         CALL MODAIBC(ICORE,MAXCOR,IUHF,ONE)
      ENDIF
#ifdef _DEBUG_LVL0
      Call check_target_lists(icore,maxcor,Irrepx,Iuhf,Iside)
#endif
C
      IF(IUHF.EQ.0 .AND. ISIDE .EQ. 1)CALL ZEROLIST(ICORE,MAXCOR,351)
      IF(IUHF.EQ.0 .AND. ISIDE .EQ. 2)CALL ZEROLIST(ICORE,MAXCOR,361)
      DISCO =ISIDE.EQ.2
      INCREM=.TRUE.
      WSPIN =IUHF.EQ.0
      IF(IUHF.EQ.0)THEN

       IF (.NOT. ADC2) THEN
          CALL TDEE_DRHFRNG(ICORE,MAXCOR,IUHF,IRREPX,DISCO,INCREM,
     &                      WSPIN,SORRNG,LISTAIBJ,TARGET(ISIDE),
     &                      TARRNG(ISIDE),ISIDE)
       ELSE 
          CALL TDEE_DRHFRNG_MODF(ICORE,MAXCOR,IUHF,IRREPX,DISCO,
     &                           INCREM,WSPIN,SORRNG,TARGET(ISIDE),
     &                           TARRNG(ISIDE),ISIDE)
       ENDIF 

      ELSE

       IF (.NOT. ADC2) THEN
          CALL TDEE_DUHFRNG(ICORE,MAXCOR,IUHF,IRREPX,DISCO,INCREM,
     &                      SORRNG,LISTAIBJ,TARGET(ISIDE),
     &                      TARRNG(ISIDE),ISIDE)
       ELSE
          CALL TDEE_DUHFRNG_MODF(ICORE,MAXCOR,IUHF,IRREPX,DISCO,
     &                           INCREM,SORRNG,TARGET(ISIDE),
     &                           TARRNG(ISIDE),ISIDE)
       ENDIF
      ENDIF
#ifdef _DEBUG_LVL0
      Call check_target_lists(icore,maxcor,Irrepx,Iuhf,Iside)
#endif
C
      RETURN
      END
