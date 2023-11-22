      Subroutine Form_genrlizd_tdens(Scr,Maxcor,Iuhf)

      Implicit Double Precision(A-H,O-Z)

      DIMENSION SCR(MAXCOR),LENVV(2),LENOO(2),LENVO(2),MAP(6)
      DIMENSION SECMOM(6),GSECMOM(6),RSECMOM(6),TM(3),GM(3),RM(3)
      DOUBLE PRECISION LNORM(200),RNORM(200)

      LOGICAL VPROP,MBPT2,CC,CCD,RCCD,DRCCD,LOGICAL LCCD,LCCSD,CC2
      CHARACTER*2 SUFFIX(32)
      CHARACTER*8 LABEL 

#include "flags.com"
#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"
#include "machsp.com"

      COMMON /CALCINFO/ NROOT(8)
      COMMON/REFTYPE/MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
      COMMON/EXTINF3/IROOT,LOCROOT,ITROOT

      DATA SUFFIX /"01","02","03","04","05","06","07","08","09","10",
     +             "11","12","13","14","15","16","17","18","19","20",
     +             "21","22","23","24","25","26","27","28","29","30",
     +             "31","32"/ 
      DATA MAP /1,4,2,5,6,3/
      DATA ONE,ONEM,ZILCH,HALF /1.0D0,-1.0D0,0.0D0,0.5D0/

      I00=100
      CALL GETREC(20,"JOABRC","RNORM   ",I00*IINTFP,RNORM)
      CALL GETREC(20,"JOABRC","LNORM   ",I00*IINTFP,LNORM)

      DO ISPIN=1,1+IUHF
         LENVV(ISPIN)=IRPDPD(1,18+ISPIN)
         LENOO(ISPIN)=IRPDPD(1,20+ISPIN)
         LENVO(ISPIN)=IRPDPD(1, 8+ISPIN)
      ENDDO
C
      IDOO=1
      IDVV =IDOO +LENOO(1)+IUHF*LENOO(2)
      IDVO =IDVV +LENVV(1)+IUHF*LENVV(2)
      IDOV =IDVO +LENVO(1)+IUHF*LENVO(2)
      IDOOG=IDOV +LENVO(1)+IUHF*LENVO(2)
      IDVVG=IDOOG+LENOO(1)+IUHF*LENOO(2)
      IDVOG=IDVVG+LENVV(1)+IUHF*LENVV(2)
      IDOVG=IDVOG+LENVO(1)+IUHF*LENVO(2)
      IDOOR=IDOVG+LENVO(1)+IUHF*LENVO(2)
      IDVVR=IDOOR+LENOO(1)+IUHF*LENOO(2)
      IDVOR=IDVVR+LENVV(1)+IUHF*LENVV(2)
      IDOVR=IDVOR+LENVO(1)+IUHF*LENVO(2)
      ITOP=IDOVR +LENVO(1)+IUHF*LENVO(2)

      MXCOR=MAXCOR-ITOP+1

      IONE=1
      CALL GETREC(20,'JOBARC','NBASTOT ',IONE,NAO)
      NAO2 = NAO* NAO
      NMO=NOCCO(1)+NVRTO(1)
      IDFLMO =ITOP
      IDFLMOG=ITOP   +NAO*NAO
      IDFLMOR=IDFLMOG+NAO*NAO
      IDPROP =IDFLMOR+(IUHF+1)*NAO*NAO
      ITMP1  =IDPROP +NAO*NAO
      ITMP2  =ITMP1  +NAO*NAO
      ITMP3  =ITMP2  +NAO*MAX(NAO,6)
      ITMP4  =ITMP3  +NAO*MAX(NAO,6)
      ITMP5  =ITMP4  +NAO*MAX(NAO,6)

      MAXLEN=0
      DO IRREPX = 1, NIRREP
         LEN = 0
         DO ISPIN = 1, IUHF+1
            LEN=LEN+IRPDPD(IRREPX,8+ISPIN)
         ENDDO
         LEN=LEN+IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
         IF(IUHF.NE.0)THEN
           LEN=LEN+IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
           LEN=LEN+IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
         ENDIF
         MAXLEN=MAX(MAXLEN,LEN)
      ENDDO

      IF (ITOP + MAXLEN .GT. MAXCOR) CALL
     +   INSMEM("@-form_genrlizd_tdens",ITOP+MAXLEN,MAXCOR)

      NTOTAL = 0
      DO IRREP = 1, NIRREP
         NTOTAL = NTOTAL + NROOT(IRREP)
      ENDDO 

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The R and L norms "
      Write(6,"(5(1x,F15.10))") (RNORM(I),I=1,NTOTAL*2)
      Write(6,*)
      Write(6,"(5(1x,F15.10))") (LNORM(I),I=1,NTOTAL*2)
#endif

      IDONE = 0
      NDONE = 0
      DO IRREPX = 1, NIRREP
         LENSZX = 0
         DO ISPIN = 1, IUHF+1
            LENSZX=LENSZX+IRPDPD(IRREPX,8+ISPIN)
         ENDDO
         LENSZX=LENSZX+IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
         IF(IUHF.NE.0)THEN
           LENSZX=LENSZX+IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
           LENSZX=LENSZX+IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
         ENDIF

         DO IROOT = 1, NROOT(IRREPX)
            IDONE = IDONE + 1
            IGET = 17 + (2*IDONE-1) + 1

#ifdef _OBSOLETED 
            IF (IDONE .EQ. NTOTAL) THEN
               Z0 = RNORM(1)
               Z0 = RNORM(IGET-18)
            ELSE
               Z0 = RNORM(IGET-18)
            ENDIF 
#endif 
            Z0 = RNORM(IGET-18)

            CALL GETLST(SCR(ITOP),IGET,1,1,3,472)
#ifdef _DEBUG_LVL0
               Write(6,*)
               write(6,"(2a,(1x,i2))") " Reading the right vector",
     +                                  " of state: ",iget
               Call checksum("list-472      :",scr(itop),lenszx,s)
#endif 
            CALL UPDATES(IRREPX,SCR(ITOP),461,2,490,IUHF)

            JDONE = 0
            DO IRREPY = 1, NIRREP
               LENSZY = 0
               DO ISPIN = 1, IUHF+1
                  LENSZY=LENSZY+IRPDPD(IRREPY,8+ISPIN)
               ENDDO
               LENSZY=LENSZY+IDSYMSZ(IRREPY,ISYTYP(1,46),ISYTYP(2,46))
               IF(IUHF.NE.0)THEN
                 LENSZY=LENSZY+IDSYMSZ(IRREPY,ISYTYP(1,44),ISYTYP(2,44))
                 LENSZY=LENSZY+IDSYMSZ(IRREPY,ISYTYP(1,45),ISYTYP(2,45))
               ENDIF

               DO JROOT = 1, NROOT(IRREPY)
                  JDONE = JDONE + 1
                  JGET = 17 + (2*JDONE-1)  + 2  
#ifdef _OBSOLETED 
                  IF (JDONE .EQ. NTOTAL) THEN
                     Z = LNORM(1)
                     Z = LNORM(JGET-19)
                  ELSE 
                     Z = LNORM(JGET-19)
                  ENDIF 
#endif 
                  Z = LNORM(JGET-19)
                  CALL GETLST(SCR(ITOP),JGET,1,1,3,472)
                  CALL DSCAL(LENSZY,Z,SCR(ITOP),1)
                  
#ifdef _DEBUG_LVL0
               write(6,"(2a,(1x,i2))") " Reading the left vector",
     +                                  " of state : ", Jget 
               Call checksum("list-472      :",scr(itop),lenszy,s)
               Write(6,*)
#endif 
                  CALL UPDATES(IRREPY,SCR(ITOP),444,0,490,IUHF)

                  CALL RESORT(SCR(ITOP),IINTFP*MAXCOR,IUHF,
     +                        IRREPX,461,454)
                  CALL RESORT(SCR(ITOP),IINTFP*MAXCOR,IUHF,
     +                        IRREPY,444,434)
C
C CALCULATE THE EXCITED STATE DENSITY MATRIX
C
                  LSTR1    =490
                  LSTL1    =490
                  LSTR1OFF =2
                  LSTL1OFF =0
                  LISTL2   =444
                  LISTR2   =461
                  LSTL2RS  =434
                  LSTR2RS  =454
                  LSTGRL   =0
                  LSTGTL   =400
                  LSTGRLOF =2
                  LSTGTLOF =0
                  LSTTMP   =90
                  LSTTMPOF =2
                  LSTT1    =90
                  LSTT1OFF =0
                  FACT     =Z*(1.0D0+DFLOAT(1-IUHF))
                  R0       = Z0
#ifdef _DEBUG_LVL0
                   Write(6, "(a,2F12.8)") "The <L|R>;R0: ", Z,R0
#endif 
                  IF(MBPT2 .OR. LCCD .OR. LCCSD .OR. RCCD .OR. DRCCD 
     &                    .OR. CC2) R0=0.0d0
                  CALL GDENS(IRREPY,IRREPX,1,SCR(IDOO),SCR(IDVV),
     &                    SCR(IDVO),SCR(IDOV),SCR(ITOP),MXCOR,IUHF,
     &                    FACT,R0,ZILCH,
     &                    ONE,LSTGRL,LSTGTL,LSTGRLOF,LSTGTLOF,LSTTMP,
     &                    LSTTMPOF,LSTR1,LSTL1,LSTR1OFF,LSTL1OFF,
     &                    LISTR2,LISTL2,LSTR2RS,LSTL2RS,LSTT1,LSTT1OFF)


                  FACT=2.0D0-DFLOAT(IUHF)
                  DO ISPIN=1,1+IUHF
                     IOFFOO =IDOO+(ISPIN-1)*LENOO(1)
                     IOFFVV =IDVV+(ISPIN-1)*LENVV(1)
                     IOFFVO =IDVO+(ISPIN-1)*LENVO(1)
                     IOFFOV =IDOV+(ISPIN-1)*LENVO(1)
                     IOFFOOG=IDOOG+(ISPIN-1)*LENOO(1)
                     IOFFVVG=IDVVG+(ISPIN-1)*LENVV(1)
                     IOFFVOG=IDVOG+(ISPIN-1)*LENVO(1)
                     IOFFOVG=IDOVG+(ISPIN-1)*LENVO(1)
                     IOFFOOR=IDOOR+(ISPIN-1)*LENOO(1)
                     IOFFVVR=IDVVR+(ISPIN-1)*LENVV(1)
                     IOFFVOR=IDVOR+(ISPIN-1)*LENVO(1)
                     IOFFOVR=IDOVR+(ISPIN-1)*LENVO(1)
C
C CALCULATE FULL UNPACKED EXCITED STATE DENSITY IN AO BASIS
C
                     CALL EXPDEN(SCR(IOFFOO),SCR(IOFFVV),SCR(IOFFVO),
     &                           SCR(IOFFOV), SCR(IDFLMO),NMO,1,
     &                           ISPIN,.FALSE.)

                     NDONE =  NDONE +  1
                     LABEL = "GTDENS"//SUFFIX(NDONE)
                     CALL PUTREC(20,"JOBARC",LABEL,NAO*NAO*IINTFP,SCR(
     &                           IDFLMO))
#ifdef _DEBUG_LVL0
                   Write(6,*)
                   Write(6,"(a,2(1x,i2))") "Irrepx,Irrepy : ",
     +                                      Irrepx,Irrepy
                   Write(6,"(a,2(1x,i2))") "Rootx,Rooty   : ",
     +                                      Iroot,Jroot
                   Call Checksum("TDENS-WRITE",SCR(idflmo),Nao2,s)
                   Write(6,*)
#endif

#ifdef _DEBUG_LVLM
                  WRITE(6,5100)
5100              FORMAT(T3,'@EDENS-I, Excited state density matrix in',
     &                      ' MO basis: ')
                     CALL SSCAL(NMO*NMO,FACT,SCR(IDFLMO),1)
                     CALL OUTPUT(SCR(IDFLMO),1,NMO,1,NMO,NMO,NMO,1)
                     CALL SSCAL(NMO*NMO,ONE/FACT,SCR(IDFLMO),1)
#endif 
                  ENDDO 
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END












