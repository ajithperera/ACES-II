      Subroutine Gden_debug(Work,Maxcor,Iuhf)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

#include "flags.com"
#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"
#include "machsp.com"
      COMMON /CALCINFO/ NROOT(8) 

      I000 = 1

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

      Print*, ITOP,MAXLEN,MAXCOR

      IF (ITOP + MAXLEN .GT. MAXCOR) CALL
     +   INSMEM("@-gtden_debug",ITOP+MAXLEN,MAXCOR)

      NTOTAL = 0
      DO IRREP = 1, NIRREP
         NTOTAL = NTOTAL + NROOT(IRREP)
      ENDDO

      IDONE = 0
      DO IRREPX = 1, 1

         LENSZ = 0
         DO ISPIN = 1, IUHF+1
            LENSZ=LENSZ+IRPDPD(IRREPX,8+ISPIN)
         ENDDO
         LENSZ=LENSZ+IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
         IF(IUHF.NE.0)THEN
           LENSZ=LENSZ+IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
           LENSZ=LENSZ+IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
         ENDIF

         Print*, IRREPX, NROOT(IRREPX),LENSZ
         DO IROOT = 1, NROOT(IRREPX)
            IDONE = IDONE + 1
            IF (IDONE .EQ. NTOTAL) THEN
               IGET = 2
               CALL GETLST(WORK(I000),IGET,1,1,1,472)
            ELSE
               IGET = 17 + (2*IDONE-1) + 1
               CALL GETLST(WOrk(I000),IGET,1,1,3,472)

#ifdef _DEBUG_LVL0
               Write(6,*) "From Tdens_debug"
               write(6,"(2a,3(1x,i2))") " Reading the right vector",
     +                                  " of state: ",
     +                             Idone,iget,iroot
               Call checksum("list-472      :",work(i000),lensz,s)
#endif
            ENDIF

            JDONE = 0
            DO JROOT = 1, IROOT
               JDONE = JDONE + 1
               JGET = 17 + (2*JDONE-1)  + 2
               IF (NDONE .EQ. NTOTAL) THEN
                  JGET = 3
                  CALL GETLST(WORK(I000),JGET,1,1,1,472)
               ELSE
                  CALL GETLST(Work(I000),JGET,1,1,3,472)
               ENDIF
#ifdef _DEBUG_LVL0
               write(6,"(2a,3(1x,i2))") " Reading the left vector",
     +                                  " of state : ",
     +                                   jdone,jget,jroot
               Call checksum("list-472      :",work(i000),lensz,s)
               Write(6,*)
#endif
            ENDDO
         ENDDO
      ENDDO

      Return
      End

           
      
