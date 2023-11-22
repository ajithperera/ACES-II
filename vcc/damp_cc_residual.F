      Subroutine damp_cc_residual(SCR,MAXCOR,IUHF,SING,
     &                            LISTT20,LISTT21,  
     &                            LISTT10,LISTT11,
     &                            DAMP_PARAMETER)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
       
      LOGICAL SING
      DIMENSION SCR(MAXCOR)

      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/FLAGS/IFLAGS(100)

      DATA ONE, ONEM / 1.0D0, -1.0D0/

      I000   = 1
      IRREPX = 1
      LENVEC = 0

      IF (SING)THEN
         LENVEC =IRPDPD(IRREPX,9)+IUHF*IRPDPD(IRREPX,10)
      ELSE
         LENVEC = 0
      ENDIF 

      DO ISPIN=3,1,IUHF-2
        LENVEC = LENVEC + IDSYMSZ(IRREPX,
     &           ISYTYP(1,60+ISPIN),ISYTYP(2,60+ISPIN))
      ENDDO

      I010 = I000 + LENVEC
      IEND = I010
      IF(IEND.GE.MAXCOR) CALL INSMEM('DAMP_CC_RESIDUAL',IEND,MAXCOR)

      IOFFT  = I000
      IOFFHT = I010

      IF (IUHF .EQ. 0) THEN
         CALL SYMETRIZE_T2(SCR, MAXCOR, 1, LISTT20+2)
         CALL SYMETRIZE_T2(SCR, MAXCOR, 1, LISTT21+2)
      ENDIF 

      IF (SING) THEN

         DO ISPIN=1,1+IUHF

            CALL GETLST(SCR(IOFFT) ,1,1,1,ISPIN,  LISTT10)
            CALL GETLST(SCR(IOFFHT),1,1,1,2+ISPIN,LISTT11)
#ifdef _DEBUG_LVLM
       call checksum("T1-DIIS1",scr(iofft), IRPDPD(IRREPX,8+ISPIN))
       call checksum("T1-DIIS1",scr(ioffht),IRPDPD(IRREPX,8+ISPIN))
#endif

            IOFFT  = IOFFT+IRPDPD(IRREPX,8+ISPIN) 
            IOFFHT = IOFFHT+IRPDPD(IRREPX,8+ISPIN)

         ENDDO

      ENDIF

      DO ISPIN=3,1,IUHF-2

         LISTT=LISTT20-1+ISPIN

         NSIZE=IDSYMSZ(IRREPX,ISYTYP(1,LISTT),ISYTYP(2,LISTT))

         CALL GETALL(SCR(IOFFT) ,NSIZE,IRREPX,LISTT20-1+ISPIN)
         CALL GETALL(SCR(IOFFHT),NSIZE,IRREPX,LISTT21-1+ISPIN)
#ifdef _DEBUG_LVLM
       call checksum("T2-DIIS1",scr(iofft), nsize)
       call checksum("T2-DIIS1",scr(ioffht), nsize)
#endif
         IOFFT  = IOFFT+NSIZE
         IOFFHT = IOFFHT+NSIZE

      ENDDO 

C
C This is T^(new) - T^(old)
C
#ifdef _DEBUG_LVLM
      Write(*,*)
      call checksum("T-NEW", SCR(I010), LENVEC)
#endif
 
      CALL DAXPY (LENVEC,-1.0D0,SCR(I000),1,SCR(I010),1)

#ifdef _DEBUG_LVLM
      RESMITER=SCR(I010-1+ISAMAX(LENVEC,SCR(I010),1))
      Write(*,*) "The largest T-Diff", RESMITER, DAMP_PARAMETER
      call checksum("T-NEW - T-OLD", SCR(I010), LENVEC)
#endif
      CALL DSCAL (LENVEC, DAMP_PARAMETER, SCR(I010), 1)
 
C For the scaled T^(new) and put it back on the same place.
C
      CALL DAXPY (LENVEC,1.0D0,SCR(I000),1,SCR(I010),1)

#ifdef _DEBUG_LVLM
      RESMITER=SCR(I010-1+ISAMAX(LENVEC,SCR(I010),1))
      Write(*,*) "The largest T-Diff", RESMITER
      call checksum("Scaled T-NEW", SCR(I010), LENVEC)
#endif
      IOFFT=I010

      IF (SING) THEN

          DO ISPIN=1,1+IUHF

             CALL PUTLST(SCR(IOFFT),1,1,1,2+ISPIN,LISTT11)
             IOFFT=IOFFT+IRPDPD(IRREPX,8+ISPIN)

          ENDDO 

      ENDIF

      DO ISPIN=3,1,IUHF-2

         LISTT=LISTT20-1+ISPIN

        NSIZE=IDSYMSZ(IRREPX,ISYTYP(1,LISTT),ISYTYP(2,LISTT))

         CALL PUTALL(SCR(IOFFT),NSIZE,IRREPX,LISTT21-1+ISPIN)

         IOFFT=IOFFT+NSIZE

      ENDDO 

      RETURN
      END

