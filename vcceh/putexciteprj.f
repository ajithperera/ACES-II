C
      SUBROUTINE PUTEXCITEPRJ(SCR, MAXCOR, IRREPX, LISTPTRN, NSIZEC,
     &                        IOFFSET, IUHF)
C
C This subroutine puts calcualted excitation pattern on the appropriate
C columns of LISTPTRN (375).
C
      IMPLICIT INTEGER (A - Z)
      DOUBLE PRECISION SCR(MAXCOR)
      CHARACTER*8 SWITCH
C
      DATA ZILCH, ONE /0.0D0, 1.0D0/
C
      CALL LOADVEC2(IRREPX, SCR, MAXCOR, IUHF, 490, 0, 443, NSIZEC,
     &              .FALSE.)
C
C Put projection vectors to appropriate lists
C    
      CALL PUTLST(SCR, IOFFSET, 1, 1, IRREPX, LISTPTRN)
C
      RETURN
      END










