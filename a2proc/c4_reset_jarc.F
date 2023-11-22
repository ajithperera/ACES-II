      Subroutine c4_reset_jarc()

      Implicit Double Precision (A-H, O-Z)
C

      Call Putrec_c4(1,'JOBARC','HAVEGEOM', 1, 1)
      Call Getrec_c4(0,'JOBARC','PASS1   ', Length, 1)
C
      Call aces_ja_truncate('JODAOUT ', 1)
CSSS      If (Length .gt. 0) Call Putrec(0,'JOBARC','PASS1   ', 1, 0)
C
      Call Putrec_c4(1,'JOBARC','PES_SCAN', 1, 1)
C
      Call Dumpja_c4("O")
      Call Dumpja_c4("I")
c
      Return 
      End
