      Subroutine Inilam(Iuhf)

      Implicit Integer(A-Z)

      Imode = 0
      Ioffw = 100

      Call Inipck(1,9,9,134,IMODE,0,1)
      Call Inipck(1,9,10,137,IMODE,0,1)
      Call Inipck(1,11,12,139,IMODE,0,1)

      If (IUHF.NE.0) Then
          Call Inipck(1,10,9,136,IMODE,0,1)
          Call Inipck(1,12,11,138,IMODE,0,1)
          Call Inipck(1,10,10,135,IMODE,0,1)
      Endif 

      Call Inipck(1,14,18,110+IOFFW,IMODE,0,1)

      IF (IUHF.NE.0)THEN
         Call Inipck(1,3,16,107+IOFFW,IMODE,0,1)
         Call Inipck(1,4,17,108+IOFFW,IMODE,0,1)
         Call Inipck(1,14,11,109+IOFFW,IMODE,0,1)
      Endif 
C
      IF (IUHF.NE.0) Call Inipck(1,10,9,117,IMODE,0,1)
      Call Inipck(1,9,10,118,IMODE,0,1)
      Call Inipck(1,9,9,123,IMODE,0,1)
      IF(IUHF.NE.0) Call Inipck(1,10,10,124,IMODE,0,1)
      Call Inipck(1,11,11,125,IMODE,0,1)
      IF(IUHF.NE.0)Call Inipck(1,12,12,126,IMODE,0,1)

      Call Inipck(1,13,11,130+IOFFW,IMODE,0,1)
      If (IUHF.NE.0) Then
         Call Inipck(1,1,9,127+IOFFW,IMODE,0,1)
         Call Inipck(1,2,10,128+IOFFW,IMODE,0,1)
         Call Inipck(1,13,18,129+IOFFW,IMODE,0,1)
      Endif
      
      Return
      End


