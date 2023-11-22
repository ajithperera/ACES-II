      SUBROUTINE UNPACK4(ABCD,A,B,C,D)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A, B, C, D, ABCD, FW

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
     
      If (iintfp.eq.2) Then
         fw=8
      End If
      If (iintfp.eq.1) Then
         fw=16
      End IF

      A=and(ishft(ABCD,-3*fw),2**fw-1)
      B=and(ishft(ABCD,-2*fw),2**fw-1)
      C=and(ishft(ABCD,-1*fw),2**fw-1)
      D=and(ABCD,2**fw-1)

      RETURN
      END
