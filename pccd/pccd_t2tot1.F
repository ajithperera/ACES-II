      Subroutine Pccd_t2tot1(Work,Maxcor,Iuhf)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

C Evaluate Sum_e T(im,ae)F(m,e) contribution to T1.
C
      CALL PCCD_FMECONT(WORK,MAXCOR,IUHF,1)
      IF (IUHF.NE.0) CALL PCCD_FMECONT(MAXCOR,IUHF,2)

C Evaluate the following three contributions to T1:
C
C   -1/2 Sum_mef T(im,ef)<ma||ef>  (T2T1AA1, T2T1AB1)
C   -1/2 Sum_men T(nm,ei)<nm||ei>  (T2T1AA2, T2T1AB2)
     
      CALL PCCD_E4S(WORK,MAXCOR,IUHF,EDUMMY)

      Return
      End 
      
