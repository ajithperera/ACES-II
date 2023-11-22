      Subroutine Do_T1crs_2RPA(ICORE, MAXCOR, IUHF)

      Implicit Double Precision (A-H, O-Z)
 
      Dimension Icore(Maxcor)

CEvaluate the following three contributions to T1:
C
C   -1/2 Sum_mef T(im,ef)<ma||ef>  (T2T1AA1, T2T1AB1)
C   -1/2 Sum_men T(nm,ei)<nm||ei>  (T2T1AA2, T2T1AB2)
C
C Also, Do the T1 = T1/{f(i,i) - f(a,a)} to get a new T1
C

      CALL T2INT1(Icore,Maxcor,IUHF,Edummy)
         
      Return
      End
      
         
