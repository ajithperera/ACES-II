      Subroutine Eompccd_ee(Work,Maxcor,Iuhf)

      Implicit Double Precision(A-H,O-Z)
    
      Dimension Work(Maxcor)

C At this point, we have the Hbar for OO-CCD,OO-LCCD and pCCD (orbital
C optimization is implied for pCCD). We can define EOM methods. 
C My feeling is that the best choice is EOM scheme that has both singles
C and doubles (R1 and R2 ..) regardless of whether we have singles in the
C grount state. It may be possible to define other variant that are better
C in terms of the accuracy and the scalling.
      
      Call Eompccd_ee_main(Work,Maxcor,Iuhf)

      Return
      End
