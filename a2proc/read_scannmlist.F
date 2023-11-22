      subroutine Read_scannmlist_init(Streaming,Stationary,Drive_IRC,
     +                                Prop)
C
      Implicit Double Precision (a-h,o-z)
C
      Logical Printdef,Streaming,Stationary,Drive_IRC 
      Character*80 Feature
      Character*10 Prop

c Call to nl_init find the *PES_SCAN namelist if it is in the ZMAT

      printdef = .True.
      call nl_init('PES_SCAN', ierr, printdef)

C The types of scan: currently supported options are 
C Streaming : This is reserved for direct dynamics runs.
C Stationary: This is primarily reserved for running single points calculations
C on a grid
C Dirve_IRC : As the name indicates, this invokes an IRC search.
 
      call nl_str('SCAN_TYPE', "Drive_IRC", Feature)

      If (Feature(1:9)  .EQ. "STREAMING")  Streaming = .True.
      If (Feature(1:10) .EQ. "STATIONARY") Stationary= .True.
      If (Feature(1:9)  .EQ. "Drive_IRC")  Drive_IRC = .True.

      call nl_str('PROP_TYPE', "Drive_IRC", Feature)
      
      If (Feature(1:9)  .EQ. "JCOUPLING")  PROP = "NMR_SPN_CC"
      If (Feature(1:6)  .EQ. "SHIFTS")     PROP = "NMR_SHIFTS"
      If (Feature(1:6)  .EQ. "ENERGY")     PROP = "SP_ENEEGY "

      call nl_term

      RETURN
      END
