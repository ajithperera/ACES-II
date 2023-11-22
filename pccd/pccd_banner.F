C
      Subroutine Pccd_banner(LU,Status)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: banner.f,v 1.3 2008/11/19 15:48:32 perera Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     banner -- write Pgmgram title/version banner
C
C SYNOPSIS
      Integer LU
C
C ARGUMENTS
C     LU      Unit to write banner to (input)
C
C DESCRIPTION
C     Write Pgmgram title/version/author banner.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
      Character*(*) PgmVer
      Parameter (PgmVer = '$Version: experimental$')
      Character*80 String
      Character*6 Status
      If (Status .Eq. "Begin ") Then
      Write (LU, '(11X,2A)') 'Starting the coupled cluster orbital',
     +                       ' optimization'
      Elseif (Status .Eq. "Finish") Then
 
      Write (LU, '(11X,2A)') 'Results from the coupled cluster orbital',
     +                       ' optimization'
      Endif 
      Write(LU,"(11x,2a)") "------------------------------------------",
     +                    "-----------"
      Write (LU, '(18X,A)') 'Authors: Ajith Perera'
      Write (LU, '(18X,A)') '         Zack Windom '
C
      Write (LU, *)
      Write (LU, '(18X,A)') '         Quantum Theory Project'
      Write (LU, '(18X,A)') '         Williamson Hall'
      Write (LU, '(18X,A)') '         University of Florida'
      Write (LU, '(18X,A)') '         Gainesville, FL 32611'
      Write (LU, *)
      String = PgmVer
      Write (LU, '(18X,A)') String(2:Len(PgmVer)-1)
      Write (LU, *)
      Write(LU,"(11x,2a)") "------------------------------------------",
     +                    "----------"
C
      Return
      End
