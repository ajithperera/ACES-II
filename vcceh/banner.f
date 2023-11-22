C
      Subroutine Banner(LU)
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
C     These parameters should hold the RCSfreeze configuration info
C
      Character*(*) PgmVer
      Parameter (PgmVer = '$Version: experimental$')
      Character*80 String
C
C      Write (LU, '(1X,A)')
C     $      'Quadratic Contribution to the EOM-CCSD Properties'
C      Write (LU, *)
C
C     Convention:  Put author names in alphabetical order.
C     Who is an author?  Anyone who had made non-trivial
C     contributions to the code.  Its up to the maintainer, really.
C
      Write (LU, '(1X,A)') 'Authors: S. Ajith Perera'
      Write (LU, '(1X,A)') '         Marcel Nooijen '
C
      Write (LU, *)
      Write (LU, '(1X,A)') '         Quantum Theory Project'
      Write (LU, '(1X,A)') '         Williamson Hall'
      Write (LU, '(1X,A)') '         University of Florida'
      Write (LU, '(1X,A)') '         Gainesville, FL 32611'
      Write (LU, *)
      String = PgmVer
      Write (LU, '(1X,A)') String(2:Len(PgmVer)-1)
      Write (LU, *)
C
      Return
      End
