C
      Subroutine Bannerjs(LU)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: bannerjs.f,v 1.3 2008/11/19 15:48:32 perera Exp $
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
      logical ccsd, mbpt, parteom, nodavid
C
      COMMON/EOMINFO/CCSD, MBPT, PARTEOM, NODAVID
C
      if (ccsd) then
        if (parteom) then
          string ='$Partitioned EOM-CCSD properties will be calculated.'
        else
          string ='$Full EOM-CCSD properties will be calculated.'
        endif
      else
        if (parteom) then
          string ='$Partitioned EOM-MBPT properties will be calculated.'
        else
          string ='$Full EOM-MBPT properties will be calculated.'
        endif
      endif
      Write (LU, '(1X,A)') String(2:Len(string)-1)
      Write (LU, *)
C
C     Convention:  Put author names in alphabetical order.
C     Who is an author?  Anyone who had made non-trivial
C     contributions to the code.  Its up to the maintainer, really.
C
      Write (LU, '(1X,A)') 'Authors: John F. Stanton'
      Write (LU, '(1X,A)') '         S. Ajith Perera'
      Write (LU, '(1X,A)') '         Marcel Nooijen'
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
