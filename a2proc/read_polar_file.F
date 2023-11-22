       Subroutine Read_polar_file(NATOMS, Orient)
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include "flags.h"
#include "machsp.com"
#include "jodaflags.com"                   
#include "mxatms.par"
C#include "istart.com"                      
C    
      Character*80 junk
      Character*1  A1
      Logical POLAR_EXSIST
      Dimension  Orient(3,3), Pols(3,3), Reoriented_Pols(3,3),
     &           Tmp(3,3)
C
      Iunit = 4
      INQUIRE(FILE='POLAR',EXIST=POLAR_EXSIST)
      IF(POLAR_EXSIST)THEN
        OPEN(UNIT=Iunit,FILE='POLAR',FORM='FORMATTED',STATUS='OLD')
        REWIND(Iunit)
      ELSE
        WRITE(6,20)
 20      format(T3,'@-READ_POLAR_file, POLAR file not present.')
        Call Errex
      ENDIF

#ifdef _DEBUG_LVL0
      Write(6,"(a,a)") "The transformation matrix from input to",
     &                 " principle axis orientation."
      call output(Orient, 1, 3, 1, 3, 3, 3, 1)
      Write(*,*)
      call xgemm("T","N",3,3,3,1.0D0,orient,3,orient,3,0.0d0,Tmp,
     &            3)
      Write(6,"(a,a)") "T(t)T; Is the transformation unitary?"
      call output(Tmp, 1, 3, 1, 3, 3, 3, 1)
      Write(*,*)
#endif
C
C Read the formatted file POLAR which contain the polarizability 
C tensor (output from a ACESII/III run). During ACESIII run, 
C construct the POLAR file from the summary,out file.  The format
C is as follows.
C---------------------------------------------------------------------------
C
C                          CCSD Polarizability Results (a.u.)
C
C---------------------------------------------------------------------------
C                  X                        Y                      Z
C  X        -272.5231287              58.2676907              15.3189671
C  Y          58.0817560            -164.2989686             -52.3853835
C  Z          15.2188241             -52.3713813             -85.5767566
C

       READ(4,*) Junk
       READ(4,*) 
       READ(4,*) Junk
       READ(4,*) 
       READ(4,*) Junk
       READ(4,*) Junk
       
       READ(4,1000) A1,(Pols(1,I), I=1, 3)
       READ(4,1000) A1,(Pols(2,I), I=1, 3)
       READ(4,1000) A1,(Pols(3,I), I=1, 3)

 1000  FORMAT(2x,A,F20.7,4x,F20.7,4x,F20.7)

#ifdef _DEBUG_LVL0
      Write(6,"(a,a)") "The Polarizablity tensor read from the ",
     &                 " POLAR input file"
      call output(POls, 1, 3, 1, 3, 3, 3, 1)
      Write(*,*)
#endif
C
C Do the transformation from the input orientation to the principle
C axix
c
      Call Xgemm("T","N",3,3,3,1.0D0,Orient,3,Pols,3,0.0D0,
     &            Tmp,3)

      Call Xgemm("N","N",3,3,3,1.0D0,Tmp,3,Orient,3,0.0D0,
     &            Reoriented_Pols,3)

#ifdef _DEBUG_LVL0
      Write(6,"(a,a)") "The Polarizablity tensor in principle axix ",
     &                  " orientation (in a.u.)"
      call output(Reoriented_Pols, 1, 3, 1, 3, 3, 3, 1)
      Write(*,*)
#endif
C    
      Tensor_avg  = (Reoriented_Pols(1,1) + Reoriented_Pols(2,2) +
     &              Reoriented_Pols(3,3))/3.0D0
      Tensor_avg  = -Tensor_avg

      Call symmet2(Reoriented_Pols,3)

      Call Eig(Reoriented_Pols,tmp,0,3,1)

      Avg_sqr     = (Reoriented_Pols(1,1)**2 + Reoriented_Pols(2,2)**2 +
     &               Reoriented_Pols(3,3)**2)
      Avg         = (Reoriented_Pols(1,1) + Reoriented_Pols(2,2) +
     &              Reoriented_Pols(3,3))

      Tensor_anis =  Dsqrt((3.0D0*Avg_sqr - Avg**2)*0.50D0)

      Write(6,"(a,F20.7,a)") "The Average Polarizability:    ", 
     &                        Tensor_avg, " (a.u.)."
      Write(6,"(a,F20.7,a)") "The Polarizability anisotropy: ",
     &                        Tensor_anis,  " (a.u.)."

      Return
      End
