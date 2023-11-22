      Subroutine Libint_get_basics(Fname,Max_atoms,Max_shells,Max_prims,
     &                            Ncenters,Nshells,Nprims,Nunique,
     &                            Nsymgen,Kka)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      Character*80 Xline,Fname 
      Character*6 Dumstring
      Logical Finish 
      Character*1 Kka(8,3)
   
      Dimension Idum(8),Scr(Max_prims),Ishl(Max_shells)
C
C Open MOL file for basis set information. The MOL file is created by
C joda by processing user input file and basis set library.
C
      Open(Unit=10,File=Fname,Form='FORMATTED',STATUS='OLD')
      Rewind(10)

C Read the first five lines. Information on those five lines 
C are not relevent in the presnt context. 
C
      read (10,'(a6,4x,9i5)') Dumstring,(Idum(i),i=1,8),Icart

      READ(10,'(A)') Xline
      READ(10,'(A)') Xline

C Read number of nuclear types (NONTYP), symmetry information
C (Nsymgen, Kka), tolerances (TLA, TLC), and coordinates are Iunits=0
C bohr,Iunits=1 Angstorm IRSTRT is not used anywhere.

      READ(10,'(2I5,9A1,1X,2E10.2,I5,I5)') Nunique, Nsymgen,
     &                                   ((Kka(J,I),I=1,3),J=1,3),
     &                                   Tol1,Tol2,Iunits,Irstart 
      READ(10,'(A)') Xline
C     
      LNP1   = 0
      LNPO   = 0
      LTNP1  = 0
      LTNPO  = 0
      Nprims = 1
C
      Ncenters  = 0
      Nshells   = 0
      Finish    = .False.  

      Do While (.Not. Finish)

         Read (10,'(A80)') Xline 
         If (Xline(1:6) .eq. 'FINISH') Finish = .True.
         If (Finish) Go to 99
C
         Ncenters = Ncenters + 1
         Read(Xline,1110) Znuc, Ijun, Nshl,(Ishl(I),I=1,NSHL)
 1110    Format(F20.1,8I5)
         READ(10,'(A)') Xline

         Do I = 1, Nshl
            Npt  = 0
            Naot = 0
            Do  I1 = 1, Ishl(I)
               Nshells = Nshells + 1
               Read(10,1120) Np1, Nao
 1120          Format (2I5)
C
               Nprims = Max(Nao,Nprims) 
               DO J=1,Np1
                  Read(10,*) A, (Scr(k),K=1,Nao)
               Enddo

            Enddo
         Enddo 
      Enddo 
C     
99    Continue 
      Close(10)

      Return
      End
      

