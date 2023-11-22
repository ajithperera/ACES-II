      Subroutine Hproject_main(P,Coords,Mass,Natoms)

      Implicit None

      Double Precision Coords(Natoms,3)
      Double Precision Mass(Natoms)
      Double Precision Threshold 
      Double Precision P

      Dimension P(3*Natoms,3*Natoms)
 
      Logical Move_2cmass
      Logical Project_rotas
      Logical Project_trans
    
      Integer Natoms
      Integer I,J
      Integer Ione 

      Data Ione /1/
C
C This routine project the rotational and translational motions 
C of Cartesian Hessian using the procedure outlined in 
C Miller, Handy and Adams, J. Chem. Phys. 92, 99, (1980).
C
C The output variables
C P         : The projecter 
C The input variables
C Natoms    : Number of real atoms.
C Mass      : Atomic massess 
C Work      : Work is a size of length scratch array. 
C Length    : Length of the scracth array. Length should be at least
C             9*Natoms*Natoms long.

C The important internal variables 
C Nx            : 3*Natoms 
C Nx2           : 9*Natoms*Natoms 
C Move-3cmass   : Move the center of the coordiantes to center of mass
C Project_rotas : Project rotations (True/False)
C Project_trans : Project translations (True/Falss)
C Threshold     : Cutoff for the length of vectors to avoid division
C                 by nearly zero during normalization. Set to 1.0D-09

      Project_rotas = .True.
      Project_trans = .True.
      Move_2cmass   = .True.
      Threshold     = 1.0D-09

C Build the projector. It is returned in Pj

      Call Project_fc(Coords,Mass,P,Natoms,Threshold,Move_2Cmass,
     +                Project_trans,Project_rotas)

      Return
      End

