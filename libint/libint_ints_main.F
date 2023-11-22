      Subroutine Libint_ints_main(Nshells,Ncenters,Nalpha,Npcoef,
     &                            Nbasis,Nprims,Spherical,Cartesian,
     &                            Angmom,Ncfps,Npfps,Ioff_alpha,
     &                            Ioff_pcoef,Nfps,Index_ccf,End_nfps,
     &                            CC_begin,CC_end,Index_int,
     &                            Atom_4shell,Alpha,Pcoef,
     &                            Shell_coords,Atom_Coords,Scale,
     &                            Alpha_pack,Pcoef_pack,Cc_begin_pack,
     &                            Cc_end_pack,Dwork,Iwork,Dleft,Ileft,
     &                            Max_cfps,Max_pfps)
      Implicit Integer(A-Z)

      Logical Spherical,Cartesian
      Logical Converged

      Double Precision Alpha,Pcoef,Atom_coords,Scale
      Double Precision Dwork
      Double Precision Shell_coords 
      Double Precision Ijij_norm 
      Double Precision Dmax,Den_Tol,Cholesky_thres

      Dimension Dwork(Dleft)
      Dimension Iwork(Ileft)

      Dimension Angmom(Nshells),Ncfps(Nshells),Npfps(Nshells)
      Dimension Ioff_pcoef(Nshells),Ioff_alpha(Nshells)
      Dimension Nfps(Nshells)
      Dimension Index_ccf(Nshells),End_nfps(Nshells)
      Dimension Cc_begin(Nprims)
      Dimension Cc_end(Nprims)
      Dimension Cc_begin_pack(4*Max_cfps)
      Dimension Cc_end_Pack(4*Max_cfps)
      Dimension Index_int(Nbasis)
      Dimension Atom_4shell(Nshells)

      Dimension Alpha(Nalpha),Pcoef(Npcoef)
      Dimension Alpha_pack(4*Max_pfps)
      Dimension Pcoef_pack(4*Max_cfps*Max_pfps)
      Dimension Atom_coords(3,Ncenters)
      Dimension Shell_coords(3,Nshells)
      Dimension Scale(Nbasis)

      Data Den_Tol,Cholesky_thres/1.0D-09,1.0D-09/
      Data Izero,Ione/0,1/

C#ifdef _DEBUG
      Print*, "Nshells,Npcoef,Cartesian,Spherical"
      Print*, Nshells,Npcoef,Cartesian,Spherical
      Print*, "Angmom"
      Print*, Angmom
      Print*, "Ncfps"
      Print*, Ncfps
      Print*, "Npfps"
      Print*, Npfps
      Print*, "Nfps"
      Print*, Nfps 
      Print*, "Ioff_pcoef"
      Print*, Ioff_pcoef
      Print*, "Ioff_alpha"
      Print*, Ioff_alpha
      Print*, "Index_ccf"
      Print*, Index_ccf
      Print*, "Index_int"
      Print*, Index_int
      Print*, "Atom_4shell"
      Print*, Atom_4shell 
      Print*, "End_nfps"
      Print*, End_nfps 
      Print*, "Cc_begin"
      Print*, Cc_begin
      Print*, "Cc_end"
      Print*, Cc_end
      Write(6,*) 
      Print*, "Alpha"
      Print*, Alpha
      Print*, "Pcoef"
      Print*, Pcoef
      Print*, "Atom Coords"
      Print*, Atom_Coords 
      Print*, "Shell_Coords"
      Print*, Shell_Coords
      Print*, "Scale"
      Print*, Scale 
C#endif 

      Call Erd_scratch_mem_calc(Nshells,Nprims,Nalpha,Npcoef,
     &                          Angmom,Ncfps,Npfps,
     &                          Atom_4shell,Shell_coords,Alpha,
     &                          Alpha_pack,Pcoef,Pcoef_pack,
     &                          Ioff_alpha,Ioff_pcoef,Cc_begin,
     &                          CC_begin_pack,Cc_end,Cc_end_pack,
     &                          Index_ccf,Spherical,Ncenters,
     &                          Max_cfps,Max_pfps,Intmax,Zmax)


      Call Erd__memory_hrr_correction(AngMom,Nshells,Spherical,
     &                                Ihrr,Zhrr)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " ERD/OED internal memory requirements"
      Write(6,"(a,4(2x,I4))") " Intmax,Zmax,Ihrr,Zhrr:",Intmax,Zmax,
     &                                            Ihrr,Zhrr
#endif 
      Intmax  = Intmax + Ihrr + 10
      Zmax    = Zmax   + Zhrr
      Nbasis2 = Nbasis*Nbasis 
      Val1_b  = 1
      Val2_b  = 1
      Val1_e  = Nbasis 
      Val2_e  = Nbasis 


C Get the maximum number of contrates function in shells and allocate
C memmory based on it for the residual and the ket integral blocks.

      Max_ncfps = Izero
      Do Ishell = 1, Nshells
         Max_ncfps = Max(Max_ncfps,Ncfps(Ishell))
      Enddo
      Max_ncfps2 = Max_ncfps*Max_ncfps

      I000 = Ione
      J000 = Ione
      J010 = J000 + Intmax
      J020 = J010 + Nbasis
      Jend = J020 + Nbasis

      I010 = I000 + Zmax 
CSSS      Iend = I010 + Nbasis*Nbasis*Nbasis*Nbasis 
      Iend = I010 + Max_ncfps2*Max_ncfps2

      If (Jend .Gt. Ileft) Call Insmem("Libint_ints_main",Jend,
     &                                  Ileft)
      If (Iend .Gt. Dleft) Call Insmem("Libint_ints_main",Iend,
     &                                  Dleft)

C Compute the diagonal block of integrals ((IJ|IJ))

      Call Compute_4cent_ijij_ints(
     &                             Val1_b,Val1_e,Val2_b,Val2_e,
     &                             Dwork(I000),Iwork(J000),
     &                             Atom_coords,Pcoef,Pcoef_pack,
     &                             Alpha,Alpha_pack,Cc_begin,
     &                             Cc_begin_pack,Cc_end,Cc_end_pack,
     &                             End_nfps,Nshells,Nbasis,
     &                             Max_cfps,Max_pfps,Zmax,Intmax,
     &                             Npfps,Ncfps,Angmom,Ioff_alpha,
     &                             Index_ccf,Atom_4shell,Ioff_pcoef,
     &                             Spherical,Nalpha,Npcoef,Index_int,
     &                             Scale,Ncenters,Nprims,DWork(I010))

      I020 = I010
CSSS      I030 = I020 + Nbasis*Nbasis*Max_ncfps2
      I030 = I020 + Nbasis2*Nbasis2 
      I040 = I030 + Nbasis2
      I050 = I040 + Nbasis2
      Iend = I050 + Nbasis2
      If (Iend .Gt. Dleft) Call Insmem("Libint_ints_main",Iend,
     &                                  Dleft)

C Save the diagonals (IJ|IJ) as the first residual (D_1) and compute the
C bra block of integrals to ket shell pair that correspond to the
C largest bra residual element.

      Call Compute_0th_residuals(DWork(I020),Dwork(I030),Iwork(J010),
     &                           Iwork(J020),Nbasis,Imax,Jmax,Ijij_norm)
C 
C Now compute the integral block (IJ|L1L2) where Ls are fixed to the
C coulmn and row ndex of the maximum diagonal element. Since integral
C codes work with shell, we fist compute the entire bock of integrals
C for shells that include the shells of row and coulmn indices of 
C the maximum diagonal element.

      Call Lookup_shell(End_nfps,Nshells,Imax,Shell_a)
      Call Lookup_shell(End_nfps,Nshells,Jmax,Shell_b)
      Ncfps_a = Ncfps(Shell_a)
      Ncfps_b = Ncfps(Shell_b)

      If (shell_a .Eq. 1) Then
         Aa1 = Ione
      Else
         Aa1 = End_nfps(Shell_a-1) + Ione
      Endif
      Aa2  = End_nfps(Shell_a)
      If (shell_b .Eq. 1) Then
          Bb1 = Ione
      Else
          Bb1 = End_nfps(Shell_b-1) + Ione
      Endif
      Bb2 = End_nfps(Shell_b)

      Val1_b = Ione
      Val2_b = Ione
      Val3_b = Aa1
      Val4_b = Bb1
      Val1_e = Nbasis
      Val2_e = Nbasis
      Val3_e = Aa2
      Val4_e = Bb2

      Call Compute_4cent_ij11_ints(
     &                             Val1_b,Val1_e,Val2_b,Val2_e,
     &                             val3_b,val3_e,Val4_b,Val4_e,
     &                             Dwork(I000),Iwork(J000),
     &                             Atom_coords,Pcoef,Pcoef_pack,
     &                             Alpha,Alpha_pack,Cc_begin,
     &                             Cc_begin_pack,Cc_end,Cc_end_pack,
     &                             End_nfps,Nshells,Nbasis,
     &                             Max_cfps,Max_pfps,Zmax,Intmax,
     &                             Npfps,Ncfps,Angmom,Ioff_alpha,
     &                             Index_ccf,Atom_4shell,Ioff_pcoef,
     &                             Spherical,Nalpha,Npcoef,Index_int,
     &                             Scale,Ncenters,Nprims,DWork(I020),
     &                             Shell_a,Shell_b,Ncfps_a,Ncfps_b)
C This is the first cholesky vector.

      Call Compute_1st_cholesky_vec(Dwork(I020),Dwork(I050),Imax,Jmax,
     &                              Aa1,Aa2,Bb1,Bb2,Ijij_norm,Nbasis,
     &                              Ncfps_a,Ncfps_b,Den_tol)

C Cholesky iterations.

      Converged  = .False. 
      Iter_count = Ione

      Do while (.Not. Converged) 

         Iter_count = Iter_count+1
         Write(6,"(2a,I3)") " Starting the Cholesky Iterations.",
     &                      " Iterations number:", Iter_count

         Call Update_diag_residuals(Dwork(I030),Dwork(I040),Dwork(I050),
     &                              Iwork(J010),Iwork(J020),Nbasis,Aa1,
     &                              Aa2,Bb1,Bb2,Kmax,Lmax,Ijij_norm,
     &                              Dmax,Iter_count-1)

         Print*, "Dmax =", Dmax
         If (Dmax .Lt. Cholesky_thres .Or. Iter_count .Eq. 20)
     &       Converged = .True.
         If (converged) Go to 10 
C Find the shell pair that Imax and Jmax belongs. Lmax and Kmax are
C indices of the largest diagonal residual Dmax

         Call Lookup_shell(End_nfps,Nshells,Lmax,Shell_a)
         Call Lookup_shell(End_nfps,Nshells,Kmax,Shell_b)
         Ncfps_a = Ncfps(Shell_a)
         Ncfps_b = Ncfps(Shell_b)

         If (shell_a .Eq. 1) Then
            Aa1 = 1
         Else
            Aa1 = End_nfps(Shell_a-1) + 1
         Endif 
         Aa2  = End_nfps(Shell_a)
         If (shell_b .Eq. 1) Then
            Bb1 = 1
         Else
            Bb1 = End_nfps(Shell_b-1) + 1
         Endif 
         Bb2 = End_nfps(Shell_b)

C Compute all the bra integrals for the ket shell pair of the indices that
C correspond to the largest diagonal residual.

         Val1_b = 1
         Val2_b = 1
         Val3_b = Aa1
         Val4_b = Bb1
         Val1_e = Nbasis
         Val2_e = Nbasis
         Val3_e = Aa2
         Val4_e = Bb2

         Call Compute_4cent_ij11_ints(Val1_b,Val1_e,Val2_b,Val2_e,
     &                                val3_b,Val3_e,Val4_b,Val4_e,
     &                                Dwork(I000),Iwork(J000),
     &                                Atom_coords,Pcoef,Pcoef_pack,
     &                                Alpha,Alpha_pack,Cc_begin,
     &                                Cc_begin_pack,Cc_end,
     &                                Cc_end_pack,End_nfps,Nshells,
     &                                Nbasis,Max_cfps,Max_pfps,
     &                                Zmax,Intmax,Npfps,Ncfps,
     &                                Angmom,Ioff_alpha,Index_ccf,
     &                                Atom_4shell,Ioff_pcoef,
     &                                Spherical,Nalpha,Npcoef,
     &                                Index_int,Scale,Ncenters,
     &                                Nprims,DWork(I020),
     &                                Shell_a,Shell_b,Ncfps_a,Ncfps_b)
    
         Iend = I050 + Nbasis2*Iter_count
         If (Iend .Gt. Dleft) Call Insmem("Libint_ints_main",Iend,
     &                                     Dleft)
         Call Compute_kth_residuals(Dwork(I020),Dwork(I040),
     &                              Dwork(I050),Iwork(J010),
     &                              Iwork(J020),Nbasis,Ncfps_a,
     &                              Ncfps_b,Aa1,Aa2,Bb1,Bb2,
     &                              Lmax,Kmax,Iter_count-1)

         Call Compute_kp1_cholesky_vec(Dwork(I040),Dwork(I050),
     &                                 Nbasis,Nbasis2,Ncfps_a,
     &                                 Ncfps_b,Aa1,Aa2,Bb1,Bb2,
     &                                 Lmax,Kmax,Ijij_norm,
     &                                 Den_tol,Iter_count) 

      Enddo

  10  Continue 

      Write(6,"(a,I3,a)") " Chelosky Iterations have converged in ",
     &                     Iter_count," iterations."
#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,"(a)") " The Chelosky vectors"
      Call output(Dwork(I050),1,Nbasis2,1,Iter_count,Nbasis2,
     &            Iter_count,1)
      Call check(Dwork(I050),Dwork(I020),Nbasis,Iter_count)
#endif 

      Return
      End
