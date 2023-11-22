c
c global variables in the program vibron are set.
c
c maximum dimensions
c
      integer max_mode, max_state, max_quanta,max_irrep,
     $    max_tile
      parameter (max_mode=50)
      parameter (max_state=30)
      parameter (max_quanta=10000)
      parameter (max_irrep=20)
      parameter (max_tile=30)
c
      M_REAL linear_couplings(max_mode, max_state, max_state)
      M_REAL quadratic_couplings(max_mode, max_mode,
     $    max_state, max_state)
      M_REAL cubic_couplings(max_mode, max_state, max_state)
      M_REAL quartic_couplings(max_mode, max_state, max_state)
      M_REAL frequencies(max_mode, max_state)
      M_REAL e_electronic(max_state, max_state)
      M_REAL e_weight(max_state)
      logical diag_quadratic, off_quadratic, linear, auto, vecfol,
     $    use_rfa, do_franck_condon, interpret, tmom_derivative,
     $    unitary, sequential, symmetric, complex, cubic, quartic,
     $    lcdspectrum
c
c The number of irreps in the full symmetry group and the number of
c modes per block.
c
      integer nirrepf, numirpf(max_irrep)
c
c The number of normal modes considered
c
      integer nmode
c
c The number of electronic states
c
      integer nstate, nstate_orig
c
      integer nsym
c
c The number of symmetric modes
c
      integer ncoupling
c
c The number of coupling modes
c
c     The number of quanta per mode.
c
      integer nquanta(max_mode)
c
c The number of electronic absorbing states
c
      integer nstate_abs
c
c The number of quanta in the absorbing states
c
      integer nquanta_abs(max_mode)
c
c the maximum number of setup state in diagonalization algorithms.
c
      integer max_setup
c
c transition moments, and their gradients along normal modes. 
c They couple initial and final states.
c
      M_REAL
     $    tmom0(3, max_state, max_state),
     $    grad_tmom(3, max_mode, max_state, max_state),
     $    tmag0(3, max_state, max_state),
     $    grad_tmag(3, max_mode, max_state, max_state)
c
c indicates either real or complex solutions in transformation
c and subsequent lanczos calculation
c
      integer n_complex
c
c To perform more effective H C multiplications and to reduce the 
c dimensions of the vector, we partition the nomal modes into tiles. For each tile
c we assign a maximum excitation level
c
      integer
     $    dim_tile(max_tile), level_tile(max_tile), tile_low(max_tile),
     $    tile_high(max_tile), ntiles, off_annihilation(max_tile),
     $    off_creation(max_tile), off_quanta(max_tile), n_vibron_stub,
     $    n_graph_stub
c
      character*30 vibron_filename, vibron_stub,
     $    graph_stub, graph_filename
c
      common /vibron_int/ nmode, nstate, nsym, ncoupling, nquanta,
     $     nirrepf, numirpf, nstate_abs, nquanta_abs, max_setup,
     $    n_complex, nstate_orig, n_vibron_stub, n_graph_stub
      common /vibron_log/ diag_quadratic, off_quadratic, linear, auto,
     $    vecfol, use_rfa, do_franck_condon, interpret, tmom_derivative,
     $    unitary, sequential, symmetric, complex, cubic, quartic,
     $    lcdspectrum
      common /vibron_real/ linear_couplings, frequencies, e_electronic,
     $    quadratic_couplings, tmom0, grad_tmom, e_weight,
     $    cubic_couplings, quartic_couplings, tmag0, grad_tmag
      common /vibron_tile/
     $    dim_tile, level_tile, tile_low,
     $    tile_high, ntiles, off_annihilation,
     $    off_creation, off_quanta
      common /vibron_char/ vibron_stub, graph_stub
c
c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End: c
