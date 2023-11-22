
              ****************************************************
              * ACES : Advanced Concepts in Electronic Structure *
              *                    Ver. 2.6.0                    *
              *               Release Candidate: 2               *
              ****************************************************

                             Quantum Theory Project
                             University of Florida
                             Gainesville, FL  32611
 @ACES_INIT_RTE: Build-specific parameters:
                 bytes per integer =  4
                 bytes per double  =  8
                 recls per integer =  4


                      ACES STATE VARIABLE REGISTRATION LOG
     ----------------------------------------------------------------------
     The ACES State Variables were initialized.
     Updating (PRINT,1): 0 -> 1
     Updating (COORDINATES,CARTESIAN): 'CARTESIAN' converts to 'CARTESIAN' -> 1
     Updating (UNITS,BOHR): 'BOHR' converts to 'BOHR' -> 1
     Updating (BASIS,STO-3G): String copy succeeded.
     Updating (SPHERICAL,ON): 'ON' converts to 'ON' -> 1
     Updating (SCF_TYPE,KS): 'KS' converts to 'KS' -> 1
     Updating (SCF_MAXCYC,30): 150 -> 30
     Updating (SCF_CONV,7): 7 -> 7
     ----------------------------------------------------------------------


                       ACES STATE VARIABLE VALIDATION LOG
     ----------------------------------------------------------------------
     Updating (PROGRAM,aces2): 'aces2' converts to 'ACES2' -> 2
     Updating (GRAD_CALC,analytical): 'analytical' converts to 'ANALYTICAL' -> 0
     Updating (DERIV_LEV,zero): 'zero' converts to 'OFF' -> 0
     Updating (CC_EXTRAPOL,DIIS): 'DIIS' converts to 'DIIS' -> 1
     Updating (ORBITALS,standard): 'standard' converts to 'STANDARD' -> 0
     Updating (ESTATE_TOL,5): -1 -> 5
     Updating (RESTART,0): '0' converts to 'OFF' -> 0
     Updating (HBARABCD,OFF): 'OFF' converts to 'OFF' -> 1
     Updating (HBARABCI,OFF): 'OFF' converts to 'OFF' -> 1
     Updating (ABCDFULL,OFF): 'OFF' converts to 'OFF' -> 2
     Updating (CACHE_RECS,45): -1 -> 45
     Updating (FILE_RECSIZ,16384): -1 -> 16384
     ----------------------------------------------------------------------

          ASV#   ASV KEY DEFINITION =    CURRENT [   DEFAULT] UNITS
          ------------------------------------------------------------
            1:                PRINT =          1 [         0] 
            2:            CALCLEVEL =          0 [         0] 
            3:            DERIV_LEV =          0 [        -1] 
            4:              CC_CONV =          7 [         7] (tol)

            5:             SCF_CONV =          7 [         7] (tol)
            6:            XFORM_TOL =         11 [        11] (tol)
            7:            CC_MAXCYC =          0 [         0] cycles
            8:           LINDEP_TOL =          5 [         5] 

            9:                  RDO =         -1 [        -1] 
           10:           SCF_EXTRAP =          1 [         1] 
           11:            REFERENCE =          0 [         0] 
           12:          CC_EXPORDER =          0 [         0] 

           13:             TAMP_SUM =          0 [         0] 
           14:            NTOP_TAMP =         15 [        15] 
           15:              DAMPSCF =         20 [        20] x 0.01
           16:           SCF_MAXCYC =         30 [       150] cycles

           18:                PROPS =          0 [         0] 
           19:              DENSITY =          0 [         0] 
           20:          SCF_EXPORDE =          6 [         6] 
           21:          CC_EXTRAPOL =          1 [         1] 

           22:            BRUECKNER =          0 [         0] 
           23:               XFIELD =          0 [         0] x 10-6
           24:               YFIELD =          0 [         0] x 10-6
           25:               ZFIELD =          0 [         0] x 10-6

           26:            SAVE_INTS =          0 [         0] 
           28:               CHARGE =          0 [         0] 
           29:          MULTIPLICTY =          1 [         1] 
           30:          CPHF_CONVER =         12 [        12] (tol)

           31:          CPHF_MAXCYC =         64 [        64] cycles
           35:               INCORE =          0 [         0] 
           36:          MEMORY_SIZE =   15000000 [  15000000] Words
           37:          FILE_RECSIZ =      16384 [        -1] Words

           38:                NONHF =          0 [         0] 
           39:             ORBITALS =          0 [        -1] 
           40:          SCF_EXPSTAR =          8 [         8] 
           41:          LOCK_ORBOCC =          0 [         0] 

           42:          FILE_STRIPE =          0 [         0] 
           43:               DOHBAR =          0 [         0] 
           44:           CACHE_RECS =         45 [        -1] 
           45:                GUESS =          0 [         0] 

           46:           JODA_PRINT =          0 [         0] 
           47:           OPT_METHOD =          0 [         0] 
           48:          CONVERGENCE =          4 [         4] H/bohr
           49:          EIGENVECTOR =          1 [         1] 

           50:              NEGEVAL =          2 [         2] 
           51:          CURVILINEAR =          0 [         0] 
           52:          STP_SIZ_CTL =          0 [         0] 
           53:             MAX_STEP =        300 [       300] millibohr

           54:            VIBRATION =          0 [         0] 
           55:            EVAL_HESS =          0 [         0] # of cyc.
           56:            INTEGRALS =          1 [         1] 
           57:          FD_STEPSIZE =          0 [         0] 10-4 bohr

           58:               POINTS =          0 [         0] 
           59:          CONTRACTION =          1 [         1] 
           60:             SYMMETRY =          0 [         0] 
           62:            SPHERICAL =          1 [         0] 

           63:          RESET_FLAGS =          0 [         0] 
           64:             PERT_ORB =          2 [         2] 
           65:             GENBAS_1 =          0 [         0] 
           66:             GENBAS_2 =          0 [         0] 

           67:             GENBAS_3 =          0 [         0] 
           68:          COORDINATES =          1 [         3] 
           69:            CHECK_SYM =          1 [         1] 
           70:            SCF_PRINT =          0 [         0] 

           71:                  ECP =          0 [         0] 
           72:              RESTART =          0 [         1] 
           73:            TRANS_INV =          0 [         0] 
           74:          HFSTABILITY =          0 [         0] 

           75:             ROT_EVEC =          0 [         0] 
           76:           BRUCK_CONV =          4 [         4] (tol)
           78:                UNITS =          1 [         0] 
           79:          FD_USEGROUP =          0 [         0] 

           80:           FD_PROJECT =          0 [         0] 
           83:                VTRAN =          0 [         0] 
           84:             HF2_FILE =          1 [         1] 
           85:             SUBGROUP =          0 [         0] 

           86:           SUBGRPAXIS =          0 [         0] 
           87:               EXCITE =          0 [         0] 
           88:            ZETA_CONV =         12 [        12] (tol)
           90:           TREAT_PERT =          0 [         0] 

           91:          ESTATE_PROP =          0 [         0] 
           92:           OPT_MAXCYC =         50 [        50] 
           93:             ABCDTYPE =          0 [         0] 
           95:           AO_LADDERS =          1 [         1] 

           96:                 FOCK =          0 [         0] 
           97:          ESTATE_MAXC =         20 [        20] 
           98:           ESTATE_TOL =          5 [        -1] (tol)
           99:            TURBOMOLE =          0 [         0] 

          100:           GAMMA_ABCD =          0 [         0] 
          101:            ZETA_TYPE =          1 [         1] 
          102:          ZETA_MAXCYC =         50 [        50] 
          103:             RESRAMAN =          0 [         0] 

          104:                  PSI =          0 [         0] 
          105:             GEOM_OPT =          0 [         0] 
          106:             EXTERNAL =          0 [         0] 
          107:          HESS_UPDATE =          0 [         0] 

          108:         INIT_HESSIAN =          0 [         0] 
          109:          EXTRAPOLATE =          0 [         0] 
          201:              EA_CALC =          0 [         0] 
          203:                 TDHF =          0 [         0] 

          204:           FUNCTIONAL =          4 [         4] 
          205:           EOM_MAXCYC =         50 [        50] cycles
          206:              EOMPROP =          0 [         0] 
          207:             ABCDFULL =          2 [         0] 

          208:           INTGRL_TOL =         14 [        14] (tol)
          209:             DAMP_TYP =          0 [         0] 
          210:             DAMP_TOL =         10 [        10] x 0.01
          211:              LSHF_A1 =          0 [         0] x 0.01

          212:              LSHF_B1 =          0 [         0] x 0.01
          213:             POLYRATE =          0 [         0] 
          214:              IP_CALC =          0 [         0] 
          216:            IP_SEARCH =          0 [         0] 

          217:               EOMREF =          0 [         0] 
          218:              SOLVENT =          0 [         0] 
          219:            EE_SEARCH =          0 [         0] 
          220:            EOM_PRJCT =          0 [         0] 

          221:               NEWVRT =          0 [         0] 
          222:             HBARABCD =          1 [         0] 
          223:             HBARABCI =          1 [         0] 
          224:             NT3EOMEE =          0 [         0] 

          225:              NOREORI =          0 [         0] 
          227:               KS_POT =          0 [         0] 
          228:             DIP_CALC =          0 [         0] 
          230:             DEA_CALC =          0 [         0] 

          232:              PROGRAM =          2 [         0] 
          233:                CCR12 =          0 [         0] 
          234:            EOMXFIELD =          0 [         0] x 10-6
          235:            EOMYFIELD =          0 [         0] x 10-6

          236:            EOMZFIELD =          0 [         0] x 10-6
          237:              INSERTF =          0 [         0] 
          238:            GRAD_CALC =          0 [         2] 
          239:            IMEM_SIZE =    3000000 [   3000000] Words

          240:              MAKERHF =          0 [         0] 
          241:           GLOBAL_MEM =          0 [         0] Words
          242:             PRP_INTS =          0 [         0] 
          243:           TRUNC_ORBS =          0 [         0] 

          244:             FNO_KEEP =          0 [         0] percent
          245:             FNO_BSSE =          0 [         0] 
          246:              NATURAL =          0 [         0] 
          248:              UNO_REF =          0 [         0] 

          249:           UNO_CHARGE =          0 [         0] 
          250:             UNO_MULT =          1 [         1] 
          251:                RAMAN =          0 [         0] 
          252:            KUCHARSKI =          0 [         0] 

          253:             SCF_TYPE =          1 [         0] 
          254:               DIRECT =          0 [         0] 
          255:         SINGLE_STORE =          0 [         0] 
          ------------------------------------------------------------

                         ACES STATE VARIABLES (STRINGS)
          ------------------------------------------------------------
          BASIS = STO-3G
          OCCUPATION = [ESTIMATED BY SCF]
          ------------------------------------------------------------
  @symmetry-i, Coordinates after  COM shift 
      0.000000000000     -1.434516490000      0.986237261223
      0.000000000000      0.000000000000     -0.124283858777
      0.000000000000      1.434516490000      0.986237261223
   Rotational constants (in cm-1): 
     9.47180       14.51296       27.26834
   Principal axis orientation for molecule:
       -1.434516490000    0.986237261223    0.000000000000
        0.000000000000   -0.124283858777    0.000000000000
        1.434516490000    0.986237261223    0.000000000000
********************************************************************************
   The full molecular point group is C2v .
   The largest Abelian subgroup of the full molecular point group is C2v .
   The computational point group is C2v .
********************************************************************************
         -----------------------------------------------
         Cartesian coordinates corresponding to internal 
                 coordinate input (Bohr) 
 ----------------------------------------------------------------
 Z-matrix   Atomic            C o o r d i n a t e s
  Symbol    Number           X              Y              Z
 ----------------------------------------------------------------
     H         1         0.00000000    -1.43451649     0.98623726
     O         8         0.00000000     0.00000000    -0.12428386
     H         1         0.00000000     1.43451649     0.98623726
 ----------------------------------------------------------------
   Interatomic distance matrix (Angstroms) 

                 H             O             H    
                [ 1]        [ 2]        [ 3]
  H    [ 1]     0.00000
  O    [ 2]     0.96000     0.00000
  H    [ 3]     1.51823     0.96000     0.00000
 @GEOPT-W, Archive file not created for single-point calculation.
@ACES2: Executing "rm -f FILES"
@ACES2: Executing "xjoda"
 @ACES2: INITIALIZATION DUMP
   mopac_guess =  F
   direct =  F
   plain_scf =  T
   nddo_guess =  F
   hf_scf =  F
   dirmp2 =  F
   fno =  F
   geomopt =  F
   raman =  F
   vibfreq =  F
   analytical_gradient =  T
   mrcc =  F
   integral command is "xvmol"
   derivative-integral command is "xvdint && xvksdint"
 @ACES2: END INITIALIZATION DUMP
@ACES2: Executing "xvmol"
 One- and two-electron integrals over symmetry-adapted AOs are calculated.

 @READIN: Spherical harmonics are used.
  
   2 types of atoms
   2 symmetry operations

   Reflection in the YZ-plane
   Reflection in the XZ-plane
  
  Integrals less than  0.10E-13 are neglected.

  Atomic type number  1 :
  -----------------------

    Nuclear charge:                        1
    number of symmetry independent atoms:  1
    highest orbital type:                  s

      1 groups of CGTOs of s type

    nuclear coordinates (in a.u.) are :

NAMN: H #1  0.000000000 -1.434516490  0.986237261
  
  Atomic type number  2 :
  -----------------------

    Nuclear charge:                        8
    number of symmetry independent atoms:  1
    highest orbital type:                  p

      1 groups of CGTOs of s type
      1 groups of CGTOs of p type

    nuclear coordinates (in a.u.) are :

NAMN: O #2  0.000000000  0.000000000 -0.124283859
  
   Internuclear distances (A) :

     for atom H #1 (coordinates :   0.00000  -1.43452   0.98624)

     H #1   3     distance is  1.518
     O #2   1     distance is  0.960

     for atom O #2 (coordinates :   0.00000   0.00000  -0.12428)


   Group multiplication table :

         1   2   3   4
         2   1   4   3
         3   4   1   2
         4   3   2   1

   Gaussian basis information :

   atom                 exponent      coefficients

 H #1  1    S   
+                 1       3.425251   0.1543290
                  2       0.623914   0.5353281
                  3       0.168855   0.4446345

 H #1  3    S   
+                 4       3.425251   0.1543290
                  5       0.623914   0.5353281
                  6       0.168855   0.4446345

 O #2  1    S   
+                 7     130.709320   0.1543290   0.0000000
                  8      23.808861   0.5353281   0.0000000
                  9       6.443608   0.4446345   0.0000000
                 10       5.033151   0.0000000  -0.0999672
                 11       1.169596   0.0000000   0.3995128
                 12       0.380389   0.0000000   0.7001155

 O #2  1    X   
+                13       5.033151   0.1559163
                 14       1.169596   0.6076837
                 15       0.380389   0.3919574

 O #2  1    Y   
+                16       5.033151   0.1559163
                 17       1.169596   0.6076837
                 18       0.380389   0.3919574

 O #2  1    Z   
+                19       5.033151   0.1559163
                 20       1.169596   0.6076837
                 21       0.380389   0.3919574
  7  7
  1  2  3  4  5  6  7
1

             SYMMETRY TRANSFORMATION INFO
              7 SYMMETRY ADAPTED BASIS FUNCTIONS 


          irreducible representation number 1

    1    1    H #1   1s     2   1   1.0   2   1.0
    2    2    O #2   1s     1   3   1.0
    3    3    O #2   1s     1   4   1.0
    4    4    O #2   2pz    1   7   1.0


          irreducible representation number 2

    5    1    O #2   2px    1   5   1.0


          irreducible representation number 3

    6    1    H #1   1s     2   1   1.0   2  -1.0
    7    2    O #2   2py    1   6   1.0
  @READIN-I, Nuclear repulsion energy :    9.1681678365 a.u.
  required memory for a1 array  8903380 words 
  required memory for a2 array  18312 words 
@ACES_MALLOC: allocated >= 34 MB of core memory
  @TWOEL-I,         54 integrals of symmetry type  I I I I
  @TWOEL-I,         45 integrals of symmetry type  I J I J
  @TWOEL-I,         39 integrals of symmetry type  I I J J
  @TWOEL-I, Total number of 2-e integrals        138.


ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.04 (      0.0)   0: 0: 0.33 (      0.3)   0: 0: 0.5 (      0.5)
@ACES2: Executing "xvmol2ja"
@ACES_MALLOC: allocated >= 57 MB of core memory
 @ACES_INIT: heap range 144875520 204875775
 @ACES_CACHE_INIT: allocated >= 2 MB of cache from core memory


ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.02 (      0.0)   0: 0: 0.01 (      0.0)   0: 0: 0.0 (      0.0)
@ACES2: Executing "xvscf_ks"
  @GETREC: record length mismatch for record COMPCLSS
           Actual:      8  Requested:      4
 @LIBRA3: Entered a3getbas.F
 
============================================================================
NAMELIST: SCF
             KEYWORD        TYPE                 VALUE               DEFAULT
--------------------  ----------  --------------------  --------------------
            print_nl     logical                  true
               print     integer                     1
           printevec     element               occ4vir
         printevec_i     element                  none
           printeval     element                   all
         printeval_i     element                  none
             convtol        real            0.00000010
               maxit     integer                    30
               aofil     logical                 false
                  ks     logical                  true
             damptyp     element                  none
             davdamp     integer                    10
             stadamp     integer                    20
              alpha1     integer                     0
               beta1     integer                     0
                 rpp     logical                  true
              rppord     integer                     6
              rppbeg     integer                     8
============================================================================
 
 
============================================================================
NAMELIST: INTGRT
             KEYWORD        TYPE                 VALUE               DEFAULT
--------------------  ----------  --------------------  --------------------
            print_nl     logical                  true
           print_scf     element                   dft
           print_int     element                   dft
          print_atom     element                 never
          print_size     element                 never
           print_mos     element                 never
           print_occ     element                 never
           potradpts     integer                    50
              radtyp     element                 Handy
            partpoly     element                 bsrad
             radscal     element                Slater
             parttyp     element                 fuzzy
           fuzzyiter     integer                     4
            radlimit        real            3.00000000
              numacc     logical                  true
            exact_ex     logical                 false
                tdks     logical                 false
             enegrid     integer                     4
             potgrid     integer                     4
             enetype     element               lebedev
             pottype     element               lebedev
                func      string  none                
               kspot      string  LDA,VWN               hf                  
              cutoff     integer                    12
============================================================================
 
  @VSCF-I, There are 4 irreducible representations.

       Irrep        # of functions
         1                  4
         2                  1
         3                  2
         4                  0


  @VSCF-I, Parameters for SCF calculation: 
             SCF reference function:  RHF
       Maximum number of iterations:   30
          Full symmetry point group: C2v 
          Computational point group: C2v 
             Initial density matrix: MOREAD                 
          SCF convergence tolerance: 10**(- 7)
       RPP convergence acceleration:  ON
               Latest start for RPP:   8
                          RPP order:   6

          Alpha population by irrep:   0    0    0    0
           Beta population by irrep:   0    0    0    0


 @INITGES-I, Routine entered. 
  @INITGES-I, Occupancies from core Hamiltonian:

          Alpha population by irrep:   3    1    1    0
           Beta population by irrep:   3    1    1    0


  --------------------------------------------------------------------
  Iteration         Total Energy            Largest Density Difference
  --------------------------------------------------------------------
       0             9.1681678365              0.0000000000E+00
       1           -69.6294236099              0.2729795902E+01
       2           -74.2590180513              0.2239047411E+01
       3           -74.5022317198              0.1530908930E+01
       4           -74.5870546454              0.1306676647E+01
       5           -74.6647224107              0.1201059147E+01
       6           -74.6949104554              0.1092617435E+01
       7           -74.7332269523              0.1041143941E+01
       8           -74.7485117141              0.9777933036E+00
       9           -74.9626256331              0.5345016290E+00
      10           -74.9626289310              0.4113352012E-03
      11           -74.9626285397              0.9869820783E-04
      12           -74.9626466715              0.2531576688E-02
      13           -74.9626442720              0.4827833571E-03
      14           -74.9626435888              0.1005077011E-03
      15           -74.9626435739              0.2493422397E-05
      16           -74.9626435731              0.1315021918E-06

  @VSCF-I, SCF has converged.

     E(SCF)=       -74.9626435731              0.1785238624E-12

     E(SCF)=       -74.9626435731              0.3022027073E-12


  ORBITAL EIGENVECTORS

  SYMMETRY BLOCK 1 (ALPHA)

                          MO #  1       MO #  2       MO #  3       MO #  4

    BASIS/ORB E          -18.27103      -0.82977      -0.14963       0.31453
                        ----------    ----------    ----------    ----------
    1; O #2  S            -0.00830       0.18232       0.26877      -0.79350
    2; O #2  S             0.99249      -0.22559       0.13081      -0.13242
    3; O #2  S             0.03394       0.77881      -0.62886       0.86912
    4; O #2  Z             0.00536       0.21320       0.75100       0.74851


  SYMMETRY BLOCK 2 (ALPHA)

                          MO #  5

    BASIS\ORB E           -0.05721
                        ----------
    5; H #1  X             1.00000


  SYMMETRY BLOCK 3 (ALPHA)

                          MO #  6       MO #  7

    BASIS\ORB E           -0.38218       0.42345
                        ----------    ----------
    6; O #2  S            -0.43964       0.83907
    7;       Y             0.61289       0.98492


  SYMMETRY BLOCK 4 (ALPHA)

 @PUTMOS-I, Writing converged MOs to NEWMOS. 
 @PUTMOS-I, Symmetry   1 Full   1 Partial   0
 @PUTMOS-I, Symmetry   2 Full   0 Partial   1
 @PUTMOS-I, Symmetry   3 Full   0 Partial   2
 @PUTMOS-I, Symmetry   4 Full   0 Partial   0
  @GETREC: record length mismatch for record COMPPERM
           Actual:     24  Requested:     12
  @GETREC: record length mismatch for record COMPSYOP
           Actual:    144  Requested:     72
  @PRJDEN-I, Analyzing reference function density.


      Trace of projected alpha density matrix =   1.000000000

      Alpha part of wavefunction is symmetric.


  ORBITAL EIGENVALUES (ALPHA)  (1H = 27.2113957 eV)

       MO #        E(hartree)               E(eV)           FULLSYM    COMPSYM
       ----   --------------------   --------------------   -------   ---------
    1     1         -18.2710313076        -497.1802627573      A1        A1 (1)
    2     2          -0.8297696237         -22.5791895691      A1        A1 (1)
    3     6          -0.3821756940         -10.3995340374      B1        B1 (3)
    4     3          -0.1496340057          -4.0717501402      A1        A1 (1)
    5     5          -0.0572077561          -1.5567028897      B2        B2 (2)
  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    6     4           0.3145346991           8.5589281593      A1        A1 (1)
    7     7           0.4234534018          11.5227580767      B1        B1 (3)
  @CRAPSO, You need    200000 words of icore memory.
         , You need  14799998 words of dore memory.
         , You need  14999998 total words of memory.
         , You have  15000000 words of memory.
  @CHECKOUT-I, Total execution time :      22.2400 seconds.
@ACES2: Executing "xintgrt"
 
============================================================================
NAMELIST: INTGRT
             KEYWORD        TYPE                 VALUE               DEFAULT
--------------------  ----------  --------------------  --------------------
            print_nl     logical                  true
           print_scf     element                   dft
           print_int     element                   dft
          print_atom     element                 never
          print_size     element                 never
           print_mos     element                 never
           print_occ     element                 never
           potradpts     integer                    50
              radtyp     element                 Handy
            partpoly     element                 bsrad
             radscal     element                Slater
             parttyp     element                 fuzzy
           fuzzyiter     integer                     4
            radlimit        real            3.00000000
              numacc     logical                  true
            exact_ex     logical                 false
                tdks     logical                 false
             enegrid     integer                     4
             potgrid     integer                     4
             enetype     element               lebedev
             pottype     element               lebedev
                func      string  none                
               kspot      string  LDA,VWN               hf                  
              cutoff     integer                    12
============================================================================
 
  @GETREC: record length mismatch for record COMPCLSS
           Actual:      8  Requested:      4
 @LIBRA3: Entered a3getbas.F

The SCF nuclear-electron attraction energy =   -196.948462236778
The SCF kinetic energy                     =     74.586856847648
The SCF coulomb energy                     =     47.335044697159
The SCF exchange energy                    =     -9.104250717660
The SCF one electron energy                =   -122.361605389130
The SCF total energy                       =    -74.962643573113

Results using the SCF density

   Total density integrates to             :       9.999988641845 electrons
   Exch tot contrib : LDA (Slater, Xalpha) :      -8.202067621493
   Exch tot contrib : Becke                :      -9.079953573653
   Exch tot contrib : Perdew-Burke-Ernzerh :      -9.023798410890
   Exch tot contrib : Perdew-Wang 91       :      -9.054758204652
   Exch tot contrib : Exact Nonlocal Excha :      -9.104250717660
   Corr tot contrib : Vosko-Wilk-Nusair    :      -0.672164959708
   Corr tot contrib : Lee-Yang-Parr        :      -0.339422187414
   Corr tot contrib : Perdew-Burke-Ernzerh :      -0.343766899186
   Corr tot contrib : Perdew-Wang 91       :      -0.366457107550
   Corr tot contrib : Wilson-Levy          :      -0.404097572480
   Corr tot contrib : Wilson-Ivanov        :      -0.350854109375
   Hybr tot contrib : Becke III LYP        :      -9.417225440431
     
 Ehar from   KS   Calculation:         -65.858392855453 a.u.
 Ex   from   KS   Calculation:          -8.202067621493 a.u.
 Ec   from   KS   Calculation:          -0.672164959708 a.u.
 Exc  from   KS   Calculation:          -8.874232581201 a.u.
 Final Result from KS Calculation:     -74.732625436654 a.u.

  @CRAPSO, You need    200000 words of icore memory.
         , You need  14799998 words of dore memory.
         , You need  14999998 total words of memory.
         , You have  15000000 words of memory.
  @CHECKOUT-I, Total execution time :       2.4800 seconds.


 @ACES2: The ACES2 program has completed successfully in  33 seconds.


