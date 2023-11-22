
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
     Updating (BASIS,CC-PVTZ): String copy succeeded.
     Updating (SPHERICAL,ON): 'ON' converts to 'ON' -> 1
     Updating (SCF_TYPE,KS): 'KS' converts to 'KS' -> 1
     ----------------------------------------------------------------------


                       ACES STATE VARIABLE VALIDATION LOG
     ----------------------------------------------------------------------
     Updating (COORDINATES,INTERNAL): 'INTERNAL' converts to 'INTERNAL' -> 0
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
            1:                PRINT =          0 [         0] 
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
           16:           SCF_MAXCYC =        150 [       150] cycles

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
           68:          COORDINATES =          0 [         3] 
           69:            CHECK_SYM =          1 [         1] 
           70:            SCF_PRINT =          0 [         0] 

           71:                  ECP =          0 [         0] 
           72:              RESTART =          0 [         1] 
           73:            TRANS_INV =          0 [         0] 
           74:          HFSTABILITY =          0 [         0] 

           75:             ROT_EVEC =          0 [         0] 
           76:           BRUCK_CONV =          4 [         4] (tol)
           78:                UNITS =          0 [         0] 
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
          BASIS = CC-PVTZ
          OCCUPATION = [ESTIMATED BY SCF]
          ------------------------------------------------------------
   5 entries found in Z-matrix 

 Job Title : NMR SPIN-SPIN COUPLING CONSTANT

  There are  5 unique internal coordinates.
  Of these,  0 will be optimized.
   User supplied Z-matrix: 
--------------------------------------------------------------------------------
       SYMBOL    BOND      LENGTH    ANGLE     ANGLE     DIHED     ANGLE
                  TO      (ANGST)    WRT      (DEG)      WRT      (DEG)
--------------------------------------------------------------------------------
        X    
        N          1         NX   
        H          2         NH         1        A    
        H          2         NH         1        A          3        D1   
        H          2         NH         1        A          3        D2   
                  *Initial values for internal coordinates* 
                      Name             Value
                        NX              1.000000
                        NH              1.012395
                        A             112.145109
                        D1            120.000000
                        D2           -120.000000
--------------------------------------------------------------------------------
  @symmetry-i, Coordinates after  COM shift 
      0.000000000000      0.000000000000     -2.017786920356
      0.000000000000      0.000000000000     -0.128060931777
      1.772020519055      0.000000000000      0.593107663183
     -0.886010259528      1.534614785529      0.593107663183
     -0.886010259528     -1.534614785529      0.593107663183
   Rotational constants (in cm-1): 
     6.34071        9.96627        9.96627
   Principal axis orientation for molecule:
        0.000000000000    0.000000000000   -2.017786920356
        0.000000000000    0.000000000000   -0.128060931777
       -0.886010259528   -1.534614785529    0.593107663183
        1.772020519055    0.000000000000    0.593107663183
       -0.886010259528    1.534614785529    0.593107663183
********************************************************************************
   The full molecular point group is C3v .
   The largest Abelian subgroup of the full molecular point group is C s .
   The computational point group is C s .
********************************************************************************
--------------------------------------------------------------------------------
   Analysis of internal coordinates specified by Z-matrix 
--------------------------------------------------------------------------------
   *The nuclear repulsion energy is   11.95411 a.u.
   *There are   2 degrees of freedom within the tot. symm. molecular subspace.
   *Z-matrix requests optimization of   0 coordinates.
   *The optimization is constrained.
   *The following   2 parameters can have non-zero 
    derivatives within the totally symmetric subspace:
             A    [  5]  NH   [  7]
   *The following   0 parameters are to be optimized:

   *The following coordinates must be varied in an  unconstrained optimization.
             A    [  5]  NH   [  7]
--------------------------------------------------------------------------------
         -----------------------------------------------
         Cartesian coordinates corresponding to internal 
                 coordinate input (Bohr) 
 ----------------------------------------------------------------
 Z-matrix   Atomic            C o o r d i n a t e s
  Symbol    Number           X              Y              Z
 ----------------------------------------------------------------
     X         0        -2.01778692     0.00000000     0.00000000
     N         7        -0.12806093     0.00000000     0.00000000
     H         1         0.59310766    -0.88601026    -1.53461479
     H         1         0.59310766     1.77202052     0.00000000
     H         1         0.59310766    -0.88601026     1.53461479
 ----------------------------------------------------------------
   Interatomic distance matrix (Angstroms) 

                 X             N             H             H             H    
                [ 1]        [ 2]        [ 3]        [ 4]        [ 5]
  X    [ 1]     0.00000
  N    [ 2]     1.00000     0.00000
  H    [ 3]     1.66979     1.01240     0.00000
  H    [ 4]     1.66979     1.01240     1.62417     0.00000
  H    [ 5]     1.66979     1.01240     1.62417     1.62417     0.00000
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
NAMN: N #1 -0.128060932  0.000000000  0.000000000
NAMN: H #2  0.593107663 -0.886010260 -1.534614786
NAMN: H #3  0.593107663  1.772020519  0.000000000
             72 SYMMETRY ADAPTED BASIS FUNCTIONS 
  @READIN-I, Nuclear repulsion energy :   11.9541086066 a.u.
  required memory for a1 array  8903380 words 
  required memory for a2 array  141128 words 
@ACES_MALLOC: allocated >= 34 MB of core memory
  @TWOEL-I,     522547 integrals of symmetry type  I I I I
  @TWOEL-I,     693102 integrals of symmetry type  I J I J
  @TWOEL-I,     365375 integrals of symmetry type  I I J J
  @TWOEL-I, Total number of 2-e integrals    1581024.


ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.94 (      0.9)   0: 0:44.43 (     44.4)   0: 0:46.3 (     46.3)
@ACES2: Executing "xvmol2ja"
@ACES_MALLOC: allocated >= 57 MB of core memory
 @ACES_INIT: heap range 144875520 204875775
 @ACES_CACHE_INIT: allocated >= 2 MB of cache from core memory


ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.13 (      0.1)   0: 0: 0.04 (      0.0)   0: 0: 0.2 (      0.2)
@ACES2: Executing "xvscf_ks"
  @GETREC: record length mismatch for record COMPCLSS
           Actual:      8  Requested:      2
 @LIBRA3: Entered a3getbas.F
 
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
 
  @VSCF-I, There are 2 irreducible representations.

       Irrep        # of functions
         1                 44
         2                 28


  @VSCF-I, Parameters for SCF calculation: 
             SCF reference function:  RHF
       Maximum number of iterations:  150
          Full symmetry point group: C3v 
          Computational point group: C s 
             Initial density matrix: MOREAD                 
          SCF convergence tolerance: 10**(- 7)
       RPP convergence acceleration:  ON
               Latest start for RPP:   8
                          RPP order:   6

          Alpha population by irrep:   0    0
           Beta population by irrep:   0    0


 @INITGES-I, Routine entered. 
  @INITGES-I, Occupancies from core Hamiltonian:

          Alpha population by irrep:   4    1
           Beta population by irrep:   4    1


  --------------------------------------------------------------------
  Iteration         Total Energy            Largest Density Difference
  --------------------------------------------------------------------
       0            11.9541086066              0.0000000000E+00
       1           -37.7717871780              0.1360392471E+02
       2           -48.5547771257              0.4743725113E+02
       3           -50.2815394632              0.4719883486E+02
       4           -50.1268358930              0.3700535900E+02
       5           -50.7514361341              0.3699852962E+02
       6           -50.2843265594              0.3577688524E+02
       7           -50.8021371810              0.3577700217E+02
       8           -50.3033849040              0.3561855526E+02
       9           -56.1992227987              0.3585547257E+02
      10           -56.1976900396              0.2473540944E+00
      11           -56.2015206718              0.4250158457E-01
      12           -56.2029385088              0.2533049677E-01
      13           -56.2043445645              0.3103517520E-01
      14           -56.2049694630              0.1646565821E-01
      15           -56.2049429804              0.7563253824E-02
      16           -56.2049358816              0.2114320598E-02
      17           -56.2049368306              0.1386956373E-03
      18           -56.2049368727              0.1045826145E-04
      19           -56.2049367589              0.4688207159E-05
      20           -56.2049367872              0.1015587879E-05

  @VSCF-I, SCF has converged.

     E(SCF)=       -56.2049367865              0.5063021220E-07

     E(SCF)=       -56.2049367865              0.3481686606E-07

  @VSCF-I, Eigenvector printing suppressed.

 @PUTMOS-I, Writing converged MOs to NEWMOS. 
 @PUTMOS-I, Symmetry   1 Full  11 Partial   0
 @PUTMOS-I, Symmetry   2 Full   7 Partial   0
  @GETREC: record length mismatch for record FULLSTGP
           Actual:     10  Requested:      8
  @GETREC: record length mismatch for record FULLPERM
           Actual:     30  Requested:     24
  @GETREC: record length mismatch for record COMPSTGP
           Actual:     10  Requested:      8
  @GETREC: record length mismatch for record COMPPERM
           Actual:     40  Requested:      8
  @GETREC: record length mismatch for record COMPSYOP
           Actual:    144  Requested:     36
  @GETREC: record length mismatch for record FULLPERM
           Actual:     30  Requested:     24
  @PRJDEN-I, Analyzing reference function density.


      Trace of projected alpha density matrix =   0.389936553

      WARNING!!!  Alpha part of wavefunction is not symmetric!


  ORBITAL EIGENVALUES (ALPHA)  (1H = 27.2113957 eV)

       MO #        E(hartree)               E(eV)           FULLSYM    COMPSYM
       ----   --------------------   --------------------   -------   ---------
    1     1         -13.8266411486        -376.2422034967                   (1)
    2     2          -0.7662264586         -20.8500913610               A   (1)
    3    45          -0.4090088864         -11.1297026535               A   (2)
    4     3          -0.4089921027         -11.1292459461               A   (1)
    5     4          -0.2177870326          -5.9262891227               A   (1)
  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    6     5           0.0142400068           0.3874904594                   (1)
    7    46           0.0936194009           2.5475145640       E       A   (2)
    8     6           0.0936199942           2.5475307087       E       A   (1)
    9     7           0.2905432865           7.9060883375               A   (1)
   10    47           0.2905889184           7.9073300442               A   (2)
   11     8           0.3902065636          10.6180652057                   (1)
   12     9           0.4175084375          11.3609873003                   (1)
   13    48           0.4944573876          13.4548756311       E           (2)
   14    10           0.4944597530          13.4549399956       E           (1)
   15    49           0.5266793112          14.3316791451                   (2)
   16    11           0.5267064832          14.3324185332                   (1)
   17    12           0.5772214786          15.7070020608                   (1)
   18    50           0.8816151081          23.9899775622                   (2)
   19    13           0.9039300545          24.5971983981                   (1)
   20    14           0.9747447543          26.5241652147                   (1)
   21    51           0.9748067756          26.5258529019                   (2)
   22    15           1.1710823585          31.8667854557                   (1)
   23    52           1.1710948675          31.8671258428                   (2)
   24    16           1.2936866468          35.2030192588                   (1)
   25    53           1.4622615212          39.7901768696                   (2)
   26    17           1.4623253874          39.7919147579                   (1)
   27    54           1.4906187872          40.5618176575       E           (2)
   28    18           1.4906278805          40.5620650988       E           (1)
   29    19           1.5733160157          42.8121246650                   (1)
   30    55           2.2110246779          60.1650674118                   (2)
   31    20           2.2110737411          60.1664024914                   (1)
   32    21           2.2798725414          62.0385138693                   (1)
   33    22           2.9254442573          79.6054212839               A   (1)
   34    23           3.0418174045          82.7720970416                   (1)
   35    56           3.0717389041          83.5863028066                   (2)
   36    24           3.0718028516          83.5880429067                   (1)
   37    57           3.1351031635          85.3105327427                   (2)
   38    25           3.1513432039          85.7524469067                   (1)
   39    58           3.1514813123          85.7562050295                   (2)
   40    59           3.3177971150          90.2818901473                   (2)
   41    26           3.3611624329          91.4619209730       E           (1)
   42    60           3.3611657563          91.4620114072       E           (2)
   43    61           3.3935224613          92.3424825099                   (2)
   44    27           3.3936425500          92.3457502912                   (1)
   45    28           3.5763422729          97.3172647473       E           (1)
   46    62           3.5763447463          97.3173320524       E           (2)
   47    29           3.5975813676          97.8952101579                   (1)
   48    30           3.7340774895         101.6094601407                   (1)
   49    31           3.8262056927         104.1163971331                   (1)
   50    63           3.8263931719         104.1214987052                   (2)
   51    32           3.9908144705         108.5956317212       E           (1)
   52    64           3.9908203410         108.5957914653       E           (2)
   53    33           4.1256956933         112.2659380469       E           (1)
   54    65           4.1256981949         112.2660061212       E           (2)
   55    66           4.2121713357         114.6190609707                   (2)
   56    34           4.2313505718         115.1409547551                   (1)
   57    35           4.5783582486         124.5835179600                   (1)
   58    67           4.7428525537         129.0596375862                   (2)
   59    36           4.7429443589         129.0621357332                   (1)
   60    37           4.8503545285         131.9849163610                   (1)
   61    68           5.2533054418         142.9497731092               A   (2)
   62    38           5.2533523000         142.9510481882               A   (1)
   63    39           5.4065730797         147.1203994529                   (1)
   64    40           5.6757371072         154.4447283120                   (1)
   65    69           5.6757768217         154.4458089992                   (2)
   66    41           5.7763400843         157.1822757309                   (1)
   67    70           5.8463155122         159.0864047888       E           (2)
   68    42           5.8463178048         159.0864671740       E           (1)
   69    71           5.8515908535         159.2299541890                   (2)
   70    72           6.7049422851         182.4508376643                   (2)
   71    43           6.7049759725         182.4517543469                   (1)
   72    44          13.1938987696         359.0244002450                   (1)
  @CRAPSO, You need    200000 words of icore memory.
         , You need  14799998 words of dore memory.
         , You need  14999998 total words of memory.
         , You have  15000000 words of memory.
  @CHECKOUT-I, Total execution time :    1034.3101 seconds.
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
           Actual:      8  Requested:      2
 @LIBRA3: Entered a3getbas.F

The SCF nuclear-electron attraction energy =   -155.452081358023
The SCF kinetic energy                     =     55.747698246878
The SCF coulomb energy                     =     39.173401173603
The SCF exchange energy                    =     -7.628063454566
The SCF one electron energy                =    -99.704383111145
The SCF total energy                       =    -56.204936785481

Results using the SCF density

   Total density integrates to             :       9.999907663928 electrons
   Exch tot contrib : LDA (Slater, Xalpha) :      -6.897194845971
   Exch tot contrib : Becke                :      -7.658464386450
   Exch tot contrib : Perdew-Burke-Ernzerh :      -7.611576121732
   Exch tot contrib : Perdew-Wang 91       :      -7.639270772727
   Exch tot contrib : Exact Nonlocal Excha :      -7.628063454566
   Corr tot contrib : Vosko-Wilk-Nusair    :      -0.626048308588
   Corr tot contrib : Lee-Yang-Parr        :      -0.318642900785
   Corr tot contrib : Perdew-Burke-Ernzerh :      -0.315693379875
   Corr tot contrib : Perdew-Wang 91       :      -0.336724293273
   Corr tot contrib : Wilson-Levy          :      -0.369498838737
   Corr tot contrib : Wilson-Ivanov        :      -0.327288963138
   Hybr tot contrib : Becke III LYP        :      -7.968532565103
     
 Ehar from   KS   Calculation:         -48.576873330915 a.u.
 Ex   from   KS   Calculation:          -6.897194845971 a.u.
 Ec   from   KS   Calculation:          -0.626048308588 a.u.
 Exc  from   KS   Calculation:          -7.523243154559 a.u.
 Final Result from KS Calculation:     -56.100116485474 a.u.

  @CRAPSO, You need    200000 words of icore memory.
         , You need  14799998 words of dore memory.
         , You need  14999998 total words of memory.
         , You have  15000000 words of memory.
  @CHECKOUT-I, Total execution time :      43.0700 seconds.


 @ACES2: The ACES2 program has completed successfully in  932 seconds.


