
              ****************************************************
              * ACES : Advanced Concepts in Electronic Structure *
              *                    Ver. 2.7.0                    *
              ****************************************************

                             Quantum Theory Project
                             University of Florida
                             Gainesville, FL  32611
 @ACES2: Executing "rm -f FILES"
 @ACES2: Executing "xjoda"
 The value of HAVEGEOM record:                     0
 The value of PES_SCAN record:                     0
 The OPTARC is here? F F
 Do we have geom? F
 Hessian in JOBARC?                    -1


                      ACES STATE VARIABLE REGISTRATION LOG
     ----------------------------------------------------------------------
     The ACES State Variables were initialized.
     Updating (REFERENCE,RHF): 'RHF' converts to 'RHF' -> 0
     Updating (FOCK,AO): 'AO' converts to 'AO' -> 1
     Updating (BASIS,STO-3G): String copy succeeded.
     Updating (MEMORY_SIZE,2GB): 15000000 -> 268435456
     Updating (SYMMETRY,OFF): 'OFF' converts to 'OFF' -> 1
     Updating (SPHERICAL,ON): 'ON' converts to 'ON' -> 1
     Updating (SCF_TYPE,KS): 'KS' converts to 'KS' -> 1
     Updating (MULTIPLICTY,1): 1 -> 1
     ----------------------------------------------------------------------


                       ACES STATE VARIABLE VALIDATION LOG
     ----------------------------------------------------------------------
     Updating (COORDINATES,INTERNAL): 'INTERNAL' converts to 'INTERNAL' -> 0
     Updating (GEOM_OPT,partial): 'partial' converts to 'PARTIAL' -> 1
 A geometry for frequency {0,1} present:                    -1
     Updating (OPT_METHOD,manr): 'manr' converts to 'MANR' -> 3
     Updating (PROGRAM,aces2): 'aces2' converts to 'ACES2' -> 2
     Updating (GRAD_CALC,analytical): 'analytical' converts to 'ANALYTICAL' -> 1
     Updating (DERIV_LEV,first): 'first' converts to 'ON' -> 1
 The value of grad calc                     1
 The Deriv. level, havegeom?, and vib? :                     1
                     0                     1                     0
 The derivative level @exit                     1
     Updating (CC_EXTRAPOL,DIIS): 'DIIS' converts to 'DIIS' -> 1
 First order props and grad flags @enter:                     0
                     1
 First order props, deriv level and grad_calc @exit:                     0
                     1                     1
     Updating (ORBITALS,standard): 'standard' converts to 'STANDARD' -> 0
     Updating (PERT_ORB,standard): 'standard' converts to 'STANDARD' -> 0
     Updating (ESTATE_TOL,5): -1 -> 5
     Updating (HBARABCD,OFF): 'OFF' converts to 'OFF' -> 1
     Updating (HBARABCI,OFF): 'OFF' converts to 'OFF' -> 1
     Updating (ABCDFULL,ON): 'ON' converts to 'ON' -> 1
     Updating (CACHE_RECS,102): -1 -> 102
     Updating (FILE_RECSIZ,131072): -1 -> 131072
     ----------------------------------------------------------------------

          ASV#   ASV KEY DEFINITION =    CURRENT [   DEFAULT] UNITS
          ------------------------------------------------------------
            1:                PRINT =          0 [         0] 
            2:            CALCLEVEL =          0 [         0] 
            3:            DERIV_LEV =          1 [        -1] 
            4:              CC_CONV =          7 [         7] (tol)

            5:             SCF_CONV =          7 [         7] (tol)
            6:            XFORM_TOL =         11 [        11] (tol)
            7:            CC_MAXCYC =          0 [         0] cycles
            8:           LINDEP_TOL =          8 [         8] 

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
           36:          MEMORY_SIZE =  268435456 [  15000000] Words
           37:          FILE_RECSIZ =     131072 [        -1] Words

           38:                NONHF =          0 [         0] 
           39:             ORBITALS =          0 [        -1] 
           40:          SCF_EXPSTAR =          8 [         8] 
           41:          LOCK_ORBOCC =          0 [         0] 

           42:          FILE_STRIPE =          0 [         0] 
           43:               DOHBAR =          0 [         0] 
           44:           CACHE_RECS =        102 [        -1] 
           45:                GUESS =          0 [         0] 

           46:           JODA_PRINT =          0 [         0] 
           47:           OPT_METHOD =          3 [         0] 
           48:          CONVERGENCE =          4 [         4] H/bohr
           49:          EIGENVECTOR =          1 [         1] 

           50:              NEGEVAL =          2 [         2] 
           51:          CURVILINEAR =          0 [         0] 
           52:          STP_SIZ_CTL =          0 [         0] 
           53:             MAX_STEP =        300 [       300] millibohr

           54:            VIBRATION =          0 [         0] 
           55:            EVAL_HESS =         -1 [        -1] # of cyc.
           56:            INTEGRALS =          1 [         1] 
           57:          FD_STEPSIZE =          0 [         0] 10-4 bohr

           58:               POINTS =          0 [         0] 
           59:          CONTRACTION =          1 [         1] 
           60:             SYMMETRY =          1 [         2] 
           62:            SPHERICAL =          1 [         1] 

           63:          RESET_FLAGS =          0 [         0] 
           64:             PERT_ORB =          0 [         2] 
           65:             GENBAS_1 =          0 [         0] 
           66:             GENBAS_2 =          0 [         0] 

           67:             GENBAS_3 =          0 [         0] 
           68:          COORDINATES =          0 [         3] 
           69:            CHECK_SYM =          1 [         1] 
           70:            SCF_PRINT =          0 [         0] 

           71:                  ECP =          0 [         0] 
           72:              RESTART =          1 [         1] 
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

           96:                 FOCK =          1 [         1] 
           97:          ESTATE_MAXC =         20 [        20] 
           98:           ESTATE_TOL =          5 [        -1] (tol)
           99:            TURBOMOLE =          0 [         0] 

          100:           GAMMA_ABCD =          0 [         0] 
          101:            ZETA_TYPE =          1 [         1] 
          102:          ZETA_MAXCYC =         50 [        50] 
          103:             RESRAMAN =          0 [         0] 

          104:                  PSI =          0 [         0] 
          105:             GEOM_OPT =          1 [         0] 
          106:             EXTERNAL =          2 [         2] 
          107:          HESS_UPDATE =          2 [         0] 

          108:         INIT_HESSIAN =          0 [         0] 
          109:          EXTRAPOLATE =          0 [         0] 
          201:              EA_CALC =          0 [         0] 
          203:                 TDHF =          0 [         0] 

          204:           FUNCTIONAL =          4 [         4] 
          205:           EOM_MAXCYC =         50 [        50] cycles
          206:              EOMPROP =          0 [         0] 
          207:             ABCDFULL =          1 [         0] 

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

          225:              NOREORI =          2 [         2] 
          227:               KS_POT =          0 [         0] 
          228:             DIP_CALC =          0 [         0] 
          230:             DEA_CALC =          0 [         0] 

          232:              PROGRAM =          2 [         0] 
          233:                CCR12 =          0 [         0] 
          234:            EOMXFIELD =          0 [         0] x 10-6
          235:            EOMYFIELD =          0 [         0] x 10-6

          236:            EOMZFIELD =          0 [         0] x 10-6
          237:              INSERTF =          0 [         0] 
          238:            GRAD_CALC =          1 [         0] 
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
          256:             OOMP_CYC =         50 [        50] 

          257:            DKH_ORDER =          0 [         0] 
          258:           UNCONTRACT =          0 [         0] 
          ------------------------------------------------------------

                         ACES STATE VARIABLES (STRINGS)
          ------------------------------------------------------------
          BASIS = STO-3G
          OCCUPATION = [ESTIMATED BY SCF]
          ------------------------------------------------------------
 -----After call to Entry-----
 The FNDDONE in Geopt:                     1
 The gradient calcs:                     1
 Hessian calc.,Geo. opt. and iarch: F T                     0
 Entering Fetchz IARCH .NE. 1                     0
   3 entries found in Z-matrix 
 
 Job Title : ANALYTICAL GEOMETRY OPTIMIZATION H2O/3-21G*
 
 @FETCHZ: The ordering of the first two atoms in the Z-matrix has been changed.
  There are                      2  unique internal coordinates.
  Of these,                      2  will be optimized.
   User supplied Z-matrix: 
--------------------------------------------------------------------------------
       SYMBOL    BOND      LENGTH    ANGLE     ANGLE     DIHED     ANGLE
                  TO      (ANGST)    WRT      (DEG)      WRT      (DEG)
--------------------------------------------------------------------------------
        H    
        O          1         R    
        H          2         R          1        A    
                  *Initial values for internal coordinates* 
                      Name             Value
                        R               0.956900
                        A             104.517000
--------------------------------------------------------------------------------
 
 @GEOPT, PASS1,  PS1EXIST=:                     0 F
 Before the call to GMETRY, XYZIN F
 In Gmetry; ncycle, ix, and ipost_vib
                     0                     0                     0
         -----------------------------------------------
         Cartesian coordinates corresponding to internal 
                 coordinate input (Angstroms) 
 ----------------------------------------------------------------
 Z-matrix   Atomic            C o o r d i n a t e s
  Symbol    Number           X              Y              Z
 ----------------------------------------------------------------
     H         1         0.00000000     0.00000000     0.00000000
     O         8         0.00000000     0.00000000     0.95690000
     H         1         0.92634935     0.00000000     1.19676349
 ----------------------------------------------------------------
The number of atoms                      3
 
 The Cartesians before translations to CM
   0.00000   0.00000   0.00000
   0.00000   0.00000   1.80828
   1.75055   0.00000   2.26156
 
 The Cartesians in center of mass coords
  -0.0979561   0.0000000  -1.7324563
  -0.0979561   0.0000000   0.0758225
   1.6525903   0.0000000   0.5290988
 
   Rotational constants (in cm-1): 
     9.53321       14.60577       27.44965
 
 The symmetry processing begins
 The Cartesians before entering symmetry auto
  -1.42995   0.98297   0.00000
   0.00000  -0.12387   0.00000
   1.42995   0.98297   0.00000
   0.00000   0.00000   0.00000
   0.00000   0.00000   0.00000
   0.00000   0.00000   0.00000
   Principal axis orientation for molecule:
       -1.429951412226    0.982974687146    0.000000000000
        0.000000000000   -0.123872715017    0.000000000000
        1.429951412226    0.982974687146    0.000000000000
********************************************************************************
   The full molecular point group is C2v .
   The largest Abelian subgroup of the full molecular point group is C2v .
   The computational point group is C1  .
********************************************************************************
 @-SYM_AUTO The Orientation matrices: ORIEN2 and ORIENT
 

                 COLUMN   1       COLUMN   2       COLUMN   3
 ROW   1      0.61209997214    0.79078038930    0.00000000000
 ROW   2      0.00000000000    0.00000000000    1.00000000000
 ROW   3      0.79078038930   -0.61209997214    0.00000000000
 

                 COLUMN   1       COLUMN   2       COLUMN   3
 ROW   1      1.00000000000    0.00000000000    0.00000000000
 ROW   2      0.00000000000    1.00000000000    0.00000000000
 ROW   3      0.00000000000    0.00000000000    1.00000000000
 
 @-SYM_AUTO The variables in /COORD/ common block
   0.00000  -1.42995   0.98297
   0.00000   0.00000  -0.12387
   0.00000   1.42995   0.98297
 
 The NEWQ array 
  -1.42995
   0.98297
   0.00000
   0.00000
  -0.12387
   0.00000
   1.42995
   0.98297
   0.00000
 
 
 @-SYM_AUTO ORIENT2=ORIEN2xORIENT

                 COLUMN   1       COLUMN   2       COLUMN   3
 ROW   1      0.61209997214    0.79078038930    0.00000000000
 ROW   2      0.00000000000    0.00000000000    1.00000000000
 ROW   3      0.79078038930   -0.61209997214    0.00000000000
 --------------B
    0.00000    1.42995    0.98297
    0.00000    0.00000   -0.12387
    0.00000   -1.42995    0.98297
 --------------B
    0.00000    1.42995    0.98297
    0.00000    0.00000   -0.12387
    0.00000   -1.42995    0.98297
 --------------B
    0.00000   -1.42995    0.98297
    0.00000    0.00000   -0.12387
    0.00000    1.42995    0.98297
 --------------B
    0.00000   -1.42995    0.98297
    0.00000    0.00000   -0.12387
    0.00000    1.42995    0.98297
 --------------B
    0.00000   -1.42995    0.98297
    0.00000    0.00000   -0.12387
    0.00000    1.42995    0.98297
 
 The Cartesians, internals and connec. @-buildb
 The NATOMS                     3
   0.00000  -1.42995   0.98297
   0.00000   0.00000  -0.12387
   0.00000   1.42995   0.98297
 
   0.00000   0.00000   0.00000
   1.80828   0.00000   0.00000
   1.80828   1.82417   0.00000
    0    0    0    1
    0    0    2    1
    0
  B MATRIX 
[  1,  1]  0.000000 [  1,  2] -0.790780 [  1,  3]  0.612100 [  1,  4]  0.000000
[  1,  5]  0.790780 [  1,  6] -0.612100 [  1,  7]  0.000000 [  1,  8]  0.000000
[  1,  9]  0.000000 [  2,  1]  0.000000 [  2,  2]  0.000000 [  2,  3]  0.000000
[  2,  4]  0.000000 [  2,  5] -0.790780 [  2,  6] -0.612100 [  2,  7]  0.000000
[  2,  8]  0.790780 [  2,  9]  0.612100 [  3,  1]  0.000000 [  3,  2] -0.338499
[  3,  3] -0.437311 [  3,  4]  0.000000 [  3,  5]  0.000000 [  3,  6]  0.874622
[  3,  7]  0.000000 [  3,  8]  0.338499 [  3,  9] -0.437311 [
  Bt MATRIX
[  1,  1]  0.000000 [  1,  2]  0.000000 [  1,  3]  0.000000 [  2,  1] -0.790780
[  2,  2]  0.000000 [  2,  3] -0.338499 [  3,  1]  0.612100 [  3,  2]  0.000000
[  3,  3] -0.437311 [  4,  1]  0.000000 [  4,  2]  0.000000 [  4,  3]  0.000000
[  5,  1]  0.790780 [  5,  2] -0.790780 [  5,  3]  0.000000 [  6,  1] -0.612100
[  6,  2] -0.612100 [  6,  3]  0.874622 [  7,  1]  0.000000 [  7,  2]  0.000000
[  7,  3]  0.000000 [  8,  1]  0.000000 [  8,  2]  0.790780 [  8,  3]  0.338499
[  9,  1]  0.000000 [  9,  2]  0.612100 [  9,  3] -0.437311 [
  G MATRIX 
[  1,  1]  2.000000 [  2,  1] -0.250667 [  2,  2]  2.000000 [  3,  1] -0.535356
[  3,  2] -0.535356 [  3,  3]  1.376609 [
  G MATRIX DETERMINANT   0.000000000000000E+000
  INVERSE G MATRIX 
[  1,  1]  0.597268 [  2,  1]  0.152955 [  2,  2]  0.597268 [  3,  1]  0.291758
[  3,  2]  0.291758 [  3,  3]  0.953349 [
  A MATRIX
[  1,  1]  0.000000 [  1,  2]  0.000000 [  1,  3]  0.000000 [  2,  1] -0.571067
[  2,  2] -0.219713 [  2,  3] -0.553424 [  3,  1]  0.237999 [  3,  2] -0.033965
[  3,  3] -0.238325 [  4,  1]  0.000000 [  4,  2]  0.000000 [  4,  3]  0.000000
[  5,  1]  0.351354 [  5,  2] -0.351354 [  5,  3]  0.000000 [  6,  1] -0.204033
[  6,  2] -0.204033 [  6,  3]  0.476650 [  7,  1]  0.000000 [  7,  2]  0.000000
[  7,  3]  0.000000 [  8,  1]  0.219713 [  8,  2]  0.571067 [  8,  3]  0.553424
[  9,  1] -0.033965 [  9,  2]  0.237999 [  9,  3] -0.238325 [
 
  R VECTOR IN ANALYZE 
 
      0.000000000000      0.000000000000      0.000000000000
 
   1   8   1
--------------------------------------------------------------------------------
   Analysis of internal coordinates specified by Z-matrix 
--------------------------------------------------------------------------------
   *The nuclear repulsion energy is    9.19786 a.u.
  R      A      (*
   *There are   2 degrees of freedom within the tot. symm. molecular subspace.
   *Z-matrix requests optimization of   2 coordinates.
   *The following   0 parameters can have non-zeroderivatives within the 
    totally symmetric subspace:

   *The following   2 parameters are to be optimized:
             R    [  1]  A    [  3]

--------------------------------------------------------------------------------
         -----------------------------------------------
         Cartesian coordinates corresponding to internal 
                 coordinate input (Angstroms) 
 ----------------------------------------------------------------
 Z-matrix   Atomic            C o o r d i n a t e s
  Symbol    Number           X              Y              Z
 ----------------------------------------------------------------
     H         1         0.00000000    -0.75669775     0.52016784
     O         8         0.00000000     0.00000000    -0.06555062
     H         1         0.00000000     0.75669775     0.52016784
 ----------------------------------------------------------------
   Interatomic distance matrix (Angstroms) 
 
                 H             O             H    
                [ 1]        [ 2]        [ 3]
  H    [ 1]     0.00000
  O    [ 2]     0.95690     0.00000
  H    [ 3]     1.51340     0.95690     0.00000

   Interatomic angles (degrees) 


  O [ 2]-H [ 1]-H [ 3]   37.74150
  H [ 1]-O [ 2]-H [ 3]  104.51700

  H [ 1]-H [ 3]-O [ 2]   37.74150

      3 interatomic angles printed.
 Creating Archive file
 
 The COORD COMMON BLOCK AT ARCHIVE
   0.00000  -1.42995   0.98297
   0.00000   0.00000  -0.12387
   0.00000   1.42995   0.98297
 
   0.00000   0.00000   0.00000
   1.80828   0.00000   0.00000
   1.80828   1.82417   0.00000
 
 @-ARCHIVE, the number of opt. cycles                     0
   1.80828   1.80828   1.82417
 Finite diffs; ignore and geomopt vars: F T
 The vib calc. related varrs:
 iflags(h_IFLAGS_vib):                     0
 The finite diffs: F
 iflags2(h_IFLAGS2_geom_opt):                     1
 first run of popt num. frq:                     0
 @ACES2: Executing "xvmol"
 One- and two-electron integrals over symmetry-adapted AOs are calculated.\n
 @READIN: Spherical harmonics are used.
  @READIN-I, Nuclear repulsion energy :    9.1978553372 a.u.
  required memory for a1 array                4451690  words 
  required memory for a2 array                   9156  words 
  @TWOEL-I,        228 integrals of symmetry type  I I I I
  @TWOEL-I, Total number of 2-e integrals        228.


ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.01 (      0.0)   0: 0: 0.01 (      0.0)   0: 0: 0.1 (      0.1)
 @ACES2: Executing "xvmol2ja"


ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.00 (      0.0)   0: 0: 0.00 (      0.0)   0: 0: 0.0 (      0.0)
 @ACES2: Executing "xvscf"
@VSCF: There are    7 functions in the AO basis.

       There are    1 irreducible representations.


       Irrep        # of functions
         1                  7
 
  @VSCF-I, Parameters for SCF calculation: 
             SCF reference function:  RHF
       Maximum number of iterations:  150
          Full symmetry point group: C2v 
          Computational point group: C1  
             Initial density matrix: CORE                   
          SCF convergence tolerance: 10**(- 7)
           Convergence acceleration: RPP
               Latest start for RPP:   8
                          RPP order:   6

  @SYMSIZ-I, Memory information:   13374099 words required.
  @SYMSIZ-I, Fock matrices are constructed from AO integral file.
  Occupancies from core Hamiltonian:

          Alpha population by irrep:   5
           Beta population by irrep:   5


The Alpha and Beta occupation vector
                     5
 
                     5
 
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
           potradpts     integer                    94                    50
           OEPSLATER     integer                     5
             overind     integer                     3
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
               ksmem     integer             150000000               1500000
                func      string  BECKE,LYP             none                
               kspot      string  BECKE,LYP             hf                  
              cutoff     integer                    12
============================================================================
 
  --------------------------------------------------------------------
  Iteration         Total Energy            Largest Error in FDS-SDF
  --------------------------------------------------------------------
       0           -63.1863869559              0.5892286227D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       1           -70.5248345024              0.6044202737D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       2           -74.3116855408              0.4820018181D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       3           -74.4774849425              0.4765546859D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       4           -74.6367127714              0.3779099731D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       5           -74.6605245616              0.3883735433D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       6           -74.7398651152              0.3211166154D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       7           -74.7395854095              0.3357154078D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       8           -74.7913143583              0.2853338346D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       9           -74.9620985226              0.1743119359D-02
The Alpha and Beta occupation vector
                     5
 
                     5
      10           -74.9620599759              0.2941813186D-03
The Alpha and Beta occupation vector
                     5
 
                     5
      11           -74.9620504841              0.3696359860D-04
The Alpha and Beta occupation vector
                     5
 
                     5
      12           -74.9620515774              0.2260889126D-06
The Alpha and Beta occupation vector
                     5
 
                     5
      13           -74.9620515707              0.6854092294D-09
The Alpha and Beta occupation vector
                     5
 
                     5

@VSCF: SCF has converged.

     E(SCF)=       -74.9620515706

 The AOBASMOS? F
 Vib calc & AOBASMOS?                     0 F
 Creating The AOBASMOS file


      Trace of projected alpha density matrix =   1.000000000

      Alpha part of wavefunction is symmetric.


  ORBITAL EIGENVALUES (ALPHA)  (1H = 27.2113957 eV)

       MO #        E(hartree)               E(eV)           FULLSYM    COMPSYM
       ----   --------------------   --------------------   -------   ---------
    1     1         -18.4830037041        -502.9483275160      A1         A (1)
    2     2          -0.8370372749         -22.7769525038      A1         A (1)
    3     3          -0.3837734398         -10.4430109307      B1         A (1)
    4     4          -0.1581613318          -4.3037905850      A1         A (1)
    5     5          -0.0658839844          -1.7927951693      B2         A (1)
  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    6     6           0.3089721553           8.4075635779      A1         A (1)
    7     7           0.4211599688          11.4603505634      B1         A (1)




ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.58 (      0.6)   0: 0: 0.89 (      0.9)   0: 0: 1.5 (      1.5)
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
           potradpts     integer                    94                    50
           OEPSLATER     integer                     5
             overind     integer                     3
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
               ksmem     integer             150000000               1500000
                func      string  BECKE,LYP             none                
               kspot      string  BECKE,LYP             hf                  
              cutoff     integer                    12
============================================================================
 

The SCF nuclear-electron attraction energy =   -196.930510980041
The SCF kinetic energy                     =     74.590152573698
The SCF coulomb energy                     =     47.281321119874
The SCF exchange energy                    =     -9.100869621340
The SCF one electron energy                =   -122.340358406343
The SCF total energy                       =    -74.962051570617

Results using the SCF density

   Total density integrates to             :       9.999986457588 electrons
   Exch tot contrib : LDA (Slater, Xalpha) :      -8.198293100036
   Exch tot contrib : Becke                :      -9.076775894800
   Exch tot contrib : Perdew-Burke-Ernzerh :      -9.020606846504
   Exch tot contrib : Perdew-Wang 91       :      -9.051555493794
   Exch tot contrib : Exact Nonlocal Excha :      -9.100869621340
   Corr tot contrib : Vosko-Wilk-Nusair    :      -0.671955068971
   Corr tot contrib : Lee-Yang-Parr        :      -0.339192217896
   Corr tot contrib : Perdew-Burke-Ernzerh :      -0.343479099470
   Corr tot contrib : Perdew-Wang 91       :      -0.366168915107
   Corr tot contrib : Wilson-Levy          :      -0.403855709071
   Corr tot contrib : Wilson-Ivanov        :      -0.350538807739
   Hybr tot contrib : Becke III LYP        :      -9.413733176127
     
 Ehar from   KS   Calculation:         -65.861181949277 a.u.
 Ex   from   KS   Calculation:          -9.076775894800 a.u.
 Ec   from   KS   Calculation:          -0.339192217896 a.u.
 Exc  from   KS   Calculation:          -9.415968112696 a.u.
 Final Result from KS Calculation:     -75.277150061973 a.u.

 Ehar with the SCF Dens and Func:      -65.861181949277 a.u.
 Ex   with the SCF Dens and Func:       -9.076775894800 a.u.
 Ec   with the SCF Dens and Func:       -0.339192217896 a.u.
 Exc  with the SCF Dens and Func:       -9.415968112696 a.u.
 Final Result with the SCF Dens and Func:    -75.277150061973 a.u.

  @CHECKOUT-I, Total execution time :       0.3409 seconds.
ACES2: Total elapsed time is        2.7 seconds
 @ACES2: Executing "xvdint"
 One- and two-electron integral derivatives are calculated
 for RHF gradients and dipole moments.
 Spherical gaussians are used.

  Cartesian Coordinates
  ---------------------

  Total number of coordinates:  9


   1   H #1     x      0.0000000000
   2            y     -1.4299514122
   3            z      0.9829746871

   4   O #2     x      0.0000000000
   5            y      0.0000000000
   6            z     -0.1238727150

   7   H #3     x      0.0000000000
   8            y      1.4299514122
   9            z      0.9829746871

 Translational invariance is used.
 
 Entering the ECP_DGRAD

                    Nuclear attraction energy gradient
                    ----------------------------------

 H #1       0.0000000000           -4.3156407555            3.2838711091
 O #2       0.0000000000            0.0000000000           -6.5677422182
 H #3       0.0000000000            4.3156407555            3.2838711091



                       Two electron energy gradient
                       ----------------------------

 H #1       0.0000000000            1.9768031331           -1.5626393205
 O #2       0.0000000000            0.0000000000            3.1252786411
 H #3       0.0000000000           -1.9768031331           -1.5626393206



                    Nuclear repulsion energy gradient
                    ---------------------------------

 H #1       0.0000000000            2.0569702027           -1.4975508088
 O #2       0.0000000000            0.0000000000            2.9951016175
 H #3       0.0000000000           -2.0569702027           -1.4975508088



                      Kinetic energy energy gradient
                      ------------------------------

 H #1       0.0000000000            0.4334950193           -0.3255012092
 O #2       0.0000000000            0.0000000000            0.6510024183
 H #3       0.0000000000           -0.4334950193           -0.3255012092



                     Renormalization energy gradient
                     -------------------------------

 H #1       0.0000000000            0.0930184862           -0.1086835184
 O #2       0.0000000000            0.0000000000            0.2173670369
 H #3       0.0000000000           -0.0930184862           -0.1086835184



                  One elctron Molecular energy gradient
                  -------------------------------------

 H #1       0.0000000000           -1.7321570472            1.3521355727
 O #2       0.0000000000            0.0000000000           -2.7042711455
 H #3       0.0000000000            1.7321570472            1.3521355728



                        Molecular energy gradient
                        -------------------------

 H #1       0.0000000000            0.2446460859           -0.2105037478
 O #2       0.0000000000            0.0000000000            0.4210074956
 H #3       0.0000000000           -0.2446460859           -0.2105037478




ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.01 (      0.0)   0: 0: 0.00 (      0.0)   0: 0: 0.0 (      0.0)
 @ACES2: Executing "xvksdint"
 fffg
 
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
           potradpts     integer                    94                    50
           OEPSLATER     integer                     5
             overind     integer                     3
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
               ksmem     integer             150000000               1500000
                func      string  BECKE,LYP             none                
               kspot      string  BECKE,LYP             hf                  
              cutoff     integer                    12
============================================================================
 
  idd=                     1
 The new gradient:
 atm    x               y               z
  1  1    0.0000000000    0.2446460859   -0.2105037478
  8  2    0.0000000000    0.0000000000    0.4210074956
  1  3    0.0000000000   -0.2446460859   -0.2105037478
                     1                     1                     1
                     2                     8                     2
                     3                     1                     3
 t=                     1                     1
  2.792392415792886E-016 -0.186946496352083       0.147368629198233     
 t=                     2                     2
 -2.012105167614231E-015 -1.630029444754655E-012 -0.294737258399256     
 t=                     3                     3
  1.727013555775450E-015  0.186946496353714       0.147368629201020     
 rotin= -2.796515395910919E-015  1.911193425740976E-012  8.561669335253221E-015
 tran= -5.852370259492889E-018 -5.852370259492889E-018 -2.081668171172169E-015
  1  1    0.0000000000    0.0576995896   -0.0631351186
  8  2    0.0000000000    0.0000000000    0.1262702372
  1  3    0.0000000000   -0.0576995896   -0.0631351186
  @CHECKOUT-I, Total execution time :       0.5699 seconds.


ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.43 (      0.4)   0: 0: 0.57 (      0.6)   0: 0: 1.0 (      1.0)
ACES2: Total elapsed time is        3.7 seconds
 @ACES2: Executing "xjoda"
 Finite diffs; ignore and geomopt vars: T T
 The vib calc. related varrs:
 iflags(h_IFLAGS_vib):                     0
 The finite diffs: F
 iflags2(h_IFLAGS2_geom_opt):                     1
 first run of popt num. frq:                     0
 The value of HAVEGEOM record:                     0
 The value of PES_SCAN record:                     0
 The OPTARC is here? T F
 Do we have geom? F
 Hessian in JOBARC?                    -1
Files copied to SAVEDIR/CURRENT: ZMAT.BAS JOBARC JAINDX OPTARC !OPTARCBK !DIPDER
 A geometry for frequency {0,1} present:                     1
 The value of grad calc                     1
 The Deriv. level, havegeom?, and vib? :                     1
                     0                     1                     0
 The derivative level @exit                     1
 First order props and grad flags @enter:                     0
                     1
 First order props, deriv level and grad_calc @exit:                     0
                     1                     1
 @-Entry, the optimization cycle                     0
   JODA beginning optimization cycle #  1.
 -----After call to Entry-----
 The FNDDONE in Geopt:                     1
 The gradient calcs:                     1
 Hessian calc.,Geo. opt. and iarch: F T                     1
 Entering Fetchz IARCH .NE. 1                     1
 
 @GEOPT, PASS1,  PS1EXIST=:                     0 F
 
 Start Reading Archive file NXM6:                      3
 The COORD COMMON BLOCK AT RETRIEVE BEFORE DECOMPRESS
   0.00000  -1.42995   0.98297
   0.00000   0.00000  -0.12387
   0.00000   1.42995   0.98297
 
   1.80828   1.80828   1.82417
 The COORD COMMON BLOCK AT RETRIEVE
   0.00000  -1.42995   0.98297
   0.00000   0.00000  -0.12387
   0.00000   1.42995   0.98297
 
   0.00000   0.00000   0.00000
 Data read from OPTARC file
 
 The COORD COMMON BLOCK/AFTER RETRIVE
   0.00000  -1.42995   0.98297
   0.00000   0.00000  -0.12387
   0.00000   1.42995   0.98297
 
   0.00000   0.00000   0.00000
 
 
 The Cartesians, internals and connec. @-buildb
 The NATOMS                     3
   0.00000  -1.42995   0.98297
   0.00000   0.00000  -0.12387
   0.00000   1.42995   0.98297
 
   0.00000   0.00000   0.00000
   1.80828   0.00000   0.00000
   1.80828   1.82417   0.00000
    0    0    0    1
    0    0    2    1
    0
  B MATRIX 
[  1,  1]  0.000000 [  1,  2] -0.790780 [  1,  3]  0.612100 [  1,  4]  0.000000
[  1,  5]  0.790780 [  1,  6] -0.612100 [  1,  7]  0.000000 [  1,  8]  0.000000
[  1,  9]  0.000000 [  2,  1]  0.000000 [  2,  2]  0.000000 [  2,  3]  0.000000
[  2,  4]  0.000000 [  2,  5] -0.790780 [  2,  6] -0.612100 [  2,  7]  0.000000
[  2,  8]  0.790780 [  2,  9]  0.612100 [  3,  1]  0.000000 [  3,  2] -0.338499
[  3,  3] -0.437311 [  3,  4]  0.000000 [  3,  5]  0.000000 [  3,  6]  0.874622
[  3,  7]  0.000000 [  3,  8]  0.338499 [  3,  9] -0.437311 [
  Bt MATRIX
[  1,  1]  0.000000 [  1,  2]  0.000000 [  1,  3]  0.000000 [  2,  1] -0.790780
[  2,  2]  0.000000 [  2,  3] -0.338499 [  3,  1]  0.612100 [  3,  2]  0.000000
[  3,  3] -0.437311 [  4,  1]  0.000000 [  4,  2]  0.000000 [  4,  3]  0.000000
[  5,  1]  0.790780 [  5,  2] -0.790780 [  5,  3]  0.000000 [  6,  1] -0.612100
[  6,  2] -0.612100 [  6,  3]  0.874622 [  7,  1]  0.000000 [  7,  2]  0.000000
[  7,  3]  0.000000 [  8,  1]  0.000000 [  8,  2]  0.790780 [  8,  3]  0.338499
[  9,  1]  0.000000 [  9,  2]  0.612100 [  9,  3] -0.437311 [
  G MATRIX 
[  1,  1]  2.000000 [  2,  1] -0.250667 [  2,  2]  2.000000 [  3,  1] -0.535356
[  3,  2] -0.535356 [  3,  3]  1.376609 [
  G MATRIX DETERMINANT   0.000000000000000E+000
  INVERSE G MATRIX 
[  1,  1]  0.597268 [  2,  1]  0.152955 [  2,  2]  0.597268 [  3,  1]  0.291758
[  3,  2]  0.291758 [  3,  3]  0.953349 [
  A MATRIX
[  1,  1]  0.000000 [  1,  2]  0.000000 [  1,  3]  0.000000 [  2,  1] -0.571067
[  2,  2] -0.219713 [  2,  3] -0.553424 [  3,  1]  0.237999 [  3,  2] -0.033965
[  3,  3] -0.238325 [  4,  1]  0.000000 [  4,  2]  0.000000 [  4,  3]  0.000000
[  5,  1]  0.351354 [  5,  2] -0.351354 [  5,  3]  0.000000 [  6,  1] -0.204033
[  6,  2] -0.204033 [  6,  3]  0.476650 [  7,  1]  0.000000 [  7,  2]  0.000000
[  7,  3]  0.000000 [  8,  1]  0.219713 [  8,  2]  0.571067 [  8,  3]  0.553424
[  9,  1] -0.033965 [  9,  2]  0.237999 [  9,  3] -0.238325 [
 
 The COORD COMMON BLOCK/at opt. start:
 
   1.80828   1.80828   1.82417
 
 The Cart. Hessian in after reading from READGH
 IAVHES flags: Cartesian exact Hess. read                     0
                     0
 If the Hessian is from VIB=EXACT, the trans/rot
 contaminants are included. If VIB=FINDIF it is
 projected.

                 COLUMN   1       COLUMN   2       COLUMN   3       COLUMN   4
 ROW   1      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   2      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   3      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   4      0.00000000000   -0.79078038930    0.61209997214    0.00000000000
 ROW   5      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   6      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   7      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   8      0.00000000000   -0.33849867214   -0.43731109936    0.00000000000
 ROW   9      0.00000000000    0.00000000000    0.00000000000    0.00000000000

                 COLUMN   5       COLUMN   6       COLUMN   7       COLUMN   8
 ROW   1      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   2      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   3      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   4      0.79078038930   -0.61209997214    0.00000000000    0.00000000000
 ROW   5      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   6      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   7     -0.79078038930   -0.61209997214    0.00000000000    0.79078038930
 ROW   8      0.00000000000    0.87462219872    0.00000000000    0.33849867214
 ROW   9      0.00000000000    0.00000000000    0.00000000000    0.00000000000

                 COLUMN   9
 ROW   1      0.00000000000
 ROW   2      0.00000000000
 ROW   3      0.00000000000
 ROW   4      0.00000000000
 ROW   5      0.00000000000
 ROW   6      0.00000000000
 ROW   7      0.61209997214
 ROW   8     -0.43731109936
 ROW   9      0.00000000000
 The Cartesian Hessian after reading from GETICFCM
 
 The COORD COMMON BLOCK/at opt. start:
 
   1.80828   1.80828   1.82417
 
 
 @-CONVQF,  The outgoing internal gradients
   -0.0842727   -0.0842727    0.0264155
 
   Internal coordinate forces and energy gradients (atomic units): 
           R       dV/dR             R       dV/dR             R       dV/dR
 
 [R    ]   1.80828  -0.08427 [R    ]   1.80828  -0.08427 [A    ]   1.82417   0.02642
 
 
 Entering convhess:form an educated guess4hess
 
 Importing Cartesian Hessian F                     0
 The Partialy trans. Cart. Hessian in geopt-before
 TWIDLE

                 COLUMN   1       COLUMN   2       COLUMN   3
 ROW   1      1.00000000000    0.00000000000    0.00000000000
 ROW   2      0.00000000000    1.00000000000    0.00000000000
 ROW   3      0.00000000000    0.00000000000    0.25000000000

                 COLUMN   1       COLUMN   2
 ROW   1      1.00000000000    0.00000000000
 ROW   2      0.00000000000    0.00000000000
 The number of opt. cycles:                     0
 
   The eigenvectors of the Hessian matrix: 

                 COLUMN   1       COLUMN   2
 ROW   1      0.00000000000    1.00000000000
 ROW   2      1.00000000000    0.00000000000
   The eigenvalues of the Hessian matrix: 
 
     0.25000    1.00000
 
 There are  0 Negative Eigenvalues.
 
 The number of degs. of freed.           at NR:    2
 
 The unmolested step
    0.1192   -0.1057
   MANR scale factor for NR step is 1.166 for R    .
   Summary of Optimization Cycle: 
   The maximum unscaled step is:    0.17457.
   Scale factor set to:  1.00000.
   Forces are in hartree/bohr and hartree/radian.
   Parameter values are in Angstroms and degrees.
--------------------------------------------------------------------------
  Parameter     dV/dR           Step          Rold            Rnew
--------------------------------------------------------------------------
    R      -0.0842727082    0.0519974002    0.9569000000    1.0088974002
    A       0.0264155112   -6.0539892178  104.5170000000   98.4630107822
--------------------------------------------------------------------------
  Minimum force:  0.026415511 / RMS force:  0.062448653
 GMETRY starting with R
  1     1.906540
  2     1.906540
  3     1.718504
  4     0.000000
  5     0.000000
  6     0.000000
  7     0.000000
  8     0.000000
  9     0.000000
 In Gmetry; ncycle, ix, and ipost_vib
                     1                     0                     0
@GMETRY-I, Decompressing R.
 Updating structure...
 GMETRY using R vector
  1     0.000000
  2     0.000000
  3     0.000000
  4     1.906540
  5     0.000000
  6     0.000000
  7     1.906540
  8     1.718504
  9     0.000000
  @GMETRY-I, Cartesian coordinates before scaling:
      1     0.000000     0.000000     0.000000
      2     1.906540     0.000000     0.000000
      3     2.187127     1.885779     0.000000
  @GMETRY-I, Cartesian coordinates after scaling:
      1     0.000000     0.000000     0.000000
      2     1.906540     0.000000     0.000000
      3     2.187127     1.885779     0.000000
  @GMETRY-I, Cartesian coordinates from Z-matrix:
      1     0.000000     0.000000     0.000000
      2     0.000000     0.000000     1.906540
      3     1.885779     0.000000     2.187127
         -----------------------------------------------
         Cartesian coordinates corresponding to internal 
                 coordinate input (Angstroms) 
 ----------------------------------------------------------------
 Z-matrix   Atomic            C o o r d i n a t e s
  Symbol    Number           X              Y              Z
 ----------------------------------------------------------------
     H         1         0.00000000     0.00000000     0.00000000
     O         8         0.00000000     0.00000000     1.00889740
     H         1         0.99791160     0.00000000     1.15737773
 ----------------------------------------------------------------
    1.0078250000
   15.9949100000
    1.0078250000
The number of atoms                      3
 
 The Cartesians before translations to CM
   0.00000   0.00000   0.00000
   0.00000   0.00000   1.90654
   1.88578   0.00000   2.18713
 
 The Cartesians in center of mass coords
  -0.1055234   0.0000000  -1.8155555
  -0.1055234   0.0000000   0.0909842
   1.7802561   0.0000000   0.3715713
 
 After translation to center of mass coordinates 
       -0.105523409878    0.000000000000   -1.815555481950
       -0.105523409878    0.000000000000    0.090984155083
        1.780256070700    0.000000000000    0.371571286306
 
  @symmetry-i, Coordinates after  COM shift 
     -0.105523409878      0.000000000000     -1.815555481950
     -0.105523409878      0.000000000000      0.090984155083
      1.780256070700      0.000000000000      0.371571286306
 
  Inertia tensor
    3.59359    0.00000   -0.70618
    0.00000    6.97703    0.00000
   -0.70618    0.00000    3.38344
 Inertia tensor:
        3.593588108568    0.000000000000   -0.706184496582
        0.000000000000    6.977028334998    0.000000000000
       -0.706184496582    0.000000000000    3.383440226431
    Diagonalized inertia tensor:
        2.774555431784    0.000000000000    0.000000000000
        0.000000000000    4.202472903215    0.000000000000
        0.000000000000    0.000000000000    6.977028334998
    Eigenvectors of inertia tensor: 
        0.653004254231    0.757354239413    0.000000000000
        0.000000000000    0.000000000000    1.000000000000
        0.757354239413   -0.653004254231    0.000000000000
    Principal axis orientation for molecular system: 
       -1.443925876715    1.105646851677    0.000000000000
        0.000000000000   -0.139331642165    0.000000000000
        1.443925876715    1.105646851677    0.000000000000
   Rotational constants (in cm-1): 
     8.62803       14.32442       21.69645
  @SYMMETRY-I, Handedness of inertial frame: 1.00000
  @SYMMETRY-I, The symmetry group is 1-fold degenerate.
 
 The symmetry processing begins
 The Cartesians before entering symmetry auto
  -1.44393   1.10565   0.00000
   0.00000  -0.13933   0.00000
   1.44393   1.10565   0.00000
   0.00000   0.00000   0.00000
   0.00000   0.00000   0.00000
   0.00000   0.00000   0.00000
   Principal axis orientation for molecule:
       -1.443925876715    1.105646851677    0.000000000000
        0.000000000000   -0.139331642165    0.000000000000
        1.443925876715    1.105646851677    0.000000000000
   Reflection in plane  3 is a valid symmetry operation.
       -1.000000000000    0.000000000000    0.000000000000
        0.000000000000   -1.000000000000    0.000000000000
        0.000000000000    0.000000000000    1.000000000000
       -1.000000000000    0.000000000000    0.000000000000
        0.000000000000    1.000000000000    0.000000000000
        0.000000000000    0.000000000000   -1.000000000000
   Rotation about  2 is a valid symmetry operation 
   Reflection in plane  1 is a valid symmetry operation.
        1.000000000000    0.000000000000    0.000000000000
        0.000000000000   -1.000000000000    0.000000000000
        0.000000000000    0.000000000000   -1.000000000000
  Symmetry bits:    5   2   0
 @-SYM_AUTO The Orientation matrices: ORIEN2 and ORIENT
 

                 COLUMN   1       COLUMN   2       COLUMN   3
 ROW   1      0.65300425423    0.75735423941    0.00000000000
 ROW   2      0.00000000000    0.00000000000    1.00000000000
 ROW   3      0.75735423941   -0.65300425423    0.00000000000
 

                 COLUMN   1       COLUMN   2       COLUMN   3
 ROW   1      1.00000000000    0.00000000000    0.00000000000
 ROW   2      0.00000000000    1.00000000000    0.00000000000
 ROW   3      0.00000000000    0.00000000000    1.00000000000
 
 @-SYM_AUTO The variables in /COORD/ common block
   0.00000  -1.44393   1.10565
   0.00000   0.00000  -0.13933
   0.00000   1.44393   1.10565
 
 The NEWQ array 
  -1.44393
   1.10565
   0.00000
   0.00000
  -0.13933
   0.00000
   1.44393
   1.10565
   0.00000
 
 
 @-SYM_AUTO ORIENT2=ORIEN2xORIENT

                 COLUMN   1       COLUMN   2       COLUMN   3
 ROW   1      0.65300425423    0.75735423941    0.00000000000
 ROW   2      0.00000000000    0.00000000000    1.00000000000
 ROW   3      0.75735423941   -0.65300425423    0.00000000000
  @CHRTAB-I, Generated transformation matrices for   1 symmetry operations.
  @TSTOPS: Results of symmetry operation testing:
--------------------------------------------------------------------------------
   Sym. Op.         Class     Trace     Result 
--------------------------------------------------------------------------------
 --------------B
    0.00000   -1.44393    1.10565
    0.00000    0.00000   -0.13933
    0.00000    1.44393    1.10565
       1               1      3.00000   Passed
--------------------------------------------------------------------------------
 @TSTOPS:                      1  operations verified.
  @SYMUNQ-I, There are   3 orbits in the C1   point group.
------------------------------------------------------------------------
    Set    Site Group                       Members of set
------------------------------------------------------------------------
      1      C1                1
      2      C1                2
      3      C1                3
------------------------------------------------------------------------
  @CHRTAB-I, Generated transformation matrices for   4 symmetry operations.
  @TSTOPS: Results of symmetry operation testing:
--------------------------------------------------------------------------------
   Sym. Op.         Class     Trace     Result 
--------------------------------------------------------------------------------
 --------------B
    0.00000    1.44393    1.10565
    0.00000    0.00000   -0.13933
    0.00000   -1.44393    1.10565
       1               1     -1.00000   Passed
 --------------B
    0.00000    1.44393    1.10565
    0.00000    0.00000   -0.13933
    0.00000   -1.44393    1.10565
       2               2      1.00000   Passed
 --------------B
    0.00000   -1.44393    1.10565
    0.00000    0.00000   -0.13933
    0.00000    1.44393    1.10565
       3               3      1.00000   Passed
 --------------B
    0.00000   -1.44393    1.10565
    0.00000    0.00000   -0.13933
    0.00000    1.44393    1.10565
       4               4      3.00000   Passed
--------------------------------------------------------------------------------
 @TSTOPS:                      4  operations verified.
  @SYMUNQ-I, There are   2 orbits in the C2v  point group.
------------------------------------------------------------------------
    Set    Site Group                       Members of set
------------------------------------------------------------------------
      1      C s               1  3
      2      C2v               2
------------------------------------------------------------------------
         -----------------------------------------------
         Cartesian coordinates corresponding to internal 
                 coordinate input (Angstroms) 
 ----------------------------------------------------------------
 Z-matrix   Atomic            C o o r d i n a t e s
  Symbol    Number           X              Y              Z
 ----------------------------------------------------------------
     H         1         0.00000000    -0.76409272     0.58508316
     O         8         0.00000000     0.00000000    -0.07373114
     H         1         0.00000000     0.76409272     0.58508316
 ----------------------------------------------------------------
   Interatomic distance matrix (Angstroms) 
 
                 H             O             H    
                [ 1]        [ 2]        [ 3]
  H    [ 1]     0.00000
  O    [ 2]     1.00890     0.00000
  H    [ 3]     1.52819     1.00890     0.00000

   Interatomic angles (degrees) 


  O [ 2]-H [ 1]-H [ 3]   40.76849
  H [ 1]-O [ 2]-H [ 3]   98.46301

  H [ 1]-H [ 3]-O [ 2]   40.76849

      3 interatomic angles printed.
 Enetring Arhcive
 
 The COORD COMMON BLOCK AT ARCHIVE
   0.00000  -1.44393   1.10565
   0.00000   0.00000  -0.13933
   0.00000   1.44393   1.10565
 
   0.00000   0.00000   0.00000
   1.90654   0.00000   0.00000
   1.90654   1.71850   0.00000
 
 @-ARCHIVE, the number of opt. cycles                     1
   1.90654   1.90654   1.71850
 @ACES2: Executing "xvmol"
 One- and two-electron integrals over symmetry-adapted AOs are calculated.\n
 @READIN: Spherical harmonics are used.
  @READIN-I, Nuclear repulsion energy :    8.7384456738 a.u.
  required memory for a1 array                4451690  words 
  required memory for a2 array                   9156  words 
  @TWOEL-I,        228 integrals of symmetry type  I I I I
  @TWOEL-I, Total number of 2-e integrals        228.


ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.01 (      0.0)   0: 0: 0.01 (      0.0)   0: 0: 0.0 (      0.0)
 @ACES2: Executing "xvmol2ja"


ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.00 (      0.0)   0: 0: 0.00 (      0.0)   0: 0: 0.0 (      0.0)
 @ACES2: Executing "xvscf"
          Alpha population by irrep:   5
           Beta population by irrep:   5


 @SORTHO-I, Orthonormalizing initial guess. 
The Alpha and Beta occupation vector
                     5
 
                     5
 
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
           potradpts     integer                    94                    50
           OEPSLATER     integer                     5
             overind     integer                     3
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
               ksmem     integer             150000000               1500000
                func      string  BECKE,LYP             none                
               kspot      string  BECKE,LYP             hf                  
              cutoff     integer                    12
============================================================================
 
  --------------------------------------------------------------------
  Iteration         Total Energy            Largest Error in FDS-SDF
  --------------------------------------------------------------------
       0           -65.8952163610              0.2260327463D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       1           -71.6460140442              0.4377512891D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       2           -74.6353075115              0.3787298344D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       3           -74.5540500511              0.4439457797D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       4           -74.5966829013              0.3976171136D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       5           -74.5199010657              0.4568480831D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       6           -74.5750827283              0.4061528015D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       7           -74.5025192835              0.4634006330D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       8           -74.5643428316              0.4101883416D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       9           -74.9519790432              0.8705921587D-01
The Alpha and Beta occupation vector
                     5
 
                     5
      10           -74.9638322141              0.1643019720D-02
The Alpha and Beta occupation vector
                     5
 
                     5
      11           -74.9637807737              0.5503144011D-04
The Alpha and Beta occupation vector
                     5
 
                     5
      12           -74.9637825791              0.1590016272D-05
The Alpha and Beta occupation vector
                     5
 
                     5
      13           -74.9637825257              0.7179763173D-08
The Alpha and Beta occupation vector
                     5
 
                     5

@VSCF: SCF has converged.

     E(SCF)=       -74.9637825255

 The AOBASMOS? T
 Vib calc & AOBASMOS?                     0 T
 Creating The AOBASMOS file


      Trace of projected alpha density matrix =   1.000000000

      Alpha part of wavefunction is symmetric.




ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.60 (      0.6)   0: 0: 0.87 (      0.9)   0: 0: 1.5 (      1.5)
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
           potradpts     integer                    94                    50
           OEPSLATER     integer                     5
             overind     integer                     3
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
               ksmem     integer             150000000               1500000
                func      string  BECKE,LYP             none                
               kspot      string  BECKE,LYP             hf                  
              cutoff     integer                    12
============================================================================
 

The SCF nuclear-electron attraction energy =   -195.985299355076
The SCF kinetic energy                     =     74.497549780674
The SCF coulomb energy                     =     46.835330808590
The SCF exchange energy                    =     -9.049809433405
The SCF one electron energy                =   -121.487749574402
The SCF total energy                       =    -74.963782525462

Results using the SCF density

   Total density integrates to             :       9.999990344360 electrons
   Exch tot contrib : LDA (Slater, Xalpha) :      -8.158153138202
   Exch tot contrib : Becke                :      -9.038863250685
   Exch tot contrib : Perdew-Burke-Ernzerh :      -8.982590869380
   Exch tot contrib : Perdew-Wang 91       :      -9.013476954712
   Exch tot contrib : Exact Nonlocal Excha :      -9.049809433405
   Corr tot contrib : Vosko-Wilk-Nusair    :      -0.669766194804
   Corr tot contrib : Lee-Yang-Parr        :      -0.337054845152
   Corr tot contrib : Perdew-Burke-Ernzerh :      -0.340521794503
   Corr tot contrib : Perdew-Wang 91       :      -0.363220316956
   Corr tot contrib : Wilson-Levy          :      -0.399146418501
   Corr tot contrib : Wilson-Ivanov        :      -0.347919294702
   Hybr tot contrib : Becke III LYP        :      -9.370865679816
     
 Ehar from   KS   Calculation:         -65.913973092057 a.u.
 Ex   from   KS   Calculation:          -9.038863250685 a.u.
 Ec   from   KS   Calculation:          -0.337054845152 a.u.
 Exc  from   KS   Calculation:          -9.375918095836 a.u.
 Final Result from KS Calculation:     -75.289891187893 a.u.

 Ehar with the SCF Dens and Func:      -65.913973092057 a.u.
 Ex   with the SCF Dens and Func:       -9.038863250685 a.u.
 Ec   with the SCF Dens and Func:       -0.337054845152 a.u.
 Exc  with the SCF Dens and Func:       -9.375918095836 a.u.
 Final Result with the SCF Dens and Func:    -75.289891187893 a.u.

  @CHECKOUT-I, Total execution time :       0.3160 seconds.
ACES2: Total elapsed time is        6.1 seconds
 @ACES2: Executing "xvdint"
 One- and two-electron integral derivatives are calculated
 for RHF gradients and dipole moments.
 Spherical gaussians are used.

  Cartesian Coordinates
  ---------------------

  Total number of coordinates:  9


   1   H #1     x      0.0000000000
   2            y     -1.4439258767
   3            z      1.1056468517

   4   O #2     x      0.0000000000
   5            y      0.0000000000
   6            z     -0.1393316422

   7   H #3     x      0.0000000000
   8            y      1.4439258767
   9            z      1.1056468517

 Translational invariance is used.
 
 Entering the ECP_DGRAD

                    Nuclear attraction energy gradient
                    ----------------------------------

 H #1       0.0000000000           -3.9248594887            3.4352938973
 O #2       0.0000000000            0.0000000000           -6.8705877946
 H #3       0.0000000000            3.9248594887            3.4352938973



                       Two electron energy gradient
                       ----------------------------

 H #1       0.0000000000            1.8608521669           -1.7394046366
 O #2       0.0000000000            0.0000000000            3.4788092733
 H #3       0.0000000000           -1.8608521669           -1.7394046366



                    Nuclear repulsion energy gradient
                    ---------------------------------

 H #1       0.0000000000            1.7867618323           -1.4371904419
 O #2       0.0000000000            0.0000000000            2.8743808837
 H #3       0.0000000000           -1.7867618323           -1.4371904419



                      Kinetic energy energy gradient
                      ------------------------------

 H #1       0.0000000000            0.3790818275           -0.3312475708
 O #2       0.0000000000            0.0000000000            0.6624951417
 H #3       0.0000000000           -0.3790818275           -0.3312475708



                     Renormalization energy gradient
                     -------------------------------

 H #1       0.0000000000            0.0897078009           -0.1106020666
 O #2       0.0000000000            0.0000000000            0.2212041333
 H #3       0.0000000000           -0.0897078009           -0.1106020666



                  One elctron Molecular energy gradient
                  -------------------------------------

 H #1       0.0000000000           -1.6693080279            1.5562538179
 O #2       0.0000000000            0.0000000000           -3.1125076359
 H #3       0.0000000000            1.6693080279            1.5562538179



                        Molecular energy gradient
                        -------------------------

 H #1       0.0000000000            0.1915441390           -0.1831508187
 O #2       0.0000000000            0.0000000000            0.3663016374
 H #3       0.0000000000           -0.1915441390           -0.1831508187




ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.01 (      0.0)   0: 0: 0.00 (      0.0)   0: 0: 0.0 (      0.0)
 @ACES2: Executing "xvksdint"
 fffg
 
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
           potradpts     integer                    94                    50
           OEPSLATER     integer                     5
             overind     integer                     3
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
               ksmem     integer             150000000               1500000
                func      string  BECKE,LYP             none                
               kspot      string  BECKE,LYP             hf                  
              cutoff     integer                    12
============================================================================
 
  idd=                     1
 The new gradient:
 atm    x               y               z
  1  1    0.0000000000    0.1915441390   -0.1831508187
  8  2    0.0000000000    0.0000000000    0.3663016374
  1  3    0.0000000000   -0.1915441390   -0.1831508187
                     1                     1                     1
                     2                     8                     2
                     3                     1                     3
 t=                     1                     1
 -5.193418255117713E-017 -0.171125333660881       0.161002089327827     
 t=                     2                     2
  4.177809538570924E-016  5.601075159233915E-014 -0.322004178655545     
 t=                     3                     3
 -3.629905028023274E-016  0.171125333660824       0.161002089327715     
 rotin=  2.745869798517517E-015 -9.061154604417254E-014 -4.554108619621406E-015
 tran=  2.856268503587752E-018  2.856268503587752E-018 -3.691491556878645E-015
  1  1    0.0000000000    0.0204188054   -0.0221487294
  8  2    0.0000000000    0.0000000000    0.0442974587
  1  3    0.0000000000   -0.0204188054   -0.0221487294
  @CHECKOUT-I, Total execution time :       0.5749 seconds.


ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.43 (      0.4)   0: 0: 0.58 (      0.6)   0: 0: 1.0 (      1.0)
ACES2: Total elapsed time is        7.1 seconds
 @ACES2: Executing "xjoda"
 Finite diffs; ignore and geomopt vars: T T
 The vib calc. related varrs:
 iflags(h_IFLAGS_vib):                     0
 The finite diffs: F
 iflags2(h_IFLAGS2_geom_opt):                     1
 first run of popt num. frq:                     0
 The value of HAVEGEOM record:                     0
 The value of PES_SCAN record:                     0
 The OPTARC is here? T F
 Do we have geom? F
 Hessian in JOBARC?                    -1
Renaming SAVEDIR/CURRENT to SAVEDIR/OLD . . . done
Files copied to SAVEDIR/CURRENT: ZMAT.BAS JOBARC JAINDX OPTARC !OPTARCBK !DIPDER
Removing old back up directory . . . done
 A geometry for frequency {0,1} present:                     1
 The value of grad calc                     1
 The Deriv. level, havegeom?, and vib? :                     1
                     0                     1                     0
 The derivative level @exit                     1
 First order props and grad flags @enter:                     0
                     1
 First order props, deriv level and grad_calc @exit:                     0
                     1                     1
 @-Entry, the optimization cycle                     1
   JODA beginning optimization cycle #  2.
 -----After call to Entry-----
 The FNDDONE in Geopt:                     1
 The gradient calcs:                     1
 Hessian calc.,Geo. opt. and iarch: F T                     1
 Entering Fetchz IARCH .NE. 1                     1
 
 @GEOPT, PASS1,  PS1EXIST=:                     0 F
   Retrieving information from last optimization cycle.
 
 Start Reading Archive file NXM6:                      3
 The COORD COMMON BLOCK AT RETRIEVE BEFORE DECOMPRESS
   0.00000  -1.44393   1.10565
   0.00000   0.00000  -0.13933
   0.00000   1.44393   1.10565
 
   1.90654   1.90654   1.71850
 The COORD COMMON BLOCK AT RETRIEVE
   0.00000  -1.44393   1.10565
   0.00000   0.00000  -0.13933
   0.00000   1.44393   1.10565
 
   0.00000   0.00000   0.00000
 Data read from OPTARC file
 
 The COORD COMMON BLOCK/AFTER RETRIVE
   0.00000  -1.44393   1.10565
   0.00000   0.00000  -0.13933
   0.00000   1.44393   1.10565
 
   0.00000   0.00000   0.00000
 
 
 The Cartesians, internals and connec. @-buildb
 The NATOMS                     3
   0.00000  -1.44393   1.10565
   0.00000   0.00000  -0.13933
   0.00000   1.44393   1.10565
 
   0.00000   0.00000   0.00000
   1.90654   0.00000   0.00000
   1.90654   1.71850   0.00000
    0    0    0    1
    0    0    2    1
    0
  B MATRIX 
[  1,  1]  0.000000 [  1,  2] -0.757354 [  1,  3]  0.653004 [  1,  4]  0.000000
[  1,  5]  0.757354 [  1,  6] -0.653004 [  1,  7]  0.000000 [  1,  8]  0.000000
[  1,  9]  0.000000 [  2,  1]  0.000000 [  2,  2]  0.000000 [  2,  3]  0.000000
[  2,  4]  0.000000 [  2,  5] -0.757354 [  2,  6] -0.653004 [  2,  7]  0.000000
[  2,  8]  0.757354 [  2,  9]  0.653004 [  3,  1]  0.000000 [  3,  2] -0.342508
[  3,  3] -0.397240 [  3,  4]  0.000000 [  3,  5]  0.000000 [  3,  6]  0.794480
[  3,  7]  0.000000 [  3,  8]  0.342508 [  3,  9] -0.397240 [
  Bt MATRIX
[  1,  1]  0.000000 [  1,  2]  0.000000 [  1,  3]  0.000000 [  2,  1] -0.757354
[  2,  2]  0.000000 [  2,  3] -0.342508 [  3,  1]  0.653004 [  3,  2]  0.000000
[  3,  3] -0.397240 [  4,  1]  0.000000 [  4,  2]  0.000000 [  4,  3]  0.000000
[  5,  1]  0.757354 [  5,  2] -0.757354 [  5,  3]  0.000000 [  6,  1] -0.653004
[  6,  2] -0.653004 [  6,  3]  0.794480 [  7,  1]  0.000000 [  7,  2]  0.000000
[  7,  3]  0.000000 [  8,  1]  0.000000 [  8,  2]  0.757354 [  8,  3]  0.342508
[  9,  1]  0.000000 [  9,  2]  0.653004 [  9,  3] -0.397240 [
  G MATRIX 
[  1,  1]  2.000000 [  2,  1] -0.147171 [  2,  2]  2.000000 [  3,  1] -0.518799
[  3,  2] -0.518799 [  3,  3]  1.181422 [
  G MATRIX DETERMINANT   0.000000000000000E+000
  INVERSE G MATRIX 
[  1,  1]  0.590726 [  2,  1]  0.124997 [  2,  2]  0.590726 [  3,  1]  0.314297
[  3,  2]  0.314297 [  3,  3]  1.122473 [
  A MATRIX
[  1,  1]  0.000000 [  1,  2]  0.000000 [  1,  3]  0.000000 [  2,  1] -0.555038
[  2,  2] -0.202316 [  2,  3] -0.622489 [  3,  1]  0.260896 [  3,  2] -0.043228
[  3,  3] -0.240654 [  4,  1]  0.000000 [  4,  2]  0.000000 [  4,  3]  0.000000
[  5,  1]  0.352722 [  5,  2] -0.352722 [  5,  3]  0.000000 [  6,  1] -0.217668
[  6,  2] -0.217668 [  6,  3]  0.481309 [  7,  1]  0.000000 [  7,  2]  0.000000
[  7,  3]  0.000000 [  8,  1]  0.202316 [  8,  2]  0.555038 [  8,  3]  0.622489
[  9,  1] -0.043228 [  9,  2]  0.260896 [  9,  3] -0.240654 [
 
 The COORD COMMON BLOCK/at opt. start:
 
   1.90654   1.90654   1.71850
 
 The Cart. Hessian in after reading from READGH
 IAVHES flags: Cartesian exact Hess. read                     0
                     0
 If the Hessian is from VIB=EXACT, the trans/rot
 contaminants are included. If VIB=FINDIF it is
 projected.

                 COLUMN   1       COLUMN   2       COLUMN   3       COLUMN   4
 ROW   1      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   2      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   3      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   4      0.00000000000   -0.75735423941    0.65300425423    0.00000000000
 ROW   5      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   6      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   7      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   8      0.00000000000   -0.34250756792   -0.39724022764    0.00000000000
 ROW   9      0.00000000000    0.00000000000    0.00000000000    0.00000000000

                 COLUMN   5       COLUMN   6       COLUMN   7       COLUMN   8
 ROW   1      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   2      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   3      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   4      0.75735423941   -0.65300425423    0.00000000000    0.00000000000
 ROW   5      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   6      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   7     -0.75735423941   -0.65300425423    0.00000000000    0.75735423941
 ROW   8      0.00000000000    0.79448045527    0.00000000000    0.34250756792
 ROW   9      0.00000000000    0.00000000000    0.00000000000    0.00000000000

                 COLUMN   9
 ROW   1      0.00000000000
 ROW   2      0.00000000000
 ROW   3      0.00000000000
 ROW   4      0.00000000000
 ROW   5      0.00000000000
 ROW   6      0.00000000000
 ROW   7      0.65300425423
 ROW   8     -0.39724022764
 ROW   9      0.00000000000
 The Cartesian Hessian after reading from GETICFCM
 
 The COORD COMMON BLOCK/at opt. start:
 
   1.90654   1.90654   1.71850
 
 
 @-CONVQF,  The outgoing internal gradients
   -0.0299275   -0.0299275    0.0065601
 
   Internal coordinate forces and energy gradients (atomic units): 
           R       dV/dR             R       dV/dR             R       dV/dR
 
 [R    ]   1.90654  -0.02993 [R    ]   1.90654  -0.02993 [A    ]   1.71850   0.00656
 
   Hessian from cycle  1 read.
   BFGS update using last two gradients and previous step.
 The number of opt. cycles:                     1
 
   The eigenvectors of the Hessian matrix: 

                 COLUMN   1       COLUMN   2
 ROW   1      0.13381363996    0.99100651348
 ROW   2     -0.99100651348    0.13381363996
   The eigenvalues of the Hessian matrix: 
 
     0.24298    0.59485
 
 There are  0 Negative Eigenvalues.
 
 The number of degs. of freed.           at NR:    2
 
 The unmolested step
    0.0751   -0.0404
   MANR scale factor for NR step is 1.100 for R    .
 
 The initial trust rad:   0.30000
 The predicted and actual energy change:   -0.00828  -0.01274
 Ratio of change:   0.64976
 The new trust rad and TAU:   0.42426   0.09193
   Summary of Optimization Cycle: 
   The maximum unscaled step is:    0.09193.
   Scale factor set to:  1.00000.
   Forces are in hartree/bohr and hartree/radian.
   Parameter values are in Angstroms and degrees.
--------------------------------------------------------------------------
  Parameter     dV/dR           Step          Rold            Rnew
--------------------------------------------------------------------------
    R      -0.0299274833    0.0309048740    1.0088974002    1.0398022743
    A       0.0065601499   -2.3133923500   98.4630107822   96.1496184323
--------------------------------------------------------------------------
  Minimum force:  0.006560150 / RMS force:  0.021664370
 GMETRY starting with R
  1     1.964941
  2     1.964941
  3     1.678127
  4     0.000000
  5     0.000000
  6     0.000000
  7     0.000000
  8     0.000000
  9     0.000000
 In Gmetry; ncycle, ix, and ipost_vib
                     2                     0                     0
@GMETRY-I, Decompressing R.
 Updating structure...
 GMETRY using R vector
  1     0.000000
  2     0.000000
  3     0.000000
  4     1.964941
  5     0.000000
  6     0.000000
  7     1.964941
  8     1.678127
  9     0.000000
  @GMETRY-I, Cartesian coordinates before scaling:
      1     0.000000     0.000000     0.000000
      2     1.964941     0.000000     0.000000
      3     2.175436     1.953634     0.000000
  @GMETRY-I, Cartesian coordinates after scaling:
      1     0.000000     0.000000     0.000000
      2     1.964941     0.000000     0.000000
      3     2.175436     1.953634     0.000000
  @GMETRY-I, Cartesian coordinates from Z-matrix:
      1     0.000000     0.000000     0.000000
      2     0.000000     0.000000     1.964941
      3     1.953634     0.000000     2.175436
         -----------------------------------------------
         Cartesian coordinates corresponding to internal 
                 coordinate input (Angstroms) 
 ----------------------------------------------------------------
 Z-matrix   Atomic            C o o r d i n a t e s
  Symbol    Number           X              Y              Z
 ----------------------------------------------------------------
     H         1         0.00000000     0.00000000     0.00000000
     O         8         0.00000000     0.00000000     1.03980227
     H         1         1.03381878     0.00000000     1.15119123
 ----------------------------------------------------------------
    1.0078250000
   15.9949100000
    1.0078250000
The number of atoms                      3
 
 The Cartesians before translations to CM
   0.00000   0.00000   0.00000
   0.00000   0.00000   1.96494
   1.95363   0.00000   2.17544
 
 The Cartesians in center of mass coords
  -0.1093204   0.0000000  -1.8667670
  -0.1093204   0.0000000   0.0981744
   1.8443138   0.0000000   0.3086690
 
 After translation to center of mass coordinates 
       -0.109320387803    0.000000000000   -1.866767014077
       -0.109320387803    0.000000000000    0.098174366587
        1.844313828206    0.000000000000    0.308668973378
 
  @symmetry-i, Coordinates after  COM shift 
     -0.109320387803      0.000000000000     -1.866767014077
     -0.109320387803      0.000000000000      0.098174366587
      1.844313828206      0.000000000000      0.308668973378
 
  Inertia tensor
    3.76227    0.00000   -0.60774
    0.00000    7.39358    0.00000
   -0.60774    0.00000    3.63131
 Inertia tensor:
        3.762272109806    0.000000000000   -0.607744948356
        0.000000000000    7.393581081389    0.000000000000
       -0.607744948356    0.000000000000    3.631308971584
    Diagonalized inertia tensor:
        3.085528110845    0.000000000000    0.000000000000
        0.000000000000    4.308052970545    0.000000000000
        0.000000000000    0.000000000000    7.393581081389
    Eigenvectors of inertia tensor: 
        0.668159736506    0.744017853625    0.000000000000
        0.000000000000    0.000000000000    1.000000000000
        0.744017853625   -0.668159736506    0.000000000000
    Principal axis orientation for molecular system: 
       -1.461951468541    1.165962235953    0.000000000000
        0.000000000000   -0.146932479201    0.000000000000
        1.461951468541    1.165962235953    0.000000000000
   Rotational constants (in cm-1): 
     8.14193       13.97337       19.50979
  @SYMMETRY-I, Handedness of inertial frame: 1.00000
  @SYMMETRY-I, The symmetry group is 1-fold degenerate.
 
 The symmetry processing begins
 The Cartesians before entering symmetry auto
  -1.46195   1.16596   0.00000
   0.00000  -0.14693   0.00000
   1.46195   1.16596   0.00000
   0.00000   0.00000   0.00000
   0.00000   0.00000   0.00000
   0.00000   0.00000   0.00000
   Principal axis orientation for molecule:
       -1.461951468541    1.165962235953    0.000000000000
        0.000000000000   -0.146932479201    0.000000000000
        1.461951468541    1.165962235953    0.000000000000
   Reflection in plane  3 is a valid symmetry operation.
       -1.000000000000    0.000000000000    0.000000000000
        0.000000000000   -1.000000000000    0.000000000000
        0.000000000000    0.000000000000    1.000000000000
       -1.000000000000    0.000000000000    0.000000000000
        0.000000000000    1.000000000000    0.000000000000
        0.000000000000    0.000000000000   -1.000000000000
   Rotation about  2 is a valid symmetry operation 
   Reflection in plane  1 is a valid symmetry operation.
        1.000000000000    0.000000000000    0.000000000000
        0.000000000000   -1.000000000000    0.000000000000
        0.000000000000    0.000000000000   -1.000000000000
  Symmetry bits:    5   2   0
 @-SYM_AUTO The Orientation matrices: ORIEN2 and ORIENT
 

                 COLUMN   1       COLUMN   2       COLUMN   3
 ROW   1      0.66815973651    0.74401785362    0.00000000000
 ROW   2      0.00000000000    0.00000000000    1.00000000000
 ROW   3      0.74401785362   -0.66815973651    0.00000000000
 

                 COLUMN   1       COLUMN   2       COLUMN   3
 ROW   1      1.00000000000    0.00000000000    0.00000000000
 ROW   2      0.00000000000    1.00000000000    0.00000000000
 ROW   3      0.00000000000    0.00000000000    1.00000000000
 
 @-SYM_AUTO The variables in /COORD/ common block
   0.00000  -1.46195   1.16596
   0.00000   0.00000  -0.14693
   0.00000   1.46195   1.16596
 
 The NEWQ array 
  -1.46195
   1.16596
   0.00000
   0.00000
  -0.14693
   0.00000
   1.46195
   1.16596
   0.00000
 
 
 @-SYM_AUTO ORIENT2=ORIEN2xORIENT

                 COLUMN   1       COLUMN   2       COLUMN   3
 ROW   1      0.66815973651    0.74401785362    0.00000000000
 ROW   2      0.00000000000    0.00000000000    1.00000000000
 ROW   3      0.74401785362   -0.66815973651    0.00000000000
  @CHRTAB-I, Generated transformation matrices for   1 symmetry operations.
  @TSTOPS: Results of symmetry operation testing:
--------------------------------------------------------------------------------
   Sym. Op.         Class     Trace     Result 
--------------------------------------------------------------------------------
 --------------B
    0.00000   -1.46195    1.16596
    0.00000    0.00000   -0.14693
    0.00000    1.46195    1.16596
       1               1      3.00000   Passed
--------------------------------------------------------------------------------
 @TSTOPS:                      1  operations verified.
  @SYMUNQ-I, There are   3 orbits in the C1   point group.
------------------------------------------------------------------------
    Set    Site Group                       Members of set
------------------------------------------------------------------------
      1      C1                1
      2      C1                2
      3      C1                3
------------------------------------------------------------------------
  @CHRTAB-I, Generated transformation matrices for   4 symmetry operations.
  @TSTOPS: Results of symmetry operation testing:
--------------------------------------------------------------------------------
   Sym. Op.         Class     Trace     Result 
--------------------------------------------------------------------------------
 --------------B
    0.00000    1.46195    1.16596
    0.00000    0.00000   -0.14693
    0.00000   -1.46195    1.16596
       1               1     -1.00000   Passed
 --------------B
    0.00000    1.46195    1.16596
    0.00000    0.00000   -0.14693
    0.00000   -1.46195    1.16596
       2               2      1.00000   Passed
 --------------B
    0.00000   -1.46195    1.16596
    0.00000    0.00000   -0.14693
    0.00000    1.46195    1.16596
       3               3      1.00000   Passed
 --------------B
    0.00000   -1.46195    1.16596
    0.00000    0.00000   -0.14693
    0.00000    1.46195    1.16596
       4               4      3.00000   Passed
--------------------------------------------------------------------------------
 @TSTOPS:                      4  operations verified.
  @SYMUNQ-I, There are   2 orbits in the C2v  point group.
------------------------------------------------------------------------
    Set    Site Group                       Members of set
------------------------------------------------------------------------
      1      C s               1  3
      2      C2v               2
------------------------------------------------------------------------
         -----------------------------------------------
         Cartesian coordinates corresponding to internal 
                 coordinate input (Angstroms) 
 ----------------------------------------------------------------
 Z-matrix   Atomic            C o o r d i n a t e s
  Symbol    Number           X              Y              Z
 ----------------------------------------------------------------
     H         1         0.00000000    -0.77363146     0.61700069
     O         8         0.00000000     0.00000000    -0.07775333
     H         1         0.00000000     0.77363146     0.61700069
 ----------------------------------------------------------------
   Interatomic distance matrix (Angstroms) 
 
                 H             O             H    
                [ 1]        [ 2]        [ 3]
  H    [ 1]     0.00000
  O    [ 2]     1.03980     0.00000
  H    [ 3]     1.54726     1.03980     0.00000

   Interatomic angles (degrees) 


  O [ 2]-H [ 1]-H [ 3]   41.92519
  H [ 1]-O [ 2]-H [ 3]   96.14962

  H [ 1]-H [ 3]-O [ 2]   41.92519

      3 interatomic angles printed.
 Enetring Arhcive
 
 The COORD COMMON BLOCK AT ARCHIVE
   0.00000  -1.46195   1.16596
   0.00000   0.00000  -0.14693
   0.00000   1.46195   1.16596
 
   0.00000   0.00000   0.00000
   1.96494   0.00000   0.00000
   1.96494   1.67813   0.00000
 
 @-ARCHIVE, the number of opt. cycles                     2
   1.96494   1.96494   1.67813
 @ACES2: Executing "xvmol"
 One- and two-electron integrals over symmetry-adapted AOs are calculated.\n
 @READIN: Spherical harmonics are used.
  @READIN-I, Nuclear repulsion energy :    8.4847451676 a.u.
  required memory for a1 array                4451690  words 
  required memory for a2 array                   9156  words 
  @TWOEL-I,        228 integrals of symmetry type  I I I I
  @TWOEL-I, Total number of 2-e integrals        228.


ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.01 (      0.0)   0: 0: 0.01 (      0.0)   0: 0: 0.0 (      0.0)
 @ACES2: Executing "xvmol2ja"


ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.00 (      0.0)   0: 0: 0.00 (      0.0)   0: 0: 0.0 (      0.0)
 @ACES2: Executing "xvscf"
          Alpha population by irrep:   5
           Beta population by irrep:   5


 @SORTHO-I, Orthonormalizing initial guess. 
The Alpha and Beta occupation vector
                     5
 
                     5
 
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
           potradpts     integer                    94                    50
           OEPSLATER     integer                     5
             overind     integer                     3
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
               ksmem     integer             150000000               1500000
                func      string  BECKE,LYP             none                
               kspot      string  BECKE,LYP             hf                  
              cutoff     integer                    12
============================================================================
 
  --------------------------------------------------------------------
  Iteration         Total Energy            Largest Error in FDS-SDF
  --------------------------------------------------------------------
       0           -65.9299673043              0.2195237908D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       1           -71.6595679540              0.4374820810D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       2           -74.5888864417              0.3977489518D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       3           -74.4443250017              0.4804121976D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       4           -74.4733102776              0.4382724725D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       5           -74.3510352160              0.5034961145D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       6           -74.4207105550              0.4519074378D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       7           -74.3159551574              0.5107551234D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       8           -74.4024512861              0.4561396585D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       9           -74.9548041611              0.5745401073D-01
The Alpha and Beta occupation vector
                     5
 
                     5
      10           -74.9590836296              0.1618550573D-02
The Alpha and Beta occupation vector
                     5
 
                     5
      11           -74.9590308452              0.1515935130D-03
The Alpha and Beta occupation vector
                     5
 
                     5
      12           -74.9590360910              0.1819764081D-05
The Alpha and Beta occupation vector
                     5
 
                     5
      13           -74.9590360240              0.8678663049D-07
The Alpha and Beta occupation vector
                     5
 
                     5

@VSCF: SCF has converged.

     E(SCF)=       -74.9590360271

 The AOBASMOS? T
 Vib calc & AOBASMOS?                     0 T
 Creating The AOBASMOS file


      Trace of projected alpha density matrix =   1.000000000

      Alpha part of wavefunction is symmetric.




ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.56 (      0.6)   0: 0: 0.91 (      0.9)   0: 0: 1.5 (      1.5)
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
           potradpts     integer                    94                    50
           OEPSLATER     integer                     5
             overind     integer                     3
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
               ksmem     integer             150000000               1500000
                func      string  BECKE,LYP             none                
               kspot      string  BECKE,LYP             hf                  
              cutoff     integer                    12
============================================================================
 

The SCF nuclear-electron attraction energy =   -195.463170538491
The SCF kinetic energy                     =     74.449083496522
The SCF coulomb energy                     =     46.591768443248
The SCF exchange energy                    =     -9.021462596002
The SCF one electron energy                =   -121.014087041969
The SCF total energy                       =    -74.959036027101

Results using the SCF density

   Total density integrates to             :       9.999987075108 electrons
   Exch tot contrib : LDA (Slater, Xalpha) :      -8.136475542733
   Exch tot contrib : Becke                :      -9.018498119852
   Exch tot contrib : Perdew-Burke-Ernzerh :      -8.962158224433
   Exch tot contrib : Perdew-Wang 91       :      -8.993019044347
   Exch tot contrib : Exact Nonlocal Excha :      -9.021462596002
   Corr tot contrib : Vosko-Wilk-Nusair    :      -0.668563113285
   Corr tot contrib : Lee-Yang-Parr        :      -0.335847760867
   Corr tot contrib : Perdew-Burke-Ernzerh :      -0.338783651250
   Corr tot contrib : Perdew-Wang 91       :      -0.361495645767
   Corr tot contrib : Wilson-Levy          :      -0.396282376464
   Corr tot contrib : Wilson-Ivanov        :      -0.346396069659
   Hybr tot contrib : Becke III LYP        :      -9.347592886739
     
 Ehar from   KS   Calculation:         -65.937573431100 a.u.
 Ex   from   KS   Calculation:          -9.018498119852 a.u.
 Ec   from   KS   Calculation:          -0.335847760867 a.u.
 Exc  from   KS   Calculation:          -9.354345880719 a.u.
 Final Result from KS Calculation:     -75.291919311819 a.u.

 Ehar with the SCF Dens and Func:      -65.937573431100 a.u.
 Ex   with the SCF Dens and Func:       -9.018498119852 a.u.
 Ec   with the SCF Dens and Func:       -0.335847760867 a.u.
 Exc  with the SCF Dens and Func:       -9.354345880719 a.u.
 Final Result with the SCF Dens and Func:    -75.291919311819 a.u.

  @CHECKOUT-I, Total execution time :       0.3399 seconds.
ACES2: Total elapsed time is        9.5 seconds
 @ACES2: Executing "xvdint"
 One- and two-electron integral derivatives are calculated
 for RHF gradients and dipole moments.
 Spherical gaussians are used.

  Cartesian Coordinates
  ---------------------

  Total number of coordinates:  9


   1   H #1     x      0.0000000000
   2            y     -1.4619514685
   3            z      1.1659622360

   4   O #2     x      0.0000000000
   5            y      0.0000000000
   6            z     -0.1469324792

   7   H #3     x      0.0000000000
   8            y      1.4619514685
   9            z      1.1659622360

 Translational invariance is used.
 
 Entering the ECP_DGRAD

                    Nuclear attraction energy gradient
                    ----------------------------------

 H #1       0.0000000000           -3.7611636001            3.4543666823
 O #2       0.0000000000            0.0000000000           -6.9087333645
 H #3       0.0000000000            3.7611636001            3.4543666823



                       Two electron energy gradient
                       ----------------------------

 H #1       0.0000000000            1.8293016389           -1.8019704916
 O #2       0.0000000000            0.0000000000            3.6039409833
 H #3       0.0000000000           -1.8293016389           -1.8019704916



                    Nuclear repulsion energy gradient
                    ---------------------------------

 H #1       0.0000000000            1.6585785728           -1.3844302842
 O #2       0.0000000000            0.0000000000            2.7688605683
 H #3       0.0000000000           -1.6585785728           -1.3844302842



                      Kinetic energy energy gradient
                      ------------------------------

 H #1       0.0000000000            0.3528634552           -0.3252356995
 O #2       0.0000000000            0.0000000000            0.6504713991
 H #3       0.0000000000           -0.3528634552           -0.3252356995



                     Renormalization energy gradient
                     -------------------------------

 H #1       0.0000000000            0.0887357349           -0.1097926705
 O #2       0.0000000000            0.0000000000            0.2195853409
 H #3       0.0000000000           -0.0887357349           -0.1097926705



                  One elctron Molecular energy gradient
                  -------------------------------------

 H #1       0.0000000000           -1.6609858371            1.6349080281
 O #2       0.0000000000            0.0000000000           -3.2698160562
 H #3       0.0000000000            1.6609858371            1.6349080281



                        Molecular energy gradient
                        -------------------------

 H #1       0.0000000000            0.1683158018           -0.1670624635
 O #2       0.0000000000            0.0000000000            0.3341249271
 H #3       0.0000000000           -0.1683158018           -0.1670624635




ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.01 (      0.0)   0: 0: 0.00 (      0.0)   0: 0: 0.0 (      0.0)
 @ACES2: Executing "xvksdint"
 fffg
 
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
           potradpts     integer                    94                    50
           OEPSLATER     integer                     5
             overind     integer                     3
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
               ksmem     integer             150000000               1500000
                func      string  BECKE,LYP             none                
               kspot      string  BECKE,LYP             hf                  
              cutoff     integer                    12
============================================================================
 
  idd=                     1
 The new gradient:
 atm    x               y               z
  1  1    0.0000000000    0.1683158018   -0.1670624635
  8  2    0.0000000000    0.0000000000    0.3341249271
  1  3    0.0000000000   -0.1683158018   -0.1670624635
                     1                     1                     1
                     2                     8                     2
                     3                     1                     3
 t=                     1                     1
  4.191133548961389E-018 -0.165401121477254       0.164907856208818     
 t=                     2                     2
  1.011411200040186E-016 -2.842170943040401E-014 -0.329815712417653     
 t=                     3                     3
 -1.043363366463110E-016  0.165401121477283       0.164907856208843     
 rotin=  1.380766693430448E-015  3.412363707572252E-015  4.552679272457821E-016
 tran=  9.959169066690301E-019  9.959169066690301E-019  9.048317650695026E-015
  1  1    0.0000000000    0.0029146804   -0.0021546073
  8  2    0.0000000000    0.0000000000    0.0043092146
  1  3    0.0000000000   -0.0029146804   -0.0021546073
  @CHECKOUT-I, Total execution time :       0.5769 seconds.


ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.45 (      0.5)   0: 0: 0.58 (      0.6)   0: 0: 1.0 (      1.0)
ACES2: Total elapsed time is       10.6 seconds
 @ACES2: Executing "xjoda"
 Finite diffs; ignore and geomopt vars: T T
 The vib calc. related varrs:
 iflags(h_IFLAGS_vib):                     0
 The finite diffs: F
 iflags2(h_IFLAGS2_geom_opt):                     1
 first run of popt num. frq:                     0
 The value of HAVEGEOM record:                     0
 The value of PES_SCAN record:                     0
 The OPTARC is here? T F
 Do we have geom? F
 Hessian in JOBARC?                    -1
Renaming SAVEDIR/CURRENT to SAVEDIR/OLD . . . done
Files copied to SAVEDIR/CURRENT: ZMAT.BAS JOBARC JAINDX OPTARC !OPTARCBK !DIPDER
Removing old back up directory . . . done
 A geometry for frequency {0,1} present:                     1
 The value of grad calc                     1
 The Deriv. level, havegeom?, and vib? :                     1
                     0                     1                     0
 The derivative level @exit                     1
 First order props and grad flags @enter:                     0
                     1
 First order props, deriv level and grad_calc @exit:                     0
                     1                     1
 @-Entry, the optimization cycle                     2
   JODA beginning optimization cycle #  3.
 -----After call to Entry-----
 The FNDDONE in Geopt:                     1
 The gradient calcs:                     1
 Hessian calc.,Geo. opt. and iarch: F T                     1
 Entering Fetchz IARCH .NE. 1                     1
 
 @GEOPT, PASS1,  PS1EXIST=:                     0 F
   Retrieving information from last optimization cycle.
 
 Start Reading Archive file NXM6:                      3
 The COORD COMMON BLOCK AT RETRIEVE BEFORE DECOMPRESS
   0.00000  -1.46195   1.16596
   0.00000   0.00000  -0.14693
   0.00000   1.46195   1.16596
 
   1.96494   1.96494   1.67813
 The COORD COMMON BLOCK AT RETRIEVE
   0.00000  -1.46195   1.16596
   0.00000   0.00000  -0.14693
   0.00000   1.46195   1.16596
 
   0.00000   0.00000   0.00000
 Data read from OPTARC file
 
 The COORD COMMON BLOCK/AFTER RETRIVE
   0.00000  -1.46195   1.16596
   0.00000   0.00000  -0.14693
   0.00000   1.46195   1.16596
 
   0.00000   0.00000   0.00000
 
 
 The Cartesians, internals and connec. @-buildb
 The NATOMS                     3
   0.00000  -1.46195   1.16596
   0.00000   0.00000  -0.14693
   0.00000   1.46195   1.16596
 
   0.00000   0.00000   0.00000
   1.96494   0.00000   0.00000
   1.96494   1.67813   0.00000
    0    0    0    1
    0    0    2    1
    0
  B MATRIX 
[  1,  1]  0.000000 [  1,  2] -0.744018 [  1,  3]  0.668160 [  1,  4]  0.000000
[  1,  5]  0.744018 [  1,  6] -0.668160 [  1,  7]  0.000000 [  1,  8]  0.000000
[  1,  9]  0.000000 [  2,  1]  0.000000 [  2,  2]  0.000000 [  2,  3]  0.000000
[  2,  4]  0.000000 [  2,  5] -0.744018 [  2,  6] -0.668160 [  2,  7]  0.000000
[  2,  8]  0.744018 [  2,  9]  0.668160 [  3,  1]  0.000000 [  3,  2] -0.340041
[  3,  3] -0.378646 [  3,  4]  0.000000 [  3,  5]  0.000000 [  3,  6]  0.757293
[  3,  7]  0.000000 [  3,  8]  0.340041 [  3,  9] -0.378646 [
  Bt MATRIX
[  1,  1]  0.000000 [  1,  2]  0.000000 [  1,  3]  0.000000 [  2,  1] -0.744018
[  2,  2]  0.000000 [  2,  3] -0.340041 [  3,  1]  0.668160 [  3,  2]  0.000000
[  3,  3] -0.378646 [  4,  1]  0.000000 [  4,  2]  0.000000 [  4,  3]  0.000000
[  5,  1]  0.744018 [  5,  2] -0.744018 [  5,  3]  0.000000 [  6,  1] -0.668160
[  6,  2] -0.668160 [  6,  3]  0.757293 [  7,  1]  0.000000 [  7,  2]  0.000000
[  7,  3]  0.000000 [  8,  1]  0.000000 [  8,  2]  0.744018 [  8,  3]  0.340041
[  9,  1]  0.000000 [  9,  2]  0.668160 [  9,  3] -0.378646 [
  G MATRIX 
[  1,  1]  2.000000 [  2,  1] -0.107125 [  2,  2]  2.000000 [  3,  1] -0.505992
[  3,  2] -0.505992 [  3,  3]  1.091493 [
  G MATRIX DETERMINANT   0.000000000000000E+000
  INVERSE G MATRIX 
[  1,  1]  0.588478 [  2,  1]  0.113897 [  2,  2]  0.588478 [  3,  1]  0.325606
[  3,  2]  0.325606 [  3,  3]  1.218063 [
  A MATRIX
[  1,  1]  0.000000 [  1,  2]  0.000000 [  1,  3]  0.000000 [  2,  1] -0.548557
[  2,  2] -0.195461 [  2,  3] -0.656447 [  3,  1]  0.269908 [  3,  2] -0.047188
[  3,  3] -0.243659 [  4,  1]  0.000000 [  4,  2]  0.000000 [  4,  3]  0.000000
[  5,  1]  0.353096 [  5,  2] -0.353096 [  5,  3]  0.000000 [  6,  1] -0.222720
[  6,  2] -0.222720 [  6,  3]  0.487317 [  7,  1]  0.000000 [  7,  2]  0.000000
[  7,  3]  0.000000 [  8,  1]  0.195461 [  8,  2]  0.548557 [  8,  3]  0.656447
[  9,  1] -0.047188 [  9,  2]  0.269908 [  9,  3] -0.243659 [
 
 The COORD COMMON BLOCK/at opt. start:
 
   1.96494   1.96494   1.67813
 
 The Cart. Hessian in after reading from READGH
 IAVHES flags: Cartesian exact Hess. read                     0
                     0
 If the Hessian is from VIB=EXACT, the trans/rot
 contaminants are included. If VIB=FINDIF it is
 projected.

                 COLUMN   1       COLUMN   2       COLUMN   3       COLUMN   4
 ROW   1      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   2      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   3      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   4      0.00000000000   -0.74401785362    0.66815973651    0.00000000000
 ROW   5      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   6      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   7      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   8      0.00000000000   -0.34004054425   -0.37864633569    0.00000000000
 ROW   9      0.00000000000    0.00000000000    0.00000000000    0.00000000000

                 COLUMN   5       COLUMN   6       COLUMN   7       COLUMN   8
 ROW   1      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   2      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   3      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   4      0.74401785362   -0.66815973651    0.00000000000    0.00000000000
 ROW   5      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   6      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   7     -0.74401785362   -0.66815973651    0.00000000000    0.74401785362
 ROW   8      0.00000000000    0.75729267137    0.00000000000    0.34004054425
 ROW   9      0.00000000000    0.00000000000    0.00000000000    0.00000000000

                 COLUMN   9
 ROW   1      0.00000000000
 ROW   2      0.00000000000
 ROW   3      0.00000000000
 ROW   4      0.00000000000
 ROW   5      0.00000000000
 ROW   6      0.00000000000
 ROW   7      0.66815973651
 ROW   8     -0.37864633569
 ROW   9      0.00000000000
 The Cartesian Hessian after reading from GETICFCM
 
 The COORD COMMON BLOCK/at opt. start:
 
   1.96494   1.96494   1.67813
 
 
 @-CONVQF,  The outgoing internal gradients
   -0.0036082   -0.0036082   -0.0006767
 
   Internal coordinate forces and energy gradients (atomic units): 
           R       dV/dR             R       dV/dR             R       dV/dR
 
 [R    ]   1.96494  -0.00361 [R    ]   1.96494  -0.00361 [A    ]   1.67813  -0.00068
 
   Hessian from cycle  2 read.
   BFGS update using last two gradients and previous step.
 The number of opt. cycles:                     2
 
   The eigenvectors of the Hessian matrix: 

                 COLUMN   1       COLUMN   2
 ROW   1      0.16717465372    0.98592729709
 ROW   2     -0.98592729709    0.16717465372
   The eigenvalues of the Hessian matrix: 
 
     0.24913    0.47519
 
 There are  0 Negative Eigenvalues.
 
 The number of degs. of freed.           at NR:    2
 
 The unmolested step
    0.0108    0.0011
   MANR scale factor for NR step is 1.013 for R    .
 
 The initial trust rad:   0.42426
 The predicted and actual energy change:   -0.00164  -0.00203
 Ratio of change:   0.80853
 The new trust rad and TAU:   0.42426   0.01099
   Summary of Optimization Cycle: 
   The maximum unscaled step is:    0.01099.
   Scale factor set to:  1.00000.
   Forces are in hartree/bohr and hartree/radian.
   Parameter values are in Angstroms and degrees.
--------------------------------------------------------------------------
  Parameter     dV/dR           Step          Rold            Rnew
--------------------------------------------------------------------------
    R      -0.0036081961    0.0040938790    1.0398022743    1.0438961532
    A      -0.0006767371    0.0615513594   96.1496184323   96.2111697916
--------------------------------------------------------------------------
  Minimum force:  0.000676737 / RMS force:  0.002595867
 GMETRY starting with R
  1     1.972678
  2     1.972678
  3     1.679202
  4     0.000000
  5     0.000000
  6     0.000000
  7     0.000000
  8     0.000000
  9     0.000000
 In Gmetry; ncycle, ix, and ipost_vib
                     3                     0                     0
@GMETRY-I, Decompressing R.
 Updating structure...
 GMETRY using R vector
  1     0.000000
  2     0.000000
  3     0.000000
  4     1.972678
  5     0.000000
  6     0.000000
  7     1.972678
  8     1.679202
  9     0.000000
  @GMETRY-I, Cartesian coordinates before scaling:
      1     0.000000     0.000000     0.000000
      2     1.972678     0.000000     0.000000
      3     2.186108     1.961098     0.000000
  @GMETRY-I, Cartesian coordinates after scaling:
      1     0.000000     0.000000     0.000000
      2     1.972678     0.000000     0.000000
      3     2.186108     1.961098     0.000000
  @GMETRY-I, Cartesian coordinates from Z-matrix:
      1     0.000000     0.000000     0.000000
      2     0.000000     0.000000     1.972678
      3     1.961098     0.000000     2.186108
         -----------------------------------------------
         Cartesian coordinates corresponding to internal 
                 coordinate input (Angstroms) 
 ----------------------------------------------------------------
 Z-matrix   Atomic            C o o r d i n a t e s
  Symbol    Number           X              Y              Z
 ----------------------------------------------------------------
     H         1         0.00000000     0.00000000     0.00000000
     O         8         0.00000000     0.00000000     1.04389615
     H         1         1.03776837     0.00000000     1.15683858
 ----------------------------------------------------------------
    1.0078250000
   15.9949100000
    1.0078250000
The number of atoms                      3
 
 The Cartesians before translations to CM
   0.00000   0.00000   0.00000
   0.00000   0.00000   1.97268
   1.96110   0.00000   2.18611
 
 The Cartesians in center of mass coords
  -0.1097380   0.0000000  -1.8742347
  -0.1097380   0.0000000   0.0984430
   1.8513598   0.0000000   0.3118732
 
 After translation to center of mass coordinates 
       -0.109738034080    0.000000000000   -1.874234689875
       -0.109738034080    0.000000000000    0.098443000245
        1.851359822277    0.000000000000    0.311873239170
 
  @symmetry-i, Coordinates after  COM shift 
     -0.109738034080      0.000000000000     -1.874234689875
     -0.109738034080      0.000000000000      0.098443000245
      1.851359822277      0.000000000000      0.311873239170
 
  Inertia tensor
    3.79328    0.00000   -0.61640
    0.00000    7.45238    0.00000
   -0.61640    0.00000    3.65911
 Inertia tensor:
        3.793276062600    0.000000000000   -0.616399819878
        0.000000000000    7.452384051534    0.000000000000
       -0.616399819878    0.000000000000    3.659107988934
    Diagonalized inertia tensor:
        3.106152506105    0.000000000000    0.000000000000
        0.000000000000    4.346231545429    0.000000000000
        0.000000000000    0.000000000000    7.452384051534
    Eigenvectors of inertia tensor: 
        0.667760000664    0.744376639553    0.000000000000
        0.000000000000    0.000000000000    1.000000000000
        0.744376639553   -0.667760000664    0.000000000000
    Principal axis orientation for molecular system: 
       -1.468415189894    1.169852528716    0.000000000000
        0.000000000000   -0.147422726949    0.000000000000
        1.468415189894    1.169852528716    0.000000000000
   Rotational constants (in cm-1): 
     8.07768       13.85062       19.38025
  @SYMMETRY-I, Handedness of inertial frame: 1.00000
  @SYMMETRY-I, The symmetry group is 1-fold degenerate.
 
 The symmetry processing begins
 The Cartesians before entering symmetry auto
  -1.46842   1.16985   0.00000
   0.00000  -0.14742   0.00000
   1.46842   1.16985   0.00000
   0.00000   0.00000   0.00000
   0.00000   0.00000   0.00000
   0.00000   0.00000   0.00000
   Principal axis orientation for molecule:
       -1.468415189894    1.169852528716    0.000000000000
        0.000000000000   -0.147422726949    0.000000000000
        1.468415189894    1.169852528716    0.000000000000
   Reflection in plane  3 is a valid symmetry operation.
       -1.000000000000    0.000000000000    0.000000000000
        0.000000000000   -1.000000000000    0.000000000000
        0.000000000000    0.000000000000    1.000000000000
       -1.000000000000    0.000000000000    0.000000000000
        0.000000000000    1.000000000000    0.000000000000
        0.000000000000    0.000000000000   -1.000000000000
   Rotation about  2 is a valid symmetry operation 
   Reflection in plane  1 is a valid symmetry operation.
        1.000000000000    0.000000000000    0.000000000000
        0.000000000000   -1.000000000000    0.000000000000
        0.000000000000    0.000000000000   -1.000000000000
  Symmetry bits:    5   2   0
 @-SYM_AUTO The Orientation matrices: ORIEN2 and ORIENT
 

                 COLUMN   1       COLUMN   2       COLUMN   3
 ROW   1      0.66776000066    0.74437663955    0.00000000000
 ROW   2      0.00000000000    0.00000000000    1.00000000000
 ROW   3      0.74437663955   -0.66776000066    0.00000000000
 

                 COLUMN   1       COLUMN   2       COLUMN   3
 ROW   1      1.00000000000    0.00000000000    0.00000000000
 ROW   2      0.00000000000    1.00000000000    0.00000000000
 ROW   3      0.00000000000    0.00000000000    1.00000000000
 
 @-SYM_AUTO The variables in /COORD/ common block
   0.00000  -1.46842   1.16985
   0.00000   0.00000  -0.14742
   0.00000   1.46842   1.16985
 
 The NEWQ array 
  -1.46842
   1.16985
   0.00000
   0.00000
  -0.14742
   0.00000
   1.46842
   1.16985
   0.00000
 
 
 @-SYM_AUTO ORIENT2=ORIEN2xORIENT

                 COLUMN   1       COLUMN   2       COLUMN   3
 ROW   1      0.66776000066    0.74437663955    0.00000000000
 ROW   2      0.00000000000    0.00000000000    1.00000000000
 ROW   3      0.74437663955   -0.66776000066    0.00000000000
  @CHRTAB-I, Generated transformation matrices for   1 symmetry operations.
  @TSTOPS: Results of symmetry operation testing:
--------------------------------------------------------------------------------
   Sym. Op.         Class     Trace     Result 
--------------------------------------------------------------------------------
 --------------B
    0.00000   -1.46842    1.16985
    0.00000    0.00000   -0.14742
    0.00000    1.46842    1.16985
       1               1      3.00000   Passed
--------------------------------------------------------------------------------
 @TSTOPS:                      1  operations verified.
  @SYMUNQ-I, There are   3 orbits in the C1   point group.
------------------------------------------------------------------------
    Set    Site Group                       Members of set
------------------------------------------------------------------------
      1      C1                1
      2      C1                2
      3      C1                3
------------------------------------------------------------------------
  @CHRTAB-I, Generated transformation matrices for   4 symmetry operations.
  @TSTOPS: Results of symmetry operation testing:
--------------------------------------------------------------------------------
   Sym. Op.         Class     Trace     Result 
--------------------------------------------------------------------------------
 --------------B
    0.00000    1.46842    1.16985
    0.00000    0.00000   -0.14742
    0.00000   -1.46842    1.16985
       1               1     -1.00000   Passed
 --------------B
    0.00000    1.46842    1.16985
    0.00000    0.00000   -0.14742
    0.00000   -1.46842    1.16985
       2               2      1.00000   Passed
 --------------B
    0.00000   -1.46842    1.16985
    0.00000    0.00000   -0.14742
    0.00000    1.46842    1.16985
       3               3      1.00000   Passed
 --------------B
    0.00000   -1.46842    1.16985
    0.00000    0.00000   -0.14742
    0.00000    1.46842    1.16985
       4               4      3.00000   Passed
--------------------------------------------------------------------------------
 @TSTOPS:                      4  operations verified.
  @SYMUNQ-I, There are   2 orbits in the C2v  point group.
------------------------------------------------------------------------
    Set    Site Group                       Members of set
------------------------------------------------------------------------
      1      C s               1  3
      2      C2v               2
------------------------------------------------------------------------
         -----------------------------------------------
         Cartesian coordinates corresponding to internal 
                 coordinate input (Angstroms) 
 ----------------------------------------------------------------
 Z-matrix   Atomic            C o o r d i n a t e s
  Symbol    Number           X              Y              Z
 ----------------------------------------------------------------
     H         1         0.00000000    -0.77705191     0.61905934
     O         8         0.00000000     0.00000000    -0.07801275
     H         1         0.00000000     0.77705191     0.61905934
 ----------------------------------------------------------------
   Interatomic distance matrix (Angstroms) 
 
                 H             O             H    
                [ 1]        [ 2]        [ 3]
  H    [ 1]     0.00000
  O    [ 2]     1.04390     0.00000
  H    [ 3]     1.55410     1.04390     0.00000

   Interatomic angles (degrees) 


  O [ 2]-H [ 1]-H [ 3]   41.89442
  H [ 1]-O [ 2]-H [ 3]   96.21117

  H [ 1]-H [ 3]-O [ 2]   41.89442

      3 interatomic angles printed.
 Enetring Arhcive
 
 The COORD COMMON BLOCK AT ARCHIVE
   0.00000  -1.46842   1.16985
   0.00000   0.00000  -0.14742
   0.00000   1.46842   1.16985
 
   0.00000   0.00000   0.00000
   1.97268   0.00000   0.00000
   1.97268   1.67920   0.00000
 
 @-ARCHIVE, the number of opt. cycles                     3
   1.97268   1.97268   1.67920
 @ACES2: Executing "xvmol"
 One- and two-electron integrals over symmetry-adapted AOs are calculated.\n
 @READIN: Spherical harmonics are used.
  @READIN-I, Nuclear repulsion energy :    8.4513060872 a.u.
  required memory for a1 array                4451690  words 
  required memory for a2 array                   9156  words 
  @TWOEL-I,        228 integrals of symmetry type  I I I I
  @TWOEL-I, Total number of 2-e integrals        228.


ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.01 (      0.0)   0: 0: 0.01 (      0.0)   0: 0: 0.0 (      0.0)
 @ACES2: Executing "xvmol2ja"


ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.00 (      0.0)   0: 0: 0.00 (      0.0)   0: 0: 0.0 (      0.0)
 @ACES2: Executing "xvscf"
          Alpha population by irrep:   5
           Beta population by irrep:   5


 @SORTHO-I, Orthonormalizing initial guess. 
The Alpha and Beta occupation vector
                     5
 
                     5
 
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
           potradpts     integer                    94                    50
           OEPSLATER     integer                     5
             overind     integer                     3
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
               ksmem     integer             150000000               1500000
                func      string  BECKE,LYP             none                
               kspot      string  BECKE,LYP             hf                  
              cutoff     integer                    12
============================================================================
 
  --------------------------------------------------------------------
  Iteration         Total Energy            Largest Error in FDS-SDF
  --------------------------------------------------------------------
       0           -65.9395750879              0.2155539976D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       1           -71.6807309739              0.4268172348D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       2           -74.6006377106              0.3931781642D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       3           -74.4467353914              0.4792179045D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       4           -74.4681049006              0.4400847871D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       5           -74.3381655513              0.5055197361D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       6           -74.4069344614              0.4554866086D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       7           -74.2979074275              0.5133014599D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       8           -74.3862998472              0.4600435479D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       9           -74.9551475797              0.5094150522D-01
The Alpha and Beta occupation vector
                     5
 
                     5
      10           -74.9582762059              0.1560439410D-02
The Alpha and Beta occupation vector
                     5
 
                     5
      11           -74.9582250860              0.1570383067D-03
The Alpha and Beta occupation vector
                     5
 
                     5
      12           -74.9582305587              0.1843212589D-05
The Alpha and Beta occupation vector
                     5
 
                     5
      13           -74.9582304899              0.1000392206D-06
The Alpha and Beta occupation vector
                     5
 
                     5
      14           -74.9582304935              0.1293705865D-08
The Alpha and Beta occupation vector
                     5
 
                     5

@VSCF: SCF has converged.

     E(SCF)=       -74.9582304935

 The AOBASMOS? T
 Vib calc & AOBASMOS?                     0 T
 Creating The AOBASMOS file


      Trace of projected alpha density matrix =   1.000000000

      Alpha part of wavefunction is symmetric.




ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.66 (      0.7)   0: 0: 0.97 (      1.0)   0: 0: 1.7 (      1.7)
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
           potradpts     integer                    94                    50
           OEPSLATER     integer                     5
             overind     integer                     3
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
               ksmem     integer             150000000               1500000
                func      string  BECKE,LYP             none                
               kspot      string  BECKE,LYP             hf                  
              cutoff     integer                    12
============================================================================
 

The SCF nuclear-electron attraction energy =   -195.395530796005
The SCF kinetic energy                     =     74.442549083819
The SCF coulomb energy                     =     46.561383460498
The SCF exchange energy                    =     -9.017938328962
The SCF one electron energy                =   -120.952981712186
The SCF total energy                       =    -74.958230493472

Results using the SCF density

   Total density integrates to             :       9.999986619639 electrons
   Exch tot contrib : LDA (Slater, Xalpha) :      -8.133769994078
   Exch tot contrib : Becke                :      -9.015964467082
   Exch tot contrib : Perdew-Burke-Ernzerh :      -8.959615369905
   Exch tot contrib : Perdew-Wang 91       :      -8.990473984394
   Exch tot contrib : Exact Nonlocal Excha :      -9.017938328962
   Corr tot contrib : Vosko-Wilk-Nusair    :      -0.668410603779
   Corr tot contrib : Lee-Yang-Parr        :      -0.335692400200
   Corr tot contrib : Perdew-Burke-Ernzerh :      -0.338550817676
   Corr tot contrib : Perdew-Wang 91       :      -0.361265492963
   Corr tot contrib : Wilson-Levy          :      -0.395889288696
   Corr tot contrib : Wilson-Ivanov        :      -0.346192275446
   Hybr tot contrib : Becke III LYP        :      -9.344692540498
     
 Ehar from   KS   Calculation:         -65.940292164510 a.u.
 Ex   from   KS   Calculation:          -9.015964467082 a.u.
 Ec   from   KS   Calculation:          -0.335692400200 a.u.
 Exc  from   KS   Calculation:          -9.351656867282 a.u.
 Final Result from KS Calculation:     -75.291949031792 a.u.

 Ehar with the SCF Dens and Func:      -65.940292164510 a.u.
 Ex   with the SCF Dens and Func:       -9.015964467082 a.u.
 Ec   with the SCF Dens and Func:       -0.335692400200 a.u.
 Exc  with the SCF Dens and Func:       -9.351656867282 a.u.
 Final Result with the SCF Dens and Func:    -75.291949031792 a.u.

  @CHECKOUT-I, Total execution time :       0.3489 seconds.
ACES2: Total elapsed time is       13.1 seconds
 @ACES2: Executing "xvdint"
 One- and two-electron integral derivatives are calculated
 for RHF gradients and dipole moments.
 Spherical gaussians are used.

  Cartesian Coordinates
  ---------------------

  Total number of coordinates:  9


   1   H #1     x      0.0000000000
   2            y     -1.4684151899
   3            z      1.1698525287

   4   O #2     x      0.0000000000
   5            y      0.0000000000
   6            z     -0.1474227269

   7   H #3     x      0.0000000000
   8            y      1.4684151899
   9            z      1.1698525287

 Translational invariance is used.
 
 Entering the ECP_DGRAD

                    Nuclear attraction energy gradient
                    ----------------------------------

 H #1       0.0000000000           -3.7581005248            3.4396925092
 O #2       0.0000000000            0.0000000000           -6.8793850183
 H #3       0.0000000000            3.7581005248            3.4396925092



                       Two electron energy gradient
                       ----------------------------

 H #1       0.0000000000            1.8376488905           -1.7995033505
 O #2       0.0000000000            0.0000000000            3.5990067009
 H #3       0.0000000000           -1.8376488905           -1.7995033505



                    Nuclear repulsion energy gradient
                    ---------------------------------

 H #1       0.0000000000            1.6462208253           -1.3727710819
 O #2       0.0000000000            0.0000000000            2.7455421637
 H #3       0.0000000000           -1.6462208253           -1.3727710819



                      Kinetic energy energy gradient
                      ------------------------------

 H #1       0.0000000000            0.3511352310           -0.3226720123
 O #2       0.0000000000            0.0000000000            0.6453440246
 H #3       0.0000000000           -0.3511352310           -0.3226720123



                     Renormalization energy gradient
                     -------------------------------

 H #1       0.0000000000            0.0890860603           -0.1093746590
 O #2       0.0000000000            0.0000000000            0.2187493179
 H #3       0.0000000000           -0.0890860603           -0.1093746590



                  One elctron Molecular energy gradient
                  -------------------------------------

 H #1       0.0000000000           -1.6716584081            1.6348747561
 O #2       0.0000000000            0.0000000000           -3.2697495121
 H #3       0.0000000000            1.6716584081            1.6348747561



                        Molecular energy gradient
                        -------------------------

 H #1       0.0000000000            0.1659904823           -0.1646285944
 O #2       0.0000000000            0.0000000000            0.3292571888
 H #3       0.0000000000           -0.1659904823           -0.1646285944




ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.01 (      0.0)   0: 0: 0.00 (      0.0)   0: 0: 0.0 (      0.0)
 @ACES2: Executing "xvksdint"
 fffg
 
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
           potradpts     integer                    94                    50
           OEPSLATER     integer                     5
             overind     integer                     3
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
               ksmem     integer             150000000               1500000
                func      string  BECKE,LYP             none                
               kspot      string  BECKE,LYP             hf                  
              cutoff     integer                    12
============================================================================
 
  idd=                     1
 The new gradient:
 atm    x               y               z
  1  1    0.0000000000    0.1659904823   -0.1646285944
  8  2    0.0000000000    0.0000000000    0.3292571888
  1  3    0.0000000000   -0.1659904823   -0.1646285944
                     1                     1                     1
                     2                     8                     2
                     3                     1                     3
 t=                     1                     1
  1.523908443147903E-016 -0.165830024260978       0.164418486510923     
 t=                     2                     2
 -2.016409741624356E-016 -2.109423746787797E-015 -0.328836973021835     
 t=                     3                     3
  5.092642421820840E-017  0.165830024260981       0.164418486510917     
 rotin=  9.539320721533561E-016 -1.265807391629542E-014  1.742239192858941E-015
 tran=  1.676294370562644E-018  1.676294370562644E-018  4.468647674116255E-015
  1  1    0.0000000000    0.0001604581   -0.0002101079
  8  2    0.0000000000    0.0000000000    0.0004202158
  1  3    0.0000000000   -0.0001604581   -0.0002101079
  @CHECKOUT-I, Total execution time :       0.5849 seconds.


ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.42 (      0.4)   0: 0: 0.59 (      0.6)   0: 0: 1.0 (      1.0)
ACES2: Total elapsed time is       14.2 seconds
 @ACES2: Executing "xjoda"
 Finite diffs; ignore and geomopt vars: T T
 The vib calc. related varrs:
 iflags(h_IFLAGS_vib):                     0
 The finite diffs: F
 iflags2(h_IFLAGS2_geom_opt):                     1
 first run of popt num. frq:                     0
 The value of HAVEGEOM record:                     0
 The value of PES_SCAN record:                     0
 The OPTARC is here? T F
 Do we have geom? F
 Hessian in JOBARC?                    -1
Renaming SAVEDIR/CURRENT to SAVEDIR/OLD . . . done
Files copied to SAVEDIR/CURRENT: ZMAT.BAS JOBARC JAINDX OPTARC !OPTARCBK !DIPDER
Removing old back up directory . . . done
 A geometry for frequency {0,1} present:                     1
 The value of grad calc                     1
 The Deriv. level, havegeom?, and vib? :                     1
                     0                     1                     0
 The derivative level @exit                     1
 First order props and grad flags @enter:                     0
                     1
 First order props, deriv level and grad_calc @exit:                     0
                     1                     1
 @-Entry, the optimization cycle                     3
   JODA beginning optimization cycle #  4.
 -----After call to Entry-----
 The FNDDONE in Geopt:                     1
 The gradient calcs:                     1
 Hessian calc.,Geo. opt. and iarch: F T                     1
 Entering Fetchz IARCH .NE. 1                     1
 
 @GEOPT, PASS1,  PS1EXIST=:                     0 F
   Retrieving information from last optimization cycle.
 
 Start Reading Archive file NXM6:                      3
 The COORD COMMON BLOCK AT RETRIEVE BEFORE DECOMPRESS
   0.00000  -1.46842   1.16985
   0.00000   0.00000  -0.14742
   0.00000   1.46842   1.16985
 
   1.97268   1.97268   1.67920
 The COORD COMMON BLOCK AT RETRIEVE
   0.00000  -1.46842   1.16985
   0.00000   0.00000  -0.14742
   0.00000   1.46842   1.16985
 
   0.00000   0.00000   0.00000
 Data read from OPTARC file
 
 The COORD COMMON BLOCK/AFTER RETRIVE
   0.00000  -1.46842   1.16985
   0.00000   0.00000  -0.14742
   0.00000   1.46842   1.16985
 
   0.00000   0.00000   0.00000
 
 
 The Cartesians, internals and connec. @-buildb
 The NATOMS                     3
   0.00000  -1.46842   1.16985
   0.00000   0.00000  -0.14742
   0.00000   1.46842   1.16985
 
   0.00000   0.00000   0.00000
   1.97268   0.00000   0.00000
   1.97268   1.67920   0.00000
    0    0    0    1
    0    0    2    1
    0
  B MATRIX 
[  1,  1]  0.000000 [  1,  2] -0.744377 [  1,  3]  0.667760 [  1,  4]  0.000000
[  1,  5]  0.744377 [  1,  6] -0.667760 [  1,  7]  0.000000 [  1,  8]  0.000000
[  1,  9]  0.000000 [  2,  1]  0.000000 [  2,  2]  0.000000 [  2,  3]  0.000000
[  2,  4]  0.000000 [  2,  5] -0.744377 [  2,  6] -0.667760 [  2,  7]  0.000000
[  2,  8]  0.744377 [  2,  9]  0.667760 [  3,  1]  0.000000 [  3,  2] -0.338504
[  3,  3] -0.377343 [  3,  4]  0.000000 [  3,  5]  0.000000 [  3,  6]  0.754687
[  3,  7]  0.000000 [  3,  8]  0.338504 [  3,  9] -0.377343 [
  Bt MATRIX
[  1,  1]  0.000000 [  1,  2]  0.000000 [  1,  3]  0.000000 [  2,  1] -0.744377
[  2,  2]  0.000000 [  2,  3] -0.338504 [  3,  1]  0.667760 [  3,  2]  0.000000
[  3,  3] -0.377343 [  4,  1]  0.000000 [  4,  2]  0.000000 [  4,  3]  0.000000
[  5,  1]  0.744377 [  5,  2] -0.744377 [  5,  3]  0.000000 [  6,  1] -0.667760
[  6,  2] -0.667760 [  6,  3]  0.754687 [  7,  1]  0.000000 [  7,  2]  0.000000
[  7,  3]  0.000000 [  8,  1]  0.000000 [  8,  2]  0.744377 [  8,  3]  0.338504
[  9,  1]  0.000000 [  9,  2]  0.667760 [  9,  3] -0.377343 [
  G MATRIX 
[  1,  1]  2.000000 [  2,  1] -0.108193 [  2,  2]  2.000000 [  3,  1] -0.503949
[  3,  2] -0.503949 [  3,  3]  1.083498 [
  G MATRIX DETERMINANT   0.000000000000000E+000
  INVERSE G MATRIX 
[  1,  1]  0.588535 [  2,  1]  0.114196 [  2,  2]  0.588535 [  3,  1]  0.326850
[  3,  2]  0.326850 [  3,  3]  1.226981 [
  A MATRIX
[  1,  1]  0.000000 [  1,  2]  0.000000 [  1,  3]  0.000000 [  2,  1] -0.548732
[  2,  2] -0.195645 [  2,  3] -0.658638 [  3,  1]  0.269666 [  3,  2] -0.047079
[  3,  3] -0.244736 [  4,  1]  0.000000 [  4,  2]  0.000000 [  4,  3]  0.000000
[  5,  1]  0.353087 [  5,  2] -0.353087 [  5,  3]  0.000000 [  6,  1] -0.222587
[  6,  2] -0.222587 [  6,  3]  0.489472 [  7,  1]  0.000000 [  7,  2]  0.000000
[  7,  3]  0.000000 [  8,  1]  0.195645 [  8,  2]  0.548732 [  8,  3]  0.658638
[  9,  1] -0.047079 [  9,  2]  0.269666 [  9,  3] -0.244736 [
 
 The COORD COMMON BLOCK/at opt. start:
 
   1.97268   1.97268   1.67920
 
 The Cart. Hessian in after reading from READGH
 IAVHES flags: Cartesian exact Hess. read                     0
                     0
 If the Hessian is from VIB=EXACT, the trans/rot
 contaminants are included. If VIB=FINDIF it is
 projected.

                 COLUMN   1       COLUMN   2       COLUMN   3       COLUMN   4
 ROW   1      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   2      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   3      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   4      0.00000000000   -0.74437663955    0.66776000066    0.00000000000
 ROW   5      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   6      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   7      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   8      0.00000000000   -0.33850436085   -0.37734326458    0.00000000000
 ROW   9      0.00000000000    0.00000000000    0.00000000000    0.00000000000

                 COLUMN   5       COLUMN   6       COLUMN   7       COLUMN   8
 ROW   1      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   2      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   3      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   4      0.74437663955   -0.66776000066    0.00000000000    0.00000000000
 ROW   5      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   6      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   7     -0.74437663955   -0.66776000066    0.00000000000    0.74437663955
 ROW   8      0.00000000000    0.75468652916    0.00000000000    0.33850436085
 ROW   9      0.00000000000    0.00000000000    0.00000000000    0.00000000000

                 COLUMN   9
 ROW   1      0.00000000000
 ROW   2      0.00000000000
 ROW   3      0.00000000000
 ROW   4      0.00000000000
 ROW   5      0.00000000000
 ROW   6      0.00000000000
 ROW   7      0.66776000066
 ROW   8     -0.37734326458
 ROW   9      0.00000000000
 The Cartesian Hessian after reading from GETICFCM
 
 The COORD COMMON BLOCK/at opt. start:
 
   1.97268   1.97268   1.67920
 
 
 @-CONVQF,  The outgoing internal gradients
   -0.0002597   -0.0002597    0.0000972
 
   Internal coordinate forces and energy gradients (atomic units): 
           R       dV/dR             R       dV/dR             R       dV/dR
 
 [R    ]   1.97268  -0.00026 [R    ]   1.97268  -0.00026 [A    ]   1.67920   0.00010
 
   Hessian from cycle  3 read.
   BFGS update using last two gradients and previous step.
 The number of opt. cycles:                     3
 
   The eigenvectors of the Hessian matrix: 

                 COLUMN   1       COLUMN   2
 ROW   1      0.24293151127    0.97004344276
 ROW   2     -0.97004344276    0.24293151127
   The eigenvalues of the Hessian matrix: 
 
     0.24736    0.43972
 
 There are  0 Negative Eigenvalues.
 
 The number of degs. of freed.           at NR:    2
 
 The unmolested step
    0.0009   -0.0005
   MANR scale factor for NR step is 1.001 for R    .
 
 The initial trust rad:   0.42426
 The predicted and actual energy change:   -0.00003  -0.00003
 Ratio of change:   0.86525
 The new trust rad and TAU:   0.42426   0.00106
   Summary of Optimization Cycle: 
   The maximum unscaled step is:    0.00106.
   Scale factor set to:  1.00000.
   Forces are in hartree/bohr and hartree/radian.
   Parameter values are in Angstroms and degrees.
--------------------------------------------------------------------------
  Parameter     dV/dR           Step          Rold            Rnew
--------------------------------------------------------------------------
    R      -0.0002597429    0.0003424595    1.0438961532    1.0442386127
    A       0.0000971582   -0.0306955669   96.2111697916   96.1804742248
--------------------------------------------------------------------------
  Minimum force:  0.000097158 / RMS force:  0.000196094
 GMETRY starting with R
  1     1.973325
  2     1.973325
  3     1.678666
  4     0.000000
  5     0.000000
  6     0.000000
  7     0.000000
  8     0.000000
  9     0.000000
 In Gmetry; ncycle, ix, and ipost_vib
                     4                     0                     0
@GMETRY-I, Decompressing R.
 Updating structure...
 GMETRY using R vector
  1     0.000000
  2     0.000000
  3     0.000000
  4     1.973325
  5     0.000000
  6     0.000000
  7     1.973325
  8     1.678666
  9     0.000000
  @GMETRY-I, Cartesian coordinates before scaling:
      1     0.000000     0.000000     0.000000
      2     1.973325     0.000000     0.000000
      3     2.185774     1.961855     0.000000
  @GMETRY-I, Cartesian coordinates after scaling:
      1     0.000000     0.000000     0.000000
      2     1.973325     0.000000     0.000000
      3     2.185774     1.961855     0.000000
  @GMETRY-I, Cartesian coordinates from Z-matrix:
      1     0.000000     0.000000     0.000000
      2     0.000000     0.000000     1.973325
      3     1.961855     0.000000     2.185774
         -----------------------------------------------
         Cartesian coordinates corresponding to internal 
                 coordinate input (Angstroms) 
 ----------------------------------------------------------------
 Z-matrix   Atomic            C o o r d i n a t e s
  Symbol    Number           X              Y              Z
 ----------------------------------------------------------------
     H         1         0.00000000     0.00000000     0.00000000
     O         8         0.00000000     0.00000000     1.04423861
     H         1         1.03816920     0.00000000     1.15666192
 ----------------------------------------------------------------
    1.0078250000
   15.9949100000
    1.0078250000
The number of atoms                      3
 
 The Cartesians before translations to CM
   0.00000   0.00000   0.00000
   0.00000   0.00000   1.97332
   1.96186   0.00000   2.18577
 
 The Cartesians in center of mass coords
  -0.1097804   0.0000000  -1.8747907
  -0.1097804   0.0000000   0.0985341
   1.8520749   0.0000000   0.3109834
 
 After translation to center of mass coordinates 
       -0.109780419306    0.000000000000   -1.874790737469
       -0.109780419306    0.000000000000    0.098534107266
        1.852074891619    0.000000000000    0.310983352605
 
  @symmetry-i, Coordinates after  COM shift 
     -0.109780419306      0.000000000000     -1.874790737469
     -0.109780419306      0.000000000000      0.098534107266
      1.852074891619      0.000000000000      0.310983352605
 
  Inertia tensor
    3.79511    0.00000   -0.61488
    0.00000    7.45704    0.00000
   -0.61488    0.00000    3.66194
 Inertia tensor:
        3.795105447177    0.000000000000   -0.614878408392
        0.000000000000    7.457040570211    0.000000000000
       -0.614878408392    0.000000000000    3.661935123034
    Diagonalized inertia tensor:
        3.110047132160    0.000000000000    0.000000000000
        0.000000000000    4.346993438051    0.000000000000
        0.000000000000    0.000000000000    7.457040570211
    Eigenvectors of inertia tensor: 
        0.667959372394    0.744197740410    0.000000000000
        0.000000000000    0.000000000000    1.000000000000
        0.744197740410   -0.667959372394    0.000000000000
    Principal axis orientation for molecular system: 
       -1.468543890546    1.170585704382    0.000000000000
        0.000000000000   -0.147515120437    0.000000000000
        1.468543890546    1.170585704382    0.000000000000
   Rotational constants (in cm-1): 
     8.07264       13.84819       19.35598
  @SYMMETRY-I, Handedness of inertial frame: 1.00000
  @SYMMETRY-I, The symmetry group is 1-fold degenerate.
 
 The symmetry processing begins
 The Cartesians before entering symmetry auto
  -1.46854   1.17059   0.00000
   0.00000  -0.14752   0.00000
   1.46854   1.17059   0.00000
   0.00000   0.00000   0.00000
   0.00000   0.00000   0.00000
   0.00000   0.00000   0.00000
   Principal axis orientation for molecule:
       -1.468543890546    1.170585704382    0.000000000000
        0.000000000000   -0.147515120437    0.000000000000
        1.468543890546    1.170585704382    0.000000000000
   Reflection in plane  3 is a valid symmetry operation.
       -1.000000000000    0.000000000000    0.000000000000
        0.000000000000   -1.000000000000    0.000000000000
        0.000000000000    0.000000000000    1.000000000000
       -1.000000000000    0.000000000000    0.000000000000
        0.000000000000    1.000000000000    0.000000000000
        0.000000000000    0.000000000000   -1.000000000000
   Rotation about  2 is a valid symmetry operation 
   Reflection in plane  1 is a valid symmetry operation.
        1.000000000000    0.000000000000    0.000000000000
        0.000000000000   -1.000000000000    0.000000000000
        0.000000000000    0.000000000000   -1.000000000000
  Symmetry bits:    5   2   0
 @-SYM_AUTO The Orientation matrices: ORIEN2 and ORIENT
 

                 COLUMN   1       COLUMN   2       COLUMN   3
 ROW   1      0.66795937239    0.74419774041    0.00000000000
 ROW   2      0.00000000000    0.00000000000    1.00000000000
 ROW   3      0.74419774041   -0.66795937239    0.00000000000
 

                 COLUMN   1       COLUMN   2       COLUMN   3
 ROW   1      1.00000000000    0.00000000000    0.00000000000
 ROW   2      0.00000000000    1.00000000000    0.00000000000
 ROW   3      0.00000000000    0.00000000000    1.00000000000
 
 @-SYM_AUTO The variables in /COORD/ common block
   0.00000  -1.46854   1.17059
   0.00000   0.00000  -0.14752
   0.00000   1.46854   1.17059
 
 The NEWQ array 
  -1.46854
   1.17059
   0.00000
   0.00000
  -0.14752
   0.00000
   1.46854
   1.17059
   0.00000
 
 
 @-SYM_AUTO ORIENT2=ORIEN2xORIENT

                 COLUMN   1       COLUMN   2       COLUMN   3
 ROW   1      0.66795937239    0.74419774041    0.00000000000
 ROW   2      0.00000000000    0.00000000000    1.00000000000
 ROW   3      0.74419774041   -0.66795937239    0.00000000000
  @CHRTAB-I, Generated transformation matrices for   1 symmetry operations.
  @TSTOPS: Results of symmetry operation testing:
--------------------------------------------------------------------------------
   Sym. Op.         Class     Trace     Result 
--------------------------------------------------------------------------------
 --------------B
    0.00000   -1.46854    1.17059
    0.00000    0.00000   -0.14752
    0.00000    1.46854    1.17059
       1               1      3.00000   Passed
--------------------------------------------------------------------------------
 @TSTOPS:                      1  operations verified.
  @SYMUNQ-I, There are   3 orbits in the C1   point group.
------------------------------------------------------------------------
    Set    Site Group                       Members of set
------------------------------------------------------------------------
      1      C1                1
      2      C1                2
      3      C1                3
------------------------------------------------------------------------
  @CHRTAB-I, Generated transformation matrices for   4 symmetry operations.
  @TSTOPS: Results of symmetry operation testing:
--------------------------------------------------------------------------------
   Sym. Op.         Class     Trace     Result 
--------------------------------------------------------------------------------
 --------------B
    0.00000    1.46854    1.17059
    0.00000    0.00000   -0.14752
    0.00000   -1.46854    1.17059
       1               1     -1.00000   Passed
 --------------B
    0.00000    1.46854    1.17059
    0.00000    0.00000   -0.14752
    0.00000   -1.46854    1.17059
       2               2      1.00000   Passed
 --------------B
    0.00000   -1.46854    1.17059
    0.00000    0.00000   -0.14752
    0.00000    1.46854    1.17059
       3               3      1.00000   Passed
 --------------B
    0.00000   -1.46854    1.17059
    0.00000    0.00000   -0.14752
    0.00000    1.46854    1.17059
       4               4      3.00000   Passed
--------------------------------------------------------------------------------
 @TSTOPS:                      4  operations verified.
  @SYMUNQ-I, There are   2 orbits in the C2v  point group.
------------------------------------------------------------------------
    Set    Site Group                       Members of set
------------------------------------------------------------------------
      1      C s               1  3
      2      C2v               2
------------------------------------------------------------------------
         -----------------------------------------------
         Cartesian coordinates corresponding to internal 
                 coordinate input (Angstroms) 
 ----------------------------------------------------------------
 Z-matrix   Atomic            C o o r d i n a t e s
  Symbol    Number           X              Y              Z
 ----------------------------------------------------------------
     H         1         0.00000000    -0.77712002     0.61944732
     O         8         0.00000000     0.00000000    -0.07806165
     H         1         0.00000000     0.77712002     0.61944732
 ----------------------------------------------------------------
   Interatomic distance matrix (Angstroms) 
 
                 H             O             H    
                [ 1]        [ 2]        [ 3]
  H    [ 1]     0.00000
  O    [ 2]     1.04424     0.00000
  H    [ 3]     1.55424     1.04424     0.00000

   Interatomic angles (degrees) 


  O [ 2]-H [ 1]-H [ 3]   41.90976
  H [ 1]-O [ 2]-H [ 3]   96.18047

  H [ 1]-H [ 3]-O [ 2]   41.90976

      3 interatomic angles printed.
 Enetring Arhcive
 
 The COORD COMMON BLOCK AT ARCHIVE
   0.00000  -1.46854   1.17059
   0.00000   0.00000  -0.14752
   0.00000   1.46854   1.17059
 
   0.00000   0.00000   0.00000
   1.97332   0.00000   0.00000
   1.97332   1.67867   0.00000
 
 @-ARCHIVE, the number of opt. cycles                     4
   1.97332   1.97332   1.67867
 @ACES2: Executing "xvmol"
 One- and two-electron integrals over symmetry-adapted AOs are calculated.\n
 @READIN: Spherical harmonics are used.
  @READIN-I, Nuclear repulsion energy :    8.4486162970 a.u.
  required memory for a1 array                4451690  words 
  required memory for a2 array                   9156  words 
  @TWOEL-I,        228 integrals of symmetry type  I I I I
  @TWOEL-I, Total number of 2-e integrals        228.


ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.01 (      0.0)   0: 0: 0.01 (      0.0)   0: 0: 0.0 (      0.0)
 @ACES2: Executing "xvmol2ja"


ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.00 (      0.0)   0: 0: 0.00 (      0.0)   0: 0: 0.0 (      0.0)
 @ACES2: Executing "xvscf"
          Alpha population by irrep:   5
           Beta population by irrep:   5


 @SORTHO-I, Orthonormalizing initial guess. 
The Alpha and Beta occupation vector
                     5
 
                     5
 
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
           potradpts     integer                    94                    50
           OEPSLATER     integer                     5
             overind     integer                     3
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
               ksmem     integer             150000000               1500000
                func      string  BECKE,LYP             none                
               kspot      string  BECKE,LYP             hf                  
              cutoff     integer                    12
============================================================================
 
  --------------------------------------------------------------------
  Iteration         Total Energy            Largest Error in FDS-SDF
  --------------------------------------------------------------------
       0           -65.9404579788              0.2157661546D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       1           -71.6831650289              0.4256095834D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       2           -74.6022061013              0.3924844006D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       3           -74.4474031863              0.4789834400D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       4           -74.4678481877              0.4401155479D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       5           -74.3370891827              0.5056883664D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       6           -74.4056726574              0.4557279144D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       7           -74.2962053236              0.5135443155D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       8           -74.3847453981              0.4603263596D+00
The Alpha and Beta occupation vector
                     5
 
                     5
       9           -74.9551900695              0.5018598589D-01
The Alpha and Beta occupation vector
                     5
 
                     5
      10           -74.9581976549              0.1542989746D-02
The Alpha and Beta occupation vector
                     5
 
                     5
      11           -74.9581470763              0.1564741368D-03
The Alpha and Beta occupation vector
                     5
 
                     5
      12           -74.9581525322              0.1841840726D-05
The Alpha and Beta occupation vector
                     5
 
                     5
      13           -74.9581524634              0.1014550579D-06
The Alpha and Beta occupation vector
                     5
 
                     5
      14           -74.9581524671              0.1263500082D-08
The Alpha and Beta occupation vector
                     5
 
                     5

@VSCF: SCF has converged.

     E(SCF)=       -74.9581524671

 The AOBASMOS? T
 Vib calc & AOBASMOS?                     0 T
 Creating The AOBASMOS file


      Trace of projected alpha density matrix =   1.000000000

      Alpha part of wavefunction is symmetric.




ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.61 (      0.6)   0: 0: 0.94 (      0.9)   0: 0: 1.6 (      1.6)
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
           potradpts     integer                    94                    50
           OEPSLATER     integer                     5
             overind     integer                     3
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
               ksmem     integer             150000000               1500000
                func      string  BECKE,LYP             none                
               kspot      string  BECKE,LYP             hf                  
              cutoff     integer                    12
============================================================================
 

The SCF nuclear-electron attraction energy =   -195.389945883638
The SCF kinetic energy                     =     74.442061245670
The SCF coulomb energy                     =     46.558745500774
The SCF exchange energy                    =     -9.017629626891
The SCF one electron energy                =   -120.947884637968
The SCF total energy                       =    -74.958152467066

Results using the SCF density

   Total density integrates to             :       9.999986564848 electrons
   Exch tot contrib : LDA (Slater, Xalpha) :      -8.133537788633
   Exch tot contrib : Becke                :      -9.015747129111
   Exch tot contrib : Perdew-Burke-Ernzerh :      -8.959397253668
   Exch tot contrib : Perdew-Wang 91       :      -8.990255601631
   Exch tot contrib : Exact Nonlocal Excha :      -9.017629626891
   Corr tot contrib : Vosko-Wilk-Nusair    :      -0.668397637010
   Corr tot contrib : Lee-Yang-Parr        :      -0.335679192171
   Corr tot contrib : Perdew-Burke-Ernzerh :      -0.338531702108
   Corr tot contrib : Perdew-Wang 91       :      -0.361246551035
   Corr tot contrib : Wilson-Levy          :      -0.395857822506
   Corr tot contrib : Wilson-Ivanov        :      -0.346175493857
   Hybr tot contrib : Becke III LYP        :      -9.344442578119
     
 Ehar from   KS   Calculation:         -65.940522840176 a.u.
 Ex   from   KS   Calculation:          -9.015747129111 a.u.
 Ec   from   KS   Calculation:          -0.335679192171 a.u.
 Exc  from   KS   Calculation:          -9.351426321282 a.u.
 Final Result from KS Calculation:     -75.291949161457 a.u.

 Ehar with the SCF Dens and Func:      -65.940522840176 a.u.
 Ex   with the SCF Dens and Func:       -9.015747129111 a.u.
 Ec   with the SCF Dens and Func:       -0.335679192171 a.u.
 Exc  with the SCF Dens and Func:       -9.351426321282 a.u.
 Final Result with the SCF Dens and Func:    -75.291949161457 a.u.

  @CHECKOUT-I, Total execution time :       0.3459 seconds.
ACES2: Total elapsed time is       16.6 seconds
 @ACES2: Executing "xvdint"
 One- and two-electron integral derivatives are calculated
 for RHF gradients and dipole moments.
 Spherical gaussians are used.

  Cartesian Coordinates
  ---------------------

  Total number of coordinates:  9


   1   H #1     x      0.0000000000
   2            y     -1.4685438905
   3            z      1.1705857044

   4   O #2     x      0.0000000000
   5            y      0.0000000000
   6            z     -0.1475151204

   7   H #3     x      0.0000000000
   8            y      1.4685438905
   9            z      1.1705857044

 Translational invariance is used.
 
 Entering the ECP_DGRAD

                    Nuclear attraction energy gradient
                    ----------------------------------

 H #1       0.0000000000           -3.7560311518            3.4401057770
 O #2       0.0000000000            0.0000000000           -6.8802115540
 H #3       0.0000000000            3.7560311518            3.4401057770



                       Two electron energy gradient
                       ----------------------------

 H #1       0.0000000000            1.8370693128           -1.8003096123
 O #2       0.0000000000            0.0000000000            3.6006192246
 H #3       0.0000000000           -1.8370693128           -1.8003096123



                    Nuclear repulsion energy gradient
                    ---------------------------------

 H #1       0.0000000000            1.6448294179           -1.3722804218
 O #2       0.0000000000            0.0000000000            2.7445608436
 H #3       0.0000000000           -1.6448294179           -1.3722804218



                      Kinetic energy energy gradient
                      ------------------------------

 H #1       0.0000000000            0.3508216460           -0.3226189929
 O #2       0.0000000000            0.0000000000            0.6452379858
 H #3       0.0000000000           -0.3508216460           -0.3226189929



                     Renormalization energy gradient
                     -------------------------------

 H #1       0.0000000000            0.0890622979           -0.1093657633
 O #2       0.0000000000            0.0000000000            0.2187315266
 H #3       0.0000000000           -0.0890622979           -0.1093657633



                  One elctron Molecular energy gradient
                  -------------------------------------

 H #1       0.0000000000           -1.6713177900            1.6358405990
 O #2       0.0000000000            0.0000000000           -3.2716811980
 H #3       0.0000000000            1.6713177900            1.6358405990



                        Molecular energy gradient
                        -------------------------

 H #1       0.0000000000            0.1657515229           -0.1644690133
 O #2       0.0000000000            0.0000000000            0.3289380266
 H #3       0.0000000000           -0.1657515229           -0.1644690133




ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.01 (      0.0)   0: 0: 0.00 (      0.0)   0: 0: 0.0 (      0.0)
 @ACES2: Executing "xvksdint"
 fffg
 
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
           potradpts     integer                    94                    50
           OEPSLATER     integer                     5
             overind     integer                     3
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
               ksmem     integer             150000000               1500000
                func      string  BECKE,LYP             none                
               kspot      string  BECKE,LYP             hf                  
              cutoff     integer                    12
============================================================================
 
  idd=                     1
 The new gradient:
 atm    x               y               z
  1  1    0.0000000000    0.1657515229   -0.1644690133
  8  2    0.0000000000    0.0000000000    0.3289380266
  1  3    0.0000000000   -0.1657515229   -0.1644690133
                     1                     1                     1
                     2                     8                     2
                     3                     1                     3
 t=                     1                     1
 -6.918401059355994E-017 -0.165745996115890       0.164471698954875     
 t=                     2                     2
  1.697655370781300E-016  3.264055692397960E-014 -0.328943397909707     
 t=                     3                     3
 -1.039365424120033E-016  0.165745996115855       0.164471698954827     
 rotin= -3.956317739595681E-016 -3.300373585264204E-014 -1.148722883658748E-015
 tran= -3.355015927433424E-018 -3.355015927433424E-018 -4.357625371653739E-015
  1  1    0.0000000000    0.0000055268    0.0000026856
  8  2    0.0000000000    0.0000000000   -0.0000053713
  1  3    0.0000000000   -0.0000055268    0.0000026856
  @CHECKOUT-I, Total execution time :       0.5699 seconds.


ACES ---SYSTEM--- (total sec) ----USER---- (total sec) -WALLCLOCK- (total sec)
TIME   0: 0: 0.44 (      0.4)   0: 0: 0.57 (      0.6)   0: 0: 1.0 (      1.0)
ACES2: Total elapsed time is       17.7 seconds
 @ACES2: Executing "xjoda"
 Finite diffs; ignore and geomopt vars: T T
 The vib calc. related varrs:
 iflags(h_IFLAGS_vib):                     0
 The finite diffs: F
 iflags2(h_IFLAGS2_geom_opt):                     1
 first run of popt num. frq:                     0
 The value of HAVEGEOM record:                     0
 The value of PES_SCAN record:                     0
 The OPTARC is here? T F
 Do we have geom? F
 Hessian in JOBARC?                    -1
Renaming SAVEDIR/CURRENT to SAVEDIR/OLD . . . done
Files copied to SAVEDIR/CURRENT: ZMAT.BAS JOBARC JAINDX OPTARC !OPTARCBK !DIPDER
Removing old back up directory . . . done
 A geometry for frequency {0,1} present:                     1
 The value of grad calc                     1
 The Deriv. level, havegeom?, and vib? :                     1
                     0                     1                     0
 The derivative level @exit                     1
 First order props and grad flags @enter:                     0
                     1
 First order props, deriv level and grad_calc @exit:                     0
                     1                     1
 @-Entry, the optimization cycle                     4
   JODA beginning optimization cycle #  5.
 -----After call to Entry-----
 The FNDDONE in Geopt:                     1
 The gradient calcs:                     1
 Hessian calc.,Geo. opt. and iarch: F T                     1
 Entering Fetchz IARCH .NE. 1                     1
 
 @GEOPT, PASS1,  PS1EXIST=:                     0 F
   Retrieving information from last optimization cycle.
 
 Start Reading Archive file NXM6:                      3
 The COORD COMMON BLOCK AT RETRIEVE BEFORE DECOMPRESS
   0.00000  -1.46854   1.17059
   0.00000   0.00000  -0.14752
   0.00000   1.46854   1.17059
 
   1.97332   1.97332   1.67867
 The COORD COMMON BLOCK AT RETRIEVE
   0.00000  -1.46854   1.17059
   0.00000   0.00000  -0.14752
   0.00000   1.46854   1.17059
 
   0.00000   0.00000   0.00000
 Data read from OPTARC file
 
 The COORD COMMON BLOCK/AFTER RETRIVE
   0.00000  -1.46854   1.17059
   0.00000   0.00000  -0.14752
   0.00000   1.46854   1.17059
 
   0.00000   0.00000   0.00000
 
 
 The Cartesians, internals and connec. @-buildb
 The NATOMS                     3
   0.00000  -1.46854   1.17059
   0.00000   0.00000  -0.14752
   0.00000   1.46854   1.17059
 
   0.00000   0.00000   0.00000
   1.97332   0.00000   0.00000
   1.97332   1.67867   0.00000
    0    0    0    1
    0    0    2    1
    0
  B MATRIX 
[  1,  1]  0.000000 [  1,  2] -0.744198 [  1,  3]  0.667959 [  1,  4]  0.000000
[  1,  5]  0.744198 [  1,  6] -0.667959 [  1,  7]  0.000000 [  1,  8]  0.000000
[  1,  9]  0.000000 [  2,  1]  0.000000 [  2,  2]  0.000000 [  2,  3]  0.000000
[  2,  4]  0.000000 [  2,  5] -0.744198 [  2,  6] -0.667959 [  2,  7]  0.000000
[  2,  8]  0.744198 [  2,  9]  0.667959 [  3,  1]  0.000000 [  3,  2] -0.338494
[  3,  3] -0.377129 [  3,  4]  0.000000 [  3,  5]  0.000000 [  3,  6]  0.754258
[  3,  7]  0.000000 [  3,  8]  0.338494 [  3,  9] -0.377129 [
  Bt MATRIX
[  1,  1]  0.000000 [  1,  2]  0.000000 [  1,  3]  0.000000 [  2,  1] -0.744198
[  2,  2]  0.000000 [  2,  3] -0.338494 [  3,  1]  0.667959 [  3,  2]  0.000000
[  3,  3] -0.377129 [  4,  1]  0.000000 [  4,  2]  0.000000 [  4,  3]  0.000000
[  5,  1]  0.744198 [  5,  2] -0.744198 [  5,  3]  0.000000 [  6,  1] -0.667959
[  6,  2] -0.667959 [  6,  3]  0.754258 [  7,  1]  0.000000 [  7,  2]  0.000000
[  7,  3]  0.000000 [  8,  1]  0.000000 [  8,  2]  0.744198 [  8,  3]  0.338494
[  9,  1]  0.000000 [  9,  2]  0.667959 [  9,  3] -0.377129 [
  G MATRIX 
[  1,  1]  2.000000 [  2,  1] -0.107661 [  2,  2]  2.000000 [  3,  1] -0.503814
[  3,  2] -0.503814 [  3,  3]  1.082514 [
  G MATRIX DETERMINANT   0.000000000000000E+000
  INVERSE G MATRIX 
[  1,  1]  0.588507 [  2,  1]  0.114047 [  2,  2]  0.588507 [  3,  1]  0.326976
[  3,  2]  0.326976 [  3,  3]  1.228132 [
  A MATRIX
[  1,  1]  0.000000 [  1,  2]  0.000000 [  1,  3]  0.000000 [  2,  1] -0.548645
[  2,  2] -0.195553 [  2,  3] -0.659050 [  3,  1]  0.269786 [  3,  2] -0.047133
[  3,  3] -0.244757 [  4,  1]  0.000000 [  4,  2]  0.000000 [  4,  3]  0.000000
[  5,  1]  0.353092 [  5,  2] -0.353092 [  5,  3]  0.000000 [  6,  1] -0.222653
[  6,  2] -0.222653 [  6,  3]  0.489515 [  7,  1]  0.000000 [  7,  2]  0.000000
[  7,  3]  0.000000 [  8,  1]  0.195553 [  8,  2]  0.548645 [  8,  3]  0.659050
[  9,  1] -0.047133 [  9,  2]  0.269786 [  9,  3] -0.244757 [
 
 The COORD COMMON BLOCK/at opt. start:
 
   1.97332   1.97332   1.67867
 
 The Cart. Hessian in after reading from READGH
 IAVHES flags: Cartesian exact Hess. read                     0
                     0
 If the Hessian is from VIB=EXACT, the trans/rot
 contaminants are included. If VIB=FINDIF it is
 projected.

                 COLUMN   1       COLUMN   2       COLUMN   3       COLUMN   4
 ROW   1      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   2      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   3      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   4      0.00000000000   -0.74419774041    0.66795937239    0.00000000000
 ROW   5      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   6      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   7      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   8      0.00000000000   -0.33849438129   -0.37712885559    0.00000000000
 ROW   9      0.00000000000    0.00000000000    0.00000000000    0.00000000000

                 COLUMN   5       COLUMN   6       COLUMN   7       COLUMN   8
 ROW   1      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   2      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   3      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   4      0.74419774041   -0.66795937239    0.00000000000    0.00000000000
 ROW   5      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   6      0.00000000000    0.00000000000    0.00000000000    0.00000000000
 ROW   7     -0.74419774041   -0.66795937239    0.00000000000    0.74419774041
 ROW   8      0.00000000000    0.75425771119    0.00000000000    0.33849438129
 ROW   9      0.00000000000    0.00000000000    0.00000000000    0.00000000000

                 COLUMN   9
 ROW   1      0.00000000000
 ROW   2      0.00000000000
 ROW   3      0.00000000000
 ROW   4      0.00000000000
 ROW   5      0.00000000000
 ROW   6      0.00000000000
 ROW   7      0.66795937239
 ROW   8     -0.37712885559
 ROW   9      0.00000000000
 The Cartesian Hessian after reading from GETICFCM
 
 The COORD COMMON BLOCK/at opt. start:
 
   1.97332   1.97332   1.67867
 
 
 @-CONVQF,  The outgoing internal gradients
   -0.0000023   -0.0000023   -0.0000112
 
   Internal coordinate forces and energy gradients (atomic units): 
           R       dV/dR             R       dV/dR             R       dV/dR
 
 [R    ]   1.97332   0.00000 [R    ]   1.97332   0.00000 [A    ]   1.67867  -0.00001
 
   Hessian from cycle  4 read.
   BFGS update using last two gradients and previous step.
 The number of opt. cycles:                     4
 
   The eigenvectors of the Hessian matrix: 

                 COLUMN   1       COLUMN   2
 ROW   1      0.21817974805    0.97590860102
 ROW   2     -0.97590860102    0.21817974805
   The eigenvalues of the Hessian matrix: 
 
     0.25633    0.42721
 
 There are  0 Negative Eigenvalues.
 
 The number of degs. of freed.           at NR:    2
 
 The unmolested step
    0.0000    0.0000
   MANR scale factor for NR step is 1.000 for R    .
 
 The initial trust rad:   0.42426
 The predicted and actual energy change:    0.00000   0.00000
 Ratio of change:   1.37785
 The new trust rad and TAU:   0.42426   0.00004
   Summary of Optimization Cycle: 
   The maximum unscaled step is:    0.00004.
   Scale factor set to:  1.00000.
   Forces are in hartree/bohr and hartree/radian.
   Parameter values are in Angstroms and degrees.
--------------------------------------------------------------------------
  Parameter     dV/dR           Step          Rold            Rnew
--------------------------------------------------------------------------
    R      -0.0000023191    0.0000015678    1.0442386127    1.0442401805
    A      -0.0000112288    0.0023996470   96.1804742248   96.1828738718
--------------------------------------------------------------------------
  Minimum force:  0.000002319 / RMS force:  0.000008108
--------------------------------------------------------------------------------
   RMS gradient is below .10000E-03.
   Convergence criterion satisfied.  Optimization completed.
--------------------------------------------------------------------------------
   Interatomic distance matrix (Angstroms) 
 
                 H             O             H    
                [ 1]        [ 2]        [ 3]
  H    [ 1]     0.00000
  O    [ 2]     1.04424     0.00000
  H    [ 3]     1.55424     1.04424     0.00000

   Interatomic angles (degrees) 


  O [ 2]-H [ 1]-H [ 3]   41.90976
  H [ 1]-O [ 2]-H [ 3]   96.18047

  H [ 1]-H [ 3]-O [ 2]   41.90976

      3 interatomic angles printed.
 
   Summary of optimized internal coordinates
   (Angstroms and degrees)
R    =     1.0442386127
A    =    96.1804742248
 
 The internal coord:
   1.973325   1.973325   1.678666
 
 The step size:
 
   0.000003   0.000003   0.000042
 
 The optimization cycle:                     5
 Joda flag, grad_calc and geoom_opt                     1                     1
 Internal flags; Hessian & geo. optimization: F F
   Frequencies of the updated Hessian at convergence:
 
    0.i    0.i    0.     0.     0.     0.  1940.  3401.  5296. 
 
   Warning: Frequencies based on the updated Hessian are not very accurate.
            A separate frequency calculation is required to get correct values.
 
 Removing SAVEDIR/CURRENT . . .
 The backup directory was successfully removed.
 @ACES2: Executing "xa2proc molden"
 @MOLDEN_MAIN: successfully created MOLDEN.INPUT                    


 @ACES2: The ACES2 program has completed successfully in                     19 
  seconds.


