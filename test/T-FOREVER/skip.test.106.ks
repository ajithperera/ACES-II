NMR SPIN-SPIN COUPLING CONSTANT
X
N 1 NX
H 2 NH* 1 A*
H 2 NH* 1 A* 3 D1
H 2 NH* 1 A* 3 D2

NX=1.000
NH=1.213
A=115.45
D1=120.
D2=-120.

*ACES2
BASIS=DZP,SPHERICAL=ON
SCF_TYPE=KS

*SCF
 ks
 print_nl=false

*INTGRT
 kspot=lda,vwn
*end

TEST.DAT
d TOTENERG
 -0.560843665644E+02
d COORD
    -2.017082434944
     0.000000000000
     0.000000000000
    -0.127356446365
     0.000000000000
     0.000000000000
     0.589844874986
    -0.900614245245
    -1.559909630785
     0.589844874986
     1.801228490490
     0.000000000000
     0.589844874986
    -0.900614245245
     1.559909630785
d NUCREP
  0.117932432361E+02
d SCFENEG
 -0.561998471086E+02
d KSSCFENG
 -0.561998471083E+02
d KSTOTELE
  0.000000000000E+00
