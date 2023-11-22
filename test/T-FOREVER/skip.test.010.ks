
# This test optimizes the geometry of H2O in a B3LYP functional.
# The potential used to generate the density is the standard HF
# coulomb repulsion and exchange.

TEST.010 NUMERICAL HFDFT OPTIMIZATION OF H2O
O
H 1 R*
H 1 R* 2 A*

R=0.943056
A=105.96879

*ACES2(CALC=SCF,SCF_TYPE=KS,GRAD_CALC=NUMERICAL
BASIS=DZP,MEMORY=4000000)

*SCF
 ks

*INTGRT
 kspot=hf
 func=b3lyp





