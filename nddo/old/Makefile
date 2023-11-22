objects=scratch.o indices.o constantsmodule.o sphericalmodule.o eig.o  tqli.o tred2.o tablesmodule.o   controlmodule.o cpusec.o dblas1.o dblas2.o dblas3.o driver.o eigsort.o default.o f1.o f2.o guess.o  hcore.o  insert.o jacobi.o local.o   overlap.o pack1.o pack2.o pm3_parameters.o pp.o ps.o \
     quadrant.o readinput.o rhf.o repel.o repam1.o rotate.o  ss.o trace.o twoe_electron.o  sto.o matprt.o atoi.o geometry.o rotcor.o  ptchg.o chgsep.o limit.o rtbis.o func.o parout.o driver1.o scfopt.o routines.o driver2.o scfq.o BDATA4.o REPPD.o ROTD.o ROTMAT.o factrl.o dparameters.o rearrange.o inighd.o rsc.o fbinom.o wstore.o ds.o dp.o dd.o doverlap.o loadparam.o matprt2.o  nodummy.o  cgfam.o cgsearch.o fcn.o main.o    am1_parameters.o tred3.o tql3.o square.o mkvector.o uhf.o dipole.o f1uhf.o f2uhf.o gradient_fast.o hcoregrad.o f1grad.o f2grad.o twoepair.o twoe_sparkle.o repelsparkle.o hcoregradsp.o gradientuhf.o f2graduhf.o hfdft.o mkd.o b3.o b88.o b3lyp.o fuzzy2.o lyp.o vwn.o pbe.o hfexc.o hessian.o initguess.o normalize.o buildgto.o ssinitial.o psinitial.o dsinitial.o ppinitial.o dpinitial.o ddinitial.o buildt.o invert.o initguessold.o

opt=  
compile90=f90
compile77=f77    
go: $(objects) ./LAPACK/*.o    
	$(compile90)  -o go $(opt) -xlang=f77  $(objects) ./LAPACK/*.o

constantsmodule.o: constantsmodule.f90
	$(compile90) -c  $(opt)  constantsmodule.f90

sphericalmodule.o: sphericalmodule.f90
	$(compile90) -c  $(opt)  sphericalmodule.f90

indices.o: indices.f90
	$(compile90) -c  $(opt)  indices.f90

scratch.o: scratch.f90
	$(compile90) -c  $(opt) scratch.f90

tablesmodule.o: tablesmodule.f90
	$(compile90) -c  $(opt)  tablesmodule.f90

controlmodule.o: controlmodule.f90
	$(compile90) -c  $(opt)  controlmodule.f90

eigsort.o: eigsort.f
	$(compile77) -c $(opt) eigsort.f

tred2.o: tred2.f
	$(compile77) -c $(opt)  tred2.f

tqli.o: tqli.f
	$(compile77) -c $(opt)  tqli.f

jacobi.o: jacobi.f
	$(compile77) -c $(opt) jacobi.f


cpusec.o: cpusec.f
	$(compile77) -c $(opt)  cpusec.f

driver.o: driver.f90 
	$(compile90) -c  $(opt)  driver.f90


eig.o: eig.f
	$(compile77) -c $(opt) eig.f

f1.o: f1.f90 indices.f90 constantsmodule.f90 tablesmodule.f90
	$(compile90) -c  $(opt) f1.f90

f1uhf.o: f1uhf.f90 indices.f90 constantsmodule.f90 tablesmodule.f90
	$(compile90) -c  $(opt) f1uhf.f90

f1grad.o: f1grad.f90 indices.f90 constantsmodule.f90
	$(compile90) -c  $(opt) f1grad.f90

f2.o: f2.f90 indices.f90 constantsmodule.f90 tablesmodule.f90
	$(compile90) -c  $(opt) f2.f90

f2uhf.o: f2uhf.f90 indices.f90 constantsmodule.f90 tablesmodule.f90
	$(compile90) -c  $(opt) f2uhf.f90

f2grad.o: f2grad.f90 indices.f90 constantsmodule.f90
	$(compile90) -c  $(opt) f2grad.f90

guess.o: guess.f90 indices.f90 constantsmodule.f90 tablesmodule.f90
	$(compile90) -c  $(opt) guess.f90


gradient_fast.o: gradient_fast.f90 constantsmodule.f90 indices.f90 controlmodule.f90 scratch.f90 tablesmodule.f90
	$(compile90) -c  $(opt) gradient_fast.f90

hcore.o: default.f90 hcore.f90 constantsmodule.f90 tablesmodule.f90 indices.f90
	$(compile90) -c  $(opt)  hcore.f90 

hcoregrad.o: hcoregrad.f90 constantsmodule.f90 tablesmodule.f90 indices.f90
	$(compile90) -c  $(opt)  hcoregrad.f90 

insert.o: insert.f90
	$(compile90) -c  $(opt)  insert.f90

local.o: local.f90 constantsmodule.f90
	$(compile90) -c  $(opt)  local.f90

overlap.o: overlap.f90 tablesmodule.f90 constantsmodule.f90 
	$(compile90) -c  $(opt)  overlap.f90

pack.o: pack.f90
	$(compile90) -c  $(opt)  pack.f90

pm3_parameters.o: pm3_parameters.f90
	$(compile90) -c  $(opt)  pm3_parameters.f90

pp.o: pp.f90 constantsmodule.f90
	$(compile90) -c  $(opt)  pp.f90

ps.o: ps.f90 constantsmodule.f90
	$(compile90) -c  $(opt)  ps.f90

quadrant.o: quadrant.f90 constantsmodule.f90
	$(compile90) -c  $(opt)  quadrant.f90

readinput.o: readinput.f90 constantsmodule.f90 tablesmodule.f90 scratch.f90 controlmodule.f90 indices.f90
	$(compile90) -c  $(opt)  readinput.f90

repel.o: repel.f90 constantsmodule.f90 tablesmodule.f90 controlmodule.f90
	$(compile90) -c  $(opt) repel.f90

repam1.o: repam1.f90 tablesmodule.f90 constantsmodule.f90
	$(compile90) -c  $(opt) repam1.f90

rotate.o: rotate.f90 constantsmodule.f90
	$(compile90) -c  $(opt)  rotate.f90

rhf.o: rhf.f90 constantsmodule.f90 indices.f90 scratch.f90 controlmodule.f90
	$(compile90) -c  $(opt)  rhf.f90

uhf.o: uhf.f90 constantsmodule.f90 indices.f90 scratch.f90 controlmodule.f90
	$(compile90) -c  $(opt)  uhf.f90

rhfoverlap.o: rhfoverlap.f90 constantsmodule.f90 indices.f90 scratch.f90 controlmodule.f90 tablesmodule.f90
	$(compile90) -c  $(opt)  rhfoverlap.f90

ss.o: ss.f90 constantsmodule.f90
	$(compile90) -c  $(opt)  ss.f90


trace.o: trace.f90
	$(compile90) -c  $(opt) trace.f90

twoe_electron.o: twoe_electron.f90 constantsmodule.f90 indices.f90 controlmodule.f90 scratch.f90 tablesmodule.f90
	$(compile90) -c  $(opt)  twoe_electron.f90

twoepair.o: twoepair.f90 constantsmodule.f90 indices.f90 controlmodule.f90 scratch.f90 tablesmodule.f90
	$(compile90) -c  $(opt)  twoepair.f90

sto.o: sto.f
	$(compile77) -c $(opt)  sto.f

matprt.o: matprt.f
	$(compile77) -c $(opt)  matprt.f

matprt2.o: matprt2.f
	$(compile77) -c $(opt)  matprt2.f

dblas1.o: dblas1.f
	$(compile77) -c  $(opt) dblas1.f

dblas2.o: dblas2.f
	$(compile77) -c  $(opt) dblas2.f

dblas3.o: dblas3.f
	$(compile77) -c  $(opt) dblas3.f

atoi.o: atoi.f90
	$(compile90) -c  $(opt) atoi.f90

geometry.o: geometry.f90 constantsmodule.f90 tablesmodule.f90
	$(compile90) -c  $(opt) geometry.f90

rotcor.o: constantsmodule.f90 
	$(compile90) -c  $(opt) rotcor.f90

ptchg.o: ptchg.f90 constantsmodule.f90 tablesmodule.f90 controlmodule.f90
	$(compile90) -c  $(opt) ptchg.f90

chgsep.o: chgsep.f90 constantsmodule.f90 
	$(compile90) -c  $(opt) chgsep.f90

limit.o: limit.f90 constantsmodule.f90 
	$(compile90) -c  $(opt) limit.f90

rtbis.o: rtbis.f
	$(compile77) -c $(opt) rtbis.f

func.o: func.f90 constantsmodule.f90 
	$(compile90) -c  $(opt) func.f90

parout.o: parout.f90 constantsmodule.f90 tablesmodule.f90 controlmodule.f90
	$(compile90) -c  $(opt) parout.f90

driver1.o: driver1.f
	$(compile77) -c $(opt) driver1.f

driver2.o: driver2.f
	$(compile77) -c $(opt) driver2.f

routines.o: routines.f
	$(compile77) -c $(opt) routines.f

scfopt.o: scfopt.f90 tablesmodule.f90 controlmodule.f90 constantsmodule.f90 indices.f90
	$(compile90) -c  $(opt) scfopt.f90

scfq.o: scfq.f90 tablesmodule.f90 controlmodule.f90 constantsmodule.f90 indices.f90
	$(compile90) -c  $(opt) scfq.f90

BDATA4.o: BDATA4.f
	$(compile77) -c $(opt) BDATA4.f

REPPD.o: REPPD.f
	$(compile77) -c $(opt) REPPD.f

ROTD.o: ROTD.f
	$(compile77) -c $(opt) ROTD.f

ROTMAT.o: ROTMAT.f
	$(compile77) -c $(opt) ROTMAT.f

factrl.o: factrl.f90
	$(compile90) -c $(opt) factrl.f90

dparameters.o: dparameters.f90
	$(compile90) -c  $(opt)  dparameters.f90

rearrange.o: rearrange.f90
	$(compile90) -c $(opt) rearrange.f90

inighd.o: inighd.f90 tablesmodule.f90
	$(compile90) -c $(opt) inighd.f90

rsc.o: rsc.f90 tablesmodule.f90
	$(compile90) -c $(opt) rsc.f90

fbinom.o: fbinom.f90 tablesmodule.f90
	$(compile90) -c $(opt) fbinom.f90

wstore.o: wstore.f90 tablesmodule.f90
	$(compile90) -c $(opt) wstore.f90

ds.o: ds.f90 constantsmodule.f90
	$(compile90) -c $(opt) ds.f90

dp.o: dp.f90 constantsmodule.f90
	$(compile90) -c $(opt) dp.f90

dd.o: dd.f90 constantsmodule.f90
	$(compile90) -c $(opt) dd.f90

doverlap.o: doverlap.f90 tablesmodule.f90 constantsmodule.f90 sphericalmodule.f90
	$(compile90) -c $(opt) doverlap.f90

loadparam.o: loadparam.f90 
	$(compile90) -c $(opt) loadparam.f90

three.o: three.f90 tablesmodule.f90 constantsmodule.f90 indices.f90
	$(compile90) -c $(opt) three.f90

om1.o: om1.f90
	$(compile90) -c  $(opt)  om1.f90

dipole.o: dipole.f90 tablesmodule.f90 constantsmodule.f90 indices.f90 scratch.f90 controlmodule.f90
	$(compile90) -c $(opt) dipole.f90

nodummy.o: nodummy.f90 tablesmodule.f90 constantsmodule.f90
	$(compile90) -c $(opt) nodummy.f90

diff_three.o: diff_three.f90 tablesmodule.f90 constantsmodule.f90
	$(compile90) -c $(opt) diff_three.f90

three_core.o: three_core.f90 tablesmodule.f90 constantsmodule.f90
	$(compile90) -c $(opt) three_core.f90

cgfam.o: cgfam.f
	$(compile77) -c $(opt) cgfam.f

cgsearch.o: cgsearch.f
	$(compile77) -c $(opt) cgsearch.f

fcn.o: fcn.f
	$(compile77) -c $(opt) fcn.f

main.o: main.f
	$(compile77) -c $(opt) main.f

gradientuhf.o: gradientuhf.f90 constantsmodule.f90 indices.f90 controlmodule.f90 scratch.f90 tablesmodule.f90
	$(compile90) -c  $(opt) gradientuhf.f90

f2graduhf.o: f2graduhf.f90 indices.f90 constantsmodule.f90 tablesmodule.f90
	$(compile90) -c  $(opt) f2graduhf.f90

hfdft.o: hfdft.f90 constantsmodule.f90 tablesmodule.f90 scratch.f90
	$(compile90) -c  $(opt) hfdft.f90

b3.o: b3.f
	$(compile77) -c $(opt) b3.f

b88.o: b88.f
	$(compile77) -c $(opt) b88.f

b3lyp.o: b3lyp.f
	$(compile77) -c $(opt) b3lyp.f

fuzzy2.o: fuzzy2.f
	$(compile77) -c $(opt) fuzzy2.f

lyp.o: lyp.f
	$(compile77) -c $(opt) lyp.f

mkd.o: mkd.f90 indices.f90 constantsmodule.f90 tablesmodule.f90 scratch.f90
	$(compile90) -c $(opt) mkd.f90

vwn.o: vwn.f
	$(compile77) -c $(opt) vwn.f

pbe.o: pbe.f
	$(compile77) -c $(opt) pbe.f

hfexc.o: hfexc.f90 indices.f90 constantsmodule.f90 tablesmodule.f90 scratch.f90
	$(compile90) -c $(opt) hfexc.f90

hfdftuhf.o: hfdftuhf.f90 constantsmodule.f90 tablesmodule.f90 scratch.f90
	$(compile90) -c  $(opt) hfdftuhf.f90

mkduhf.o: mkduhf.f90 indices.f90 constantsmodule.f90 tablesmodule.f90 scratch.f90
	$(compile90) -c $(opt) mkduhf.f90

hfexcuhf.o: hfexcuhf.f90 indices.f90 constantsmodule.f90 tablesmodule.f90 scratch.f90
	$(compile90) -c $(opt) hfexcuhf.f90

am1_parameters.o: am1_parameters.f90
	$(compile90) -c  $(opt)  am1_parameters.f90

default.o: tablesmodule.f90 controlmodule.f90
	$(compile90) -c $(opt) default.f90

tql3.o: tql3.f
	$(compile77) -c   -fast tql3.f

tred3.o: tred3.f
	$(compile77) -c  -fast tred3.f

pack1.o: pack1.f90 indices.f90
	$(compile90) -c  $(opt)  pack1.f90

pack2.o: pack2.f90 indices.f90
	$(compile90) -c  $(opt)  pack2.f90

square.o: square.f90 indices.f90
	$(compile90) -c  $(opt)  square.f90

mkvector.o: mkvector.f90 indices.f90
	$(compile90) -c  $(opt)  mkvector.f90

twoe_sparkle.o: twoe_sparkle.f90 constantsmodule.f90 indices.f90 controlmodule.f90 scratch.f90 tablesmodule.f90
	$(compile90) -c  $(opt)  twoe_sparkle.f90

repelsparkle.o: repelsparkle.f90 constantsmodule.f90 tablesmodule.f90 controlmodule.f90
	$(compile90) -c  $(opt) repelsparkle.f90

hcoregradsp.o: hcoregradsp.f90 constantsmodule.f90 tablesmodule.f90 indices.f90
	$(compile90) -c  $(opt)  hcoregradsp.f90 

hessian.o: hessian.f90 constantsmodule.f90 tablesmodule.f90  scratch.f90 controlmodule.f90
	$(compile90) -c  $(opt)  hessian.f90 

initguess.o: initguess.f90 tablesmodule.f90 constantsmodule.f90 indices.f90 scratch.f90
	$(compile90) -c $(opt) initguess.f90

normalize.o: normalize.f90 constantsmodule.f90
	$(compile90) -c $(opt) normalize.f90

buildgto.o: buildgto.f90 constantsmodule.f90
	$(compile90) -c $(opt) buildgto.f90

ssinitial.o: ssinitial.f90 constantsmodule.f90
	$(compile90) -c $(opt) ssinitial.f90

psinitial.o: psinitial.f90 constantsmodule.f90
	$(compile90) -c $(opt) psinitial.f90

dsinitial.o: dsinitial.f90 constantsmodule.f90
	$(compile90) -c $(opt) dsinitial.f90

ppinitial.o: ppinitial.f90 constantsmodule.f90
	$(compile90) -c $(opt) ppinitial.f90

dpinitial.o: dpinitial.f90 constantsmodule.f90
	$(compile90) -c $(opt) dpinitial.f90

ddinitial.o: ddinitial.f90 constantsmodule.f90
	$(compile90) -c $(opt) ddinitial.f90

buildt.o: buildt.f90 tablesmodule.f90 constantsmodule.f90 sphericalmodule.f90
	$(compile90) -c $(opt) buildt.f90

invert.o: invert.f
	$(compile77) -c $(opt) invert.f

initguessold.o: initguessold.f90 tablesmodule.f90 constantsmodule.f90 indices.f90 scratch.f90
	$(compile90) -c $(opt) initguessold.f90



