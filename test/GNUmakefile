
# to get a list of failed tests, filter the output of gmake with:
#    sed -n '/FAILED/ {s/\(.*\):.*/out\.\1/;p;}'

#SHELL=/bin/bsh
pwd:=$(shell pwd)

GENBAS:=${pwd}/GENBAS
ECPDATA:=${pwd}/ECPDATA
AUTODIR=TESTDIR
HANDDIR=WORKDIR
JUNK=tmp tmp.* *.tmp out.* rec.* core

GETRECS=0
recs  = TOTENERG COORD FORCECON
recs += E_AVG E_CLOSED E_CORE E_OPEN
recs += NUCDIPOL NUCENEG NUCREP
recs += PARENERG REFDIPOL
recs += SCFENEG TOTENER2 VACENEG
recs += E0 E1 E2 KSSCFENG KSTOTELE S2SCF

zmats:=$(wildcard zmat.*)
scripts:=$(wildcard script.*)

SKIPOUT=0
ifneq (${SKIPOUT},0)
   outs:=$(wildcard out.*)
   skips:=$(patsubst out.%, %, ${outs})
   zmats:=$(filter-out ${skips}, ${zmats})
   scripts:=$(filter-out ${skips}, ${scripts})
endif

workdir:=${HANDDIR}
auto=0
ifeq (${GETRECS},0)
   define runzmat
      rm -rf ${workdir} && mkdir ${workdir}
      cp ${pwd}/$@ ${workdir}/ZMAT && ln -s ${GENBAS} ${workdir}/GENBAS
      out=${pwd}/out.$@; i=0; \
         printf "%s:\t" $@; \
         cd ${workdir} && xaces2 >$${out} 2>&1 || i=1; \
         if test $${i} -eq 0; \
         then sed '1,/TEST\.DAT/d' ZMAT >test.dat \
               && xa2proc test test.dat >>$${out} 2>&1 \
               && echo PASSED \
               || echo "FAILED in test module"; \
         else echo "FAILED in xaces2"; \
         fi
      test ${auto} -ne 0 && rm -rf ${workdir} || true
   endef
   define runscript
      rm -rf ${workdir} && mkdir ${workdir}
      chmod u+x $@
      out=${pwd}/out.$@; i=0; \
         printf "%s:\t" $@; \
         cd ${workdir} && GENBAS=${GENBAS} ECPDATA=${ECPDATA} USERDIR=${pwd} \
            ${pwd}/$@ >$${out} 2>&1 || i=1; \
         if test $${i} -eq 0; \
         then sed '1,/TEST\.DAT/d' ${pwd}/$@ >test.dat \
               && xa2proc test test.dat >>$${out} 2>&1 \
               && echo PASSED \
               || echo "FAILED in test module"; \
         else echo "FAILED in script"; \
         fi
      test ${auto} -ne 0 && rm -rf ${workdir} || true
   endef
else
   define runzmat
      rm -rf ${workdir} && mkdir ${workdir}
      cp ${pwd}/$@ ${workdir}/ZMAT && ln -s ${GENBAS} ${workdir}/GENBAS
      out=${pwd}/out.$@; i=0; recf=${pwd}/rec.$@; rm -f $${recf}; \
         printf "%s records:\t" $@; \
         cd ${workdir} && xaces2 > $${out} || i=1; \
         test $${i} -eq 0 && for rec in ${recs}; \
         do idim=`xa2proc jasum \
                   | egrep "[^A-Z0-9]$${rec} *[^A-Z0-9]" \
                   | cut -b 60-80`; \
            test -n "$$idim" && idim=$$((idim/2)) || idim=0; \
            if test $$idim -ne 0; \
            then echo d $$rec >> $$recf; \
                 xa2proc jareq d $$rec $$idim \
                  | egrep -v '(\@|^ *$$)' | tail -$$idim >> $$recf || i=2; \
            fi; \
         done; \
         (test $${i} -eq 0 && echo DONE) || \
         (test $${i} -eq 1 && echo FAILED in xaces2) || \
         (test $${i} -eq 2 && echo FAILED in record module)
      test ${auto} -ne 0 && rm -rf ${workdir} || true
   endef
   define runscript
      rm -rf ${workdir} && mkdir ${workdir}
      chmod u+x $@
      out=${pwd}/out.$@; i=0; recf=${pwd}/rec.$@; rm -f $${recf}; \
         printf "%s records:\t" $@; \
         cd ${workdir} && GENBAS=${GENBAS} ECPDATA=${ECPDATA} USERDIR=${pwd} \
            ${pwd}/$@ > $${out} || i=1; \
         test $${i} -eq 0 && for rec in ${recs}; \
         do idim=`xa2proc jasum \
                   | egrep "[^A-Z0-9]$${rec} *[^A-Z0-9]" \
                   | cut -b 60-80`; \
            test -n "$$idim" && idim=$$((idim/2)) || idim=0; \
            if test $$idim -ne 0; \
            then echo d $$rec >> $$recf; \
                 xa2proc jareq d $$rec $$idim \
                  | egrep -v '(\@|^ *$$)' | tail -$$idim >> $$recf || i=2; \
            fi; \
         done; \
         (test $${i} -eq 0 && echo DONE) || \
         (test $${i} -eq 1 && echo FAILED in script) || \
         (test $${i} -eq 2 && echo FAILED in record module)
      test ${auto} -ne 0 && rm -rf ${workdir} || true
   endef
endif

################################################################################

SUFFIXES :=
.SUFFIXES:
.PHONY: all which FORCE clean

all : workdir:=${AUTODIR}
all : auto=1
all : which ${zmats} ${scripts}
	@echo; echo Finished running all tests.

FORCE: ;

which : FORCE
	@echo; echo ACESII binary: $$(which xaces2); echo

${zmats} : FORCE
	@$(runzmat)

${scripts} : FORCE
	@$(runscript)

clean : FORCE
	rm -rf ${AUTODIR} ${HANDDIR} ${JUNK}
	@echo; ls -F; echo

