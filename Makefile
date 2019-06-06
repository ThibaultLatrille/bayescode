# ==============================================================================================================
#  COMPILATION
# ==============================================================================================================
.PHONY: all # Requires: cmake 3.1.0 or better
all: _build
	@cd _build ; make --no-print-directory -j8

_build: CMakeLists.txt # default mode is release
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake ..

.PHONY: coverage
coverage:
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake -DCOVERAGE_MODE=ON ..
	@make --no-print-directory test

.PHONY: debug
debug:
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake -DDEBUG_MODE=ON ..
	@make --no-print-directory

.PHONY: release
release:
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake ..
	@make --no-print-directory

.PHONY: clean
clean:
	@rm -rf _build
	@rm -rf _build_coverage
	@rm -rf _test
	@rm -rf _aamutsel
	@rm -rf _dated

# ==============================================================================================================
#  CODE QUALITY
# ==============================================================================================================
.PHONY: format # Requires: clang-format
format:
	@clang-format -i `find -name *.*pp`

# ==============================================================================================================
#  TESTING
# ==============================================================================================================
POINTS=2

.PHONY: run-unit-tests
run-unit-tests: all
	@echo "\n\e[35m\e[1m== Tree test ==================================================================\e[0m"
	_build/tree_test
	@echo "\n\n\e[35m\e[1m== All sequential tests =======================================================\e[0m"
	_build/all_tests
	@echo "\n\n\e[35m\e[1m== MPI par test ===============================================================\e[0m"
	mpirun -np 3 _build/mpi_par_test

.PHONY: run-app-tests
run-app-tests: all
	@rm -rf _test
	@mkdir _test
	@echo "\n\e[35m\e[1m== Globom run ===============================================================\e[0m"
	@_build/globom -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u ${POINTS} _test/globom_gal4
	@echo "\n\e[35m\e[1m== Globom restart ===========================================================\e[0m"
	_build/globom _test/globom_gal4
	@echo "\n\e[35m\e[1m== Globom read ==============================================================\e[0m"
	_build/readglobom _test/globom_gal4
	@echo "\n\e[35m\e[1m== CodonM2a run =============================================================\e[0m"
	_build/codonm2a -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u ${POINTS} _test/codonM2a_gal4
	@echo "\n\e[35m\e[1m== CodonM2a restart =========================================================\e[0m"
	_build/codonm2a _test/codonM2a_gal4
	@echo "\n\e[35m\e[1m== CodonM2a read ============================================================\e[0m"
	_build/readcodonm2a _test/codonM2a_gal4
	@echo "\n\e[35m\e[1m== MutSel run ===============================================================\e[0m"
	_build/aamutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u ${POINTS} _test/aamutsel_gal4
	@echo "\n\e[35m\e[1m== MutSel restart ===========================================================\e[0m"
	_build/aamutsel _test/aamutsel_gal4
	@echo "\n\e[35m\e[1m== MutSel read ==============================================================\e[0m"
	_build/readaamutsel _test/aamutsel_gal4
	@echo "\n\e[35m\e[1m== MutSel read site-profiles ================================================\e[0m"
	_build/readaamutsel --ss _test/aamutsel_gal4
	@echo "\n\e[35m\e[1m== MutSel with polymorphism run =============================================\e[0m"
	_build/aamutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -p -u ${POINTS} _test/aamutsel_gal4_poly
	@echo "\n\e[35m\e[1m== MutSel Multiple omega run ================================================\e[0m"
	_build/mutselomega -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick --freeomega --omegancat 3 -u ${POINTS} _test/mutselomega_gal4
	@echo "\n\e[35m\e[1m== MutSel Multiple omega restart ============================================\e[0m"
	_build/mutselomega _test/mutselomega_gal4
	@echo "\n\e[35m\e[1m== MutSel Multiple omega read ===============================================\e[0m"
	_build/readmutselomega _test/mutselomega_gal4
	@echo "\n\e[35m\e[1m== Diffsel double sparse ====================================================\e[0m"
	@make --no-print-directory diffseldsparse
	@echo "\n\e[35m\e[1m== Dated MutSel run =========================================================\e[0m"
	_build/nodemutsel --ncat 3 -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u ${POINTS} _test/nodemutsel_gal4
	@echo "\n\e[35m\e[1m== Dated MutSel restart =====================================================\e[0m"
	_build/nodemutsel _test/nodemutsel_gal4
	@echo "\n\e[35m\e[1m== Dated MutSel read ========================================================\e[0m"
	_build/readnodemutsel --ss _test/nodemutsel_gal4

# @make --no-print-directory run-multigeneglobom-test
.PHONY: run-multigeneglobom-test
run-multigeneglobom-test: all
	@echo "\n\e[35m\e[1m== Multigene Single Omega ===================================================\e[0m"
	cd data/small_multigene && mpirun -np 2 ../../_build/multigeneglobom -t tree.nwk -a verysmall.list  -u ${POINTS} tmp

.PHONY: testpr
test:
	@make --no-print-directory run-unit-tests
	@make --no-print-directory run-app-tests

.PHONY: aamutsel
aamutsel: _build
	@cd _build ; make --no-print-directory -j8 aamutsel readaamutsel
	@rm -rf _aamutsel
	@mkdir _aamutsel
	_build/aamutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u 10 _aamutsel/gal4
	_build/aamutsel _aamutsel/gal4
	_build/readaamutsel _aamutsel/gal4
	_build/readaamutsel --ss _aamutsel/gal4
	_build/aamutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u 10 -p _aamutsel/gal4_poly
	_build/aamutsel _aamutsel/gal4_poly
	_build/readaamutsel _aamutsel/gal4_poly
	_build/readaamutsel --ss _aamutsel/gal4_poly

.PHONY: mutselomega
mutselomega: _build
	@cd _build ; make --no-print-directory -j8 mutselomega readmutselomega
	@rm -rf _mutselomega
	@mkdir _mutselomega
	_build/mutselomega -a data/samhd1/samhd1.ali -t data/samhd1/samhd1.tree --omegashift 0.0 --freeomega --omegancat 1 -u 20 _mutselomega/samhd1
	_build/mutselomega _mutselomega/samhd1
	_build/readmutselomega _mutselomega/samhd1
	_build/readmutselomega --ss _mutselomega/samhd1

.PHONY: dated
dated: _build
	@cd _build ; make --no-print-directory -j8 nodemutsel
	@rm -rf _dated
	@mkdir _dated
	_build/nodemutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick --ncat 3 -u ${POINTS} _dated/node_gal4
	_build/nodemutsel _dated/node_gal4
	_build/nodemutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick --ncat 3 -u ${POINTS} -p _dated/node_poly_gal4
	_build/nodemutsel _dated/node_poly_gal4

.PHONY: diffseldsparse
diffseldsparse: all
		@rm -f delme*.*
		_build/diffseldsparse -a data/besnard/cyp_small.phy -t data/besnard/cyp_coding.Chrysithr_root.nhx -e 1 -u 3 tmp
