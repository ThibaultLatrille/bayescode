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
	@rm -rf _mutseldm5
	@rm -rf _mutselomega
	@rm -rf _dated
	@rm -rf _branchOmegaNe

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
	@echo "\n\e[35m\e[1m== Node Omega run ===========================================================\e[0m"
	_build/nodeomegamutsel -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u ${POINTS} _test/nodeomega_gal4
	@echo "\n\e[35m\e[1m== Node Omega restart =======================================================\e[0m"
	_build/nodeomegamutsel _test/nodeomega_gal4
	@echo "\n\e[35m\e[1m== Node MutSel run ==========================================================\e[0m"
	_build/nodemutsel --ncat 3 -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u ${POINTS} _test/nodemutsel_gal4
	@echo "\n\e[35m\e[1m== Node MutSel restart ======================================================\e[0m"
	_build/nodemutsel _test/nodemutsel_gal4
	@echo "\n\e[35m\e[1m== Node MutSel read =========================================================\e[0m"
	_build/readnodemutsel --ss _test/nodemutsel_gal4
	@echo "\n\e[35m\e[1m== branch Omega, Ne - Site MutSel run =======================================\e[0m"
	_build/branchOmegaNeSiteMutsel --ncat 3 -a data/polymorphism/gal4.ali -t data/polymorphism/gal4.newick -u ${POINTS} _test/branchomegamutsel_gal4
	@echo "\n\e[35m\e[1m== branch Omega, Ne - Site MutSel read ======================================\e[0m"
	_build/readbranchOmegaNeSiteMutsel --ss _test/branchomegamutsel_gal4

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
	_build/mutselomega -a data/bglobin/bglobin.phy -t data/bglobin/bglobin.tre --omegashift 0.0 --freeomega --omegancat 3 -u 30 --flatfitness _mutselomega/flat_bglobin
	_build/readmutselomega -b 10 _mutselomega/flat_bglobin
	_build/mutselomega _mutselomega/flat_bglobin
	_build/mutselomega -a data/bglobin/bglobin.phy -t data/bglobin/bglobin.tre --omegashift 0.0 --fitness_profiles data/bglobin/bglobin.prefs --deltaomegaarray data/bglobin/deltaomegaarray.csv -u 30 _mutselomega/clamped_bglobin
	_build/readmutselomega -b 10 _mutselomega/clamped_bglobin
	_build/mutselomega _mutselomega/clamped_bglobin
	_build/mutselomega -a data/bglobin/bglobin.phy -t data/bglobin/bglobin.tre --omegashift 0.0 --deltaomegaarray data/bglobin/deltaomegaarray.csv --ncat 10 -u 30 _mutselomega/clamped_bglobin
	_build/readmutselomega -b 10 _mutselomega/clamped_bglobin
	_build/mutselomega _mutselomega/clamped_bglobin
	_build/mutselomega -a data/bglobin/bglobin.phy -t data/bglobin/bglobin.tre --omegashift 0.0 --freeomega --omegancat 3 --ncat 30 -u 30 _mutselomega/bglobin
	_build/readmutselomega -b 10 _mutselomega/bglobin
	_build/mutselomega _mutselomega/bglobin

.PHONY: DM5
DM5: _build
	@cd _build ; make --no-print-directory -j8 mutseldm5 readmutseldm5
	@rm -rf _mutseldm5
	@mkdir _mutseldm5
	_build/mutseldm5 -a data/bglobin/bglobin.phy -t data/bglobin/bglobin.tre  --omegashift 0.0 --freeomega --omegancat 10 --flatfitness --fixp0 --p0 0.0 -u 30 --hypermean_threshold 0.0 --hyperinvshape_threshold 10.0 _mutseldm5/flat_bglobin
	_build/readmutseldm5 _mutseldm5/flat_bglobin
	_build/mutseldm5 _mutseldm5/flat_bglobin
	_build/readmutseldm5 --ss _mutseldm5/flat_bglobin
	_build/mutseldm5 -a data/bglobin/bglobin.phy -t data/bglobin/bglobin.tre  --omegashift 1.0 --freeomega --omegancat 10 --ncat 30 -u 30 _mutseldm5/bglobin
	_build/readmutseldm5 _mutseldm5/bglobin
	_build/mutseldm5 _mutseldm5/bglobin

.PHONY: branchOmegaNe
branchOmegaNe: _build
	@cd _build ; make --no-print-directory -j8 branchOmegaNeSiteMutsel readbranchOmegaNeSiteMutsel
	@rm -rf _branchOmegaNe
	@mkdir _branchOmegaNe
	@echo "\n\e[35m\e[1m== branch Omega, Ne - Site MutSel run =======================================\e[0m"
	_build/branchOmegaNeSiteMutsel --ncat 3 -a data/bglobin/bglobin.phy -t data/bglobin/bglobin.PopSize.nhx --node_popsize_tag "PopulationSize" -u 20 _branchOmegaNe/branchomegamutsel_bglobin
	@echo "\n\e[35m\e[1m== branch Omega, Ne - Site MutSel read ======================================\e[0m"
	_build/readbranchOmegaNeSiteMutsel --newick _branchOmegaNe/branchomegamutsel_bglobin
	@echo "\n\e[35m\e[1m== branch Omega, Ne - Site MutSel read ======================================\e[0m"
	_build/readbranchOmegaNeSiteMutsel --ss _branchOmegaNe/branchomegamutsel_bglobin

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
