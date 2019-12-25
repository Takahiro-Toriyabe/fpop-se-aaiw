clear
set more off
set obs 1000
set seed 1120002

// Generate data
gen z = rnormal()
gen y2 = -0.1 * runiformint(-1, 3) + 0.5 * z + rnormal()
gen y1 = 0.2 * runiformint(-1, 3) + 0.2 * z + rnormal()
gen y0 = rnormal()

gen d = rnormal(z) < 0
gen d1 = d * rnormal(z) < 0
gen d2 = d * (1 - d1)

gen sample = runiformint(1, 100) <= 80

gen y = d2 * y2 + d1 * y1 + (1- d) * y0

// Usual analysis
reg y d1 d2 z if sample, robust

// Finite-population-adjusted SE
fpopse y d1 d2 if sample, control(z) rho(0.8)
local se_d1 = _se[d1]
local se_d2 = _se[d2]

// If you want re-calculate the SE with different value of rho
fpopse_recalc, rho(0.9)

// Check if "noresidualize" option works
foreach var in d1 d2 {
	qui reg `var' z if sample
	qui predict `var'_resid, residual
}

fpopse y d1_resid d2_resid if sample, control(z) rho(0.8) noresidualize
assert abs(_se[d1_resid] - `se_d1') / `se_d1' < 0.1^8
assert abs(_se[d2_resid] - `se_d2') / `se_d2' < 0.1^8
