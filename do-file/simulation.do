clear
set more off
set matsize 11000
set seed 1120002

local rho = 0.01
local n = 100000
local psi = 2
local sigma =1

local nsim = 10000

// Generate population
set obs `n'

forvalues j = 1(1)10 {
	gen z`j' = rnormal()
}

gen xi = rnormal()
gen theta = rnormal(z1 * `psi', `sigma') 

// Simulation
forvalues j = 1(1)7 {
	matrix B`j' = J(`nsim', 1, .)
	matrix SE`j' = J(`nsim', 2, .)
}

forvalues r = 1(1)`nsim' {
	
	// Re-draw random variables
	capture drop flag_sample u y y_alt?
	gen flag_sample = runiform() < `rho'
	gen u = rnormal()
	
	// Update observed outcome
	gen y = u * theta + xi
	gen y_alt1 = u * rnormal(z1 * 0, `sigma') + xi
	gen y_alt2 = u * rnormal(z1 * `psi', 0) + xi
	gen y_alt3 = u * rnormal(z1 * 0, 0) + xi

	// Column 1
	qui fpopse y u if flag_sample, control(z1) rho(`rho')

	matrix B1[`r', 1] = _b[u]
	matrix V1 = e(V), e(V_conv)
	matrix SE1[`r', 1] = sqrt(V1[1, 1])
	matrix SE1[`r', 2] = sqrt(V1[1, 2])
	
	// Column 2
	qui fpopse y u if flag_sample, control(z1-z10) rho(`rho')
	
	matrix B2[`r', 1] = _b[u]
	matrix V2 = e(V), e(V_conv)
	matrix SE2[`r', 1] = sqrt(V2[1, 1])
	matrix SE2[`r', 2] = sqrt(V2[1, 2])
	
	// Column 3
	qui fpopse y u if flag_sample & _n <= `n' / 10, control(z1) rho(`rho')
	
	matrix B3[`r', 1] = _b[u]
	matrix V3 = e(V), e(V_conv)
	matrix SE3[`r', 1] = sqrt(V3[1, 1])
	matrix SE3[`r', 2] = sqrt(V3[1, 2])
	
	// Column 4
	qui fpopse y u if _n <= 1000, control(z1) rho(1)
	
	matrix B4[`r', 1] = _b[u]
	matrix V4 = e(V), e(V_conv)
	matrix SE4[`r', 1] = sqrt(V4[1, 1])
	matrix SE4[`r', 2] = sqrt(V4[1, 2])
	
	// Column 5
	qui fpopse y_alt1 u if flag_sample, control(z1) rho(`rho')

	matrix B5[`r', 1] = _b[u]
	matrix V5 = e(V), e(V_conv)
	matrix SE5[`r', 1] = sqrt(V5[1, 1])
	matrix SE5[`r', 2] = sqrt(V5[1, 2])
	
	// Column 6
	qui fpopse y_alt2 u if flag_sample, control(z1) rho(`rho')

	matrix B6[`r', 1] = _b[u]
	matrix V6 = e(V), e(V_conv)
	matrix SE6[`r', 1] = sqrt(V6[1, 1])
	matrix SE6[`r', 2] = sqrt(V6[1, 2])
	
	// Column 7
	qui fpopse y_alt3 u if flag_sample, control(z1) rho(`rho')

	matrix B7[`r', 1] = _b[u]
	matrix V7 = e(V), e(V_conv)
	matrix SE7[`r', 1] = sqrt(V7[1, 1])
	matrix SE7[`r', 2] = sqrt(V7[1, 2])
}

// Display results
clear
set obs `nsim'
forvalues j = 1(1)7 {
	svmat B`j'
	svmat SE`j'
	
	foreach i in 1 2 {
		qui gen b_u95`j'`i' = B`j' + 1.96 * SE`j'`i'
		qui gen b_l95`j'`i' = B`j' - 1.96 * SE`j'`i'
		qui gen coverage`j'`i' = inrange(0, b_l95`j'`i', b_u95`j'`i')
	}

	tabstat SE`j'1 coverage`j'1 SE`j'2 coverage`j'2, format(%04.3f)
}

