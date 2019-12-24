// Main command to calculate the finite-population-adjusted SE
capture program drop fpopse
program define fpopse, eclass
	syntax varlist(min=2 numeric) [if] [aw], ///
		control(varlist min=1 numeric fv) rho(numlist min=1 max=1)
		
	marksample touse
	
	// Check the value of rho
	_fpopse_check_rho, rho(`rho')

	// Parse varlist
	tempname y treatvars
	local `y': word 1 of `varlist'
	local `treatvars' = subinstr("`varlist'", "``y'' ", "", 1)
	
	// Residualize treatment variables
	tempname treatvars_resid eta
	local `treatvars_resid' ""
	foreach var of varlist ``treatvars'' {
		tempvar `var'_resid
		qui reg `var' `control' [`weight'`exp'] if `touse'
		qui predict ``var'_resid' if e(sample), residual
		local `treatvars_resid' ``treatvars_resid'' ``var'_resid'
	}
	qui matrix accum `eta' = ``treatvars_resid'' if `touse', nocons

	// Estimate delta_ehw
	tempname b V N df_r delta_ehw_vars delta_ehw
	qui reg `varlist' `control' [`weight'`exp'] if `touse', robust
	foreach stat in b V N df_r {
		matrix ``stat'' = e(`stat')
	}
	
	tempvar e
	qui predict `e' if e(sample), residual
	
	local `delta_ehw_vars' ""
	foreach var of varlist ``treatvars_resid'' {
		tempvar `var'_resid_e
		qui gen ``var'_resid_e' = `var' * `e' if e(sample)
		local `delta_ehw_vars' ``delta_ehw_vars'' ``var'_resid_e'
	}
	qui matrix accum `delta_ehw' = ``delta_ehw_vars'' if e(sample), nocons
	
	// Estimation of delta_cond
	tempname delta_cond_vars delta_cond
	local `delta_cond_vars' ""
	foreach var of varlist ``delta_ehw_vars'' {
		tempvar `var'_resid
		qui reg `var' `control' [`weight'`exp'] if `touse'
		qui predict ``var'_resid' if e(sample), residual
		local `delta_cond_vars' ``delta_cond_vars'' ``var'_resid'
	}
	qui matrix accum `delta_cond' = ``delta_cond_vars'' if e(sample), nocons

	// Calculate upper bound of the variance
	tempname var_ehw var_ub
	matrix `var_ehw' = inv(`eta') * `delta_ehw' * inv(`eta')
	matrix `var_ub' = inv(`eta') * (`rho' * `delta_cond' + (1-`rho') * `delta_ehw') * inv(`eta')

	// Return estimation results
	tempname n_treatvars b_return V_return V_conv N_return df_r_return
	local `n_treatvars' = wordcount("``treatvars''")
	matrix `b_return' = `b'[1, 1..``n_treatvars'']
	matrix `V_return' = `var_ub'
	matrix `V_conv' = `V'[1..``n_treatvars'', 1..``n_treatvars'']
	
	local `N_return' = `N'[1, 1]
	local `df_r_return' = `df_r'[1, 1]

	matrix rownames `V_return' = `: rownames `V_conv''
	matrix colnames `V_return' = `: colnames `V_conv''

	ereturn post `b_return' `V_return', depname("``y''") obs(``N_return'') ///
		dof(``df_r_return'') esample(`touse')
	
	foreach mat in eta delta_ehw delta_cond V_conv {
		ereturn matrix `mat' = ``mat''
	}
	ereturn scalar rho = `rho'
	ereturn display
end

// Program to re-calculate SE for different rho
capture program drop fpopse_recalc
program define fpopse_recalc, rclass
	syntax , rho(numlist min=1 max=1)
	
	// Check the value of rho
	_fpopse_check_rho, rho(`rho')
	
	tempname eta delta_ehw delta_cond V
	foreach mat in eta delta_ehw delta_cond {
		matrix ``mat'' = e(`mat')
	}

	matrix `V' = inv(`eta') * (`rho' * `delta_cond' + (1-`rho') * `delta_ehw') * inv(`eta')
	matrix rownames `V' = `: rownames e(V)'
	matrix colnames `V' = `: colnames e(V)'

	display _skip(2)
	display "{bf:Variance-covariance matrix with rho = `rho'}"
	matlist `V'
	return matrix V = `V'
end

// Program to check if rho takes valid value
capture program drop _fpopse_check_rho
program define _fpopse_check_rho
	syntax , rho(numlist min=1 max=1)
	
	capture assert inrange(`rho', 0, 1)
	if _rc {
		display as error "Invalid value: rho should be between 0 and 1"
		exit 198
	}
end
