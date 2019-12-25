# fpop-se-aaiw

Calculate finite population adjusted SE (Stata program) based on Abadie et. al. (2019) "Sampling-based vs. Design-based Uncertainty in Regression Analysis"

## Main command

```Stata
fpopse depvar treatvars [if] [aw], rho(#) [control(varlist)] [noresidualize]
```

- depvar: dependent variable
- treatvars: treatment variable(s) generating design-based uncertainty
- rho: sample-population ratio (in [0, 1])
- control: control variable(s) treated as non-random variables (If not specified, a constant term is used)
- noresidualize: specify it if the treatment variables are already residualized

### Return values

- `e(V_conv)`: conventional (heteroskedasticity-robust) variance matrix
- `e(gamma)`: Outer product of treatment variables (residualized by the controls)
- `e(delta_ehw)`: $\delta_{ehw}$ in the paper (multiplied by the sample size)
- `e(delta_z)`: $\delta_{z}$ in the paper (multiplied by the sample size)

## After estimation...

You can re-calculate the variance matrix with a different value of rho by using

```Stata
fpopse_recalc, rho(#)
```

which returns the corresponding variance matrix in `r(V)`