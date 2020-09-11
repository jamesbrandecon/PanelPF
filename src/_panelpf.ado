capture program drop panelpf
program define panelpf
syntax varlist [if],  [start(name) gmm_options(string) iv(varlist)]

if `"`gmm_options'"' == "" & `"`start'"' == "" {
	di as error "You must specify a matrix of starting values using the start() or gmm_options() options"
	exit
}

capture mat init = `start'

// Check if xtset
capture xtset
if _rc != 0 {
	di as error "Data must be xtset before running panelpf"
	exit
}

// Extract y and x from varlist
tokenize `varlist'
local y `1'
macro shift 
local x `*'

local v "`x'"
global count: word count `v'

// Set lags of x as default IVs 
if `"`iv'"' == "" {
	local iv = "`x'"
}

local nparam = $count + 2
local neq = $count + 1

// Give default options for stata's gmm command 
if `"`gmm_options'"' == "" {
	local gmm_options = "from(init) nparameters(`nparam') nequations(`neq') twostep instruments(1:`iv', nocons) winitial(identity) conv_maxiter(5000)"
}
di "Beginning two-step GMM estimation"

gmm _panelpf2 `if', y("`y'") x("`x'") `gmm_options'
end

capture program drop _panelpf2
program define _panelpf2
syntax varlist [if], at(name) y(varlist) x(varlist)
tempvar ly fy g lg fg
tempvar e1 e2 e3 omega lomega fomega
qui gen `g' = 0
qui gen `lg' = 0
qui gen `fg' = 0

// Number of x variables
local v "`x'"
global count: word count `v'

// Make lags and forwards,
gen `fy' = f.`y'
gen `ly' = l.`y'

foreach i of numlist 1(1)$count {
	tempvar lx_`i' fx_`i' x_`i'
	local tvar: word `i' of `x'
	local t "`tvar'"
	qui gen `lx_`i'' = l.`t'
	qui gen `x_`i'' = `t'
	qui gen `fx_`i'' = f.`tvar'
}


// Initial values for parameters
foreach i of numlist 1(1)$count {
	tempname beta`i'
	scalar `beta`i'' = `at'[1,`i']
	qui replace `fg' = `fg' + `beta`i''*`fx_`i''
	qui replace `g' = `g' + `beta`i''*`x_`i''
	qui replace `lg' = `lg' + `beta`i''*`lx_`i''
}
tempname rho1 rho2
local rstart $count
scalar `rho1' = `at'[1,`rstart'+1]
scalar `rho2' = `at'[1,`rstart'+2]

// Compute omega for each period (conditional on parameter values)
qui gen double `fomega' = `fy' - `fg' `if'
qui gen double `omega' = `y' - `g' `if'
qui gen double `lomega' = `ly' - `lg' `if'

qui gen double `e1' = `fomega' - `rho1'*`omega'- `rho2'*`omega'^2 `if'
qui sum `e1' `if'
replace `e1' = `e1' - `r(mean)' `if'

qui gen double `e2' = `e1'*`lomega' `if'
qui gen double `e3' = `e1'*`lomega'^2 `if' 

qui replace `1' = `e1' `if' // Productivity innovation t+1
qui replace `2' = `e2' `if' // Enforcing first-order markov process
qui replace `3' = `e3' `if' // Enforcing first-order markov process, quadratic term
end


capture program drop _panelpf
program define _panelpf
syntax varlist [if], at(name)
tempvar e1 e2 e3 omega lomega fomega
tempvar lv fv lk fk ly fy

// Initial values for parameters
tempname betav betak rho1 rho2
scalar `betav' = `at'[1,1]
scalar `betak' = `at'[1,2]

scalar `rho1' = `at'[1,3]
scalar `rho2' = `at'[1,4]

// Compute omega for each period (conditional on parameter values)
qui gen double `fomega' = `fy' - `betav'*`fv' - `betak'*`fk' `if'
qui gen double `omega' = y - `betav'*v - `betak'*k  `if'
qui gen double `lomega' = `ly' - `betav'*`lv' - `betak'*`lk' `if'

qui gen double `e1' = `fomega' - `rho1'*`omega'- `rho2'*`omega'^2 `if'
qui sum `e1' `if'
replace `e1' = `e1' - `r(mean)' `if'

qui gen double `e2' = `e1'*`lomega' `if'
qui gen double `e3' = `e1'*`lomega'^2 `if'


qui replace `1' = `e1' `if' // Productivity innovation t+1
qui replace `2' = `e2' `if' // Enforcing first-order markov process
qui replace `3' = `e3' `if' // Enforcing first-order markov process, quadratic term
end

capture program drop _panelpfTL
program define _panelpfTL
syntax varlist [if], at(name)
tempvar e1 e2 e3 omega lomega fomega
tempvar lv fv lk fk ly fy

gen `lv' = l.v
gen `fv' = f.v
gen `lk' = l.k
gen `fk' = f.k
gen `ly' = l.y
gen `fy' = f.y

// Initial values for parameters
tempname bv bk bvk bk2 bv2 rho1 rho2
scalar `bv' = `at'[1,1]
scalar `bk' = `at'[1,2]
scalar `bvk' = `at'[1,3]
scalar `bv2' = `at'[1,4]
scalar `bk2' = `at'[1,5]

scalar `rho1' = `at'[1,6]
scalar `rho2' = `at'[1,7]

// Compute omega for each period (conditional on parameter values)
qui gen double `fomega' = `fy' - `bv'*`fv' - `bk'*`fk' - `bvk'*`fv'*`fk' - `bv2'*`fv'*`fv' - `bk2'*`fk'*`fk' `if'
qui gen double `omega' = y - `bv'*v - `bk'*k  - `bvk'*v*k - `bv2'*v*v - `bk2'*k*k `if'
qui gen double `lomega' = `ly' - `bv'*`lv' - `bk'*`lk' - `bvk'*`lv'*`lk' - `bv2'*`lv'*`lv' - `bk2'*`lk'*`lk' `if'

qui gen double `e1' = `fomega' - `rho1'*`omega'- `rho2'*`omega'^2 `if'
qui sum `e1' `if'
replace `e1' = `e1' - `r(mean)' `if'

qui gen double `e2' = `e1'*`lomega' `if'
qui gen double `e3' = `e1'*`lomega'^2 `if'

qui replace `1' = `e1' `if' // Productivity innovation t+1
qui replace `2' = `e2' `if' // Enforcing first-order markov process
qui replace `3' = `e3' `if' // Enforcing first-order markov process, quadratic term
end
