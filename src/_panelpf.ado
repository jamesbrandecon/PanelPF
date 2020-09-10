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
