// --------------------------------------------------------------
// Code to simulate panel data on firm production in which 
// (1) all inputs choices are correlated with productivity
// (2) some inputs are chosen based on current productivity 
// (3) these relationships are not invertible, and 
// (4) productivity follows a quadratic Markov process
//
// The program _panelpf.ado provides a new GMM estimator for the 
// parameters of this model. 
//
// Note: Production function is Cobb-Douglas here. Translog can be
// estimated using _panelpfTL
//
// Written by James Brand
// Last Updated: July 6, 2020
// --------------------------------------------------------------


// Initialize simulated dataset
local N = 1000 // # firms
local T = 10 // Simulating 10 periods, will only keep period 9 (latest period with lag and forward)

local betak = 0.2 // coefficient on "capital", chosen in period t-1
local betav = 0.8 // coefficient on "flexible input", chosen in period t

local rho_1 = 0.95 // Linear term in Markov process
local rho_2 = -0.025 // Quadratic term

qui do _panelpf.ado

// Loop over S simulation samples, storing estimates for each sample
local S = 100
foreach s of numlist 1(1)`S' {
clear
qui {
set obs `N'

gen id = mod(_n,`N') + 1 // Generate firm ID for xtset
expand `T' // Make panel of firms
bys id: gen time = _n

xtset id time

// Generate firm productivity
gen innovation = 0.3*rnormal()

gen omega = rnormal()
foreach t of numlist 2(1)`T' {
	replace omega  = `rho_1'*l.omega + `rho_2'*l.omega^2 + innovation if time==`t'
}

// Generate endogenous inputs which are not an invertible function of omega or l.omega
gen k = rnormal()
gen v = rnormal()

foreach t of numlist 2(1)`T' {
	replace k  = 0.6*l.k + 0.3*l.omega + 0.7*rnormal() if time==`t'
	replace v  = 0.6*l.v + 0.3*omega + 0.7*rnormal() if time==`t'
}

gen error = 0.5*rnormal()
gen y = `betak'*k + `betav'*v + omega + error

// Restrict sample, make lags and forwards
keep if time>=8

mat init = 0.5, 0.5, 1, 0
}
panelpf y k v, gmm_options(from(init) nparameters(4) nequations(3) twostep instruments(1:k v, nocons) winitial(identity) conv_maxiter(50))

mat b = e(b)
local beta1`s' = b[1,1]
local beta2`s' = b[1,2]
local beta3`s' = b[1,3]
local beta4`s' = b[1,4]

// OLS reg for comparison
qui regress y k v
mat b = e(b)
local ols1`s' = b[1,1]
local ols2`s' = b[1,2]

// IV reg for comparison
qui ivreg2 y (k v = l.k l.v)
mat b = e(b)
local iv1`s' = b[1,1]
local iv2`s' = b[1,2]
}

// Compare results to OLS
clear 
set obs `S'
gen bk_panel  = .
gen bv_panel  = .
gen rho1_panel  = .
gen rho2_panel  = .

gen bk_ols = .
gen bv_ols = .

gen bk_iv = .
gen bv_iv = .

foreach s of numlist 1(1)`S' { 
	replace bk_panel = `beta1`s'' in `s'
	replace bv_panel = `beta2`s'' in `s'
	replace rho1_panel = `beta3`s'' in `s'
	replace rho2_panel = `beta4`s'' in `s'
	
	replace bk_ols = `ols1`s'' in `s'
	replace bv_ols = `ols2`s'' in `s'
	
	replace bk_iv = `iv1`s'' in `s'
	replace bv_iv = `iv2`s'' in `s'
}
sum
