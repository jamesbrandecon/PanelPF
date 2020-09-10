# panelpf

## Description
This is a dynamic panel estimator for a production function in which
- Inputs are correlated with unobserved productivity
- Productivity evolves as a nonlinear (quadratic in current version) Markov process
- No input is 1:1 with productivity (i.e. traditional approaches to inverting out productivity are misspecified)

## Installation and Usage
For now, simply download `_panelpf.ado`. For output y, inputs x, syntax is simply `panelpf y x, gmm_options(string)` where `gmm_options` contains the set of instruments for x and other options to pass to Stata’s `gmm` command.

## Examples

- `cobbDouglasExample.do`
 This example simulates a model of production in which
- all inputs choices are correlated with productivity
- some inputs are chosen based on current productivity
- these relationships are not invertible, and
- productivity follows a quadratic Markov process

Then I estimate the production function according to my procedure, OLS, and IV regressions. For comparisons to the production function estimator by Ackerberg, Caves, and Frazer (2015), see my paper titled “Estimating Productivity and Markups Under Imperfect Competition”.
