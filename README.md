# panelpf

## Description
This is a dynamic panel estimator for a production function in which
- Inputs are correlated with unobserved productivity
- Productivity evolves as a nonlinear (quadratic in current version) Markov process
- No input is 1-to-1 with productivity, unconditionally or conditional on other inputs (i.e. traditional approaches to inverting out productivity are misspecified)

## Installation and Usage
For now, download `_panelpf.ado`. For output y, inputs x, syntax is `panelpf y x, gmm_options(string)` where `gmm_options` contains the set of instruments for x and other options to pass to Stata’s `gmm` command.

Note: as described in the paper, the estimating equation for this command is in terms of f.y and f.x (i.e. everything is shifted forward a period). Shift IVs accordingly (see example file). 

## Examples

`cobbDouglasExample.do`: This example simulates a model of production in which
- all inputs choices are correlated with productivity
- some inputs are chosen based on current productivity
- these relationships are not invertible, and
- productivity follows a quadratic Markov process

Then I estimate the production function according to my procedure, OLS, and IV regressions. For comparisons to the production function estimator by Ackerberg, Caves, and Frazer (2015), see my paper titled “Estimating Productivity and Markups Under Imperfect Competition”.

## To-do
- Next/Very soon: Make lagged explanatory variables default IVs, set default gmm options to reduce user input
- Near future: Allow for higher polynomial orders in productivity Markov process 
