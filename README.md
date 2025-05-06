# Macro_teaching_material_StochasticGrowthModel_RBCmodel

this repository collects the codes I typically use in the Foundations of Macroeconomics (1st semester MA course at WU (Vienna University of Economics and Business))

growth.m: solves the log-linearized deterministic growth model with the method of undetermined coefficients and via matrix algebra (standard eigenvalue-eigenvector decomposition) and computes transition paths

stochasticgrowth.m: solves the log-linearized deterministic growth model with the method of undetermined coefficients and via matrix algebra (standard eigenvalue-eigenvector decomposition) and computes impulse responses to a total factor productivity shock

stochasticgrowth_Klein.m: maps the log-linearized system of difference equation of the stochastic growth model into the format expected by solab.m, which solves the model via a QZ (generalized eigenvalue-eigenvector decomposition). Uses solab.m, qzdiv.m, qzswitch.m

stochasticgrowth_loglin.mod: Dynare code to solve the already log-linearized stochastic growth model, gives solution in terms of variables in percentage deviations from steady state (hat variables)

stochasticgrowth_level.mod: Dynare code to obtain a solution to the stochastic growth model in terms of absolute deviations from steady state (level variables)

stochasticgrowth_loglin2.mod: Dynare code that codes up the nonlinear system of the  stochastic growth model and uses variable transformation (x_t = exp(log(x_t))) to obtain a solution to the stochastic growth model in terms of percentage deviations from steady state (log-variables)

stochasticgrowth_loglin3.mod: Dynare code that codes up the nonlinear system of the  stochastic growth model and uses variable transformation (x_t = x_ss*exp(x_hat_t))) to obtain a solution to the stochastic growth model in terms of percentage deviations from steady state (hat-variables)

rbc_levels.mod: Dynare code for a standard RBC model with log preferences, social planner problem, gives solution in terms of level variables (absolute deviations from stst)

rbc_decentralized_levels.mod: Dynare code for a standard RBC model with log preferences, decentralized economy, gives solution in terms of level variables (absolute deviations from stst)

rbc.mod: Dynare code for a standard RBC model with log preferences, social planner problem, gives solution in terms of percentage deviation from stst variables

rbc_decentralized.mod: Dynare code for a standard RBC model with log preferences, decentralized economy, gives solution in terms of percentage deviation from stst variables
