#! python
import aga3

#set units and univseral constants (required)
units = "IP"
l_units, T_units, P_units, dP_units, rho_units, mflow_units, vflow_units, alpha_units, mu_units = aga3.set_units(units)
#inputs (required)
T_f	= 0.00
P_f = 200.0
dP = 56.2110
d_orifice = 0.2220050
alpha_orifice = 9.25 * (10**-6)
D_pipe = 0.3355317
alpha_pipe = 6.2 * (10**-6)
#inputs (optional, set to default values)
downstream_tap = False
T_r = aga3.T_r
k = 1.3198
mu = 9.0850 * (10**-6)
#fluid density at base and flowing
#NOTE: fluid density should be calculated using AGA8, to be implemented in the future
rho_f = 2.0466
rho_b = 0.116198

#a. At T_f, calculate terms that depend only upon orifice geometry, d, D, beta, E_v, orifice discharge constants
d = aga3.thermal_expansion(alpha_orifice, d_orifice, T_r, T_f)
D = aga3.thermal_expansion(alpha_pipe, D_pipe, T_r, T_f)
beta = aga3.diameter_ratio(d, D)
E_v = aga3.velocity_factor(beta)
C_d0, C_d1, C_d2, C_d3, C_d4 = aga3.discharge_constants(D, beta)
print "d =	%r %s" % (d, l_units)
print "D =	%r %s" % (D, l_units)
print "beta =	%r" % beta
print "E_v = 	%r" % E_v
print "C_d0 =	%r" % C_d0
print "C_d1 = 	%r" % C_d1
print "C_d2 = 	%r" % C_d2
print "C_d3 = 	%r" % C_d3
print "C_d4 = 	%r" % C_d4

#b. Calculate flowing pressure
if downstream_tap:
	P_f = aga3.upstream_pressure(P_f, dP)
print "P_f = 	%r %s" % (P_f, P_units)

#c. Calculate required fluid properties at T_f, P_f, and other specified conditions
#NOTE: AGA8 not implemented yet so update densities manually
rho_f = rho_f
rho_b = rho_b
print "rho_f =	%r %s" % (rho_f, rho_units)
print "rho_b =	%r %s" % (rho_b, rho_units)

#d. Calculate the appropriate fluid expansion factor
if k > 0: #check if compressible fluid
	Y = aga3.expansion_factor(beta, dP, P_f, k)
else:
	Y = 1.0
print "Y = 	%r" % Y

#e. Calculate the iteration flow factor
F_I = aga3.iteration_flow_factor(d, D, dP, E_v, mu, rho_f, Y)
print "F_I = 	%r" % F_I

#f. Determine the converged value of C_dFT, set bounds flag if outside certainty statement in standard
C_dFT, C_df = aga3.discharge_coefficient(C_d0, C_d1, C_d2, C_d3, C_d4, F_I)
print "C_dFT =	%r" % C_dFT
print "C_df = 	%r" % C_df 

#g. Calculate final values of mass flow (q_m) and volume flows (q_v, q_b) 
q_m = aga3.mass_flow(C_dFT, d, dP, E_v, rho_f, Y)
q_v = aga3.actual_flow(C_dFT, d, dP, E_v, rho_f, Y)
q_b = aga3.base_flow(C_dFT, d, dP, E_v, rho_b, rho_f, Y)
print "q_m = 	%r %s" % (q_m, mflow_units)
print "q_v = 	%r %s" % (q_v, vflow_units)
print "q_b = 	%r %s" % (q_b, vflow_units)