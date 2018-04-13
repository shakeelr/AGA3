#AGA 3 Calculation
#revision 0.1 (minor bug fixes)
#Shakeel Rajwani
#27/04/2017

#! python
import math

#initialize global variables using US units
R = 10.7316 			#universal gas constant
M_r = 28.9625 			#molecular weight of air
N_c = 323.279			#unit conversion factor (orifice flow)
N_Ic = 0.000623582		#unit conversion for iteration flow factor
N3 = 27.7070  			#unit conversion factor
N4 = 1.0 				#unit conversion factor (discharge coeffecient)
N5 = 459.67				#unit conversion factor (absolute temperature)
T_r = 68.0	 			#reference temperature
Pi = 3.14159

def round_sigfigs(num, sig_figs):
    """
	Round to specified number of sigfigs.
	Python built in round is to decimal places, while AGA requires rounding to signficant digits
	"""
    if num != 0:
        return round(num, -int(math.floor(math.log10(abs(num))) - (sig_figs - 1)))
    else:
        return 0

def set_units(units):
	"""
	Set units and universal constants used in the AGA calculations
	"""
	global R
	global N_c
	global N_Ic
	global N3
	global N4
	global N5
	global T_r
	global Pi
	#returns units and sets global constants
	if units == 'US':
		l_units = "in" 			#Orifice and meter tube internal diameter
		T_units = "deg F" 		#Temperatures
		P_units = "psia" 		#Pressures
		dP_units = "in H2O" 	#Orifice differential pressure
		rho_units = "lbm/ft^3" 	#Fluid density
		mflow_units = "lbm/hr"	#mass flows
		vflow_units = "ft^3/hr"	#volume flows
		alpha_units = "in/in-F" #thermal expansion coefficents
		mu_units = "cP"			#viscosity
		R = 10.7316 			#universal gas constant
		M_r = 28.9625 			#molecular weight of air
		N_c = 323.279			#unit conversion factor (orifice flow)
		N_Ic = 0.000623582		#unit conversion for iteration flow factor
		N3 = 27.7070  			#unit conversion factor
		N4 = 1.0 				#unit conversion factor (discharge coeffecient)
		N5 = 459.67				#unit conversion factor (absolute temperature)
		T_r = 68.0	 			#reference temperature
		Pi = 3.14159
	elif units == 'IP':
		l_units = "ft"
		T_units = "deg F"
		P_units = "psia"
		dP_units = "in H2O"
		rho_units = "lbm/ft^3"
		mflow_units = "lbm/hr"
		vflow_units = "ft^3/hr"
		alpha_units = "ft/ft-F"
		mu_units = "lbm/ft-s"
		R = 10.7316
		M_r = 28.9625
		N_c = 46552.1
		N_Ic = 0.0773327
		N3 = 27.7070
		N4 = 0.08333333
		N5 = 459.67
		T_r = 68.0
		Pi = 3.14159
	elif units == 'SI':
		l_units = "m"
		T_units = "deg K"
		P_units = "Pa"
		dP_units = "Pa"
		rho_units = "kg/s"
		mflow_units = "kg/s"
		vflow_units = "m^3/s"
		alpha_units = "m/m-K"
		mu_units = "cP"
		R = 8314.51
		M_r = 28.9625
		N_c = 1.0
		N_Ic = 1.0
		N3 = 1.0
		N4 = 0.0254
		N5 = 0.0
		T_r = 293.15
		Pi = 3.14159
	else: #default to metric
		l_units = "mm"
		T_units = "deg C"
		P_units = "bar"
		dP_units = "millibar"
		rho_units = "kg/m^3"
		mflow_units = "kg/hr"
		vflow_units = "m^3/hr"
		alpha_units = "mm/mm-C"
		mu_units = "cP"
		R = 0.0831451
		M_r = 28.9625
		N_c = 0.0360000
		N_Ic = 0.100000
		N3 =  1000.00
		N4 = 25.4
		N5 = 273.15
		T_r = 20.0
		Pi = 3.14159
	return l_units, T_units, P_units, dP_units, rho_units, mflow_units, vflow_units, alpha_units, mu_units

def thermal_expansion(alpha, d_r, T_r, T_f):
	"""
	AGA3 Procedure 4.3.2.1/4.3.2.2
	Calculation of Orifice Plate Bore/Meter Tube Diameter from Reference/Measured Diameter
	"""
	#alpha = coefficient of thermal expansion
	#d_r = orifice/pipe diameter at measured/reference temperature
	#T_r = measured/reference temperature
	#T_f = flowing temperature
	#returns d = diameter of orifice plate bore/meter tube at flowing conditions
	d = d_r * (1 + alpha * (T_f - T_r))
	return round_sigfigs(d, 6)

def diameter_ratio(d, D):
	"""
	AGA3 Procedure 4.3.2.3
	Calculation of Flowing Diameter Ratio (beta) from Meter Tube and Orifice Bore Diameters
	"""
	#d = orifice plate bore diameter at flowing conditions
	#D = meter internal diameter at flowing conditions
	#returns beta = ratio of orifice plate bore diameter to meter tube internal diameter at flowing conditions
	beta = d / D
	return round_sigfigs(beta, 6)

def velocity_factor(beta):
	"""
	AGA3 Procedure 4.3.2.4
	Calculation of Velocity of Approach Factor, E_v
	"""
	#beta = ratio of orifice plate bore diameter to meter tube internal diameter at flowing conditions
	#returns E_v = velocity of approach factor
	E_v = 1 / ((1 - beta**4)**(1/2.0))
	return round_sigfigs(E_v, 6)
	
def discharge_constants(D, beta):
	"""
	AGA3 Procedure 4.3.2.5
	Calculation of Flange-Tapped Orifice Plate Coefficient of Discharge Constants
	"""
	#D = meter tube internal diameter at flowing conditions
	#beta = ratio of orifice plate bore diameter to meter tube internal diameter
	#returns C_d0, C_d1, C_d2, C_d3, C_d4 = orifice plate coefficents of discharge
	global N4 #unit conversion factor (discharge coeffecient)
	
	#Set parameter values A1 through A6 and S1 through S8 used in discharge equation
	A0, A1, A2, A3, A4, A5, A6 = 0.5961, 0.0291, -0.229, 0.003, 2.8, 0.000511, 0.021
	S1, S2, S3, S4, S5, S6, S7, S8 = 0.0049, 0.0433, 0.0712, -0.1145, -0.2300, -0.0116, -0.5200, -0.1400
	
	#Step 1.  Calculate the dimensionless upstream tap location, L1, and dimensionless downstream tap location, L2
	L1 = N4 / D
	L1 = round_sigfigs(L1, 6)
	L2 = N4 / D
	L2 = round_sigfigs(L2, 6)
	
	#Step 2.  Calculate the dimensionless downstream dam height, M2
	M2 = (2 * L2)/ (1 - beta)
	M2 = round_sigfigs(M2, 6)
	
	#Step 3.  Calculate upstream tap correction factor, T_u
	T_u = (S2 + S3 * 2.71828**(-8.5*L1) + S4 * 2.71828**(-6.0*L2)) * (beta**4 / (1 - beta**4))
	T_u = round_sigfigs(T_u, 6)
	
	#Step 4.  Calculate the downstream tap correction factor, T_d
	T_d = S6 * (M2 + S7 * M2**1.3) * beta**1.1
	T_d = round_sigfigs(T_d, 6)
	
	#Step 5.  Calculate small pipe correction factor, T_s
	if D > (A4 * N4):
		T_s = 0.0
	else:
		T_s = A3 * (1 - beta) * (A4 - D / N4)
	T_s = round_sigfigs(T_s, 6)
	
	#Step 6.  Calculate the orifice plate coefficient of discharge constants at Reynolds number of 4000
	C_d0 = A0 + A1 * beta**2 + A2 * beta**8 + T_u + T_d + T_s
	C_d1 = A5 * beta**0.7 * 250**0.7
	C_d2 = A6 * beta**4 * 250**0.35
	C_d3 = S1 * beta**4 * beta**0.8 * 4.75**0.8 * 250**0.35
	C_d4 = (S5 * T_u + S8 * T_d) * beta**0.8 * 4.75**0.8
	return round_sigfigs(C_d0, 6), round_sigfigs(C_d1, 6), round_sigfigs(C_d2, 6), round_sigfigs(C_d3,6), round_sigfigs(C_d4,6)

def upstream_pressure(P_f, dP):
	"""
	AGA3 Procedure 4.3.2.6
	Calculation of Upstream Flowing Fluid Pressure from Downstream Static Pressure
	"""
	#P_f = flowing pressure (downstream tap)
	#dP = orifice differential pressure
	#returns P_f = flowing pressure (upstream tap)
	global N3 #unit conversion factor
	P_f = (dP / N3) + P_f
	return round_sigfigs(P_f, 6)

def expansion_factor(beta, dP, P_f, k):
	"""
	AGA3 Procedure 4.3.2.7
	Calculation of Compressible Fluid Expansion Factor
	"""
	#beta = ratio of orifice plate bore diameter to meter tube internal diameter calculated at flowing conditions
	#dP = orifice differential pressure
	#P_f = flowing pressure
	#k = isentropic exponent
	#returns Y = expansion factor
	global N3 #unit conversion factor
	
	#Step 1.  Calculate the orifice differental pressure to flowing pressure ration, x
	x = dP / (N3 * P_f)
	x = round_sigfigs(x, 6)
	
	#Step 2.  Calculate expansion factor pressure constant, Y_p
	Y_p = (0.41 + 0.35 * beta**4) / k 
	Y_p = round_sigfigs(Y_p, 6)
	
	#Step 3.  Calculate the expansion factor
	Y = 1 - (Y_p * x)
	return round_sigfigs(Y, 6)

def iteration_flow_factor(d, D, dP, E_v, mu, rho_f, Y):
	"""
	AGA3 Procedure 4.3.2.8
	Calculation of Iteration Flow Factor
	"""
	#d = orifice plate bore diameter calcluated at flowing temperature
	#D = meter tube internal diameter calculated at flowing temperature
	#dP = orifice differential pressure
	#E_v = velocity of approach factor
	#mu = absolute viscosity of fluid flowing
	#rho_f = density of fluid at flowing conditions
	#Y = expansion factor
	#returns F_I = iteration flow factor
	global N_Ic #unit conversion for iteration flow factor
	#Step 1.  Calculate iteration flow factor intermediate values
	F_Ic = (4000 * N_Ic * D * mu) / (E_v * Y * d**2)
	F_Ic = round_sigfigs(F_Ic, 6)
	F_Ip = (2 * rho_f * dP)**(1/2.0)
	F_Ip = round_sigfigs(F_Ip, 7)
	
	#Step 2.  Test for limiting value of iteration flow factor and limit accordingly
	if F_Ic < 1000 * F_Ip:
		F_I = F_Ic / F_Ip
	else:
		F_I = 1000
	return round_sigfigs(F_I, 6)

def discharge_coefficient(C_d0, C_d1, C_d2, C_d3, C_d4, F_I):
	"""
	AGA3 Procedure 4.3.2.9
	Calculation of Flange-Tapped Orifice Plate Coefficient of Discharge
	"""
	#C_d0, C_d1, C_d2, C_d3, C_d4 = orifice plate coefficents of discharge
	#F_I = iteration flow factor
	#returns C_dFT, C_df = orifice plate coefficient of discharge, bounds flag
	
	#Set constants
	X_c = 1.142139337256165 #value of X where low Reynolds number switch occurs
	A = 4.343524261523267 #correlation constant for low Reynolds number
	B = 3.764387693320165 #correlation constant for low Reynolds number
	
	#Step 1.  Initialize C_dFT to a value at infinite Reynolds number
	C_dFT = C_d0
	
	#Steps 2, 3, 4, 5: 
	#2. Calculate X, the ratio of 4000 to the assumed Reynolds number
	#3. Calculate correlation value of C_dFT, F_c, at the assumed Reynolds number, X and the 
	#derivaitive of the correlation with respect to the assumed value of C_dFT, D_c, then
	#4. Calculate the amount to change the guess for C_dFT, dC_d
	#5. Iterate until the absolute value of dC_d is less than 0.000005
	dC_d = 1 	
	while abs(dC_d) >= 0.000005:
		#Step 2
		X = F_I / C_dFT
		#Step 3
		if X < X_c:
			F_c = C_d0 + (C_d1 * X**0.35 + C_d2 + C_d3 * X**0.8) * X**0.35 + C_d4 * X**0.8
			D_c = (0.7 * C_d1 * X**0.35 + 0.35 * C_d2 + 1.15 * C_d3 * X**0.8) * X**0.35 + 0.8 * C_d4 * X**0.8
		else:
			F_c = C_d0 + C_d1 * X**0.7 + (C_d2 + C_d3 * X**0.8) * (A - B / X) + C_d4 * X**0.8
			D_c = 0.7 * C_d1 * X**0.7 + (C_d2 + C_d3 * X**0.8) * B / X + 0.8 * C_d3 * (A - B / X) * X**0.8 + 0.8 * C_d4 * X**0.8
		#Step 4
		dC_d = (C_dFT - F_c) / (1 + D_c / C_dFT)
		C_dFT = C_dFT - dC_d
	
	#Step 6.  If the value of X is greater than 1, set bounds flag
	if X > 1:
		C_df = True
	else:
		C_df = False
	
	return round_sigfigs(C_dFT, 6), C_df
	
def mass_flow(C_dFT, d, dP, E_v, rho_f, Y):
	"""
	AGA3 Procdure 4.3.2.10
	Calculation of Mass Flow Rate
	"""
	#C_dFT = converged orifice plate coefficient of discharge
	#d = orifice plate bore diameter calcluated at flowing temperature
	#dP = orifice differential pressure
	#E_v = velocity of approach factor
	#rho_f = density of the fluid at flowing conditions
	#Y = expansion factor
	#returns q_m = mass flow rate
	global N_c #unit conversion factor (orifice flow)
	global Pi
	
	#Step 1.  Calculate the mass flow factor according to the formula:
	F_mass = (Pi / 4) * N_c * E_v * d**2
	F_mass = round_sigfigs(F_mass, 6)
	
	#Step 2.  Calculate mass flow rate
	q_m = F_mass * C_dFT * Y * (2 * rho_f * dP)**(1/2.0)
	return round_sigfigs(q_m, 6)
	
def actual_flow(C_dFT, d, dP, E_v, rho_f, Y):
	"""
	AGA3 Procedure 4.3.2.11
	Calculation of Volume Flow Rate at Flowing (Actual) Conditions
	"""
	#C_dFT = converged orifice plate coefficient of discharge
	#d = orifice plate bore diameter calcluated at flowing temperature
	#dP = orifice differential pressure
	#E_v = velocity of approach factor
	#rho_f = density of the fluid at flowing conditions
	#Y = expansion factor
	#returns q_v = mass flow rate
	global N_c #unit conversion factor (orifice flow)
	global Pi

	#Step 1.  Calculate the mass flow factor according to the formula:
	F_mass = (Pi / 4) * N_c * E_v * d**2
	F_mass = round_sigfigs(F_mass, 6)	
	
	#Step 2.  Calculate volume flow rate
	q_v = (F_mass * C_dFT * Y * (2 * rho_f * dP)**(1/2.0)) / rho_f
	return round_sigfigs(q_v, 6)

def base_flow(C_dFT, d, dP, E_v, rho_b, rho_f, Y):
	"""
	AGA3 Procedure 4.3.2.12
	Calculation of Volume Flow Rate at Base (Standard) Conditions
	"""
	#C_dFT = converged orifice plate coefficient of discharge
	#d = orifice plate bore diameter calcluated at flowing temperature
	#dP = orifice differential pressure
	#E_v = velocity of approach factor
	#rho_b = density of the fluid at base conditions
	#rho_f = density of the fluid at flowing conditions
	#Y = expansion factor
	#returns q_v = mass flow rate
	global N_c #unit conversion factor (orifice flow)
	global Pi

	#Step 1.  Calculate the mass flow factor according to the formula:
	F_mass = (Pi / 4) * N_c * E_v * d**2
	F_mass = round_sigfigs(F_mass, 6)
	
	#Step 2.  Calculate volume flow rate
	q_b = (F_mass * C_dFT * Y * (2 * rho_f * dP)**(1/2.0)) / rho_b
	return round_sigfigs(q_b, 6)
	
def flowing_density_ideal():
	"""
	AGA3 Procedure 4.3.3.1
	Calculation of Natural Gas Flowing Density Using Ideal Gas Relative Density (Specific Gravity)
	"""
	pass #need to implement AGA8 before this function can be implemented
	
def  base_density_ideal():
	"""
	AGA3 Procedure 4.3.3.2
	Calculation of Natural Gas Base Density Using Ideal Gas Relative Density (Specific Gravity)
	"""
	pass #need to implement AGA8 before this function can be implemented
	
def flowing_density_real():
	"""
	AGA3 Procedure 4.3.3.3
	Calculation of Natural Gas Flowing Density Using Real Gas Relative Density (Specific Gravity)
	"""
	pass #need to implement AGA8 before this function can be implemented
	
def base_density_real():
	"""
	AGA3 Procedure 4.3.3.4
	Calculation of Natural Gas Base Density Using Real Gas Relative Density (Specific Gravity)
	"""
	pass #need to implement AGA8 before this function can be implemented
	
def isentropic_exponent():
	"""
	AGA3 Procedure 4.3.3.5
	Natural Gas Isentropic Exponent
	"""
	#returns k = isentropic exponent for calculation of expansion factor
	return 1.3198
	
def viscosity():
	"""
	AGA3 Procedure 4.3.3.6
	Natural Gas Viscosity
	"""
	#returns mu = absolute viscosity of the flowing fluid for calculation of Reynolds number or iteration flow factor
	return 0.010268 #assuming metric units (cP)