"""
Example case based on the sCO2 heat exchanger study by:
Josh Van Meter (2006) Experimental Investigation of a Printed Circuit Heat 
Exchanger using Supercritical Carbon Dioxide and Water as Heat Transfer Media,
MsC Thesis, Kansas State Univeristy

Used for validation of the micro channel heat exchanger model and for 
verification for supercritical fluid simualtions.

Author: Ingo Jahn
Last: Modified 23/03/2017
""" 

# set fluid conditions at heat exchanger inlet and outlet
F.fluidH = 'CO2' 
F.TH_in = 273.15+88.
F.mdotH = 100./3600
F.PH_in = 8.e6
F.fluidC = 'water'
F.TC_in = 273.15+36.0
F.mdotC = 700./3600
F.PC_in = 1.e5
F.T_ext = 295 # External Temperature (K) optional

F.T0 = [ ] 

# set geometry for heat exhanger - required settings depend on type
G.HXtype = 'micro-channel'
G.N_R = 21    # number of rows in HX
G.N_C = 53    # number of columns in HX matrix
G.t_1 = 1.31e-3  #
G.t_2 = 0.48e-3 # 
G.t_casing = 0.5*(44.6e-3 + 32.5e-3) #
G.HX_L = 0.5* (1.230 + 1.378)    # length of HX (m)
G.d_tube = 1.506e-3 # tube diameter
G.k_wall = 16 # thermal conductivity (W / m /K)
G.epsilonH = 0. # roughness height for H channel
G.epsilonC = 0. # roughness height for C cahannel

# Set modelling parameters
M.N_cell = 40 # number of cells
M.flag_axial = 1
M.external_loss = 0
M.Nu_CorrelationH = 2 
M.Nu_CorrelationC = 5 
M.f_CorrelationH = 1 
M.co_flow = 0

