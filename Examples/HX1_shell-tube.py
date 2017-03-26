"""
Template input file for HX_solver.py
""" 

# set fluid conditions at heat exchanger inlet and outlet
F.fluidH = 'CO2' 
F.TH_in = 756.170449065
F.mdotH = 5.
F.PH_in = 7.68e6
F.PH_out = 7.68e6
F.TC_in = 324.91481002
F.mdotC = 5.
F.PC_in = 12.e6
F.PC_out = 12.e6
F.T_ext = 295 # External Temperature (K) optional

F.T0 = [ ] 

# set geometry for heat exhanger - required settings depend on type
G.HXtype = 'shell-tube'
G.N_T = 176    # number of tubes
G.d = 8.e-3    # tube inner diameter
G.D = 10.e-3    # tube outer diameter
G.DD = 207.e-3 # shell inner diameter 
G.t_casing = 1.e-6 #
G.HX_L = 2*0.620    # length of HX (m)
G.k_wall = 30 # thermal conductivity (W / m /K)
G.epsilonH = 0. # roughness height for H channel
G.epsilonC = 0. # roughness height for C cahannel

# Set modelling parameters
M.N_cell = 20 # number of cells
M.flag_axial = 1
M.external_loss = 0
M.Nu_CorrelationH = 5 
M.Nu_CorrelationC = 6 
M.f_CorrelationH = 1 
M.co_flow = 1
