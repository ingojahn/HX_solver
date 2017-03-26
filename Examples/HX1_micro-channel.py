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
G.HXtype = 'micro-channel'
G.N_R = 100    # number of rows in HX
G.N_C = 400    # number of columns in HX matrix
G.t_1 = 2e-3  #
G.t_2 = 0.5e-3 # 
G.t_casing = 5e-3 #
G.HX_L = 1.    # length of HX (m)
G.d_tube = 1.5e-3 # tube diameter
G.k_wall = 16 # thermal conductivity (W / m /K)
G.epsilonH = 0. # roughness height for H channel
G.epsilonC = 0. # roughness height for C cahannel

# Set modelling parameters
M.N_cell = 10 # number of cells
M.flag_axial = 1
M.external_loss = 0
M.Nu_CorrelationH = 2 
M.Nu_CorrelationC = 2 
M.f_CorrelationH = 1 
M.co_flow = 0

