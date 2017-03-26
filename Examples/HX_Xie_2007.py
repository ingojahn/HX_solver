"""
Example case based on heat exchangers data from.
G.N. Xie, Q.W. Wang, M. Zeng, L.Q. Luo (2007), Heat transfer analysis for 
shell-and-tube heat exchangers with experimental data by artificial neural 
network approach, Applied Thermal Engineering 27 (2007) 1096-1104

Used for validation of the shell and tube heat exchanger model with oil and water. 
The following assumptiosn were made:
- use Therminol T66 as fluid for oil side
- working pressure on both sides is 1 bara (1.e5 Pa)
- to account for forward and backward passage the effective tube length has been doubled. 
- heat conduction in the casing has been ignored (t_casing = 0)

Author: Ingo Jahn
Last: Modified 23/03/2017
"""

#import CoolProp.CoolProp as CP 

# set geometry for heat exhanger - required settings depend on type
G.HXtype = 'shell-tube'
G.N_T = 176/2    # number of tubes (reduced as number given in paper is for both directions)
G.d = 8.e-3    # tube inner diameter
G.D = 10.e-3    # tube outer diameter
G.DD = 207.e-3 # shell inner diameter 
G.t_casing = 1.e-6 #
G.HX_L = 0.620    # length of HX (m)
G.k_wall = 30 # thermal conductivity (W / m /K)
G.epsilonH = 0. # roughness height for H channel
G.epsilonC = 0. # roughness height for C cahannel

# Boundary Conditions
water = 'water'
T_w = 29.3+273.15   # water inlet temperature
Re_w = 3094       # water inlet Reynolds number
mu_w = CP.PropsSI('V','T',T_w,'P',6.e5,water)    # viscosity
rho_w = CP.PropsSI('D','T',T_w,'P',6.e5,water)   # density
V_w = Re_w * mu_w / G.d / rho_w
mdot_w = rho_w*V_w*G.N_T * G.d**2 / 4*np.pi

oil = 'INCOMP::T66'
T_o = 59.8+273.15   # oil inlet temperature
Re_o = 1825          # oil inlet Reynolds number
mu_o = CP.PropsSI('V','T',T_o,'P',6.e5,oil)    # viscosity
rho_o = CP.PropsSI('D','T',T_o,'P',6.e5,oil)   # density
V_o_max = Re_o * mu_o / G.D / rho_o
#mdot_o = rho_o * V_o* ( G.DD**2/ 4*np.pi - G.N_T * G.D**2 / 4*np.pi )
mdot_o = rho_o * V_o_max* 50e-3**2 * np.pi


# set fluid conditions at heat exchanger inlet and outlet
F.fluidH = 'water' 
F.TH_in = T_w
F.PH_in = 6.e5
F.PH_out = 6.e5
F.mdotH = mdot_w
F.fluidC = oil
F.TC_in = T_o
F.mdotC = 5.
F.PC_in = 6.e5
F.PC_out = 6.e5
F.mdotC = mdot_o
F.T0 = [ ] 

print "Masslfows: mdot_oil=",F.mdotC, ';  mdot_water=',F.mdotH, '(kg/s)'

# Set modelling parameters
M.N_cell = 10# number of cells
M.flag_axial = 1
M.external_loss = 0
M.Nu_CorrelationH = 5 
M.Nu_CorrelationC = 6 
M.f_CorrelationH = 0 
M.co_flow = 0
