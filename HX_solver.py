#! /usr/bin/env python
"""
Python Code to evaluate on- and off-design performance of heat exchangers.

Function has two oeprating modes:
(1) Stand alone
This evaluates the heat exchanger performance and can be used to plot 
temperature traces inside the heat exchanger.

(2) imported
The function is can be called by Cycle.py to allow the quasi-steady evaluation 
of heat exchanger performance as part of Cycle off-design modelling.

Version: 1.0
Author: Ingo Jahn
Last Modified: 26/03/2017
"""

import numpy as np
import CoolProp.CoolProp as CP
import scipy as sci
from scipy import optimize
import matplotlib.pyplot as plt 
import sys as sys 
import os as os 
from getopt import getopt

###
###
class Fluid:
    def __init__(self):
        self.fluidH  = []
        self.fluidC  = []
        self.TH_in  = [] # (K)
        self.mdotH  = [] # (kg/s)
        self.PH_in  = [] # pressure (Pa)
        self.PH_out = [] # pressure (Pa)
        self.TC_in  = [] # (K)
        self.mdotC  = [] # (kg/s)
        self.PC_in  = [] # pressure (Pa)
        self.PC_out = [] # pressure (Pa)
        self.T_ext = [] # ambient surroundign temperature
        self.T0 = []
    ###
    def check(self,M):
        if not self.fluidH and not self.fluidC:
            raise MyError('Neither F.fluidH or F.fluidC specified')
        if self.fluidH and not self.fluidC:
            self.fluidC = self.fluidH
            print 'fluidC not specified and set to equal fluidH'
        if self.fluidC and not self.fluidH:
            self.fluidH = self.fluidC
            print 'fluidH not specified and set to equal fluidC'
        if not self.TH_in:
            raise MyError('F.TH_in not specified correctly')
        if not self.mdotH:
            raise MyError('F.mdot_H not specified')
        if not self.PH_in:
            raise MyError('F.PH_in not specified')
        if not self.PH_out:
            self.PH_out = self.PH_in
            print 'PH_out not specified and set to equal PH_in'
        if not self.TC_in:
            raise MyError('F.TC_in not specified')
        if not self.mdotC:
            self.mdotC = self.mdotH
            print 'mdotC not specified and set to equal mdotH'
        if self.mdotC < 0. or self.mdotH < 0.:
            raise MyError('Mass flow rates in both channels needs to be > 0.')
        if not self.PC_in:
            raise MyError('F.PC_in not specified')
        if not self.PC_out:
            self.PC_out = self.PC_in
            print 'PC_out not specified and set to equal PC_in'  
        #if M.external_loss is 1:
        #    if not self.T_ext:
        #        raise MyError('T_ext not specified') 
        # 
        if M.print_flag:
            print 'Check of Fluid input parameters: Passed' 

    ###
    def get_T0(self,M):
        """
        function to set initial/old conditions
        """    
        if len(self.T0) is 0:
            T0 = np.zeros((4*M.N_cell))


            # set TWH
            T0[M.N_cell:2*M.N_cell] = np.ones(M.N_cell) * 0.5*(self.TH_in+self.TC_in)
            # set TWC 
            T0[2*M.N_cell:3*M.N_cell] = np.ones(M.N_cell) * 0.5*(self.TH_in+self.TC_in) 
            if M.co_flow:
                # set TH
                T0[0:M.N_cell] = self.TH_in + np.arange(M.N_cell)/float(M.N_cell)*(self.TC_in-self.TH_in)*0.25
                # set TC
                T0[3*M.N_cell:4*M.N_cell] = self.TC_in - np.arange(M.N_cell)/float(M.N_cell)*(self.TC_in-self.TH_in)*0.25
            else:
                # set TH
                T0[0:M.N_cell] = self.TH_in + np.arange(M.N_cell)/float(M.N_cell)*(self.TC_in-self.TH_in)*0.5
                # set TC
                T0[3*M.N_cell:4*M.N_cell] = self.TC_in - (1.-np.arange(M.N_cell)/float(M.N_cell))*(self.TC_in-self.TH_in)*0.5
                #T0[4*M.N_cell:5*M.N_cell] = np.ones(M.N_cell) * 0.5* (self.PH_in + self.PH_out) 
            #T0[5*M.N_cell:6*M.N_cell] = np.ones(M.N_cell) * 0.5* (self.PC_in + self.PC_out)  
            return T0
        else:
            return self.T0                      
###
###
class Model:
    def __init__(self):
        self.optim ='root:hybr'
        self.N_cell = [] # number of cells
        self.co_flow = 0  # default is to analyse counterflow heat exchangers. Set to 1 for co flow
        self.flag_axial = [] # whether axial heat conduction is considered
        self.external_loss = [] # whether external heat loss is considered
        self.Nu_CorrelationH = [] # set correlation for heat transfer in H channel
        self.Nu_CorrelationC = [] # set correlation for heat transfer in C channel
        self.f_CorrelationH = [] # set correlation for friction factor in H channel
        self.f_CorrelationC = [] # set correlation for friction factor in C channel
        self.H_dP_poly = [] # polynominal coefficients for pressure drop in H channel
        self.C_dP_poly = [] # polynominal coefficients for pressure drop in H channel       
    ###
    def check(self):
        if not self.N_cell:
            raise MyError('M.N_cell not specified')   
        if not self.flag_axial:
            self.flag_axial = 0
            print 'M.flag_axial no specified defaulting to 0'   
        if not self.Nu_CorrelationH:
            raise MyError('M.Nu_CorrelationH not specified')  
        if not self.Nu_CorrelationC:
            self.Nu_CorrelationC = self.Nu_CorrelationH
            print 'Nu_CorrelationC not specified and set to equal Nu_CorrelationH' 
        if not isinstance(self.f_CorrelationH,int):
            raise MyError('M.f_CorrelationH not specified')  
        if not isinstance(self.f_CorrelationC,int):
            self.f_CorrelationC = self.f_CorrelationH
            print 'f_CorrelationC not specified and set to equal f_CorrelationH' 
        if not (self.co_flow == 0 or self.co_flow == 1):
            raise MyError('M.co_flow not defined correctly. Set to 0 or 1')  
        if len(self.H_dP_poly) > 0 and self.f_CorrelationH != 0:
            raise MyError('Pressure drop polynominal H_dP_poly can only be used in conjunction with f_CorrelationH = 0')    
        if len(self.C_dP_poly) > 0 and self.f_CorrelationC != 0:
            raise MyError('Pressure drop polynominal C_dP_poly can only be used in conjunction with f_CorrelationC = 0')                    
        # 
        if self.print_flag:
            print 'Check of Model input parameters: Passed'
###
###
    def set_poly(self,F): 
        # use polynominal to update outlet pressure
        if len(self.H_dP_poly) > 0:
            temp = []
            for a in range(len(self.H_dP_poly)):
                temp.append(self.H_dP_poly*F.mdotC**a)
            F.PC_out = sum(temp) 

        if len(self.C_dP_poly) > 0: # if polynominal is specified use this to calc PH_out
            temp = []
            for a in range(len(self.C_dP_poly)):
                temp.append(self.C_dP_poly*F.mdotC**a)
            F.PC_out = sum(temp) 
###
###
class Geometry:
    def __init__(self):
        self.HXtype = []
        self.k_wall = [] # thermal conductivity (W / mk)
        self.type_list = ['micro-channel','shell-tube']
    ###
    def check_initialise(self,M):
        if not self.HXtype:
            raise MyError('G.HXtype not specified') 
        if not any(self.HXtype in s for s in self.type_list):
            raise MyError('Specified type in G.HXtype not supported. Use:'+self.type_list) 

        # look at microchannel case
        if self.HXtype == 'micro-channel':
            self.micro_check()
            self.micro_init()
        # 
        if self.HXtype == 'shell-tube':
            self.shell_tube_check()
            self.shell_tube_init()
        #
        if M.print_flag:
            print 'Check of Geometry input parameterrs: Passed' 
    ###
    def micro_check(self):
        if not self.N_R:
            raise MyError('G.N_R not specified')  
        if not self.N_C:
            raise MyError('G.N_C not specified') 
        if not self.t_1:
            raise MyError('G.t_1 not specified')  
        if not self.t_2:
            raise MyError('G.t_2 not specified') 
        if not self.t_casing:
            raise MyError('G.t_casing not specified') 
        if not self.HX_L:
            raise MyError('G.HX_L not specified')  
        if not self.d_tube:
            raise MyError('G.d_tube not specified')
        if not self.k_wall:
            raise MyError('G.k_wall not specified') 
    ###
    def micro_init(self):
        self.Area = self.N_C * (self.d_tube + self.t_2) * (2*self.N_R-1) * self.HX_L 
        #self.Area = self.N_C * self.N_R * (self.d_tube * np.pi) * self.HX_L 

        self.AH = self.N_R*self.N_C * self.d_tube**2/4 *np.pi # total flow area (m2)
        self.Lc_H = self.d_tube # characteristic length (m)
        self.Dh_H = self.d_tube
        self.AC = self.N_R*self.N_C * self.d_tube**2/4 *np.pi # total flow area (m2) 
        self.Lc_C = self.d_tube # characteristic length (m)
        self.Dh_C = self.d_tube
    
        self.t_wall = self.t_1 # (m)
        
        L1 = self.N_C*self.d_tube  + (self.N_C-1)*self.t_2   + 2*self.t_casing
        L2 = 2*self.N_R*self.d_tube+ (2*self.N_R-1)*self.t_1 + 2*self.t_casing
        self.A_wall = L1*L2 - self.AH - self.AC

        self.L_wall = self.HX_L  
    ###
    def shell_tube_check(self):
        if not self.N_T:
            raise MyError('G.N_T not specified')  
        if not self.d:
            raise MyError('G.d not specified') 
        if not self.D:
            raise MyError('G.D not specified')  
        if not self.DD:
            raise MyError('G.DD not specified') 
        if not self.t_casing:
            raise MyError('G.t_casing not specified') 
        if not self.HX_L:
            raise MyError('G.HX_L not specified')  
        if not self.k_wall:
            raise MyError('G.k_wall not specified') 
    ###
    def shell_tube_init(self):
        self.Area = 0.5*(self.d +self.D) * np.pi * self.HX_L * self.N_T
        self.AH = self.N_T * self.d**2/4 *np.pi # total flow area (m2)
        self.Lc_H = self.d # characteristic length (m)
        self.Dh_H = self.d # hydraulic diameter 
        self.AC = self.DD**2/4*np.pi - self.N_T* self.D**2/4.*np.pi # total flow area (m2)
        self.AC = self.AC 
        self.Lc_C = self.D # characteristic length (m)
        Perimeter = self.DD * np.pi + self.N_T * self.D*np.pi
        self.Dh_C = 4.* self.AC / Perimeter # hydraulic diameter
    
        self.t_wall = (self.D-self.d)/2. # (m)
        
        self.A_wall = self.DD*np.pi*self.t_wall + self.N_T * (self.D**2 - self.d**2) / 4.*np.pi

        self.L_wall = self.HX_L  
##
##
def calc_Nu(P,Tm,Tp,Tw, mdot,A_total,L_c,fluid ,Correlation):
    """
    function to calculate local Nusselt number based on bulk flow properties
    Inputs:
    P - bulk pressure (Pa)
    Tm - bulk temperature to left (K)
    Tp - bulk temperature to right (K)
    Tw - temperature of wall (K)
    mdot - total mass flow rate (kg/s)
    A_total - total flow area (m**2)
    L_c - characteristic length (m)
    fluid - fluid type
    Correlation - select correlation to be used:
        1 - Yoon correlation for sCO2 in pipes for Tb > Tpc
        2 - S.M. Liao and T.S Zhaou correlation for micorchannels
            http://www.me.ust.hk/~mezhao/pdf/33.PDF 
            - horizontal pipes
        3 - S.M. Liao and T.S Zhaou correlation for micorchannels
            http://www.me.ust.hk/~mezhao/pdf/33.PDF 
            - vertical pipes, upwards flow
        4 - - S.M. Liao and T.S Zhaou correlation for micorchannels
            http://www.me.ust.hk/~mezhao/pdf/33.PDF 
            - vertical pipes, downwards flow

    Outputs:
    Nu - Nussel Number
    """
    Tb = 0.5*(Tm+Tp) # bulk temperature

    # calculate Prandl number
    Pr = CP.PropsSI('PRANDTL','P',P,'T', Tb ,fluid)

    # claculate Reynolds number        
    rho_b = CP.PropsSI('DMASS','P',P,'T', Tb ,fluid)
    U = abs( mdot / (rho_b * A_total))
    mu_b = CP.PropsSI('VISCOSITY','P',P,'T', Tb ,fluid)
    Re = rho_b * U * L_c / mu_b

    if Correlation is 1:
        # Yoon correlation for horizontal pipes
        Nu = 0.14 * Re**0.69 * Pr**0.66
    elif Correlation is 2:
        # Liao correlation for horizontal pipes
        rho_w = CP.PropsSI('DMASS','P',P,'T', Tw ,fluid)
        Gr = abs( 9.80665 * (rho_b-rho_w)*rho_b*L_c**3 / mu_b**2)
        Cp_b = CP.PropsSI('CPMASS','P',P,'T', Tb ,fluid)
        i_b = CP.PropsSI('UMASS','P',P,'T', Tb ,fluid)
        i_w = CP.PropsSI('UMASS','P',P,'T', Tw ,fluid)
        if Tw == Tb:
            Nu = 0.
        else:
            Cp_bar = (i_w-i_b) / (Tw-Tb)
            #print 'Gr: ', Gr, 'Cp_bar: ', Cp_bar
            Nu = 0.124 * Re**0.8 * Pr**0.4 * (Gr/Re**2)**0.203 * (rho_w/rho_b)**0.842 * (Cp_bar / Cp_b)**0.384 
    elif Correlation is 3:
        # Liao correlation for vertical pipes - upwards flow
        rho_w = CP.PropsSI('DMASS','P',P,'T', Tw ,fluid)
        rho_mid = CP.PropsSI('DMASS','P',P,'T', 0.5*(Tw+Tb) ,fluid)
        rho_m = 1/(Tw-Tb) * (Tw-Tb)/6. * (rho_b + 4*rho_mid + rho_w)  # integration using Simpsons rule
        Gr_m = abs( 9.80665 * (rho_b-rho_m)*rho_b*L_c**3 / mu_b**2)
        Cp_b = CP.PropsSI('CPMASS','P',P,'T', Tb ,fluid)
        i_b = CP.PropsSI('UMASS','P',P,'T', Tb ,fluid)
        i_w = CP.PropsSI('UMASS','P',P,'T', Tw ,fluid)
        if Tw == Tb:
            Nu = 0.
        else:
            Cp_bar = (i_w-i_b) / (Tw-Tb)
            #print 'Gr: ', Gr, 'Cp_bar: ', Cp_bar
            Nu = 0.354 * Re**0.8 * Pr**0.4 * (Gr_m/Re**2.7)**0.157 * (rho_w/rho_b)**1.297 * (Cp_bar / Cp_b)**0.296 
    elif Correlation is 4:
        # Liao correlation for vertical pipes - downwards flow
        rho_w = CP.PropsSI('DMASS','P',P,'T', Tw ,fluid)
        rho_mid = CP.PropsSI('DMASS','P',P,'T', 0.5*(Tw+Tb) ,fluid)
        rho_m = 1/(Tw-Tb) * (Tw-Tb)/6. * (rho_b + 4*rho_mid + rho_w)  # integration using Simpsons rule
        Gr_m = abs( 9.80665 * (rho_b-rho_m)*rho_b*L_c**3 / mu_b**2)
        Cp_b = CP.PropsSI('CPMASS','P',P,'T', Tb ,fluid)
        i_b = CP.PropsSI('UMASS','P',P,'T', Tb ,fluid)
        i_w = CP.PropsSI('UMASS','P',P,'T', Tw ,fluid)
        if Tw == Tb:
            Nu = 0.
        else:
            Cp_bar = (i_w-i_b) / (Tw-Tb)
            #print 'Gr: ', Gr, 'Cp_bar: ', Cp_bar
            Nu = 0.643 * Re**0.8 * Pr**0.4 * (Gr_m/Re**2.7)**0.186 * (rho_w/rho_b)**2.154 * (Cp_bar / Cp_b)**0.751 
    elif Correlation is 5:
        # For flow in circular pipes (from HLT) use for incompressible fluids (water/oil)
        # Dittus Boelter Equation
        if Tw > Tb: # heating of fluid
            n = 0.3
        else:       # cooling of fluid
            n = 0.4
        Nu = 0.023 * Re**0.8 * Pr**n
    elif Correlation is 6:
        # Correlation for shaell side as per paper by Xie et al.
        C = 0.16442; m=0.65582 
        e = np.exp( C + m * np.log(Re) )
        Nu = e * Pr**(1./3.)
    elif Correlation is 7:
        if Re == 0. or Re < 0.0001:
            # use Natural Convection relationship
            rho_w = CP.PropsSI('DMASS','P',P,'T', Tw ,fluid)
            beta = 1./Tb
            Gr = abs( 9.80665 * rho_w**2 * L_c**3 * (Tw - Tb) * beta / mu_b**2)            
            GrPr = Gr*Pr
            if GrPr <= 1e9:
                C = 0.53; n=1/4
            else:
                C = 0.126; n=1/3            
            Nu = C * GrPr **n
        if Re > 0.0001 and Re <= 0.004:
             C = 0.437; n = 0.0985       
        elif Re >0.004 and Re <= 0.09:
             C = 0.565; n = 0.136
        elif Re > 0.09 and Re <= 1.:
            C = 0.8; n = 0.280
        elif Re > 1.0 and Re <= 35.:
            C = 0.795; n= 0.384
        elif Re > 35. and Re <= 5000.:
            C = 0.583; n = 0.471
        elif Re > 5000. and Re <= 50000.:
            C = 0.148; n = 0.633
        elif Re > 50000. and Re <= 200000:
            C = 0.0208; n = 0.814
        else:
            raise MyError('Correlation outside range of valid Nusselt numbers.')

    else:
        raise MyError('Correlation option for Nusselt number calculation not implemented.')

    return Nu
##
##
def calc_friction(P, Tm, Tp, mdot, A_total, Dh, fluid, Correlation, epsilon = 0.):
    """
    function to calculate local friction factor based on local geometry and bulk flow properties
    Inputs:
    P - bulk pressure (Pa)
    Tm - bulk temperature to left (K)
    Tp - bulk temperature to right (K)
    mdot - total mass flow rate (kg/s)
    A_total - total flow area (m**2)
    Dh - Hydraulic Diameter (m)
    fluid - fluid type
    Correlation - select correlation to be used:
        1 - autmatically switches between laminar and turbulent flow
        2 - laminar flow - circular pipe
        3 - turbulent flow - rough pipes (Haaland's formula)
    epsilon - roughness height (m)
    
    Outputs:
    f - friction factor
    """
    Tb = 0.5*(Tm+Tp) # bulk temperature

    # claculate Reynolds number        
    rho_b = CP.PropsSI('DMASS','P',P,'T', Tb ,fluid)
    U = abs(mdot / (rho_b * A_total))
    mu_b = CP.PropsSI('VISCOSITY','P',P,'T', Tb ,fluid)
    Re = rho_b * U * Dh / mu_b
    if Correlation is 1:
        if Re < 2300:
            f = 64. / Re
        else:
            temp = np.log10( (epsilon/Dh / 3.7)*1.11 + (6.9/Re ))
            f = (-1.8 * temp )**-2.
    elif Correlation is 2:
        # laminar pipe flow
        f = 64. / Re
    elif Correlation is 3:
        # tubulent rough pipes
        temp = np.log10( (epsilon/Dh / 3.7)*1.11 + (6.9/Re ))
        f = (-1.8 * temp )**-2.
    
    else:
        raise MyError('Correlation option for friction factor calculation not implemented.')

    # print Gr / Re**2

    return f
###
###
def equations(T, M, G, F, flag=0):
    """
    function that evaluates the steady state energy balance for each cell
    Inputs:
    T - vector containg temperatures and pressures at various interface points
    M - class containign model parameters
    G - class containing geometry parameters
    F - class containing fluid boundary conditions
    flag - allows operation fo function to be altered
            0 - default for operation
            1 - output temperature and heat fluxes
            2 - returns pressure traces
    Outputs:
    error - vector containing miss-balance in energy equations for different lcoations
    """
    TH, TWH, TWC, TC= open_T(T,F.TH_in,F.TC_in,F.PH_in,F.PC_in,M,F)
                                        
    # print 'Temperature',TH,TWH,TC
    if flag is 1:
        Q1 = []; Q2 = []; Q3 = [] 
        Q4 = []; Q5 = []; Q6 = []
        Q7 = []; Q8 = []

    error = np.zeros(4*M.N_cell)

    # calculate Pressure distribution in both pipes based on current temperatures
    # pressure is calculated in flow direction. 
    PH = [F.PH_in]; PC=[F.PC_in]
    for i in range(M.N_cell): 
        if M.f_CorrelationH == 0: # apply linear pressure drop if no correlation is specificed
            PH.append(F.PH_in -  (F.PH_in - F.PH_out) / (M.N_cell-1.) * (i+1))
        elif M.f_CorrelationH == 1:
            # calculate pressure drop due to friction
            f= calc_friction(PH[i-1],TH[i],TH[i+1],F.mdotH,G.AH,G.Dh_H, F.fluidH , M.f_CorrelationH, G.epsilonH)
            rhoH = CP.PropsSI('D','P', PH[i-1], 'T', TH[i] ,F.fluidH)
            V = F.mdotH / rhoH / G.AH # calculate flow velocity 
            h_f = f * G.HX_L/(M.N_cell+1.) /G.Dh_H * V*V / (2.*9.81) # calculate friction head loss
            dP = h_f * rhoH * 9.81
            PH.append(PH[i] - dP)
        else:
            raise MyError('Corelation type not implemented')

        if M.f_CorrelationC == 0: # apply linear pressure drop if no correlation is specificed
            PC.append(F.PC_in - (F.PC_in - F.PC_out) / (M.N_cell-1.) * (i+1))
        elif M.f_CorrelationC == 1:
            if M.co_flow:
                f= calc_friction(PC[i-1],TC[i],TC[i+1],F.mdotC,G.AC,G.Dh_C, F.fluidC , M.f_CorrelationC, G.epsilonC)
                rhoC = CP.PropsSI('D','P', PC[i-1], 'T', TC[i] ,F.fluidC)
                V = F.mdotC / rhoC / G.AC # calculate flow velocity 
                h_f = f * G.HX_L/(M.N_cell+1.) /G.Dh_C * V*V / (2.*9.81) # calculate friction head loss
                dP = h_f * rhoC * 9.81
                PC.append(PC[i] - dP)
            else:
                # calculate pressure drop due to friction
                j = M.N_cell - i - 1  # 
                f= calc_friction(PC[i-1],TC[j],TC[j+1],F.mdotC,G.AC,G.Dh_C, F.fluidC , M.f_CorrelationC, G.epsilonC)
                rhoC = CP.PropsSI('D','P', PC[i-1], 'T', TC[j] ,F.fluidC)
                V = F.mdotC / rhoC / G.AC # calculate flow velocity 
                h_f = f * G.HX_L/(M.N_cell+1.) /G.Dh_C * V*V / (2.*9.81) # calculate friction head loss
                dP = h_f * rhoC * 9.81
                PC.append(PC[i] - dP)
        else:
            raise MyError('Corelation type not implemented')
            
   
    if not M.co_flow: # 
        # reverse direction of PC as flow is from i = -1 to i = 0
        PC = list(reversed(PC))

    #print 'PH', PH
    #print 'PC', PC

    # calculate energy balance for high pressure stream (H); low pressure stream (C) and dividing wall
    for i in range(M.N_cell):

        #print 'In HX_solver'
        #print 'TH', TH
        #print 'PH', PH
        kH = CP.PropsSI('CONDUCTIVITY','P', (0.5*(PH[i]+PH[i+1])), 'T', (0.5*(TH[i]+TH[i+1])) ,F.fluidH)
        kHm = CP.PropsSI('CONDUCTIVITY','P', PH[i],                'T', TH[i] ,F.fluidH)
        kHp = CP.PropsSI('CONDUCTIVITY','P', PH[+1],               'T', TH[i+1] ,F.fluidH)
        kC = CP.PropsSI('CONDUCTIVITY','P', (0.5*(PC[i]+PC[i+1])), 'T', (0.5*(TC[i]+TC[i+1])) ,F.fluidC)
        kCm = CP.PropsSI('CONDUCTIVITY','P', PC[i],                'T', TC[i] ,F.fluidC)
        kCp = CP.PropsSI('CONDUCTIVITY','P', PC[i+1],              'T', TC[i+1] ,F.fluidC)

        NuH = calc_Nu((0.5*(PH[i]+PH[i+1])),TH[i],TH[i+1],TWH[i], F.mdotH, G.AH, G.Lc_H, F.fluidH, M.Nu_CorrelationH)
        NuC = calc_Nu((0.5*(PC[i]+PC[i+1])),TC[i],TC[i+1],TWC[i], F.mdotC, G.AC, G.Lc_C, F.fluidC, M.Nu_CorrelationC)

        hH = NuH * kH / G.Lc_H 
        hC = NuC * kC / G.Lc_C

        # heat trasnfer in H channel         
        q1_conv = CP.PropsSI('HMASS','P',PH[i],'T',TH[i],F.fluidH) * F.mdotH
        if i == 0:
            q1_cond = - kHm * G.AH/ (G.L_wall/M.N_cell/2.) * ( 0.5*(TH[i] + TH[i+1]) - TH[i] )
        else:
            q1_cond = - kHm * G.AH/ (G.L_wall/M.N_cell)    * ( 0.5*(TH[i] + TH[i+1]) - 0.5*(TH[i] + TH[i-1]) )
        q1 = q1_conv + q1_cond
                
        q2_conv = CP.PropsSI('HMASS','P',PH[i+1],'T',TH[i+1],F.fluidH) * F.mdotH 
        if i == M.N_cell-1:
            q2_cond = - kHp * G.AH/ (G.L_wall/M.N_cell/2.) * (      TH[i+1]            - 0.5*(TH[i] + TH[i+1]) )
        else:
            q2_cond = - kHp * G.AH/ (G.L_wall/M.N_cell)    * ( 0.5*(TH[i+1] + TH[i+2]) - 0.5*(TH[i] + TH[i+1]) )
        q2 = q2_conv + q2_cond
        q3 = hH * G.Area/M.N_cell * (0.5*(TH[i]+TH[i+1]) - TWH[i])      

        # Heat Transfer in wall
        q7_temp = (G.k_wall * G.Area/M.N_cell)/G.t_wall * (TWH[i] - TWC[i]) 
        if M.flag_axial is 1:
            if i == 0:
                q7_p = - G.k_wall* G.A_wall/ (G.L_wall/M.N_cell) * ( 0.5*(TWH[i+1]+TWC[i+1]) - 0.5*(TWH[i]  +TWC[i]  ))
                q7_m = 0.
            elif i == M.N_cell-1:
                q7_p = 0.
                q7_m = - G.k_wall* G.A_wall/ (G.L_wall/M.N_cell) * ( 0.5*(TWH[i]  +TWC[i]  ) - 0.5*(TWH[i-1]+TWC[i-1]))
            else:
                q7_p = - G.k_wall* G.A_wall/ (G.L_wall/M.N_cell) * ( 0.5*(TWH[i+1]+TWC[i+1]) - 0.5*(TWH[i]  +TWC[i]  ))
                q7_m = - G.k_wall* G.A_wall/ (G.L_wall/M.N_cell) * ( 0.5*(TWH[i]  +TWC[i]  ) - 0.5*(TWH[i-1]+TWC[i-1]))
        else:
            q7_p = 0.
            q7_m = 0.
        q7h = q7_temp - 0.5 * (-q7_p +q7_m)
        q7c = q7_temp + 0.5 * (-q7_p +q7_m)

        q4 = hC * G.Area/M.N_cell * (TWC[i] - 0.5*(TC[i]+TC[i+1]))
        qC_cond = - kC * G.AC/ (G.L_wall/M.N_cell) * ( TC[i+1] - TC[i])

        # Heat Transfer in C Channel        
        q5_conv = CP.PropsSI('HMASS','P',PC[i],'T',TC[i],F.fluidC) * -F.mdotC
        if i == 0:
            q5_cond = - kCm * G.AC/ (G.L_wall/M.N_cell/2.) * ( 0.5*(TC[i] + TC[i+1]) - TC[i] )
        else:
            q5_cond = - kCm * G.AC/ (G.L_wall/M.N_cell)    * ( 0.5*(TC[i] + TC[i+1]) - 0.5*(TC[i] + TC[i-1]) )
        q5 = q5_conv + q5_cond
        q6_conv = CP.PropsSI('HMASS','P',PC[i+1],'T',TC[i+1],F.fluidC) * -F.mdotC
        if i == M.N_cell-1:
            q6_cond = - kCp * G.AC/ (G.L_wall/M.N_cell/2.) * (      TC[i+1]            - 0.5*(TC[i] + TC[i+1]) )
        else:
            q6_cond = - kCp * G.AC/ (G.L_wall/M.N_cell)    * ( 0.5*(TC[i+1] + TC[i+2]) - 0.5*(TC[i] + TC[i+1]) )
        q6 = q6_conv + q6_cond

        # calculate mis-match in energy fluxes
        #print q1, q2, q3
        error[i] = q1-q2-q3
        error[M.N_cell+i] = q3-q7h
        error[2*M.N_cell+i] = q4-q7c
        error[3*M.N_cell+i] = q4+q5-q6

        if flag is 1:
            Q1.append(q1); Q2.append(q2); Q3.append(q3)
            Q4.append(q4); Q5.append(q5); Q6.append(q6)
            Q7.append(q7h); Q8.append(q7c)

        #print q1,q2,q3,q4,q5,q6,q7
        #print 'Error', error
        #print i

    if flag is 0:  
        return error
    elif flag is 1:
        return error, T, Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8
    elif flag is 2:
        return PH, PC
    else:
        raise MyError('flag option not defined.')        
###
###
def open_T(T,TH_in,TC_in,PH_in,PC_in,M,F):
    """
    function to unpack the Temperature vector T into the 6 vectors
    TH, TWH, TWC, TC, PH, PC
    """
    N_cell = M.N_cell

    TH = np.zeros(N_cell+1)
    TWH = np.zeros(N_cell)
    TWC = np.zeros(N_cell)
    TC = np.zeros(N_cell+1)    
    TH[0] = TH_in
    TH[1:N_cell+2] = T[0:N_cell]
    TWH = T[N_cell:2*N_cell]
    TWC = T[2*N_cell:3*N_cell]
    if M.co_flow:
        TC[0] = TC_in
        TC[1:N_cell+2] = T[3*N_cell:4*N_cell] 
    else:
        TC[0:N_cell] = T[3*N_cell:4*N_cell] 
        TC[N_cell] =TC_in 



    return TH,TWH,TWC,TC #,PH,PC
###
###
def main(uoDict):
    """
    main function
    """
    # create string to collect warning messages
    warn_str = "\n"

    # main file to be executed 
    jobFileName = uoDict.get("--job", "test")

    # strip .py extension form jobName
    jobName = jobFileName.split('.')
    jobName = jobName[0]

    # create classes to store input data
    M = Model()
    F = Fluid()
    G = Geometry()

    # set print_flag (can be overwritten from jobfile)
    M.print_flag = 1
    if uoDict.has_key("--noprint"):
        M.print_flag = 0

    # Execute jobFile, this creates all the variables
    execfile(jobFileName,globals(),locals())

    if M.print_flag:
        print "Input data file read"
    
    # Check that required input data has been provided
    M.check()
    F.check(M)
    M.set_poly(F)
    G.check_initialise(M)

    # Initialise temperature vector
    T0 = F.get_T0(M)

    # change flow direction of cold channel if co_flow
    if M.co_flow:
        F.mdotC = -F.mdotC

    # print 'T0', T0
    
    # set up tuple of optional inputs for use by fsolve
    args = (M, G, F, 0)   

    if M.optim == 'fsolve':
        T , infodict, status, mesg = sci.optimize.fsolve(equations,T0,args=args, full_output=1)
    elif M.optim == 'root:hybr':
        sol = sci.optimize.root(equations,T0,args=args,method='hybr',options={'xtol':1.e-12})
        status = sol.status
        T = sol.x
        mesg = sol.message
    elif gdata.optim == 'root:lm':
        sol = sci.optimize.root(equations,A0,args=args,method='lm',options={'eps':1.e-3, 'xtol':1.e-12, 'ftol':1e-12})
        status = sol.status
        A = sol.x
        mesg = sol.message
    elif gdata.optim == 'root:Newton-CG':
        sol = sci.optimize.root(equations,A0,args=args,method='lm',options={'eps':1.e-3, 'xtol':1.e-12})
        status = sol.status
        A = sol.x
        mesg = sol.message
    elif gdata.optim == 'root:df-sane':
        sol = sci.optimize.root(equations,A0,args=args,method='df-sane',options={'ftol':1.e-12})
        status = sol.status
        A = sol.x
        mesg = sol.message
    else:
        raise MyError("gdata.optim = '' not set preoperly.")

    if status is not 1:
        print mesg    
        raise MyError('HX_solver.py: fsolve unable to converge.')

    TH, TWH, TWC, TC = open_T(T,F.TH_in,F.TC_in,F.PH_in,F.PC_in,M,F)
                                # open_T(T,F.TH_in,F.TC_in,M.N_cell)

    # create pressure traces for output    
    PH, PC= equations(T, M, G, F, 2) 
                               
    if M.print_flag:
        print "Plotting results"
        plot_HX(TH,TWH,TWC,TC,M.N_cell)
        error, T, Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8= equations(T, M, G, F,  1)
        plot_HXq(Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, M.N_cell, G.HX_L)
        plot_Hp(PH,PC,M.N_cell)

        #print Q1, '\n', Q2, '\n', Q3, '\n', Q4, '\n', Q5, '\n', Q6, '\n', Q7 , '\n', Q8           
      
        print "\n"
        print 'Temperatures:'
        print 'Hot-channel           ' , TH
        print 'Hot-channel wall temp ' , TWH
        print 'Cold-channel wall temp' , TWC
        print 'Cold-channel          ' , TC
        print '\n'
        print 'Pressures:'
        print 'Hot-channel  ' , PH
        print 'Cold-channel ' , PC
      
        print "\n \n"
        print "Power Transferred - (H) channel"
        print 'T_in    (K): %.2f ' %(TH[0])
        print 'T_out   (K): %.2f ' %(TH[-1])
        print 'Delta T (K): %.2f ' %(abs(TH[0]-TH[-1]))
        P_inH =  CP.PropsSI('HMASS','P',F.PH_in,'T',TH[0],F.fluidH) * F.mdotH   
        P_outH =  CP.PropsSI('HMASS','P',F.PH_out,'T',TH[-1],F.fluidH) * F.mdotH 
        rho_in = CP.PropsSI('D','P',F.PH_in,'T',TH[0],F.fluidH)
        rho_out = CP.PropsSI('D','P',F.PH_out,'T',TH[-1],F.fluidH)
        mu_in = CP.PropsSI('V','P',F.PH_in,'T',TH[0],F.fluidH) 
        mu_out = CP.PropsSI('V','P',F.PH_out,'T',TH[-1],F.fluidH)   
        print 'Power   (kW): %.2f ' %((P_inH - P_outH)/1e3)
        print 'Reynolds number (in) : %.2f ' %( rho_in *F.mdotH/G.AH/rho_in * G.Lc_H / mu_in ) 
        print 'Reynolds number (out): %.2f ' %( rho_out *F.mdotH/G.AH/rho_out * G.Lc_H / mu_out ) 
        print "\n"
        print "Power Transferred - (C) channel"
        print 'T_in    (K): %.2f ' %(TC[-1])
        print 'T_out   (K): %.2f ' %(TC[0])
        print 'Delta T (K): %.2f ' %(abs(TC[0]-TC[-1]))
        P_inC =  CP.PropsSI('HMASS','P',F.PC_in,'T',TC[0],F.fluidC) * F.mdotC   
        P_outC =  CP.PropsSI('HMASS','P',F.PC_out,'T',TC[-1],F.fluidC) * F.mdotC    
        rho_in = CP.PropsSI('D','P',F.PC_in,'T',TC[0],F.fluidC)
        rho_out = CP.PropsSI('D','P',F.PC_out,'T',TC[-1],F.fluidC)
        mu_in = CP.PropsSI('V','P',F.PC_in,'T',TC[0],F.fluidC) 
        mu_out = CP.PropsSI('V','P',F.PC_out,'T',TC[-1],F.fluidC)   
        print 'Power   (kW): %.2f ' %((P_inC - P_outC)/1e3)
        print 'Reynolds number (in) : %.2f ' %( rho_in *F.mdotC/G.AC/rho_in * G.Lc_C / mu_in ) 
        print 'Reynolds number (out): %.2f ' %( rho_out *F.mdotC/G.AC/rho_out * G.Lc_C / mu_out ) 
        print "\n \n"   

        print 'Heat Transfer Info:'
        DT_A = TH[0] - TC[0]; DT_B = TH[-1] - TC[-1]
        T_LM = (DT_B - DT_A) / np.log(DT_B/DT_A)
        #T_LM = abs(((TH[-1]-TC[0]) - (TH[0]-TC[-1])) / (np.log( (TH[-1]-TC[0])/(TH[0]-TC[-1]) ) ))
        print 'Delta T Log Mean (K): %.2f' %( T_LM )
        print 'HTC  (W /(m K): %.2f' %( abs(P_inC - P_outC) / G.Area / abs(T_LM) )

        print "\n \n" 
        plt.draw()


        plt.pause(1) # <-------
        print '\n \n'
        raw_input("<Hit Enter To Close Figures>")
        plt.close()

    PH_out = PH[-1]
    PC_out = PC[0]
    TH_out = TH[-1]
    TC_out = TC[0]

    return PH_out, TH_out, PC_out, TC_out, PH, TH, PC, TC, T0
###
###
def plot_HX(TH,TWH,TWC,TC,N_cell):
    fig = plt.figure()
    plt.plot(np.linspace(0,1.,num=N_cell+1),TH,'--',label="(H) Channel")
    plt.plot(np.linspace(0.5/N_cell,1.-0.5/N_cell,num=N_cell),TWH,'o--',label="Wall")
    plt.plot(np.linspace(0.5/N_cell,1.-0.5/N_cell,num=N_cell),TWC,'o-',label="Wall")
    plt.plot(np.linspace(0,1.,num=N_cell+1),TC,label = "(C) Channel")
    plt.ylabel('Temperature (K)')
    plt.xlabel('Position along Heat Exchanger (normalised)')
    plt.title('Temperature Distributions in Heat Exchanger')
    plt.legend(loc=2)
###
###
def plot_Hp(PH,PC,N_cell):
    fig, ax1  = plt.subplots()
    l1 = ax1.plot(np.linspace(0,1.,num=N_cell+1),PH,'--',label="Pressure (H) Channel")
    ax2 = ax1.twinx()
    l2 = ax2.plot(np.linspace(0,1.,num=N_cell+1),PC, label = "Pressure (C) Channel")
    ax1.set_ylabel('Pressure H channel (Pa)')
    ax2.set_ylabel('Pressure C channel (Pa)')
    plt.xlabel('Position along Heat Exchanger (normalised)')
    plt.title('Pressure Distributions in Heat Exchanger')
    lines = l1+l2
    labels = [l.get_label() for l in lines]
    plt.legend(lines,labels,loc=2)

###
###
def plot_HXq(Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, N_cell,Length):
    fig = plt.figure()
    plt.plot(np.linspace(0,1.,num=N_cell+1)[0:-1]            ,np.array(Q1)/1.e3,label="q1")
    plt.plot(np.linspace(0,1.,num=N_cell+1)[1:]              ,np.array(Q2)/1.e3,label="q2")
    plt.plot(np.linspace(0.5/N_cell,1.-0.5/N_cell,num=N_cell),np.array(Q3)/1.e3,label="q3")
    plt.plot(np.linspace(0.5/N_cell,1.-0.5/N_cell,num=N_cell),np.array(Q4)/1.e3,label="q4")
    plt.plot(np.linspace(0,1.,num=N_cell+1)[0:-1]            ,np.array(Q5)/1.e3,label="q5")
    plt.plot(np.linspace(0,1.,num=N_cell+1)[1:]               ,np.array(Q6)/1.e3,label="q6")
    plt.plot(np.linspace(0.5/N_cell,1.-0.5/N_cell,num=N_cell),np.array(Q7)/1.e3,label="q7")
    plt.plot(np.linspace(0.5/N_cell,1.-0.5/N_cell,num=N_cell),np.array(Q8)/1.e3,label="q8")
    plt.ylabel('Heat Flow (kW)')
    plt.xlabel('Position along Heat Exchanger (normalised)')
    plt.title('Energy Fluxes in Heat Exchanger')
    plt.legend(loc=2)
    ###
    ###
    fig = plt.figure()
    plt.plot(np.linspace(0.5/N_cell,1.-0.5/N_cell,num=N_cell),np.array(Q3)/1.e3/(Length/N_cell),label="q3")
    plt.plot(np.linspace(0.5/N_cell,1.-0.5/N_cell,num=N_cell),np.array(Q4)/1.e3/(Length/N_cell),label="q4")
    plt.plot(np.linspace(0.5/N_cell,1.-0.5/N_cell,num=N_cell),np.array(Q7)/1.e3/(Length/N_cell),label="q7")
    plt.plot(np.linspace(0.5/N_cell,1.-0.5/N_cell,num=N_cell),np.array(Q8)/1.e3/(Length/N_cell),label="q8")
    plt.ylabel('Flow per unit length (kW/m)')
    plt.xlabel('Position along Heat Exchanger (normalised)')
    plt.title('Energy Fluxes in Heat Exchanger')
    plt.legend(loc=4)
###
###
shortOptions = ""
longOptions = ["help", "job=", "noprint"]
###
###
def printUsage():
    print ""
    print "Usage: HX_solver.py [--help] [--job=<jobFileName>] [--noprint]"
    print "\n"
    print " --help      Dispay help"
    print "\n"
    print " --job=      Use this to specify the job file." 
    print "\n"
    print " --noprint   This suppresses on screen outputs."
    return
###
###
class MyError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
###
###
if __name__ == "__main__":
    userOptions = getopt(sys.argv[1:], shortOptions, longOptions)
    uoDict = dict(userOptions[0])
    
    if len(userOptions[0]) == 0 or uoDict.has_key("--help"):
        printUsage()
        sys.exit(1)

    # execute the code
    try:
        main(uoDict)
        print "\n \n"
        print "SUCESS."
        print "\n \n"

    except MyError as e:
        print "This run has gone bad."
        print e.value
        sys.exit(1)
    
