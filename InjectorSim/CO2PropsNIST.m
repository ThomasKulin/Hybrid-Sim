function [Props] = CO2PropsNIST(T,rho)
%#
%% INTRODUCTION
%==========================================================================
% Purpose:
% Returns a sturcture containing the thermodynamic properties for saturated

% carbon dioxide (CO2) for a given density (kg/m^3) and temperature (K).
% Valid temperature and pressure ranges are:
% 216.59 K ≤ T ≤ 303.69K 0 MPa ≤ P ≤ 800 MPa
%--------------------------------------------------------------------------
% The properties in the structure are:
%--------------------------------------------------------------------------
% P X s u cv cp h c
% rho_l s_l u_l cv_l cp_l h_l c_l
% rho_v s_v u_v cv_v cp_v h_v c_v
%--------------------------------------------------------------------------
% Inputs:
% T - Temperature, K
% rho - Density, kg/m3
% Outputs:
% P - Pressure, MPa
% X - Quality, (Vapor Mass/Total Fluid Mass)
% s - Specific Entropy, kJ/(kg*K)
% u - Specific Internal Energy, kJ/kg
% cv - Specific Heat at Constant Volume, kJ/(kg*K)
% cp - Specific Heat at Constant Pressure, kJ/(kg*K)
% h - Specific Enthalpy, kJ/kg
% c - Speed of Sound, m/s
% rho - Density, kg/m3
% state - -1 = Negative Input, 0 = Liquid, 1 = Saturated, 2 = Gas
% l - Liquid designator
% v - Vapor designator
%--------------------------------------------------------------------------
% Revision History:
% Written for a CO2 Blowdown model developed at Utah State University by
% Brian Solomon
%--------------------------------------------------------------------------
% Based upon the Helmholtz Energy based equations of state described by
% Span, R. and Wagner, W. in
% "A New Equation of State for Carbon Dioxide Covering the Fluid Region
% from the Triple-Point Temperature to 1100 K at Pressures up to 800 MPa"

% Journal of Physical and Chemical Reference Data
% Vol 25, No. 6, 1996. Pp 1509-1596
%--------------------------------------------------------------------------

%% CALCULATE THE SATURATION PROPERTIES
%==========================================================================
% CO2 Constants
%--------------------------------------------------------------------------
R = 0.1889241; % Gas constant, kJ/(kg*K)
Tt = 216.592; % Triple point temperature, K (Eq. 3.1)
Pt = 0.51795; % Triple point pressure, MPa (Eq. 3.2)
Tc = 304.1282; % Critical temperature, K (Eq. 3.3)
Pc = 7.3773; % Critical pressure, MPa (Eq. 3.4)
rhoc = 467.6; % Critical density, kg/m3 (Eq. 3.5)
%{
%% Load the CO2 properties created by the NIST Webbook
%==========================================================================
[labels,T_NIST,y] = readColData('SatCO2Props.txt',25); %doesnt work use
importdata as follows:
data = importdata('SatCO2Props.txt')
y = data.data
T_NIST = y(:,1);
P_NIST = y(:,2); % Pressure (MPa)
rho_l_NIST = y(:,3); % Density (l, kg/m3)
rho_v_NIST = y(:,15); % Density (v, kg/m3)
u_l_NIST = y(:,5); % Internal Energy (l, kJ/kg)
u_v_NIST = y(:,17); % Internal Energy (v, kJ/kg)
h_l_NIST = y(:,6); % Enthalpy (l, kJ/kg)
h_v_NIST = y(:,18); % Enthalpy (v, kJ/kg)
s_l_NIST = y(:,7); % Entropy (l, J/g*K)
s_v_NIST = y(:,19); % Entropy (v, J/g*K)
cv_l_NIST = y(:,8); % Cv (l, J/g*K)
cv_v_NIST = y(:,20); % Cv (v, J/g*K)
cp_l_NIST = y(:,9); % Cp (l, J/g*K)
cp_v_NIST = y(:,21); % Cp (v, J/g*K)
c_l_NIST = y(:,10); % Sound Speed (l, m/s)
c_v_NIST = y(:,22); % Sound Speed (v, m/s)
%}
SatCO2Props = open('SatCO2Props.mat');

T_NIST = SatCO2Props.T_NIST; % Temperature (K)
P_NIST = SatCO2Props.P_NIST; % Pressure (MPa)
rho_l_NIST = SatCO2Props.rho_l_NIST; % Density (l, kg/m3)
rho_v_NIST = SatCO2Props.rho_v_NIST; % Density (v, kg/m3)
u_l_NIST = SatCO2Props.u_l_NIST; % Internal Energy (l, kJ/kg)
u_v_NIST = SatCO2Props.u_v_NIST; % Internal Energy (v, kJ/kg)
h_l_NIST = SatCO2Props.h_l_NIST; % Enthalpy (l, kJ/kg)
h_v_NIST = SatCO2Props.h_v_NIST; % Enthalpy (v, kJ/kg)
s_l_NIST = SatCO2Props.s_l_NIST; % Entropy (l, J/g*K)
s_v_NIST = SatCO2Props.s_v_NIST; % Entropy (v, J/g*K)
cv_l_NIST = SatCO2Props.cv_l_NIST; % Cv (l, J/g*K)
cv_v_NIST = SatCO2Props.cv_v_NIST; % Cv (v, J/g*K)
cp_l_NIST = SatCO2Props.cp_l_NIST; % Cp (l, J/g*K)
cp_v_NIST = SatCO2Props.cp_v_NIST; % Cp (v, J/g*K)
c_l_NIST = SatCO2Props.c_l_NIST; % Sound Speed (l, m/s)
c_v_NIST = SatCO2Props.c_v_NIST; % Sound Speed (v, m/s)
%% Interpolate to find the liquid and vapor properties
%==========================================================================
rho_l = interp1(T_NIST,rho_l_NIST,T,'linear','extrap');
rho_v = interp1(T_NIST,rho_v_NIST,T,'linear','extrap');
P = interp1(T_NIST,P_NIST,T,'linear','extrap');
s_l = interp1(T_NIST,s_l_NIST,T,'linear','extrap');
s_v = interp1(T_NIST,s_v_NIST,T,'linear','extrap');
u_l = interp1(T_NIST,u_l_NIST,T,'linear','extrap');
u_v = interp1(T_NIST,u_v_NIST,T,'linear','extrap');
cp_l = interp1(T_NIST,cp_l_NIST,T,'linear','extrap');
cp_v = interp1(T_NIST,cp_v_NIST,T,'linear','extrap');
cv_l = interp1(T_NIST,cv_l_NIST,T,'linear','extrap');
cv_v = interp1(T_NIST,cv_v_NIST,T,'linear','extrap');
h_l = interp1(T_NIST,h_l_NIST,T,'linear','extrap');
h_v = interp1(T_NIST,h_v_NIST,T,'linear','extrap');
c_l = interp1(T_NIST,c_l_NIST,T,'linear','extrap');
c_v = interp1(T_NIST,c_v_NIST,T,'linear','extrap');

%% CALCULATE THE PROPERTIES AT THE GIVEN TEMPERATURE AND DENSITY

%==========================================================================
% Account for negative density or temperature
%--------------------------------------------------------------------------
if rho < 0 || T < 0
X = NaN;
P = NaN;
s_v = NaN; u_v = NaN; cp_v = NaN;
cv_v = NaN; h_v = NaN; c_v = NaN;
s_l = NaN; u_l = NaN; cp_l = NaN;
cv_l = NaN; h_l = NaN; c_l = NaN;
s = s_v*X + s_l*(1-X);
u = u_v*X + u_l*(1-X);
h = h_v*X + h_l*(1-X);
cp = cp_v*X + cp_l*(1-X);
cv = cv_v*X + cv_l*(1-X);
c = NaN;
state = -1;
% SATURATED
%--------------------------------------------------------------------------
else
X = (rho_v/rho)*((rho_l-rho)/(rho_l-rho_v));
s = s_v*X + s_l*(1-X);
u = u_v*X + u_l*(1-X);
h = h_v*X + h_l*(1-X);
cp = cp_v*X + cp_l*(1-X);
cv = cv_v*X + cv_l*(1-X);
c = c_v*X + c_l*(1-X);
state = 1;
end

%% CREATE THE OUTPUT STRUCTURE
%==========================================================================
Props.P = P; % Pressure
Props.X = X; % Quality
Props.s = s; Props.s_l = s_l; Props.s_v = s_v; % Entropy

Props.u = u; Props.u_l = u_l; Props.u_v = u_v; % Internal Energy
Props.cv = cv; Props.cv_l = cv_l; Props.cv_v = cv_v; % Specific Heat at Constant Volume
Props.cp = cp; Props.cp_l = cp_l; Props.cp_v = cp_v; % Specific Heat at Constant Pressure
Props.h = h; Props.h_l = h_l; Props.h_v = h_v; % Enthalpy
Props.c = c; Props.c_l = c_l; Props.c_v = c_v; % Speed of Sound
Props.rho_l = rho_l;Props.rho_v = rho_v; % Density
Props.state = state; % State
%Props.gma l = cp_l/cv_l;
%Props.gma v = cp_v/cv_v;
end