function [Props] = CO2Props(T,rho)
%#
%% INTRODUCTION
%==========================================================================
% Purpose:
% Returns a_sturcture containing the thermodynamic properties for carbon
% dioxide (CO2) for a given density (kg/m^3) and temperature (K). Valid
% temperature and pressure ranges are:
% 216 K ≤ T ≤ 1100K 0 MPa ≤ P ≤ 800 MPa
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
% cv - Specific Heat at Constant_volume, kJ/(kg*K)
% cp - Specific Heat at Constant_pressure, kJ/(kg*K)
% h - Specific Enthalpy, kJ/kg
% c - Speed of Sound, m/s
% rho - Density, kg/m3
% state - -1 = Negative Input, 0 = Liquid, 1 = Saturated, 2 = Gas
% l - Liquid designator
% v - Vapor designator
%--------------------------------------------------------------------------
% Revision History:

% Written for a CO2 Blowdown model developed at Utah State University by
% Matthew Wilson
% 4130 Old Main Hill
% Logan, UT 84322-4130
% Recommented and checked by:
% Brian Solomon
%--------------------------------------------------------------------------
% Based upon the Helmholtz Energy based equations of state described by
% Span, R. and Wagner, W. in
% "A New Equation of State for Carbon Dioxide Covering the Fluid Region
% from the Triple-Point Temperature to 1100 K at_pressures up to 800 MPa"
% Journal of Physical and Chemical Reference Data
% Vol 25, No. 6, 1996. Pp 1509-1596
%--------------------------------------------------------------------------

%% CALCULATE THE SATURATION PROPERTIES
%==========================================================================
% CO2 Constants
%--------------------------------------------------------------------------
R = 0.1889241; % Gas constant, kJ/(kg*K)
Tt = 216.592; % Triple point temperature, K (Eq. 3.1)
Pt = 0.51795; % Triple point_pressure, MPa (Eq. 3.2)
Tc = 304.1282; % Critical temperature, K (Eq. 3.3)
Pc = 7.3773; % Critical pressure, MPa (Eq. 3.4)
rhoc = 467.6; % Critical density, kg/m3 (Eq. 3.5)
% d = rho/rhoc; % Reduced density, ∆ = phi/phi critical, phi - mass density
t = Tc/T; % Inverse reduced temperature, tau = T critical/T, T - temperature

% Coefficients
%--------------------------------------------------------------------------
% Melting pressure coefficients, Section 3.3
a_m = [1955.5390 2055.4593]';
% Sublimation pressure coefficients, Section 3.4

a_s = [-14.740846 2.4327015 -5.3061778]';
% Vapor pressure coefficients, Section 3.5
a_p = [-7.0602087 1.9391218 -1.6463597 -3.2995634]';
t_p = [1.0 1.5 2.0 4.0]';
% Saturated liquid density coeffiecients, Section 3.6
a_l = [1.9245108 -0.62385555 -0.32731127 0.39245142]';
t_l = [0.34 0.5 10/6 11/6]';
% Saturated vapor density coeffiecients, Section 3.7
a_v = [-1.7074879 -0.82274670 -4.6008549 -10.111178 -29.742252]';
t_v = [0.34 0.5 1.0 7/3 14/3]';

% Phase Property Calculations
%--------------------------------------------------------------------------
% P melt = Pt* (1 + a_m(1)*(T/Tt-1) + a_m(2)*(T/Tt-1)^2); % Melting Pressure, Eq. 3.10
% P sub = Pt*exp(Tt/T*(a_s(1)*(1-T/Tt) + ... % Sublimation Presssure, Eq. 3.12
% a_s(2)*(1-T/Tt)^1.9 + ...
% a_s(3)*(1-T/Tt)^2.9));
P_sat = Pc * exp(Tc/T*sum(a_p.*(1-T/Tc).^t_p)); % Vapor Pressure, Eq. 3.13
rho_l = rhoc * exp(sum(a_l.*(1-T/Tc).^t_l)); % Saturated Liquid Density, Eq. 3.14
rho_v = rhoc * exp(sum(a_v.*(1-T/Tc).^t_v)); % Saturated Vapor Density, Eq. 3.15

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
% GAS
%--------------------------------------------------------------------------
elseif rho < rho_v || imag(rho_l)
	[P s u cp cv h c] = Helmholtz(t,rho/rhoc);
	X = 1;
	s_v = s; u_v = u; cp_v = cp;
	cv_v = cv; h_v = h; c_v = c;
	s_l = NaN; u_l = NaN; cp_l = NaN;
	cv_l = NaN; h_l = NaN; c_l = NaN;
	state = 2;
	%P = P/1E6;
% LIQUID
%--------------------------------------------------------------------------
elseif rho > rho_l
	[P s u cp cv h c] = Helmholtz(t,rho/rhoc);
	X = 0;
	s_l = s; u_l = u; cp_l = cp;
	cv_l = cv; h_l = h; c_l = c;
	s_v = NaN; u_v = NaN; cp_v = NaN;
	cv_v = NaN; h_v = NaN; c_v = NaN;
	state = 0;
	%P = P/1E6;
% MELTING???
%--------------------------------------------------------------------------
%
% SATURATED
%--------------------------------------------------------------------------

else
	X = rho_v*(rho_l - rho)/(rho*(rho_l-rho_v));
	[P s_l u_l cp_l cv_l h_l c_l] = Helmholtz(t,rho_l/rhoc);
	[P s_v u_v cp_v cv_v h_v c_v] = Helmholtz(t,rho_v/rhoc);
	s = s_v*X + s_l*(1-X);
	u = u_v*X + u_l*(1-X);
	h = h_v*X + h_l*(1-X);
	cp = cp_v*X + cp_l*(1-X);
	cv = cv_v*X + cv_l*(1-X);

	c = NaN;
	P = P_sat;
	state = 1;
end

%% CREATE THE OUTPUT STRUCTURE
%==========================================================================
Props.P = P; % Pressure
Props.X = X; % Quality
Props.s = s; Props.s_l = s_l; Props.s_v = s_v; % Entropy
Props.u = u; Props.u_l = u_l; Props.u_v = u_v; % Internal Energy
Props.cv = cv; Props.cv_l = cv_l; Props.cv_v = cv_v; % Specific Heat at Constant_volume
Props.cp = cp; Props.cp_l = cp_l; Props.cp_v = cp_v; % Specific Heat at Constant_pressure
Props.h = h; Props.h_l = h_l; Props.h_v = h_v; % Enthalpy
Props.c = c; Props.c_l = c_l; Props.c_v = c_v; % Speed of Sound
Props.rho_l = rho_l;Props.rho_v = rho_v; % Density
Props.state = state; % State
%Props.gma_l = cp_l/cv_l;
%Props.gma_v = cp_v/cv_v;
end
