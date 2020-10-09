clc; clear all;

% Create a temperature and density vector to calculate CO2 properties at
%--------------------------------------------------------------------------
T = linspace(216.59,300,50); % Temperature vector, K (-27.67 to 116.33F)
%T = [T linspace(300,304.1,15)];
T = [T linspace(300,307,50)];
%rho = linspace(200,1100); % Density vector, kg/m3
rho = 500;

% Load the CO2 properties created by the NIST Webbook
%--------------------------------------------------------------------------
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
% n = 3;
% T_NIST = downsample(y(:,1),n);
% P_NIST = downsample(y(:,2),n);
% rho_l_NIST = downsample(y(:,3),n);
% rho_v_NIST = downsample(y(:,15),n);
% 
% u_l_NIST = downsample(y(:,5),n);
% u_v_NIST = downsample(y(:,17),n);
% h_l_NIST = downsample(y(:,6),n);
% h_v_NIST = downsample(y(:,18),n);
% s_l_NIST = downsample(y(:,7),n);
% s_v_NIST = downsample(y(:,19),n);
% cv_l_NIST = downsample(y(:,8),n);
% cv_v_NIST = downsample(y(:,20),n);
% cp_l_NIST = downsample(y(:,9),n);
% cp_v_NIST = downsample(y(:,21),n);
% c_l_NIST = downsample(y(:,10),n);
% c_v_NIST = downsample(y(:,22),n);

% Loop through temperature and density vectors to populate the CO2 properties
%--------------------------------------------------------------------------
for i = 1:length(T)
	for j = 1:length(rho)
		Props = CO2Props(T(i),rho(j));
		P(i,j) = Props.P; % Pressure
		X(i,j) = Props.X; % Quality
		s(i,j) = Props.s; % Entropy
		s_l(i,j) = Props.s_l;
		s_v(i,j) = Props.s_v;
		u(i,j) = Props.u; % Internal Energy
		u_l(i,j) = Props.u_l;
		u_v(i,j) = Props.u_v;
		cv(i,j) = Props.cv; % Specific Heat at Constant Volume
		cv_l(i,j) = Props.cv_l;
		cv_v(i,j) = Props.cv_v;
		cp(i,j) = Props.cp; % Specific Heat at Constant Pressure
		cp_l(i,j) = Props.cp_l;
		cp_v(i,j) = Props.cp_v;
		h(i,j) = Props.h; % Enthalpy
		h_l(i,j) = Props.h_l;
		h_v(i,j) = Props.h_v;

		c(i,j) = Props.c; % Speed of Sound
		c_l(i,j) = Props.c_l;
		c_v(i,j) = Props.c_v;
		rho_l(i,j) = Props.rho_l; % Density
		rho_v(i,j) = Props.rho_v;
		state(i,j) = Props.state; % State
	end
end

%Shift energy terms so that the reference matches
%--------------------------------------------------------------------------
u_diff = u_l_NIST(1) - u_l(1)
s_diff = s_l_NIST(1) - s_l(1)
h_diff = h_l_NIST(1) - h_l(1)

% Change the independent variables to temperature and pressure.
%--------------------------------------------------------------------------
%[rhoArr gmaArr] = CO2tp(T_NIST, P_NIST);

% Parse mosster data CO2 density
%--------------------------------------------------------------------------
rho_parse=-0.136928567045648.*T.^2+68.960274348176327.*T-7677.810285569415;

% Plot the pressure vs temperature
%--------------------------------------------------------------------------
figure(1);
plot(T,P(:,1),'b'); hold on; grid on;
plot(T_NIST,P_NIST,'ob');
%title('CO2 Pressure');
legend('Helmholtz','NIST','Location','SE');
xlabel('Temperature, K'); ylabel('Pressure, MPa');
saveas(1,'CO2 Pressure.png')
%saveas(1,'LaTeX Plots/CO2 Pressure.eps','epsc')

% Plot the density vs temperature

%--------------------------------------------------------------------------
figure(2);
plot(T,rho_l(:,1),'b'); hold on; grid on;
plot(T_NIST,rho_l_NIST,'ob');
plot(T,rho_v(:,1),'--r')
plot(T_NIST,rho_v_NIST,'sr')
%plot(T,rho_parse);
%title('CO2 Density');
legend('Liquid-Helmholtz','Liquid-NIST','Vapor-Helmholtz','Vapor-NIST');
%legend('Liquid-Helmholtz','Liquid-NIST','Vapor-Helmholtz','Vapor-NIST','Monster');
xlabel('Temperature, K'); ylabel('Density, kg/m3');
saveas(2,'CO2 Density.png')
%saveas(2,'LaTeX Plots/CO2 Density.eps','epsc')

% Plot the internal energy vs temperature
%--------------------------------------------------------------------------
figure(3);
plot(T,u_l(:,1),'b'); hold on; grid on;
plot(T_NIST,u_l_NIST,'ob');
plot(T,u_v(:,1),'--r')
plot(T_NIST,u_v_NIST,'sr')
title('CO2 Internal Energy');
legend('Liquid-Helmholtz','Liquid-NIST','Vapor-Helmholtz','Vapor-NIST');
xlabel('Temperature, K'); ylabel('Internal Energy, kJ/kg');
saveas(3,'CO2 Internal Energy.png')
%saveas(3,'LaTeX Plots/CO2 Internal Energy.eps','epsc')

% Plot the enthalpy vs temperature
%--------------------------------------------------------------------------
figure(4);
plot(T,h_l(:,1),'b');hold on; grid on;
plot(T_NIST,h_l_NIST,'ob');
plot(T,h_v(:,1),'--r')
plot(T_NIST,h_v_NIST,'sr')
%title('CO2 Enthalpy');

legend('Liquid-Helmholtz','Liquid-NIST','Vapor-Helmholtz','Vapor-NIST','Location','SE');
%legend('Liquid-NIST','Vapor-NIST');
xlabel('Temperature, K'); ylabel('Enthalpy, kJ/kg');
saveas(4,'CO2 Enthalpy.png')
%saveas(4,'LaTeX Plots/CO2 Enthalpy.eps','epsc')

% Plot the entropy vs temperature
%--------------------------------------------------------------------------
figure(5);
plot(T,s_l(:,1),'b'); hold on; grid on;
plot(T_NIST,s_l_NIST,'ob');
plot(T,s_v(:,1),'--r')
plot(T_NIST,s_v_NIST,'sr')
%title('CO2 Entropy');
legend('Liquid-Helmholtz','Liquid-NIST','Vapor-Helmholtz','Vapor-NIST');
xlabel('Temperature, K'); ylabel('Entropy, J/g*K');
saveas(5,'CO2 Entropy.png')
%saveas(5,'LaTex Plots/CO2 Entropy.eps','epsc')

% Plot the specific heat at constant volume vs temperature
%--------------------------------------------------------------------------
figure(6);
plot(T,cv_l(:,1),'b'); hold on; grid on;
plot(T_NIST,cv_l_NIST,'ob');
plot(T,cv_v(:,1),'--r')
plot(T_NIST,cv_v_NIST,'sr')
title('CO2 Specific Heat at Constant Volume');
legend('Liquid-Helmholtz','Liquid-NIST','Vapor-Helmholtz','Vapor-NIST');
xlabel('Temperature, K'); ylabel('Specific Heat at Constant Volume, J/g*K');
saveas(6,'CO2 Specific Heat at Constant Volume.png')
%saveas(6,'LaTeX Plots/CO2 Specific Heat at Constant Volume.eps','epsc')

% Plot the specific heat at constant pressure vs temperature
%--------------------------------------------------------------------------
figure(7);

plot(T,cp_l(:,1),'b'); hold on; grid on;
plot(T_NIST,cp_l_NIST,'ob');
plot(T,cp_v(:,1),'--r')
plot(T_NIST,cp_v_NIST,'sr')
title('CO2 Specific Heat at Constant Pressure');
legend('Liquid-Helmholtz','Liquid-NIST','Vapor-Helmholtz','Vapor-NIST');
xlabel('Temperature, K'); ylabel('Specific Heat at Constant Pressure, J/g*K');
saveas(7,'CO2 Specific Heat at Constant Pressure.png')
%saveas(7,'LaTeX Plots/CO2 Specific Heat at Constant Pressure.eps','epsc')


% Plot the speed of sound vs temperature
%--------------------------------------------------------------------------
figure(8);
plot(T,c_l(:,1),'b'); hold on; grid on;
plot(T_NIST,c_l_NIST,'ob');
plot(T,c_v(:,1),'--r')
plot(T_NIST,c_v_NIST,'sr')
title('CO2 Speed of Sound');
legend('Liquid-Helmholtz','Liquid-NIST','Vapor-Helmholtz','Vapor-NIST');
xlabel('Temperature, K'); ylabel('Speed of Sound, m/s');
saveas(8,'CO2 Speed of Sound.png')
%saveas(8,'LaTeX Plots/CO2 Speed of Sound.eps','epsc')

% Plot the enthalpy vs density
%--------------------------------------------------------------------------
figure(9);
%plot(T,h_l(:,1),'b');
plot(rho_l,h_l,'--b');hold on; grid on;
%plot(T,h_v(:,1),'r')
plot(rho_v,h_v,'r')
%title('CO2 Enthalpy');
%legend('Liquid-Helmholtz','Liquid-NIST','Vapor-Helmholtz','Vapor-NIST');
legend('Liquid-Helmholtz','Vapor-Helmholtz');

xlabel('Denisty, kg/m3'); ylabel('Enthalpy, kJ/kg');
saveas(9,'CO2 Enthalpy vs Density.png')
%saveas(9,'LaTex Plots/CO2 Enthalpy vs Density.eps','epsc')