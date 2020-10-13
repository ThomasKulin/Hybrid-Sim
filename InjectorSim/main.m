clear
clc;
close all;
% Initial Conditions
%==========================================================================
% Tank Properties
%--------------------------------------------------------------------------
V = .01824; % Tank volume, m3
M_0 = 13; % Initial tank fluid mass, kg
T1_0 = 274.25; % Initial tank fluid temperature, K
R = 0.1889241; % Gas constant, kJ/(kg*K)
% Injector Properties
%--------------------------------------------------------------------------
n_inj = 20; % Number of injector orifices (UNTESTED)-THOMAS
d_inj = 0.0015; % Injector diameter, m
Ac = (pi/4)*(d_inj)^2; % Injector cross sectional area; m3
Cd = 0.8; % Injector discharge coefficient
Pamb = 85.9e-3; % Atmospheric pressure, MPa. (4500ft)
% Time Iteration
%--------------------------------------------------------------------------
tstop = 15; % Stop time, s
dt = .5; % Step time, s
species = "N2O";  % pick CO2 or N2O
NIST = true;  % pull thermo data from NIST database or from approximation eqn(only CO2)
Inj_Switch = 1; % Isentropic = 1, Adiabatic = 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[a,b,c,d,e,f,g,h,i,j] = InjectorSim(V,M_0,T1_0,n_inj,d_inj,Cd,Pamb,species,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize state
%--------------------------------------------------------------------------
rho1_0 = M_0/V;
if ~NIST
    Props1_0 = CO2Props(T1_0,rho1_0);  % equation based CO2 Model
else
    Props1_0 = PropsNIST(T1_0,rho1_0,species);  % Data table lookup
end
P1_0 = Props1_0.P;
rhol_0 = Props1_0.rho_l;
rhov_0 = Props1_0.rho_v;
X1_0 = (rhov_0/rho1_0)*(rhol_0-rho1_0)/(rhol_0-rhov_0);
h1_0 = Props1_0.h;
H1_0 = M_0*Props1_0.h;
st_0 = Props1_0.state;

t = 0; % Column 1: Time (sec)
M = M_0; % Column 2: Tank fluid mass, M (kg)
rho1 = rho1_0; % Column 3: Tank fluid density, rho1 (kg/mˆ3)
T1 = T1_0; % Column 4: Tank temperature, T1 (K)
P1 = P1_0; % Column 5: Tank pressure, P1 (MPa)
X1 = X1_0; % Column 6: Tank quality, X1 ()
h1 = h1_0; % Column 7: Tank specific enthalpy, h1 (kJ/kg)
H1 = H1_0; % Column 8: Tank total enthalpy, H1 (J) ??? Units
mdot = 0; % Column 9: Tank mass flow rate, mdot (kg/s)
P2 = Pamb; % Column 10: Injector outlet pressure, P2 (Pa)
st = st_0; % Column 11: Fluid state, -1=- Input, 0=Liq, 1=Sat, 2=Gas


State = [t, M, rho1, T1, P1, X1, h1, H1, mdot, P2, st];
State2 = [0, 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 0;
while t < tstop

    i = i + 1

    % Exit the loop when the tank is out of fluid
    %----------------------------------------------------------------------
    if M < 0.1; break; end

    % Initial Tank Fluid Properties
    %======================================================================
    %Note: this guess MUST yeild a quality less than 1 and >0
    guess=[300 300]; %[T,rho]
    
    %Set up function (need to match pressure and quality)
    if ~NIST
        pFunc = @(v) [getfield(CO2Props(v(1),v(2)),'P')-P1; getfield(CO2Props(v(1),v(2)),'X')-X1];
    else
        pFunc = @(v) [getfield(PropsNIST(v(1),v(2),species),'P')-P1; getfield(PropsNIST(v(1),v(2),species),'X')-X1];
    end

    % Solve for [T,rho] of saturated but pure liquid at P1
    % lsqnonlin tries to find a T and rho that makes pFunc = 0
    v1 = lsqnonlin(pFunc,guess,0,inf,optimset('Display','off','TolFun',1e-14));
    T1 = v1(1); rho1 = v1(2);
    if ~NIST
        Props1 = CO2Props(T1,rho1);
    else
        Props1 = PropsNIST(T1,rho1,species);
    end
    Pv1 = Props1.P; % Fluid Vapor Pressure, MPa
    rhoL1 = Props1.rho_l; % Fluid liquid density, kg/m3
    h1 = Props1.h; % Fluid specific enthalpy, kJ/kg
    H1 = M*h1; % Fluid total enthalpy, kJ
    s1 = Props1.s; % Fluid entropy, kJ/kg*K

    % Propegate properties across the injector
    %======================================================================

    %Assume an isentropic (s1=s2) or adiabatic (h1=h2) injector , solve for properties at P2
    if Inj_Switch == 1
        %Set up function to match pressure and entropy
        if ~NIST
            pFunc = @(v) [getfield(CO2Props(v(1),v(2)),'P')-P2 getfield(CO2Props(v(1),v(2)),'s')-s1];
        else
            pFunc = @(v) [getfield(PropsNIST(v(1),v(2),species),'P')-P2 getfield(PropsNIST(v(1),v(2),species),'s')-s1];
        end
    elseif Inj_Switch == 2
        %Set up function to match pressure and enthalpy
        if ~NIST
            pFunc = @(v) [getfield(CO2Props(v(1),v(2)),'P')-P2 getfield(CO2Props(v(1),v(2)),'h')-h1];
        else
            pFunc = @(v) [getfield(PropsNIST(v(1),v(2),species),'P')-P2 getfield(PropsNIST(v(1),v(2),species),'h')-h1];
        end
    end

    % Solve for T2 & rho2 downstream of the injector
    %----------------------------------------------------------------------
    guess=[300 300]; %[T,rho]
    v2 = lsqnonlin(pFunc,guess,0,inf,optimset('Display','off','TolFun',1e-14));
    T2 = v2(1); rho2 = v2(2);
    % Solve for the remaining properties downstream of the injector
    %----------------------------------------------------------------------
    if ~NIST
        Props2=CO2Props(T2,rho2); % Downstream fluid properties
    else
        Props2=PropsNIST(T2,rho2,species);
    end
    Pv2 = Props2.P; % Downstream fluid vapor pressure
    h2 = Props2.h; % Downstream fluid enthalpy

    % Calculate the Mass Flow
    %======================================================================
    % Non-Equalibrium Parameter
    %k = sqrt(abs(P1-P2)/abs(Pv2-P1)); % Whitmores Equation
    k = sqrt((P1-P2)/(Pv1-P2)); % Spencer and Stanfords Equation
    % Weighting Coefficient
    W = (1/(k+1));
    % Incompressible fluid mass flow rate
    mdot_inc = Ac*sqrt(2*rhoL1*(P1-P2)*1e6)*n_inj;
    % Homogeneous Equilibrium mass flow rate
    mdot_HEM = rho2*Ac*sqrt(2*(h1-h2))*n_inj;
    % Weighted Non-Homogeneous Equilibrium (modified Stanford) mass flow rate
    mdot = Cd*((1-W) * mdot_inc + W * mdot_HEM); % Shannon's Theory
    
    % Update the upstream fluid properties for the next step
    %======================================================================
    M = M - mdot*dt; % Update tank fluid mass
    Hdot = h1*mdot; % Enthalpy flow rate
    H1 = H1 - Hdot*dt; % Update tank total enthalpy

    % Calculate the new tank enthalpy and density
    %---------------------------------------------------------------------
    rho1 = M/V; % Update tank specific density
    h1 = H1/M; % Update tank specific enthalpy

    % Calculate the new tank temperature
    %----------------------------------------------------------------------
    % Create a function for lsqnonlin to solve T(rho,h)
    if ~NIST
        pFunc = @(T_Unknown) getfield(CO2Props(T_Unknown,rho1),'h')-h1;
    else
        pFunc = @(T_Unknown) getfield(PropsNIST(T_Unknown,rho1,species),'h')-h1;
    end
    % Since T Unknown is not pre defined in pFunc, lsqnonlin will find a ...T for rho Known and h Known
    T1 = lsqnonlin(pFunc,300,0,inf,optimset('Display','off','TolFun',1e-14));
    % [T1] = CO2 rho h 2T(rho1,h1,300);

    % Calulate the new tank pressure and quality
    %----------------------------------------------------------------------
    if ~NIST
    	Props1 = CO2Props(T1,rho1);
    else
        Props1 = PropsNIST(T1,rho1,species);
    end
    P1 = Props1.P;
    X1 = Props1.X;

    st1 = Props1.state;

    % Update the state
    %----------------------------------------------------------------------
    t = t + dt;

    State = [State; t, M, rho1, T1, P1, X1, h1, H1, mdot, Pamb, st1];
    State2 = [State2; mdot_inc, mdot_HEM];
end

 % Up-pack the state vector for ploting
 %==========================================================================
t = State(:,1); % Column 1: Time (sec)
M = State(:,2); % Column 2: Tank fluid mass, M (kg)
rho1 = State(:,3); % Column 3: Tank fluid density, rho1 (kg/mˆ3)
T1 = State(:,4); % Column 4: Tank temperature, T1 (K)

P1 = State(:,5); % Column 5: Tank pressure, P1 (Pa)
X1 = State(:,6); % Column 6: Tank quality, X1 ()
h1 = State(:,7); % Column 7: Tank specific enthalpy, h1 (kJ/kg)
H1 = State(:,8); % Column 8: Tank total enthalpy, H1 (J) ??? Units
mdot = State(:,9); % Column 9: Tank mass flow rate, mdot (kg/s)
P2 = State(:,10); % Column 10: Injector outlet pressure, P2 (Pa)
st = State(:,11); % Column 11: Fluid state, -1=- Input, 0=Liq, ...1=Sat, 2=Gas

mdot_inc = State2(:,1); % mdot_incompressible
mdot_HEM = State2(:,2); % mdot_HEM

%% Plot
%==========================================================================
% Column 2 - Tank fluid mass
%--------------------------------------------------------------------------
figure(2)
plot(t,M); hold on; grid on;
title('Tank Mass');
xlim([0 tstop]);
xlabel('Time, s'); ylabel('Tank Fluid Mass, kg');
saveas(2,'Tank Mass.png')

% Column 3 - Tank fluid specific density
%--------------------------------------------------------------------------
% figure(3);
% plot(t,rho1); hold on; grid on;
% title('rho1');
% xlabel('Time, s'); ylabel('Density, kg/m3');
% saveas(3,'Tank Density.png')

% Column 4 - Tank temperature
%--------------------------------------------------------------------------
figure(4);
plot(t,T1); hold on; grid on;
title('Temperature');
xlim([0 tstop]);
xlabel('Time, s'); ylabel('Temperature, K');
saveas(4,'Tank Temperature.png')

% Column 5 - Tank pressure
%--------------------------------------------------------------------------
figure (5)
plot(t,P1); hold on; grid on;
title('P1');
xlim([0 tstop]);
xlabel('Time, s'); ylabel('Pressure, MPa');
saveas(5,'Tank Pressure.png')

% Column 6 - Tank quality
%--------------------------------------------------------------------------
% figure (6)
% plot(t,X1); hold on; grid on;
% ylim([0 1]);
% title('X1');
% xlabel('Time, s'); ylabel('Quality');
% legend('Model Predicted');
% saveas(6,'Tank Quality.png')

% Coumn 7 - Tank specific enthalpy
%--------------------------------------------------------------------------
% figure (7)
% plot(t,h1); hold on; grid on;
% title('h1');
% xlabel('Time, s'); ylabel('Specific Enthalpy, kJ/kg-K');
% legend('Model Predicted');
% saveas(7,'Tank Specific Enthalpy.png')

% Coumn 8 - Tank total enthalpy
%--------------------------------------------------------------------------
% figure (8)
% plot(t,H1); hold on; grid on;
% title('H1'); xlabel('time, s'); ylabel('Enthalpy');
% legend('Model Predicted');
% saveas(8,'Tank Total Enthalpy.png')

% Column 9 - Mass flow rate
%--------------------------------------------------------------------------
% figure(9);
% plot(t,mdot,'b'); hold on; grid on;
% xlabel('Time, s');
% ylabel('Mass Flow Rate, kg/s');
% xlim([0 tstop]); %ylim([0 3]);
% saveas(9,'Tank Mass Flow Rate Venturi.png')

figure(90);
plot(t,mdot,'b'); hold on; grid on;
plot(t,mdot_inc,':b');
plot(t,mdot_HEM,'--bo');
xlabel('Time, s'); ylabel('Mass Flow Rate, kg/s');
xlim([0 tstop]); ylim([0 1.6]);
legend('Model Predicted','Incompressible Model','HEM Model');
saveas(90,'Tank Mass Flow Rate.png')


% Column 11 - Tank fluid state
%--------------------------------------------------------------------------
% figure(11);
% plot(t,st); grid on;
% title('Fluid State');
% xlabel('Time, s'); ylabel('Fluid State');
% ylim([-2 3]);
% legend('Model Predicted');
% saveas(11,'Tank Fluid State.png')
% %