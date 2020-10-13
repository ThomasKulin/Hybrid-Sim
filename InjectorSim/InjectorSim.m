function [M, rho1, T1, P1, X1, h1, H1, mdot, st1] = InjectorSim(V,M_0,T1_0,n_inj,d_inj,Cd,Pamb,species,dt)
   %% INTRODUCTION
    %==========================================================================
    % Purpose:
    % Returns a sturcture containing the thermodynamic properties of the Ox
    % tank and injector flow properties
    %--------------------------------------------------------------------------
    % The properties in the output structure are:
    % M: Tank fluid mass [kg]
    % rho1: Tank fluid density, rho1 (kg/mˆ3)
    % T1: Tank temperature, T1 (K)
    % P1: Tank pressure, P1 (MPa)
    % X1: Tank quality, X1 ()
    % h1: Tank specific enthalpy, h1 (kJ/kg)
    % H1: Tank total enthalpy, H1 (kJ) ??? Units
    % mdot: Tank mass flow rate, mdot (kg/s)
    % st: Fluid state, -1=- Input, 0=Liq, ...1=Sat, 2=Gas
    %--------------------------------------------------------------------------
    % Inputs:
    % V: Tank volume [M^3]
    % M_0: Initial tank fluid mass [kg]
    % T1_0: Initial tank fluid temperature [k]
    % n_inj: Number of injector orifices
    % d_inj: Injector orifice diameter [m]
    % Cd: Injector discharge coefficient
    % P2: Downstream pressure [MPa]
    % species: Species going through the injector (note, humans are not permitted to enter)
    % dt: Time step [s]
    
 
    %%Calculate the stuff
    %==========================================================================
    % Tank Properties
    %--------------------------------------------------------------------------
    %V=double(V); % Tank volume, m3
    %M_0=double(M_0); % Initial tank fluid mass, kg
    %T1_0=double(T1_0); % Initial tank fluid temperature, K
    R = 0.1889241; % Gas constant, kJ/(kg*K)
    % Injector Properties
    %--------------------------------------------------------------------------
    %n_inj % Number of injector orifices (UNTESTED)-THOMAS
    %d_inj % Injector diameter, m
    Ac = (pi/4)*(d_inj)^2; % Injector cross sectional area; m3
    %Cd=double(Cd); % Injector discharge coefficient
    %Pamb=double(Pamb); % Atmospheric pressure, MPa. (4500ft)
    % Time Iteration
    %--------------------------------------------------------------------------
    %dt=double(dt); % Step time, s
    %species  % pick CO2 or N2O
    NIST = true;  % pull thermo data from NIST database or from approximation eqn(only CO2)
    Inj_Switch = 1; % Isentropic = 1, Adiabatic = 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


end