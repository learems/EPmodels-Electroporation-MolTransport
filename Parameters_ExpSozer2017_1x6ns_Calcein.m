%% MODEL PARAMETERS

% Time stepping
par.tlist1 = 'range(0,5e-10,1e-8) 10^{range(-8,(-3+8)/50,-3)} 1e-3+10^{range(log10(1e-3),(log10(2e2)-log10(1e-3))/350, log10(2e2))}';                                 
par.tlist2 = 'range(0,5e-9,1e-8) 10^{range(-8,(-3+8)/10,-3)}  1e-3+10^{range(log10(1e-3),(log10(2e2)-log10(1e-3))/50, log10(2e2))}';

% Pulse parameters and temperature
par.pulseshape = 'trapezoidal'; % Pulse shape
par.Eapp = 20e6; % Applied electric field [V/m]
par.tpulse = 6e-9; % Pulse duration [s]
par.trise = 5e-9; % Pulse rise time [s]
par.tfall = 8e-9; % Pulse fall time [s]
par.npulse = 1; % Number of pulses
par.fp = 1e-10; % Repetition frequency [Hz]
par.T = 295.15; % Temperature [K]

% Geometrical parameters
par.rcell = 7e-6; % Cell radius [m]
par.dcell = 2*par.rcell; % Cell height [m] (for cylindrical cells)
par.del = 30*par.rcell; % Radius of the extracellular domain [m] 

% Cell electrical parameters
par.sigma_e = 1.7; % External conductivity [S/m]
par.sigma_i = 0.3; % Internal conductivity [S/m]
par.epsilon_e = 72; % External permittivity []
par.epsilon_i = par.epsilon_e; % Internal permittivity []
par.eps0 = 8.8542e-12; % Vacuum permittivity [F/m]
par.dm = 5e-9; % Membrane thickness [m]
par.Cm = par.eps0*5 / par.dm; % Membrane capacitance [F/m^2]
par.Gm = 9.5e-9 / par.dm; % Membrane conductance [S/m^2]
par.Urest = -50e-3; % Resting voltage (Um = Vi - Ve)

% Parameters for molecular transport
% Idx 1 = Calcein
par.Nsolutes = 1; % Number of solute molecules simulated; Choose either 1 or 3
par.r0 = 0.19e-9; % Radius of charge carrier
par.l0 = 2*par.r0; % Length of charge carrier
par.r1 = 0.58e-9; % Radius of solute 1
par.l1 = 1.89e-9; % Length of solute 1
par.z1 = -3.61; % Valence of solute 1
par.D1_e = 4.66e-10; % Extracellular diff. coefficient of solute 1 [m^2/s]
par.D1_i = par.D1_e/4; % Intracellular diff. coefficient of solute 1 [m^2/s]
par.c1_i0 = 0; % Initial intracellular concentration of solute 1 [mol/m^3]
par.c1_e0 = 0.2; % Initial extracellular concentration of solute 1 [mol/m^3]
%%%% Following parameters are irrelevant when Nsolutes = 1
% par.r2 = 0.61e-9; % Radius of solute 2
% par.l2 = 1.46e-9; % Length of solute 2
% par.z2 = 1; % Valence of solute 2
% par.D2_i = 0; % Intracellular diff. coefficient of solute 2 [m^2/s]
% par.D2_e = 0; % Extracellular diff. coefficient of solute 2 [m^2/s]
% par.c2_i0 = 0; % Initial intracellular concentration of solute 2 [mol/m^3]
% par.c2_e0 = 0; % Initial extracellular concentration of solute 2 [mol/m^3]
% par.r12 = 0.61e-9; % Radius of bound solutes 1&2
% par.l12 = 1.46e-9; % Length of bound solutes 1&2
% par.z12 = 1; % Valence of bound solutes 1&2
% par.D12_i = 0; % Intracellular diff. coefficient of bound solutes 1&2 [m^2/s]
% par.D12_e = 0; % Extracellular diff. coefficient of bound solutes 1&2 [m^2/s]
% par.c12_i0 = 0; % Initial intracellular concentration of bound solutes 1&2 [mol/m^3]
% par.c12_e0 = 0; % Initial extracellular concentration of bound solutes 1&2 [mol/m^3]
% par.kass = 0; % Association constant [1/(s*mol/m^3)]
% par.kdis = 0; % Dissociation constant [1/s]
% 1 M = 1 mol/l = 1e3 mol/m^3;

% Depth for convolved image
par.y_opt = 1e-6;