%% ELECTROPORATION PARAMETERS
k = 1.38065e-23;

% Electroporation parameters
par.rstar = 0.51e-9; % Pore radius at hydrophobic and hydrophilic transition [m]
par.rpmin = 0.8e-9; % Local minimum in the pore radius space at Um = 0 [m]
par.rpmax = 50e-9; % Maximum pore radius [m]
par.Wstar_kT = 45; % Pore creation energy barrier [kT]
par.a = 1e9; % Pore creation rate density [1/(m^2*s)]
par.alpha = 0; % Asymmetric pore creation constant [kT/V] 
par.beta = 1/0.258^2; % Symmetric pore creation constant [kT/V^2]
par.Dp = 5e-14; % Diffusion coefficient in the pore radius space [m^2/s]
par.gamma = 2e-11; % Edge tension [N]
par.Gamma0 = 1e-6; % Surface tension of nonelectroporated membrane [N/m]
par.Gamma1 = 2e-2; % Hydrocarbon-water interfacial tension [N/m]
par.Fmax = 6.9e-10; % Parameter in the electrical force tending to expand the pore [N/V^2]
par.rh = 0.95e-9; % Parameter in the electrical force tending to expand the pore [m]
par.rt = 0.23e-9; % Parameter in the electrical force tending to expand the pore [m]
par.nu = 0.25; % Relative length of the pore entrance
par.fprot = 0; % Areal fraction of proteins
par.taup = 1.5; % Pore resealing time constant [s]

% Options
% Expression for conductance of a single pore
par.expr_Gp = 'Li';      % Choose between: 'Li', 'Smith', 'Toroidal'
% Expression for pore surface energy
par.expr_Wsurf = 'Li';   % Choose between: 'Li', 'Smith'
% Expression for molar flux across the membrane
par.expr_molflux = 'Li'; % Choose between: 'Li', 'Smith', 'Combination'

% Discretization
par.axi2D = true;
par.Ntheta = 60; % Number of discretization angles
par.dcomp2 = 5e-9; % Separation between lines in Component 2
par.rmesh0 = 0.8e-6; % Max mesh size within the cell
par.rmesh1 = 1e-9; % Max mesh size in pore radius space on [2 nm, rpmax - 2 nm]
par.rmesh2 = 1e-10; % Max mesh size in pore radius space on [rpmin, 2 nm] and [rpmax - 2 nm, rpmax]
par.rmesh3 = 2.5e-12; % Max mesh size in pore radius space around the point rp = rstar
