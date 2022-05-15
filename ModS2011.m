% Electroporation model 
% based on Smith 2011 with modifications
% Created in Comsol 5.2, updated in Comsol 5.6
% Author: Lea Rems
% Created: Dec 13, 2021
% Version of the script: 5b

%% Some prerequisites for creating the model

% Related to pore energy function at Um = 0 V
% This part is necessary for determining some of the model parameters
kT = 1.38065e-23*par.T;
rp = linspace(par.rstar,3e-9,501)';
Fedge = -2*pi*par.gamma;
if strcmp(par.expr_Wsurf,'Li')
    Fsurf = 2*pi*rp*par.Gamma0;
    Fsurf_rpmin = 2*pi*par.rpmin*par.Gamma0;
    Wsurf_rstar = -pi*par.rstar^2*par.Gamma0;
    Wsurf_rpmin = -pi*par.rpmin^2*par.Gamma0;
elseif strcmp(par.expr_Wsurf,'Smith')
    Fsurf = (2*pi*rp + 0.5*pi*par.dm*(2-pi))*par.Gamma0;
    Fsurf_rpmin = (2*pi*par.rpmin + 0.5*pi*par.dm*(2-pi))*par.Gamma0;
    dAlp_rstar = pi*(par.rstar + 0.5*par.dm)^2 - pi^2*0.5*par.dm*(par.rstar + 0.5*par.dm) + 2*pi*(0.5*par.dm)^2;
    dAlp_rpmin = pi*(par.rpmin + 0.5*par.dm)^2 - pi^2*0.5*par.dm*(par.rpmin + 0.5*par.dm) + 2*pi*(0.5*par.dm)^2;
    Wsurf_rstar = -dAlp_rstar*par.Gamma0;
    Wsurf_rpmin = -dAlp_rpmin*par.Gamma0;
end

% Determine pore destruction barrier Wd for the chosen taup
% Solution for Wd is found where fun1 = 0. Zero crossing of fun1 is determined by interpolation. 
Wd_test= linspace(10,30,10001)*kT;
fun1 = par.taup - 0.5*sqrt(pi)*(par.rpmin-par.rstar)^2/par.Dp*(Wd_test./kT).^(-1.5).*exp(Wd_test./kT);
Wd = interp1(fun1,Wd_test,0);

% Determine B and b
% Solution for b is found whete fun1 = 0. Zero crossing of fun1 is determined by interpolation. 
btest = linspace(0.1,10,10001);
fun1 = -Wd + 2*pi*par.gamma*(par.rstar-par.rpmin) + (2*pi*par.gamma - Fsurf_rpmin)./(btest.*par.rstar.^btest.*par.rpmin.^(-btest-1)).*(1 - (par.rstar/par.rpmin).^btest) + Wsurf_rstar - Wsurf_rpmin;
par.b = interp1(fun1,btest,0);
par.B = (2*pi*par.gamma - Fsurf_rpmin)./(par.b*par.rstar.^par.b*par.rpmin^(-par.b-1));
Fsteric = par.B*par.b*(par.rstar)^par.b*(rp).^(-par.b-1);

% Determine neq_star
Fp = Fsteric + Fedge + Fsurf;
% Equilibrium pore density distribution and the total pore density are equal to:
% neq = neq_star*exp(cumtrapz(rp,Fp/kT))
% Neq = trapz(rp,neq) = neq_star*trapz(rp,exp(cumtrapz(rp,Fp/kT)))
par.Neq = 0.01*par.taup*par.a*(1-par.fprot); 
% We add the factor 0.01 such that Neq/a has negligible contribution to the resealing time. 
par.neq_star = par.Neq/trapz(rp,exp(cumtrapz(rp,Fp./kT)));

% Determine parameter rd in Fdipole
Alipedge = 0.5*pi^2*par.dm*(rp + 0.5*par.dm) - 2*pi*(0.5*par.dm)^2;
Alipedge_star = 0.5*pi^2*par.dm*(par.rstar + 0.5*par.dm) - 2*pi*(0.5*par.dm)^2;
Olipedge = (2*pi^2*rp + (pi^2 - 8)*par.dm)./(8*pi*rp + 4*(pi-2)*par.dm);
Olipedge_star = (2*pi^2*par.rstar + (pi^2 - 8)*par.dm)/(8*pi*par.rstar + 4*(pi-2)*par.dm);
expr = Alipedge.*sin(Olipedge)./(Alipedge_star.*sin(Olipedge_star));
P = polyfit(rp-par.rstar,expr,1);
par.rd = 1/P(1);

% Determine parameter Ur which is used to decribe the resting potential
if strcmp(par.expr_Gp,'Li')
    par.sigma_p = (par.sigma_e-par.sigma_i)/(log(par.sigma_e/par.sigma_i) + (par.sigma_e == par.sigma_i)) + par.sigma_e*(par.sigma_e == par.sigma_i);
    Gpmin = 2*pi*par.rpmin^2*par.sigma_p/(pi*par.rpmin + 2*par.dm);
elseif strcmp(par.expr_Gp,'Smith')
    H0 = 0.4; K0 = 0.65; % Values of H0 and K0 for charge carrier with radius of 0.19 nm
    par.sigma_p = 2*par.sigma_e*par.sigma_i/(par.sigma_e + par.sigma_i);
    Gpmin = 2*pi*par.rpmin^2*par.sigma_p*K0*H0/(pi*par.rpmin*K0*H0 + par.dm);    
elseif strcmp(par.expr_Gp,'Toroidal')
    par.sigma_p = (par.sigma_e-par.sigma_i)/(log(par.sigma_e/par.sigma_i) + (par.sigma_e == par.sigma_i)) + par.sigma_e*(par.sigma_e == par.sigma_i);
    Gpmin = 2*pi*par.rpmin^2*par.sigma_p/(pi*par.rpmin + par.dm*(1 + tanh(0.1360*(par.rpmin*1e9-1.2408))));
end
Gep_eq = par.Neq*Gpmin; % Membrane conductace due to pores at Um = 0 V
par.Ur = (par.Gm + Gep_eq)/par.Gm*par.Urest; % Parameter to account for the resting voltage [V]

% Diffusion coeficients within the membrane
if strcmp(par.expr_molflux,'Li') || strcmp(par.expr_molflux,'Combination')
    par.D1_m = (par.D1_e-par.D1_i)/(log(par.D1_e/par.D1_i) + (par.D1_e == par.D1_i)) + par.D1_e*(par.D1_e == par.D1_i);
    if par.Nsolutes == 3
        par.D2_m = (par.D2_e-par.D2_i)/(log(par.D2_e/par.D2_i) + (par.D2_e == par.D2_i)) + par.D2_e*(par.D2_e == par.D2_i);
        par.D12_m = (par.D12_e-par.D12_i)/(log(par.D12_e/par.D12_i) + (par.D12_e == par.D12_i)) + par.D12_e*(par.D12_e == par.D12_i);
    end
elseif strcmp(par.expr_molflux,'Smith')
    par.D1_m = 2*par.D1_i*par.D1_e/(par.D1_i+par.D1_e);
    if par.Nsolutes == 3
        par.D2_m = 2*par.D2_i*par.D2_e/(par.D2_i+par.D2_e);
        par.D12_m = 2*par.D12_i*par.D12_e/(par.D12_i+par.D12_e);
    end
end

%% Clear variables which are not needed anymore
clear Gpmin Gep_eq kT rp neq_star Neq btest Wd Wd_test fun1
clear Fsteric Fedge Fsurf Fp Fsurf_rpmin Wsurf_rstar Wsurf_rpmin dAlp_rstar dAlp_rpmin
clear Alipedge Alipedge_star Olipedge Olipedge_star expr P H0 K0

%% CREATE MODEL
fprintf('Creating model ...................................................\n')

clear model
import com.comsol.model.*
import com.comsol.model.util.*
model = ModelUtil.create('Model');
model.comments(['']);

% Create component 1 (the cell in the electric field)
model.modelNode.create('comp1');
model.geom.create('geom1', 2);
if par.axi2D
    model.geom('geom1').axisymmetric(true);
else
    model.component('comp1').spatialCoord({'r' 'z' 'w'});
end
model.mesh.create('mesh1', 'geom1');
% Create component 2 (pore radius space)
model.modelNode.create('comp2');
model.geom.create('geom2', 2);
model.mesh.create('mesh2', 'geom2');
model.frame('material2').coord({'rp' 'y' 'z'});

%% PARAMETERS

% Pulse parameters and temperature
model.param.set('Eapp', [num2str(par.Eapp),'[V/m]'], 'Applied electric field');
model.param.set('Uapp', 'Eapp*del', 'Applied voltage');
model.param.set('tpulse', [num2str(par.tpulse),'[s]'], 'Pulse duration');
model.param.set('trise', [num2str(par.trise),'[s]'], 'Pulse rise time');
model.param.set('tfall', [num2str(par.tfall),'[s]'], 'Pulse fall time');
model.param.set('npulse', num2str(par.npulse), 'Number of pulses');
model.param.set('fp', [num2str(par.fp),'[Hz]'], 'Pulse repetition frequency');
model.param.set('T', [num2str(par.T),'[K]'], 'Temperature');

% Geometrical parameters
model.param.set('rcell', [num2str(par.rcell),'[m]'], 'Cell radius');
model.param.set('dcell', [num2str(par.dcell),'[m]'], 'Height of the cell, relevant only in cylindrical coordinate system');
model.param.set('del', [num2str(par.del),'[m]'], 'Distance between the electrodes');

% Cell electrical parameters
model.param.set('sigma_e', [num2str(par.sigma_e),'[S/m]'], 'Medium conductivity');
model.param.set('epsilon_e', [num2str(par.epsilon_e)], 'Medium permittivity');
model.param.set('sigma_i', [num2str(par.sigma_i),'[S/m]'], 'Cytoplasmic conductivity');
model.param.set('epsilon_i', [num2str(par.epsilon_i)], 'Cytoplasmic permittivity');
model.param.set('sigma_rat', 'sigma_i/sigma_e', 'Ratio between the intra- and extracellular conductivity');
model.param.set('Cm', [num2str(par.Cm),'[F/m^2]'], 'Plasma membrane capacitance');
model.param.set('Gm', [num2str(par.Gm),'[S/m^2]'], 'Plasma membrane conductance');
model.param.set('dm', [num2str(par.dm),'[m]'], 'Plasma membrane thickness');
model.param.set('dp', 'dm/2', 'Thickness of the internal pore region');
model.param.set('sigma_m', 'Gm*dm', 'Plasma membrane conductivity');
model.param.set('epsilon_m', 'Cm*dm/eps0', 'Plasma membrane permittivity');
model.param.set('Urest', [num2str(par.Urest),'[V]'], 'Resting voltage');
model.param.set('Ur', [num2str(par.Ur),'[V]'], 'Parameter to account for the resting voltage');

%=======================%%%%%%%%%%%%%%%%%%%%%%%%%%========================%
    % Electroporation parameters
    model.param.set('rstar', [num2str(par.rstar),'[m]'], 'Pore radius at hydrophobic and hydrophilic transition [m]');
    model.param.set('rpmin', [num2str(par.rpmin),'[m]'], 'Local minimum in the pore radius space at Um = 0');
    model.param.set('rpmax', [num2str(par.rpmax),'[m]'], 'Maximum pore radius');
    model.param.set('a', [num2str(par.a),'[1/(m^2*s)]'], 'Electroporation parameter - prefactor');
    model.param.set('alpha', [num2str(par.alpha),'[1/V]'], 'Asymmetric pore creation constant');
    model.param.set('beta', [num2str(par.beta),'[1/V^2]'], 'Symmetric pore creation constant');
    model.param.set('neq_star', [num2str(par.neq_star),'[1/m^3]'], 'Equilibrium pore density at rp = rstar and Um = 0 V');
    model.param.set('Dp', [num2str(par.Dp),'[m^2/s]'], 'Diffusion coefficient in the pore radius space');
    model.param.set('gamma', [num2str(par.gamma),'[N]'], 'Edge tension');
    model.param.set('Gamma0', [num2str(par.Gamma0),'[N/m]'], 'Surface tension of nonelectroporated membrane');
    model.param.set('Gamma1', [num2str(par.Gamma1),'[N/m]'], 'Hydrocarbon-water interfacial tension');
    model.param.set('Fmax', [num2str(par.Fmax),'[N/V^2]'], 'Parameter in the electrical force tending to expand the pore');
    model.param.set('rh', [num2str(par.rh),'[m]'], 'Parameter in the electrical force tending to expand the pore');
    model.param.set('rt', [num2str(par.rt),'[m]'], 'Parameter in the electrical force tending to expand the pore');
    model.param.set('rd', [num2str(par.rd),'[m]'], 'Parameter in Fdipole');
    model.param.set('B', [num2str(par.B),'[J]'], 'Parameter in Wster pore energy');
    model.param.set('b', [num2str(par.b)], 'Parameter in Wster pore energy');    
    model.param.set('nu', [num2str(par.nu)], 'Relative length of the pore mouth');
    model.param.set('sigma_p', [num2str(par.sigma_p),'[S/m]'], 'Effective conductivity inside the pore');
    model.param.set('fprot', [num2str(par.fprot)], 'Areal fraction of proteins');
    model.param.set('taup', [num2str(par.taup),'[s]'], 'Pore destruction rate');
    
    % Discretization parameters
    par.Nnodes = par.Ntheta + 1; % Number of membrane nodes    
    model.param.set('Ntheta', num2str(par.Ntheta), 'Number of discretization angles');
    model.param.set('Nnodes', num2str(par.Nnodes), 'Number of membrane nodes');
    model.param.set('dtheta', 'pi/Ntheta', 'Discretization angle');
    model.param.set('dcomp2', [num2str(par.dcomp2),'[m]'], 'Separation between lines in Component 2');
%=========================================================================%

% Physical constants
model.param.set('eps0', '8.8542e-12[F/m]', 'Vacuum permittivity');
model.param.set('qe', '1.602e-19[A*s]', 'Elementary charge');
model.param.set('k', '1.38065e-23[J/K]', 'Boltzmann constant');
model.param.set('F', '96485[A*s/mol]', 'Faraday constant');
model.param.set('R', '8.314[J/(K*mol)]', 'Universal gas constant');
model.param.set('NA', '6.022e23[1/mol]','Avogadro constant');

% Parameters for molecular transport
model.param.set('r0', [num2str(par.r0),'[m]'], 'Radius of charge carrier');
model.param.set('l0', [num2str(par.l0),'[m]'], 'Length of charge carrier');
model.param.set('r1', [num2str(par.r1),'[m]'], 'Radius of solute 1');
model.param.set('l1', [num2str(par.l1),'[m]'], 'Length of solute 1');
model.param.set('z1', [num2str(par.z1)], 'Valence of solute 1');
model.param.set('D1_i', [num2str(par.D1_i),'[m^2/s]'], 'Intracellular diff. coefficient of solute 1');
model.param.set('D1_e', [num2str(par.D1_e),'[m^2/s]'], 'Extracellular diff. coefficient of solute 1');
model.param.set('D1_m', [num2str(par.D1_m),'[m^2/s]'], 'Diff. coefficient of solute 1 within the membrane');
model.param.set('c1_i0', [num2str(par.c1_i0),'[mol/m^3]'], 'Initial intracellular concentration of solute 1');
model.param.set('c1_e0', [num2str(par.c1_e0),'[mol/m^3]'], 'Initial extracellular concentration of solute 1');
if par.Nsolutes == 3
    model.param.set('r2', [num2str(par.r2),'[m]'], 'Radius of solute 2');
    model.param.set('l2', [num2str(par.l2),'[m]'], 'Length of solute 2');
    model.param.set('z2', [num2str(par.z2)], 'Valence of solute 2');
    model.param.set('D2_i', [num2str(par.D2_i),'[m^2/s]'], 'Intracellular diff. coefficient of solute 2');
    model.param.set('D2_e', [num2str(par.D2_e),'[m^2/s]'], 'Extracellular diff. coefficient of solute 2');
    model.param.set('D2_m', [num2str(par.D2_m),'[m^2/s]'], 'Diff. coefficient of solute 2 within the membrane');
    model.param.set('c2_i0', [num2str(par.c2_i0),'[mol/m^3]'], 'Initial intracellular concentration of solute 2');
    model.param.set('c2_e0', [num2str(par.c2_e0),'[mol/m^3]'], 'Initial extracellular concentration of solute 2');
    model.param.set('r12', [num2str(par.r12),'[m]'], 'Radius of bound solutes 1&2');
    model.param.set('l12', [num2str(par.l12),'[m]'], 'Length of bound solutes 1&2');
    model.param.set('z12', [num2str(par.z12)], 'Valence of bound solutes 1&2');
    model.param.set('D12_i', [num2str(par.D12_i),'[m^2/s]'], 'Intracellular diff. coefficient of bound solutes 1&2');
    model.param.set('D12_e', [num2str(par.D12_e),'[m^2/s]'], 'Extracellular diff. coefficient of bound solutes 1&2');
    model.param.set('D12_m', [num2str(par.D12_m),'[m^2/s]'], 'Diff. coefficient of solute 12 within the membrane');
    model.param.set('c12_i0', [num2str(par.c12_i0),'[mol/m^3]'], 'Initial intracellular concentration of bound solutes 1&2');
    model.param.set('c12_e0', [num2str(par.c12_e0),'[mol/m^3]'], 'Initial extracellular concentration of bound solutes 1&2');
    model.param.set('kass', [num2str(par.kass),'[1/(s*mol/m^3)]'], 'Association constant');
    model.param.set('kdis', [num2str(par.kdis),'[1/s]'], 'Dissociation constant');
end
model.param.set('a1', '-1.2167', 'Bungay&Brenner hindrance');
model.param.set('a2', '1.5336', 'Bungay&Brenner hindrance');
model.param.set('a3', '-22.5083', 'Bungay&Brenner hindrance');
model.param.set('a4', '-5.6117', 'Bungay&Brenner hindrance');
model.param.set('a5', '-0.3363', 'Bungay&Brenner hindrance');
model.param.set('a6', '-1.216', 'Bungay&Brenner hindrance');
model.param.set('a7', '1.647', 'Bungay&Brenner hindrance');

%% FUNCTIONS
model.func.create('an1', 'Analytic');
model.func('an1').set('funcname', 'trapezoidal');
model.func('an1').set('expr', 'flc1hs(t-trise/2,trise/2) - flc1hs(t-tpulse-tfall/2,tfall/2)');
model.func('an1').set('args', 't');
model.func('an1').set('fununit', '1');
model.func('an1').set('argunit', 's');
model.func('an1').set('periodic', true);
model.func('an1').set('periodicupper', '1/fp');

model.func.create('an2', 'Analytic');
model.func('an2').set('funcname', 'exponential');
model.func('an2').set('expr', 'exp(-t/tpulse)*flc1hs(t-trise/2,trise/2)');
model.func('an2').set('args', 't');
model.func('an2').set('fununit', '1');
model.func('an2').set('argunit', 's');
model.func('an2').set('periodic', true);
model.func('an2').set('periodicupper', '1/fp');

model.func.create('an0', 'Analytic');
model.func('an0').set('expr', '1-flc1hs(t-npulse*1/fp,0)');
model.func('an0').set('args', {'t'});
model.func('an0').set('argunit', 's');
model.func('an0').set('fununit', '1');

%% GEOMETRY
% ---------------------  Component 1 -------------------------------------
% Rectangular extracellular domain
model.geom('geom1').create('r1', 'Rectangle');
model.geom('geom1').feature('r1').set('size', {'del/2' 'del'});
model.geom('geom1').feature('r1').set('pos', {'0' '-del/2'});
% Cell membrane
model.geom('geom1').create('pol1', 'Polygon');
model.geom('geom1').feature('pol1').set('x', 'rcell*sin({range(0,dtheta,pi)})');
model.geom('geom1').feature('pol1').set('y', 'rcell*cos({range(0,dtheta,pi)})');

% ---------------------  Component 2 -------------------------------------
model.geom('geom2').create('pol1', 'Polygon');
model.geom('geom2').feature('pol1').set('x', 'rstar, rpmin, 2*rpmin, max(2*rpmin,rpmax-2[nm]), rpmax');
model.geom('geom2').feature('pol1').set('y', '0, 0, 0, 0, 0');
model.geom('geom2').feature('pol1').set('type', 'open');
model.geom('geom2').create('arr1', 'Array');
model.geom('geom2').feature('arr1').set('size', {'1' 'Nnodes'});
model.geom('geom2').feature('arr1').set('displ', {'0' 'dcomp2'});
model.geom('geom2').feature('arr1').selection('input').set({'pol1'});
model.geom('geom2').run;

%% SELECTIONS
% ---------------------  Component 1 -------------------------------------
% Indices of the domains
extdomain = 1;
intdomain = 2;

% Indices of anode/cathode
anode = min(mphselectcoords(model,'geom1',[par.del/2, -par.del/2],'boundary'));
cathode = min(mphselectcoords(model,'geom1',[par.del/2, par.del/2],'boundary'));
elinsul = max(mphselectcoords(model,'geom1',[par.del/2, par.del/2],'boundary'));

% Indices of membrane points and edges from theta = 0 to theta = pi
mempoints_raw = intersect(mphgetadj(model,'geom1','point','domain',intdomain),mphgetadj(model,'geom1','point','domain',extdomain));
mempoints_coords_raw = mphgetcoords(model,'geom1','point',mempoints_raw);
[theta_i,I] = sort(atan2(mempoints_coords_raw(1,:),mempoints_coords_raw(2,:)));
mempoints = mempoints_raw(:,I);
mempoints_coords = mempoints_coords_raw(:,I);
for i = 1:length(mempoints)-1
    memedges(1,i) = intersect(mphgetadj(model,'geom1','boundary','point',mempoints(i)),mphgetadj(model,'geom1','boundary','point',mempoints(i+1)));
end

model.selection.create('sel1', 'Explicit');
model.selection('sel1').model('comp1');
model.selection('sel1').geom('geom1', 1);
model.selection('sel1').set(memedges);
model.selection('sel1').label('Membrane');

% ---------------------  Component 2 -------------------------------------
Nnodes = par.Nnodes;
rpedges = [(1:Nnodes)' (Nnodes*1+1:Nnodes*2)' (Nnodes*2+1:Nnodes*3)' (Nnodes*3+1:Nnodes*4)'];
rppoints_coords = mphgetcoords(model,'geom2','point',1:Nnodes);

par.extdomain = extdomain;
par.intdomain = intdomain;
par.mempoints = mempoints;
par.memedges = memedges;

%% COUPLINGS
% ---------------------  Component 1 -------------------------------------
% Average over the cell membrane
model.cpl.create('aveop1', 'Average','geom1');
model.cpl('aveop1').selection.named('sel1');
model.cpl('aveop1').set('axisym', true);
% Average over the intracellular domain
model.cpl.create('aveop2', 'Average', 'geom1');
model.cpl('aveop2').selection.set(intdomain);
model.cpl('aveop2').set('axisym', true);
% Integration over the intracellular domain
model.cpl.create('intop0', 'Integration', 'geom1');
model.cpl('intop0').selection.set(intdomain);
model.cpl('intop0').set('axisym', true);

% ---------------------  Component 2 -------------------------------------
% Integration over each line in the pore radius space
for i = 1:Nnodes
    intopname = ['intop',num2str(i)];
    model.cpl.create(intopname, 'Integration', 'geom2');
    model.cpl(intopname).selection.geom('geom2', 1);
    model.cpl(intopname).selection.set(rpedges(i,:));
    model.cpl(intopname).set('intorder', '4');
end

%% VARIABLES
% ---------------------  Component 1 -------------------------------------
% Ext. medium
model.variable.create('var1');
model.variable('var1').label('Ext. medium');
model.variable('var1').model('comp1');
model.variable('var1').selection.geom('geom1', 2);
model.variable('var1').selection.set(extdomain);
model.variable('var1').set('Ve', 'V');
model.variable('var1').set('c1', 'c1_e');
if par.Nsolutes == 3
    model.variable('var1').set('c2', 'c2_e');
    model.variable('var1').set('c12', 'c12_e');
end 

% Cytoplasm
model.variable.create('var2');
model.variable('var2').label('Cytoplasm');
model.variable('var2').model('comp1');
model.variable('var2').selection.geom('geom1', 2);
model.variable('var2').selection.set(intdomain);
model.variable('var2').set('Vi', 'V2');
model.variable('var2').set('c1', 'c1_i');
if par.Nsolutes == 3
    model.variable('var2').set('c2', 'c2_i');
    model.variable('var2').set('c12', 'c12_i');
end

% Membrane
model.variable.create('var3');
model.variable('var3').label('Membrane');
model.variable('var3').model('comp1');
model.variable('var3').selection.named('sel1');
model.variable('var3').set('theta', 'atan2(r,z)');
model.variable('var3').set('Um', 'Vi-Ve +1[mV]'); % 1 mV is added to avoid numerical problems
model.variable('var3').set('sigma_tot', 'sigma_m + sigma_ep');
model.variable('var3').set('sigma_ep', 'Gep*dm');
model.variable('var3').set('Pe1', 'F/(R*T)*z1*Um');
if par.Nsolutes == 3
    model.variable('var3').set('Pe2', 'F/(R*T)*z2*Um');
    model.variable('var3').set('Pe12', 'F/(R*T)*z12*Um');
end
if strcmp(par.expr_molflux,'Li')
    % Expression from Li et al. BBA 2013 that takes into account the difference between intracellular and extracellular conductivity,
    % however it does not take into account the hindrance and partitioning factors. 
    model.variable('var3').set('Jc1', 'PAD*D1_m*(Pe1+log(sigma_rat))*(sigma_rat-1)*(c1_e-c1_i*exp(Pe1))/(dm*log(sigma_rat)*(1-sigma_rat*exp(Pe1)))');
    if par.Nsolutes == 3
        model.variable('var3').set('Jc2', 'PAD*D2_m*(Pe2+log(sigma_rat))*(sigma_rat-1)*(c2_e-c2_i*exp(Pe2))/(dm*log(sigma_rat)*(1-sigma_rat*exp(Pe2)))');
        model.variable('var3').set('Jc12', 'PAD*D12_m*(Pe12+log(sigma_rat))*(sigma_rat-1)*(c12_e-c12_i*exp(Pe12))/(dm*log(sigma_rat)*(1-sigma_rat*exp(Pe12)))');
    end
    % Note: The definition of Jc1, Jc2, and Jc12 depends on how you define the transmembrane voltage. 
    % Here we define it as Um = Vi - Ve, as was done in Li et al. BBA 2013 (http://dx.doi.org/10.1016/j.bbamem.2012.08.014). 
    % In their earlier paper from 2011, they defined Um = Ve - Vi; therefore, the expressions for Jci were different. 
elseif strcmp(par.expr_molflux,'Smith')
    % Expression from Smith, PhD thesis, 2011. This takes into account the hindrance and partitioning factors, 
    % but does not take into account the difference between intracellular and extracellular conductivity. 
    model.variable('var3').set('Jc1', 'PADa*D1_m*(Pe1/dm*(c1_e/(1-exp(Pe1) + (Pe1==0)) + c1_i/(1-exp(-Pe1) + (Pe1==0))) + (c1_i-c1_e)/dm*(Pe1==0))');
    if par.Nsolutes == 3
        model.variable('var3').set('Jc2', 'PADb*D2_m*(Pe2/dm*(c2_e/(1-exp(Pe2) + (Pe2==0)) + c2_i/(1-exp(-Pe2) + (Pe2==0))) + (c2_i-c2_e)/dm*(Pe2==0))');
        model.variable('var3').set('Jc12', 'PADc*D12_m*(Pe12/dm*(c12_e/(1-exp(Pe12) + (Pe12==0)) + c12_i/(1-exp(-Pe12) + (Pe12==0))) + (c12_i-c12_e)/dm*(Pe12==0))');
    end
elseif strcmp(par.expr_molflux,'Combination')
    % Expression that combines the two expressions above: This expression takes into account both the hindrance and partitioning factors, 
    % as well as the the difference between intracellular and extracellular conductivity.
    model.variable('var3').set('Jc1', 'PADa*D1_m*(Pe1+log(sigma_rat))*(sigma_rat-1)*(c1_e-c1_i*exp(Pe1))/(dm*log(sigma_rat)*(1-sigma_rat*exp(Pe1)))');
    if par.Nsolutes == 3
        model.variable('var3').set('Jc2', 'PADb*D2_m*(Pe2+log(sigma_rat))*(sigma_rat-1)*(c2_e-c2_i*exp(Pe2))/(dm*log(sigma_rat)*(1-sigma_rat*exp(Pe2)))');
        model.variable('var3').set('Jc12', 'PADc*D12_m*(Pe12+log(sigma_rat))*(sigma_rat-1)*(c12_e-c12_i*exp(Pe12))/(dm*log(sigma_rat)*(1-sigma_rat*exp(Pe12)))');
    end    
end

% Integration variables
model.variable.create('var4');
model.variable('var4').label('Integration variables');
model.variable('var4').model('comp1');
model.variable('var4').set('Gep_ave', 'aveop1(Gep)');
model.variable('var4').set('PADtot', 'aveop1(PAD)');
model.variable('var4').set('c1_avg', 'aveop2(c1)');
if par.Nsolutes == 3
    model.variable('var4').set('c2_avg', 'aveop2(c2)');
    model.variable('var4').set('c12_avg', 'aveop2(c12)');
end
    
% ---------------------  Component 2 -------------------------------------
% Pore creation
model.variable.create('var5');
model.variable('var5').label('Pore creation and expansion');
model.variable('var5').model('comp2');
model.variable('var5').set('Jpc', 'a*fc*exp(alpha*Um + beta*Um^2)', 'Pore creation rate');
model.variable('var5').set('Jpd', 'a*(1-fprot)/neq_star*exp((1-(rpmin/rstar)^2)*(alpha*Um + beta*Um^2))', 'Pore destruction rate');
model.variable('var5').set('fc', 'max(0,1 - fprot - fporated)', 'Fraction of membrane area available for pore creation');

if strcmp(par.expr_Wsurf,'Li')
    model.variable('var5').set('Gamma_eff', '2*Gamma1 - (2*Gamma1-Gamma0)/(1-PAD)^2');
    model.variable('var5').set('Wsurf', '-pi*rp^2*Gamma_eff');
    model.variable('var5').set('Fsurf', '-d(Wsurf,rp)');
elseif strcmp(par.expr_Wsurf,'Smith')
    model.variable('var5').set('Gamma_eff', '2*Gamma1 - (2*Gamma1-Gamma0)/(1-fAlp)^2');
    model.variable('var5').set('dAlp', 'pi*(rp + 0.5*dm)^2 - pi^2*0.5*dm*(rp + 0.5*dm) + 2*pi*(0.5*dm)^2');
    model.variable('var5').set('Wsurf', '-dAlp*Gamma_eff');
    model.variable('var5').set('Fsurf', '-d(Wsurf,rp)');
end    
model.variable('var5').set('Fsteric', 'B*b*(rstar[1/m])^b*(rp[1/m])^(-b-1)*1[1/m]');
model.variable('var5').set('Fedge', '-2*pi*gamma');
model.variable('var5').set('Fel', 'Fmax/(1+rh/(rp+rt))*Um^2');
model.variable('var5').set('Fdipole', 'alpha*k*T*Up/rd');
model.variable('var5').set('Fp', 'Fsteric + Fedge + Fsurf + Fel + Fdipole');

% Partitioning and hindrance
model.variable.create('var6');
model.variable('var6').label('Partitioning and hindrance');
model.variable('var6').model('comp2');
model.variable('var6').set('w0', '5.3643[1/F]*qe^2/(k*T)*(rp[1/m])^(-1.803)');
model.variable('var6').set('w1','(z1*min(1,dp/l1))^2*w0');
model.variable('var6').set('vm0', 'Up*F/(R*T)');
model.variable('var6').set('vm1', 'Up*z1*F/(R*T)');
model.variable('var6').set('K0', 'w0/((w0*(1-2*nu)+2*nu)*exp(w0)-2*nu)*(abs(Up)<1e-6)+(-1+exp(vm0))/((w0*exp(w0-nu*vm0)-nu*vm0)/(w0-nu*vm0)*exp(vm0)-(w0*exp(w0+nu*vm0)+nu*vm0)/(w0+nu*vm0)+(abs(Up)<1e-6))*(abs(Up)>=1e-6)');
model.variable('var6').set('K1', 'w1/((w1*(1-2*nu)+2*nu)*exp(w1)-2*nu)*(abs(Up)<1e-6)+(-1+exp(vm1))/((w1*exp(w1-nu*vm1)-nu*vm1)/(w1-nu*vm1)*exp(vm1)-(w1*exp(w1+nu*vm1)+nu*vm1)/(w1+nu*vm1)+(abs(Up)<1e-6))*(abs(Up)>=1e-6)');
model.variable('var6').set('lam0', 'r0/rp*(r0/rp < 1) + (r0/rp >= 1)');
model.variable('var6').set('fA0', '(1 - lam0)^2');
model.variable('var6').set('ft0', '9/4*pi^2*sqrt(2)*(1-lam0)^(-2.5)*(1 + a1*(1-lam0) + a2*(1-lam0)^2) +a3 + a4*lam0 + a5*lam0^2 + a6*lam0^3 + a7*lam0^4');
model.variable('var6').set('fD00', '6*pi/ft0');
model.variable('var6').set('fD0', 'fD00/(fD00 + (1-fD00)*(min(dp,l0)/(2*r0)))');
model.variable('var6').set('H0', 'fA0*fD0');
model.variable('var6').set('lam1', 'r1/rp*(r1/rp < 1) +(r1/rp >= 1)');
model.variable('var6').set('fA1', '(1 - lam1)^2');
model.variable('var6').set('ft1', '9/4*pi^2*sqrt(2)*(1-lam1 + (lam1==1))^(-2.5)*(1 + a1*(1-lam1) + a2*(1-lam1 + (lam1==1))^2) +a3 + a4*lam1 + a5*lam1^2 + a6*lam1^3 + a7*lam1^4');
model.variable('var6').set('fD10', '6*pi/ft1');
model.variable('var6').set('fD1', 'fD10/(fD10 + (1-fD10)*(min(dp,l1)/(2*r1)))');
model.variable('var6').set('H1', 'fA1*fD1');
if par.Nsolutes == 3
    model.variable('var6').set('w2','(z2*min(1,dp/l2))^2*w0');
    model.variable('var6').set('w12','(z12*min(1,dp/l12))^2*w0');
    model.variable('var6').set('vm2', 'Up*z2*F/(R*T)');
    model.variable('var6').set('vm12', 'Up*z12*F/(R*T)');
    model.variable('var6').set('K2', 'w2/((w2*(1-2*nu)+2*nu)*exp(w2)-2*nu)*(abs(Up)<1e-6)+(-1+exp(vm2))/((w2*exp(w2-nu*vm2)-nu*vm2)/(w2-nu*vm2)*exp(vm2)-(w2*exp(w2+nu*vm2)+nu*vm2)/(w2+nu*vm2)+(abs(Up)<1e-6))*(abs(Up)>=1e-6)');
    model.variable('var6').set('K12', 'w12/((w12*(1-2*nu)+2*nu)*exp(w12)-2*nu)*(abs(Up)<1e-6)+(-1+exp(vm12))/((w12*exp(w12-nu*vm12)-nu*vm12)/(w12-nu*vm12)*exp(vm12)-(w12*exp(w12+nu*vm12)+nu*vm12)/(w12+nu*vm12)+(abs(Up)<1e-6))*(abs(Up)>=1e-6)');
    model.variable('var6').set('lam2', 'r2/rp*(r2/rp < 1) + (r2/rp >= 1)');
    model.variable('var6').set('fA2', '(1 - lam2)^2');
    model.variable('var6').set('ft2', '9/4*pi^2*sqrt(2)*(1-lam2 + (lam2==1))^(-2.5)*(1 + a1*(1-lam2) + a2*(1-lam2 + (lam2==1))^2) +a3 + a4*lam2 + a5*lam2^2 + a6*lam2^3 + a7*lam2^4');
    model.variable('var6').set('fD20', '6*pi/ft2');
    model.variable('var6').set('fD2', 'fD20/(fD20 + (1-fD20)*(min(dp,l2)/(2*r2)))');
    model.variable('var6').set('H2', 'fA2*fD2');
    model.variable('var6').set('lam12', 'r12/rp*(r12/rp < 1) + (r12/rp >= 1)');
    model.variable('var6').set('fA12', '(1 - lam12)^2');
    model.variable('var6').set('ft12', '9/4*pi^2*sqrt(2)*(1-lam12 + (lam12==1))^(-2.5)*(1 + a1*(1-lam12) + a2*(1-lam12 + (lam12==1))^2) +a3 + a4*lam12 + a5*lam12^2 + a6*lam12^3 + a7*lam12^4');
    model.variable('var6').set('fD120', '6*pi/ft12');
    model.variable('var6').set('fD12', 'fD120/(fD120 + (1-fD120)*(min(dp,l12)/(2*r12)))');
    model.variable('var6').set('H12', 'fA12*fD12');
end

% Integration variables
model.variable.create('var7');
model.variable('var7').label('Integration variables');
model.variable('var7').model('comp2');
if strcmp(par.expr_Gp,'Li')
    model.variable('var7').set('Gp', '2*pi*rp^2*sigma_p/(pi*rp + 2*dm)');
elseif strcmp(par.expr_Gp,'Smith')
    model.variable('var7').set('Gp', '2*pi*rp^2*sigma_p*K0*H0/(pi*rp*K0*H0 + 2*dp)');
elseif strcmp(par.expr_Gp,'Toroidal')
    model.variable('var7').set('Gp', '2*pi*rp^2*sigma_p/(pi*rp + dm*(1 + tanh(0.1360*(rp[1/m]*1e9-1.2408))))');
end
% Transpore voltage
model.variable('var7').set('Gp0', '2*pi*rp^2*sigma_p/(pi*rp + 2*dp)');
model.variable('var7').set('Gpp', '(sigma_p*pi*rp^2)/dp');
model.variable('var7').set('Up', 'Gp0/Gpp*Um');

for i = 1:Nnodes
    model.variable('var7').set(['Gep',num2str(i,'%02.0f')], ['intop',num2str(i),'(n*Gp)']);
end
for i = 1:Nnodes
    model.variable('var7').set(['fporated',num2str(i,'%02.0f')], ['intop',num2str(i),'(n*pi*(rp + 0.5*dm)^2)']);
    model.variable('var7').set(['fAlp',num2str(i,'%02.0f')], ['intop',num2str(i),'(n*dAlp)']);
    model.variable('var7').set(['PAD',num2str(i,'%02.0f')], ['intop',num2str(i),'(n*pi*rp^2)']);
    model.variable('var7').set(['PADa',num2str(i,'%02.0f')], ['intop',num2str(i),'(n*pi*rp^2*K1*H1)']);
    if par.Nsolutes == 3
        model.variable('var7').set(['PADb',num2str(i,'%02.0f')], ['intop',num2str(i),'(n*pi*rp^2*K2*H2)']);
        model.variable('var7').set(['PADc',num2str(i,'%02.0f')], ['intop',num2str(i),'(n*pi*rp^2*K12*H12)']);
    end
end

% Common variables defined over seperate membrane edges/points
% ---------------------  Component 1 -------------------------------------
lastvaridx = 7;
r1 = mempoints_coords(1,:);
y1 = mempoints_coords(2,:);

for i = 1:Nnodes-1
    varname = ['var',num2str(lastvaridx+i)];
    model.variable.create(varname);
    model.variable(varname).model('comp1');
    ltheta = sqrt((y1(i+1)-y1(i))^2 + (r1(i+1)-r1(i))^2);
    % Values of Gep and PAD along the membrane between the nodes - linear interpolation 
    model.variable(varname).set('fporated', ['(comp2.fporated',num2str(i+1,'%02d'),' -comp2.fporated',num2str(i,'%02d'),')/(',num2str(ltheta),'*1[m])*sqrt((r - ',num2str(r1(i)),')^2 + (z - ',num2str(y1(i)),')^2) + comp2.fporated',num2str(i,'%02d')]);
    model.variable(varname).set('Gep', ['(comp2.Gep',num2str(i+1,'%02d'),' -comp2.Gep',num2str(i,'%02d'),')/(',num2str(ltheta),'*1[m])*sqrt((r - ',num2str(r1(i)),')^2 + (z - ',num2str(y1(i)),')^2) + comp2.Gep',num2str(i,'%02d')]);
    model.variable(varname).set('PAD', ['(comp2.PAD',num2str(i+1,'%02d'),' -comp2.PAD',num2str(i,'%02d'),')/(',num2str(ltheta),'*1[m])*sqrt((r - ',num2str(r1(i)),')^2 + (z - ',num2str(y1(i)),')^2) + comp2.PAD',num2str(i,'%02d')]);
    model.variable(varname).set('PADa', ['(comp2.PADa',num2str(i+1,'%02d'),' -comp2.PADa',num2str(i,'%02d'),')/(',num2str(ltheta),'*1[m])*sqrt((r - ',num2str(r1(i)),')^2 + (z - ',num2str(y1(i)),')^2) + comp2.PADa',num2str(i,'%02d')]);
    if par.Nsolutes == 3
        model.variable(varname).set('PADb', ['(comp2.PADb',num2str(i+1,'%02d'),' -comp2.PADb',num2str(i,'%02d'),')/(',num2str(ltheta),'*1[m])*sqrt((r - ',num2str(r1(i)),')^2 + (z - ',num2str(y1(i)),')^2) + comp2.PADb',num2str(i,'%02d')]);
        model.variable(varname).set('PADc', ['(comp2.PADc',num2str(i+1,'%02d'),' -comp2.PADc',num2str(i,'%02d'),')/(',num2str(ltheta),'*1[m])*sqrt((r - ',num2str(r1(i)),')^2 + (z - ',num2str(y1(i)),')^2) + comp2.PADc',num2str(i,'%02d')]);
    end
    model.variable(varname).selection.geom('geom1', 1); 
    model.variable(varname).selection.set(memedges(i));
end
clear r1 y1
lastvaridx = lastvaridx + Nnodes-1;

% ---------------------  Component 2 -------------------------------------
for i = 1:Nnodes
    varname = ['var',num2str(lastvaridx+i)];
    model.variable.create(varname);
    model.variable(varname).model('comp2');
    model.variable(varname).set('Um', ['comp1.at0(',num2str(mempoints_coords(1,i)),',',num2str(mempoints_coords(2,i)),',comp1.Um)']);
    model.variable(varname).set('Gep', ['Gep',num2str(i,'%02.0f')]); 
    model.variable(varname).set('fporated', ['fporated',num2str(i,'%02.0f')]);
    model.variable(varname).set('fAlp', ['fAlp',num2str(i,'%02.0f')]);
    model.variable(varname).set('PAD', ['PAD',num2str(i,'%02.0f')]);
    model.variable(varname).set('PADa', ['PADa',num2str(i,'%02.0f')]);
    if par.Nsolutes == 3
        model.variable(varname).set('PADb', ['PADb',num2str(i,'%02.0f')]);
        model.variable(varname).set('PADc', ['PADc',num2str(i,'%02.0f')]);
    end
    model.variable(varname).selection.geom('geom2', 1);
    model.variable(varname).selection.set(rpedges(i,:));   
end
%% MATERIALS
model.material.create('mat1', 'Common', 'comp1');
model.material.create('mat2', 'Common', 'comp1');
model.material.create('mat3', 'Common', 'comp1');
model.material('mat1').selection.set(extdomain);
model.material('mat2').selection.set(intdomain);
model.material('mat3').selection.named('sel1');
model.material('mat1').propertyGroup('def').set('electricconductivity', {'sigma_e' '0' '0' '0' 'sigma_e' '0' '0' '0' 'sigma_e'});
model.material('mat1').propertyGroup('def').set('relpermittivity', {'epsilon_e' '0' '0' '0' 'epsilon_e' '0' '0' '0' 'epsilon_e'});
model.material('mat2').propertyGroup('def').set('electricconductivity', {'sigma_i' '0' '0' '0' 'sigma_i' '0' '0' '0' 'sigma_i'});
model.material('mat2').propertyGroup('def').set('relpermittivity', {'epsilon_i' '0' '0' '0' 'epsilon_i' '0' '0' '0' 'epsilon_i'});
model.material('mat3').propertyGroup('def').set('electricconductivity', {'sigma_tot' '0' '0' '0' 'sigma_tot' '0' '0' '0' 'sigma_tot'});
model.material('mat3').propertyGroup('def').set('relpermittivity', {'epsilon_m' '0' '0' '0' 'epsilon_m' '0' '0' '0' 'epsilon_m'});

%% PHYSICS
% Electric Currents
model.physics.create('ec', 'ConductiveMedia', 'geom1');
model.physics('ec').selection.set(extdomain);
model.physics('ec').prop('ShapeProperty').set('order_electricpotential', '2');
model.physics('ec').create('pot1', 'ElectricPotential', 1);
model.physics('ec').feature('pot1').selection.set(anode);
model.physics('ec').feature('pot1').set('V0', ['Uapp*',par.pulseshape,'(t)*an0(t)']);
model.physics('ec').create('pot2', 'ElectricPotential', 1);
model.physics('ec').feature('pot2').selection.set(cathode);
model.physics('ec').create('dimp1', 'DistributedImpedance', 1);
model.physics('ec').feature('dimp1').selection.named('sel1');
model.physics('ec').feature('dimp1').set('Vref', 'V2');
model.physics('ec').feature('dimp1').set('ds', 'dm');
model.physics('ec').create('ncd1', 'NormalCurrentDensity', 1);
model.physics('ec').feature('ncd1').selection.named('sel1');
model.physics('ec').feature('ncd1').set('nJ', '-Gm*Ur');
if par.axi2D == false
    model.physics('ec').prop('d').set('d', 'dcell');
end
model.physics.create('ec2', 'ConductiveMedia', 'geom1');
model.physics('ec2').selection.set(intdomain);
model.physics('ec2').feature('init1').set('V2', 'Urest');
model.physics('ec2').create('dimp2', 'DistributedImpedance', 1);
model.physics('ec2').feature('dimp2').selection.named('sel1');
model.physics('ec2').feature('dimp2').set('Vref', 'V');
model.physics('ec2').feature('dimp2').set('ds', 'dm');
model.physics('ec2').create('ncd2', 'NormalCurrentDensity', 1);
model.physics('ec2').feature('ncd2').selection.named('sel1');
model.physics('ec2').feature('ncd2').set('nJ', 'Gm*Ur');

% CoefficientFormBoundaryPDE
model.physics.create('cb', 'CoefficientFormBoundaryPDE', 'geom2');
model.physics('cb').field('dimensionless').field('n');
model.physics('cb').field('dimensionless').component({'n'});
model.physics('cb').prop('Units').set('DependentVariableQuantity', 'none');
model.physics('cb').prop('Units').set('CustomDependentVariableUnit', '1/m^3');
model.physics('cb').prop('Units').set('CustomSourceTermUnit', '1/(m^3*s)');
model.physics('cb').feature('init1').set('n', 'neq_star*exp(-1/(k*T)*(B*(rstar/rp)^b+ 2*pi*gamma*rp) + (B+2*pi*gamma*rstar)/(k*T))');
model.physics('cb').create('flux1', 'FluxBoundary', 0);
model.physics('cb').feature('flux1').selection.set(1:Nnodes);
model.physics('cb').create('zflx1', 'ZeroFluxBoundary', 0);
model.physics('cb').feature('zflx1').selection.set(4*Nnodes+1:5*Nnodes);
model.physics('cb').feature('cfeq1').set('c', {'Dp' '0' '0' 'Dp'});
model.physics('cb').feature('cfeq1').set('f', '0');
model.physics('cb').feature('cfeq1').set('al', {'-Fp*Dp/(k*T)' '0'});
model.physics('cb').feature('flux1').set('g', 'Jpc');
model.physics('cb').feature('flux1').set('q', 'Jpd');

% Transport of diluted species #1
model.physics.create('tds', 'DilutedSpecies', 'geom1');
model.physics('tds').field('concentration').field('c1_e');
model.physics('tds').selection.set(extdomain);
model.physics('tds').prop('ShapeProperty').set('order_concentration', '2');
model.physics('tds').prop('TransportMechanism').set('Migration', '1');
model.physics('tds').feature('cdm1').set('V_src', 'root.comp1.V');
model.physics('tds').feature('cdm1').set('minput_temperature', 'T');
model.physics('tds').create('conc1', 'Concentration', 1);
model.physics('tds').feature('conc1').selection.set([anode cathode elinsul]);
model.physics('tds').create('fl1', 'Fluxes', 1);
model.physics('tds').feature('fl1').selection.named('sel1');
if par.Nsolutes == 3
    model.physics('tds').field('concentration').component({'c1_e' 'c2_e' 'c12_e'});
    model.physics('tds').feature('cdm1').set('z', {'z1'; 'z2'; 'z12'});
    model.physics('tds').feature('cdm1').set('D_c1_e', {'D1_e'; '0'; '0'; '0'; 'D1_e'; '0'; '0'; '0'; 'D1_e'});
    model.physics('tds').feature('cdm1').set('D_c2_e', {'D2_e'; '0'; '0'; '0'; 'D2_e'; '0'; '0'; '0'; 'D2_e'});
    model.physics('tds').feature('cdm1').set('D_c12_e', {'D12_e'; '0'; '0'; '0'; 'D12_e'; '0'; '0'; '0'; 'D12_e'});
    model.physics('tds').create('reac1', 'Reactions', 2);
    model.physics('tds').feature('reac1').selection.set(extdomain);
    model.physics('tds').feature('reac1').set('R_c1_e', '-kass*c2_e*c1_e + kdis*c12_e');
    model.physics('tds').feature('reac1').set('R_c2_e', '-kass*c2_e*c1_e + kdis*c12_e');
    model.physics('tds').feature('reac1').set('R_c12_e', '+kass*c2_e*c1_e - kdis*c12_e');
    model.physics('tds').feature('init1').set('initc', {'c1_e0'; 'c2_e0'; 'c12_e0'});
    model.physics('tds').feature('conc1').set('species', {'1'; '1'; '1'});
    model.physics('tds').feature('conc1').set('c0', {'c1_e0'; 'c2_e0'; 'c12_e0'});
    model.physics('tds').feature('fl1').set('species', {'1'; '1'; '1'});
    model.physics('tds').feature('fl1').set('N0', {'Jc1'; 'Jc2'; 'Jc12'});
elseif par.Nsolutes == 1
    model.physics('tds').field('concentration').component({'c1_e'});
    model.physics('tds').feature('cdm1').set('z', {'z1'});
    model.physics('tds').feature('cdm1').set('D_c1_e', {'D1_e'; '0'; '0'; '0'; 'D1_e'; '0'; '0'; '0'; 'D1_e'});
    model.physics('tds').feature('init1').set('initc', {'c1_e0'});
    model.physics('tds').feature('conc1').set('species', {'1'});
    model.physics('tds').feature('conc1').set('c0', {'c1_e0'});
    model.physics('tds').feature('fl1').set('species', {'1'});
    model.physics('tds').feature('fl1').set('N0', {'Jc1'});    
end

% Transport of diluted species #2
model.physics.create('tds2', 'DilutedSpecies', 'geom1');
model.physics('tds2').field('concentration').field('c1_i');
model.physics('tds2').selection.set(intdomain);
model.physics('tds2').prop('ShapeProperty').set('order_concentration', '2');
model.physics('tds2').prop('TransportMechanism').set('Migration', '1');
model.physics('tds2').feature('cdm1').set('V_src', 'root.comp1.V2');
model.physics('tds2').feature('cdm1').set('minput_temperature', 'T');
model.physics('tds2').create('fl2', 'Fluxes', 1);
model.physics('tds2').feature('fl2').selection.named('sel1');
if par.Nsolutes == 3
    model.physics('tds2').field('concentration').component({'c1_i' 'c2_i' 'c12_i'});
    model.physics('tds2').feature('cdm1').set('z', {'z1'; 'z2'; 'z12'});
    model.physics('tds2').feature('cdm1').set('D_c1_i', {'D1_i'; '0'; '0'; '0'; 'D1_i'; '0'; '0'; '0'; 'D1_i'});
    model.physics('tds2').feature('cdm1').set('D_c2_i', {'D2_i'; '0'; '0'; '0'; 'D2_i'; '0'; '0'; '0'; 'D2_i'});
    model.physics('tds2').feature('cdm1').set('D_c12_i', {'D12_i'; '0'; '0'; '0'; 'D12_i'; '0'; '0'; '0'; 'D12_i'});
    model.physics('tds2').feature('init1').set('initc', {'c1_i0'; 'c2_i0'; 'c12_i0'});
    model.physics('tds2').create('reac1', 'Reactions', 2);
    model.physics('tds2').feature('reac1').selection.set(intdomain);
    model.physics('tds2').feature('reac1').set('R_c1_i', '-kass*c2_i*c1_i + kdis*c12_i');
    model.physics('tds2').feature('reac1').set('R_c2_i', '-kass*c2_i*c1_i + kdis*c12_i');
    model.physics('tds2').feature('reac1').set('R_c12_i', '+kass*c2_i*c1_i - kdis*c12_i');
    model.physics('tds2').feature('fl2').set('species', {'1'; '1'; '1'});
    model.physics('tds2').feature('fl2').set('N0', {'-Jc1'; '-Jc2'; '-Jc12'});
elseif par.Nsolutes == 1
    model.physics('tds2').field('concentration').component({'c1_i'});
    model.physics('tds2').feature('cdm1').set('z', {'z1'});
    model.physics('tds2').feature('cdm1').set('D_c1_i', {'D1_i'; '0'; '0'; '0'; 'D1_i'; '0'; '0'; '0'; 'D1_i'});
    model.physics('tds2').feature('init1').set('initc', {'c1_i0'});
    model.physics('tds2').feature('fl2').set('species', {'1'});
    model.physics('tds2').feature('fl2').set('N0', {'-Jc1'});
end
%% MESH
model.mesh('mesh1').create('ftri1', 'FreeTri');
model.mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.mesh('mesh1').feature('ftri1').feature('size1').selection.geom('geom1', 2);
model.mesh('mesh1').feature('ftri1').feature('size1').selection.set(intdomain);
model.mesh('mesh1').feature('ftri1').feature('size1').set('custom', 'on');
model.mesh('mesh1').feature('ftri1').feature('size1').set('hmaxactive', true);
model.mesh('mesh1').feature('ftri1').feature('size1').set('hmax', num2str(par.rmesh0));
model.mesh('mesh1').run;

model.mesh('mesh2').create('edg1', 'Edge');
model.mesh('mesh2').feature('edg1').selection.geom('geom2');
model.mesh('mesh2').feature('edg1').create('size1', 'Size');
model.mesh('mesh2').feature('edg1').create('size2', 'Size');
model.mesh('mesh2').feature('edg1').create('size3', 'Size');
model.mesh('mesh2').feature('edg1').feature('size1').selection.geom('geom2', 1);
model.mesh('mesh2').feature('edg1').feature('size2').selection.geom('geom2', 1);
model.mesh('mesh2').feature('edg1').feature('size3').selection.geom('geom2', 0);
select = rpedges(:,3);
model.mesh('mesh2').feature('edg1').feature('size1').selection.set(select(:));
model.mesh('mesh2').feature('edg1').feature('size1').set('custom', 'on');
model.mesh('mesh2').feature('edg1').feature('size1').set('hmaxactive', true);
model.mesh('mesh2').feature('edg1').feature('size1').set('hmax', num2str(par.rmesh1));
select = rpedges(:,[2 4]);
model.mesh('mesh2').feature('edg1').feature('size2').selection.set(select(:));
model.mesh('mesh2').feature('edg1').feature('size2').set('custom', 'on');
model.mesh('mesh2').feature('edg1').feature('size2').set('hmaxactive', true);
model.mesh('mesh2').feature('edg1').feature('size2').set('hmax', num2str(par.rmesh2));
model.mesh('mesh2').feature('edg1').feature('size2').set('hgradactive', true);
model.mesh('mesh2').feature('edg1').feature('size2').set('hgrad', 1.1);
select = 1:Nnodes;
model.mesh('mesh2').feature('edg1').feature('size3').selection.set(select(:));
model.mesh('mesh2').feature('edg1').feature('size3').set('custom', 'on');
model.mesh('mesh2').feature('edg1').feature('size3').set('hmaxactive', true);
model.mesh('mesh2').feature('edg1').feature('size3').set('hmax', num2str(par.rmesh3));
model.mesh('mesh2').feature('edg1').feature('size3').set('hgradactive', true);
model.mesh('mesh2').feature('edg1').feature('size3').set('hgrad', 1.1);
model.mesh('mesh2').run;

% Study 1: Electroporation
model.study.create('std1');
model.study('std1').create('time', 'Transient');
model.study('std1').feature('time').set('activate', {'ec' 'on' 'ec2' 'on' 'cb' 'on' 'tds' 'off' 'tds2' 'off' 'frame:material1' 'on' 'frame:material2' 'on'});
model.study('std1').feature('time').set('tlist', par.tlist1);
model.study('std1').feature('time').set('usertol', true);
model.study('std1').feature('time').set('rtol', 0.001);
model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature.remove('fcDef');
model.sol('sol1').feature('t1').set('tlist', par.tlist1);
model.sol('sol1').feature('t1').set('rtol', 0.001);
model.sol('sol1').feature('t1').set('tstepsbdf', 'strict');
model.sol('sol1').feature('t1').set('initialstepbdfactive', true);
model.sol('sol1').feature('t1').set('initialstepbdf', '1e-11');

% Study 2: Molecular transport
model.study.create('std2');
model.study('std2').create('time', 'Transient');
model.study('std2').feature('time').set('activate', {'ec' 'off' 'ec2' 'off' 'cb' 'off' 'tds' 'on' 'tds2' 'on' 'frame:material1' 'on' 'frame:material2' 'on'});
model.study('std2').feature('time').set('tlist', par.tlist2);
model.study('std2').feature('time').set('useinitsol', true);
model.study('std2').feature('time').set('usesol', true);
model.study('std2').feature('time').set('notsolmethod', 'sol');
model.study('std2').feature('time').set('notstudy', 'std1');
model.study('std2').feature('time').set('notsolnum', 'all');
model.sol.create('sol2');
model.sol('sol2').study('std2');
model.sol('sol2').attach('std2');
model.sol('sol2').create('st1', 'StudyStep');
model.sol('sol2').create('v1', 'Variables');
model.sol('sol2').create('t1', 'Time');
model.sol('sol2').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol2').feature('t1').create('d1', 'Direct');
model.sol('sol2').feature('t1').feature.remove('fcDef');
model.sol('sol2').feature('v1').set('notsolmethod', 'sol');
model.sol('sol2').feature('v1').set('notsol', 'sol1');
model.sol('sol2').feature('v1').set('notsolnum', 'all');
model.sol('sol2').feature('t1').set('tlist', par.tlist2);
model.sol('sol2').feature('t1').set('tstepsbdf', 'strict');
model.sol('sol2').feature('t1').set('initialstepbdfactive', true);
model.sol('sol2').feature('t1').set('initialstepbdf', '1e-11');

%% RESULTS
model.result.dataset.create('rev1', 'Revolve2D');
model.result.dataset.create('cpl1', 'CutPlane');
model.result.dataset('rev1').set('data', 'dset3');
model.result.dataset('rev1').set('revangle', 180);
model.result.dataset('cpl1').set('quickplane', 'zx');

%%
% out = model;
fprintf('Model created ....................................................\n')
