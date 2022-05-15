%% Main script for running the simulations 
% This is an example script that shows how to simulate selected experiments. 

% This script needs to be run from Comsol with Matlab!

% It is suggested not to run this entire script at once. A single simulation
% for one pulse duration and amplitude can take >1 hour. Thus, running the 
% entire script could take days. 
% Instead, select the simulations you would like to run, make your own 
% script based on the examples below and run the simulations from your
% script. 

% Import Comsol utilities 
import com.comsol.model.*
import com.comsol.model.util.*
ModelUtil.showProgress(1);

%% Experiment by Puc et al. - Calculations with S2011 model

% Import parameters
Parameters_ExpPuc2003
Parameters_ModS2011
modname = 'ModS2011';
foldername = 'Puc_S2011';
mkdir(foldername)

EF0 = (0.4:0.2:2.0)*1e5;
tp_v = [1e-4*ones(size(EF0)), 1e-3*ones(size(EF0))];
EF_v = [EF0 EF0];

% Create models and run simulations
for ii = 1%:length(tp_v)
    tp = tp_v(ii);
    EF = EF_v(ii);
    par.Eapp = EF;
    par.tpulse = tp;
    
    filename = ['tp',num2str(tp*1e6),'us_EF',num2str(EF./1e3),'Vmm','.mph'];
    filepath = [foldername,'/',filename];
    if ~isfile(filepath)
        fprintf(filename);
        fprintf('\n');
        % Create model
        run(modname)
        
        % Run model
        model.sol('sol1').runAll;
        model.sol('sol2').runAll;
        
        % Save model
        mphsave(model,filepath);
    end
end

%% Experiment by Puc et al. - Calculations with LL2011 model
% Note that LL2011 model is equivalent to S2011 model, except that it uses
% different values of model parameters and a few different expressions. 

% Import parameters
Parameters_ExpPuc2003
Parameters_ModS2011_LL2011
modname = 'ModS2011';
foldername = 'Puc_LL2011';
mkdir(foldername)

EF0 = (0.4:0.2:2.0)*1e5;
tp_v = [1e-4*ones(size(EF0)), 1e-3*ones(size(EF0))];
EF_v = [EF0 EF0];

% Create models and run simulations
for ii = 1:length(tp_v)
    tp = tp_v(ii);
    EF = EF_v(ii);
    par.Eapp = EF;
    par.tpulse = tp;
    
    filename = ['tp',num2str(tp*1e6),'us_EF',num2str(EF./1e3),'Vmm','.mph'];
    filepath = [foldername,'/',filename];
    if ~isfile(filepath)
        fprintf(filename);
        fprintf('\n');
        % Create model
        run(modname)
        
        % Run model
        model.sol('sol1').runAll;
        model.sol('sol2').runAll;
        
        % Save model
        mphsave(model,filepath);
    end
end

%% Experiment by Canatella et al. - Calculations with S2011 model

% Import parameters
Parameters_ExpCanatella2001
Parameters_ModS2011
modname = 'ModS2011';
foldername = 'Canatella_S2011';
mkdir(foldername)

% All combinations of pulse parameters
tp_v = [0.05, 0.05, 0.05, 0.09, 0.09, 0.09, 1.10, 1.10, 1.10, 1.10, 2.80, 2.80, 2.80, 2.80, 2.80, 2.80, ...
        5.30, 5.30, 5.30, 5.30, 5.30, 5.30, 10.0, 10.0, 10.0, 10.0, 10.0, 21.0, 21.0, 21.0, 21.0, 21.0, 10.0, 21.0]*1e-3;
EF_v = [1.20, 2.03, 2.85, 1.37, 2.04, 3.16, 0.28, 0.93, 1.22, 1.81, 0.37, 0.45, 0.55, 0.77, 1.09, 1.35, ...
        0.37, 0.45, 0.55, 0.64, 0.78, 0.95, 0.37, 0.40, 0.45, 0.55, 0.69, 0.28, 0.35, 0.40, 0.44, 0.58, 0.94, 0.79]*1e5;

% Create models and run simulations
for ii = 1:length(tp_v)
    tp = tp_v(ii);
    EF = EF_v(ii);
    par.Eapp = EF;
    par.tpulse = tp;
    
    filename = ['tp',num2str(tp*1e6),'us_EF',num2str(EF./1e3),'Vmm','.mph'];
    filepath = [foldername,'/',filename];
    if ~isfile(filepath)
        fprintf(filename);
        fprintf('\n');
        % Create model
        run(modname)
        
        % Run model
        model.sol('sol1').runAll;
        model.sol('sol2').runAll;
        
        % Save model
        mphsave(model,filepath);
    end
end

%% Experiment by Sozer et al., ns pulses - Calculations with S2011 model
for ii = 1:4
    switch ii
        case 1
            Parameters_ExpSozer2017_1x6ns_Calcein
            filename = 'Calcein_1x6ns.mph';
        case 2
            Parameters_ExpSozer2017_1x6ns_Propidium
            filename = 'Propidium_1x6ns.mph';
        case 3
            Parameters_ExpSozer2017_10x6ns_Calcein
            filename = 'Calcein_10x6ns.mph';
        case 4    
            Parameters_ExpSozer2017_10x6ns_Propidium
            filename = 'Propidium_10x6ns.mph';
    end
    
    % Import parameters related to the electroporation model
    Parameters_ModS2011
    modname = 'ModS2011';
    foldername = 'Sozer_S2011';
    mkdir(foldername)
    
    % Adapt the pore destruction rate to meet the experiment
    par.taup = 25;
    % Adapt discretization parameters to improve accuracy for ns pulses 
    par.rmesh1 = 1e-10;
    par.rmesh2 = 1e-11;
    par.Ntheta = 80;

    filepath = [foldername,'/',filename];
    if ~isfile(filepath)
        fprintf(filename);
        fprintf('\n');
        % Create model
        run(modname)
        
        % Run model
        model.sol('sol1').runAll;
        model.sol('sol2').runAll;
        
        % Save model
        mphsave(model,filepath);
    end
end

%% Experiment by Gabriel and Teissie - Calculations with LL2011 model

% Import parameters 
Parameters_ExpGabriel1999
Parameters_ModS2011_LL2011
modname = 'ModS2011';
foldername = 'Gabriel_LL2011';
mkdir(foldername)

filename = 'LL2011.mph';
filepath = [foldername,'/',filename];
if ~isfile(filepath)
    % Create model
    run(modname)
    
    % Run model
    model.sol('sol1').runAll;
    model.sol('sol2').runAll;
    
    % Save model
    mphsave(model,filepath);
end