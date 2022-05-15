%% Extract data
% This is an example script that shows how to extract data from the saved mph files. 

%% Extract data from simulations of Puc et al. experiments with S2011 model. 
close all

EF0 = (0.4:0.2:2.0)*1e5;
tp_v = [1e-4*ones(size(EF0)), 1e-3*ones(size(EF0))];
EF_v = [EF0 EF0];
foldername = 'Puc_S2011';

clear t c1_avg
for ii = 1:length(tp_v)
    tp = tp_v(ii);
    EF = EF_v(ii);
    filename = ['tp',num2str(tp*1e6),'us_EF',num2str(EF./1e3),'Vmm','.mph'];
    filepath = [foldername,'/',filename];
    if isfile(filepath)
        model = mphload([filepath]);
        tt = mpheval(model,'t','edim',0,'selection',1,'dataset','dset3','solnum','all','dataonly','on');
        cc = mpheval(model,'c1_avg','edim',0,'selection',1,'dataset','dset3','solnum','all','dataonly','on');
        if (round(tt(end)) == 200)
            t(ii) = tt(end);
            c1_avg(ii) = cc(end);
        else % If the simulation failed and did not reach the end
            t(ii) = nan;
            c1_avg(ii) = nan;
        end
    else
        t(ii) = nan;
        c1_avg(ii) = nan;
    end
end

% Plot
figure; hold on; box on
h1 = plot(EF0*1e-5,c1_avg(1:9),'-','LineWidth',1.5);
h2 = plot(EF0*1e-5,c1_avg(10:18),'-','LineWidth',1.5);
set(gca,'FontSize',14)
xlim([0 2.2])
ylim([0 0.2])
xlabel('Eapp (kV/cm)')
ylabel('Concentration (mM)')
legend('100 \mus','1 ms')
title('Lucifer Yellow')

%% Extract data from simulations of Gabriel and Teisse experiment with S2011 model. 
% Plot 2D concentration profiles
close all
foldername = 'Gabriel_LL2011';
filename = 'LL2011.mph';
filepath = [foldername,'/',filename];

Parameters_ExpGabriel1999
Parameters_ModS2011_LL2011

model = mphload(filepath);
get_Gabriel_plots(model,par)
