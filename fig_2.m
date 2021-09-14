%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds the fitness gradient of a model over varying sex ratio. 
% Generates Figure 2 and figure S.1 by changing setupPars.experiment
% Tested 31/08/21 using MATLAB 2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
clear
close all

%% Declaring variables and packing into a structure
b = 25;
q = 0.5;
s = 0.5;
delta = 1;
d = 0.1;

uMJ = 0.01;
uFJ = 0.05;
uMA = 0.05;
uFA = 0.05;
tau = 0.1;

k = who;
j = cell(length(k),1);
for i = 1 : length(k)
    j{i,1} = eval(k{i});
end
pars = cell2struct(j,k,1);

clearvars -except 'pars'

%% Declaring options
% Running parameters
tspan = [0,1000];
tol = 1e-4;
searchCount = 100;

% Options
setupPars.tol = tol;
setupPars.searchCount = searchCount;
setupPars.tspan = tspan;
setupPars.solver = 'MEX';
plotting = 1;

% Select which experiment to run
% setupPars.experiment = 'uMJAgeUnions';
setupPars.experiment = 'uMAAgeUnions';
% setupPars.experiment = 'uMAAgeUnionsFrequency';

% Pre-define arrays to vary parameters over
sArray = 0.01 : 0.01 : 0.99;

% Sets up algebraic fitness values
if strcmp(setupPars.experiment,'uMJAgeUnions')
    run('symbolicSetupAgeUnions.m');
    popini = [0.2,0.2,0.2,0.2,0.1];
    dwfun = @(MJ,FJ,MA,FA,U,s,uMJ,uFJ) dw(MJ,FJ,MA,FA,U,pars.b,pars.q,s,uMJ,uFJ,pars.uMA,pars.uFA,pars.tau,pars.delta,pars.d);
    uMJArray = 0.5;
    dwArray = NaN(length(uMJArray),length(sArray));
    peak = NaN(1,length(uMJArray));
elseif strcmp(setupPars.experiment,'uMAAgeUnions')
    run('symbolicSetupAgeUnions.m');
    popini = [0.2,0.2,0.2,0.2,0.1];
    dwfun = @(MJ,FJ,MA,FA,U,s,uMA,uFA) dw(MJ,FJ,MA,FA,U,pars.b,pars.q,s,pars.uMJ,pars.uFJ,uMA,uFA,pars.tau,pars.delta,pars.d);
    uMAArray = 0.5;
    dwArray = NaN(length(uMAArray),length(sArray));
    peak = NaN(1,length(uMAArray));
elseif strcmp(setupPars.experiment,'uMAAgeUnionsFrequency')
    run('symbolicSetupAgeUnionsFrequency.m');
    popini = [0.2,0.2,0.2,0.2,0.1];
    wfun = @(MJ,FJ,MA,FA,U,s,smut,uMA,uFA) w(MJ,FJ,MA,FA,U,pars.b,s,smut,pars.uMJ,pars.uFJ,uMA,uFA,pars.tau,pars.delta,pars.d);
    dwfun = @(MJ,FJ,MA,FA,U,s,uMA,uFA) dw(MJ,FJ,MA,FA,U,pars.b,s,pars.uMJ,pars.uFJ,uMA,uFA,pars.tau,pars.delta,pars.d);
    uMAArray = 0.25;
    dwArray = NaN(length(uMAArray),length(sArray));
    peak = NaN(1,length(uMAArray));
end

% Select which experiment to run
Run_uMJ = 0;
Run_uMA = 0;
Run_uMAFrequency = 0;
if strcmp(setupPars.experiment,'uMJAgeUnions')
    Run_uMJ = 1;
elseif strcmp(setupPars.experiment,'uMAAgeUnions')
    Run_uMA = 1;
elseif strcmp(setupPars.experiment,'uMAAgeUnionsFrequency')
    Run_uMAFrequency = 1;
end

%% Finding fitness gradient
% Finds the fitness gradient whilst varying juvenile mortality
if Run_uMJ
    setupPars.tol = tol;
    setupPars.searchCount = searchCount;

    for i = 1 : length(uMJArray)
        pars.uMJ = uMJArray(i);
        
        dwArray(i,:) = calculateFitness(pars,setupPars,sArray,popini,dwfun);
        disp(100*i/length(uMJArray));
    end
end

% Finds the fitness gradient whilst varying adult mortality
if Run_uMA
    setupPars.tol = tol;
    setupPars.searchCount = searchCount;

    for i = 1 : length(uMAArray)
        pars.uMA = uMAArray(i);
        
        dwArray(i,:) = calculateFitness(pars,setupPars,sArray,popini,dwfun);
        disp(100*i/length(uMAArray));
    end
end

% Finds the fitness gradient whilst varying adult mortality in the 
% frequency base model 
if Run_uMAFrequency
    setupPars.tol = tol;
    setupPars.searchCount = searchCount;
    changeInSign = NaN(1,length(uMAArray));
    
    for i = 1 : length(uMAArray)
        pars.uMA = uMAArray(i);
        
        dwArray(i,:) = calculateFitness(pars,setupPars,sArray,popini,dwfun);
        disp(100*i/length(uMAArray));
    end
end

%% Plotting

if plotting
    figure()
    hold on
    grid on
    box on

    axis([0,1,-5,5])
%     midpoint = ceil(size(dwArray,1)/2)-5;
%     plot(sArray,dwArray(midpoint,:),'LineWidth',2);
    plot(sArray,dwArray,'LineWidth',2);
    
    scatter(0.5,0,'xk');
    
    markers1 = linspace(0,0.5,100);
    plot(markers1,zeros(length(markers1),1),'--k');
    markers2 = linspace(-5,0,100);
    plot(0.5*ones(length(markers2),1),markers2,'--k');
    
    text(0.5,0.3,'$s^{*}=0.5$','interpreter','latex','fontsize',12);
    
    xlabel('sex ratio at birth (proportion male), $s$','interpreter','latex','fontsize',12);
    ylabel('fitness gradient','interpreter','latex','fontsize',12);
    
    % Save to PDF
    if(exist('save2pdf','file'))
        save2pdf('fig_2.pdf')
    end
end

