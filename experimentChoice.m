%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script used by calculateFitness.m to assign the correct
% model solver and fitness gradient for the specific experiment.
% Tested 31/08/21 using MATLAB 2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(experiment,'uMJAgeUnions')
    if strcmp(solver,'MEX')
        runModelSolve = @(trait) modelSolveAgeUnions(tspan(end),popini(1),popini(2),popini(3),popini(4),popini(5),pars.b,pars.q,trait,pars.uMJ,pars.uFJ,pars.uMA,pars.uFA,pars.tau,pars.delta,pars.d);
        dw = @(popinput,trait) generaldw(popinput(1),popinput(2),popinput(3),popinput(4),popinput(5),trait,pars.uMJ,pars.uFJ);
    else
        disp('Unrecognised experiment');
    end
elseif strcmp(experiment,'uMAAgeUnions')
    if strcmp(solver,'MEX')
        runModelSolve = @(trait) modelSolveAgeUnions(tspan(end),popini(1),popini(2),popini(3),popini(4),popini(5),pars.b,pars.q,trait,pars.uMJ,pars.uFJ,pars.uMA,pars.uFA,pars.tau,pars.delta,pars.d);
        dw = @(popinput,trait) generaldw(popinput(1),popinput(2),popinput(3),popinput(4),popinput(5),trait,pars.uMA,pars.uFA);
    else
        disp('Unrecognised experiment');
    end
elseif strcmp(experiment,'uMAAgeUnionsFrequency')
    if strcmp(solver,'MEX')
        runModelSolve = @(trait) modelSolveAgeUnionsFrequency(tspan(end),popini(1),popini(2),popini(3),popini(4),popini(5),pars.b,trait,pars.uMJ,pars.uFJ,pars.uMA,pars.uFA,pars.tau,pars.delta,pars.d);
        dw = @(popinput,trait) generaldw(popinput(1),popinput(2),popinput(3),popinput(4),popinput(5),trait,pars.uMA,pars.uFA);
    else
        disp('Unrecognised experiment');
    end
else
    disp('Unrecognised experiment')
end
