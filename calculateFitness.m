%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function run by sexRatioSelection.m. Calculates the fitness gradient
% for each value of the sex ratio inputted and outputs that in an array
% Tested 31/08/21 using MATLAB 2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dwArray = calculateFitness(pars,setupPars,traitArray,popini,generaldw)

% Unpacking the setup parameters
tol = setupPars.tol;
searchCount = setupPars.searchCount;
tspan = setupPars.tspan;
experiment = setupPars.experiment;
solver = setupPars.solver;

% Pre-define array for efficiency
dwArray = NaN(length(traitArray),1);

% Selects the correct fitness gradient function for the experiment
run('experimentChoice.m');

% Calculates the fitness gradient for every sex ratio value
for j = 1 : length(traitArray)
    [t,pop] = runModelSolve(traitArray(j));
    dwArray(j) = dw(pop(end,:),traitArray(j));
end
end