%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sets up symbolic system to find the fitness and its derivative
% Tested 31/08/21 using MATLAB 2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define all the synbolic variables and populations
syms MJ FJ MA FA U d uMJ uFJ uMA uFA b delta smut s tau xi

%% Defining system
NJ = MJ + FJ;
Nu = MA + FA;
NA = Nu + 2*U;
N = NJ + NA;
% xi = (1/N)*dN/dt
xi_sub = b*U - uMJ*MJ - uFJ*FJ - uMA*MA - uFA*FA - (uMA + uFA)*U;
% Define matrix A
A = [-uMJ - tau - xi, 0, 0, 0, 0.5*b*((s+smut)/2), 0.5*b*((s+smut)/2); ...
    0, -uFJ - tau - xi, 0, 0, 0.5*b*(1-((s+smut)/2)), 0.5*b*(1-((s+smut)/2));...
    tau, 0, -uMA - (delta*FA)/(Nu) - xi, 0, uFA + d, 0;...
    0, tau, 0, -uFA - (delta*MA)/(Nu) - xi, 0, uMA + d;...
    0, 0, (delta*FA)/(Nu), 0, -(uMA + uFA + d + xi), 0;...
    0, 0, 0, (delta*MA)/(Nu), 0, -(uMA + uFA + d + xi)];
% Define matrix B
B = [0, 0, 0, 0, 0.5*b*((s+smut)/2), 0.5*b*((s+smut)/2); ...
    0, 0, 0, 0, 0.5*b*(1-((s+smut)/2)), 0.5*b*(1-((s+smut)/2));...
    0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0];
% Calculate matrix V
V = B - A;

%% Completing calculations 
% Find inverse of V
inverseV = inv(V);
% Calculate next generation matrix
nextGen = B*inverseV;
% Calculate eigenvalues
eigenvalues = eig(nextGen);
% Find symbollic fitness, w
sw = simplify(subs(eigenvalues(6),xi,xi_sub));
% Differentiate the invasion fitness to find the selection gradient
sdw = simplify(subs(diff(sw,smut),smut,s));

%% Convert to normal matlab functions
w = matlabFunction(sw,'vars',{'MJ','FJ','MA','FA','U','b','s','smut','uMJ','uFJ','uMA','uFA','tau','delta','d'});

dw = matlabFunction(sdw,'vars',{'MJ','FJ','MA','FA','U','b','s','uMJ','uFJ','uMA','uFA','tau','delta','d'});

%% Clear the symbollic variables
symObj = syms;
cellfun(@clear,symObj);

clear symObj
disp('Setup complete')
