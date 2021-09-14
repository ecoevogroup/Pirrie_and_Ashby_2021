%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sets up symbolic system to find the fitness and its derivative
% Tested 31/08/21 using MATLAB 2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define all the synbolic variables and populations
syms MJ FJ MA FA U d q uMJ uFJ uMA uFA b delta smut s tau 

%% Defining system
NJ = MJ + FJ;
Nu = MA + FA;
NA = Nu + 2*U;
N = NJ + NA;
% Define matrix A
A = [-uMJ - tau, 0, 0, 0, 0.5*b*(1-q*N)*((s+smut)/2), 0.5*b*(1-q*N)*((s+smut)/2); ...
    0, -uFJ - tau, 0, 0, 0.5*b*(1-q*N)*(1-((s+smut)/2)), 0.5*b*(1-q*N)*(1-((s+smut)/2));...
    tau, 0, -uMA - (delta*FA)/(Nu), 0, uFA + d, 0;...
    0, tau, 0, -uFA - (delta*MA)/(Nu), 0, uMA + d;...
    0, 0, (delta*FA)/(Nu), 0, -(uMA + uFA + d), 0;...
    0, 0, 0, (delta*MA)/(Nu), 0, -(uMA + uFA + d)];
% Define matrix B
B = [0, 0, 0, 0, 0.5*b*(1-q*N)*((s+smut)/2), 0.5*b*(1-q*N)*((s+smut)/2); ...
    0, 0, 0, 0, 0.5*b*(1-q*N)*(1-((s+smut)/2)), 0.5*b*(1-q*N)*(1-((s+smut)/2));...
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
sw = eigenvalues(6);
% Differentiate the invasion fitness to find the selection gradient
sdw = subs(diff(sw,smut),smut,s);

%% Convert to normal matlab functions
dw = matlabFunction(sdw,'vars',{'MJ','FJ','MA','FA','U','b','q','s','uMJ','uFJ','uMA','uFA','tau','delta','d'});

w = matlabFunction(sw,'vars',{'MJ','FJ','MA','FA','U','b','q','s','smut','uMJ','uFJ','uMA','uFA','tau','delta','d'});

%% Clear the symbollic variables
symObj = syms;
cellfun(@clear,symObj);

clear symObj
disp('Setup complete')
