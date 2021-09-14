function sexRatioSweep

% sexRatioSweep.m
%
% Generates data for Figure S.1 in "Does differential mortality after parental
% investment affect sex ratio evolution? No."

% Fixed parameters
b = 25;
d = 0.1;
s1 = 0.5;
delta = 1;
tau = 0.1;
uFJ = 0.05;
uFA = 0.05;
maxSteps = 1e7;
res1 = 41;
res2 = 21;
simtotal = 100;
fullOutput = 0;

% Variable parameters
S2 = linspace(0,1,res1);
Q = [0,1e-4];
cullingThreshold = [1e4,inf];
survivingProportion = [0.1,1];

% Set up data arrays
FREQ1 = NaN*zeros(length(S2),res2,length(Q),simtotal);
FREQ2 = NaN*zeros(length(S2),res2,length(Q),simtotal);

% Run paramater sweep
for simnum=1:simtotal
    % Vary juvenile male mortality
    uMA = uFA;
    UMJ = linspace(1e-3,2*uFJ,res2);
    
    for k=1:length(Q)
        q = Q(k);
        for j=1:length(UMJ)
            uMJ = UMJ(j);
            tic;
            for i=1:length(S2)
                s2 = S2(i);
                for l=1:simtotal
                    [~,pops] = sexRatioSimulation(b,d,q,s1,s2,delta,tau,uMJ,uFJ,uMA,uFA,cullingThreshold(k),survivingProportion(k),maxSteps,fullOutput);
                    FREQ1(i,j,k,l) = pops(1)/sum(pops);
                end
            end
            toc;
            PROGRESS = [j/length(UMJ),k/length(Q)]
        end
    end
    
    % Vary adult male mortality
    uMJ = uFJ;
    UMA = linspace(1e-3,2*uFA,res2);
    
    for k=1:length(Q)
        q = Q(k);
        for j=1:length(UMA)
            uMA = UMA(j);
            tic;
            for i=1:length(S2)
                s2 = S2(i);
                for l=1:simtotal
                    [~,pops] = sexRatioSimulation(b,d,q,s1,s2,delta,tau,uMJ,uFJ,uMA,uFA,cullingThreshold(k),survivingProportion(k),maxSteps,fullOutput);
                    FREQ2(i,j,k,l) = pops(1)/sum(pops);
                end
            end
            toc;
            PROGRESS = [j/length(UMA),k/length(Q)]
        end
    end
end
clear i j k l s2 q uMA uMJ PROGRESS ans
save('sexRatioSweep.mat')
