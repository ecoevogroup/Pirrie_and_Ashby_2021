function [T,pops] = sexRatioSimulation(b,d,q,s1,s2,delta,tau,uMJ,uFJ,uMA,uFA,cullingThreshold,survivingProportion,maxSteps,fullOutput)

% sexRatioSimulation.m
%
% Generates data for Figure S.1 in "Does differential mortality after parental
% investment affect sex ratio evolution? No."

% Set up initial populations
MJ1 = 1000;
FJ1 = 1000;
MA1 = 1000;
FA1 = 1000;
MJ2 = MJ1;
FJ2 = FJ1;
MA2 = MA1;
FA2 = FA1;
UM1F1 = 0;
UM1F2 = 0;
UM2F1 = 0;
UM2F2 = 0;
N1 = MJ1 + FJ1 + MA1 + FA1 + 2*UM1F1 + UM1F2 + UM2F1;
N2 = MJ2 + FJ2 + MA2 + FA2 + 2*UM2F2 + UM1F2 + UM2F1;

% Set up event population changes
popChange = zeros(32,12);
popChange(1:2,1:2) = eye(2); % Births of type 1
popChange(3:4,5:6) = eye(2); % Births of type 2
popChange(5:6,1:2) = -eye(2); % Ageing of type 1
popChange(7:8,5:6) = -eye(2); % Ageing of type 2
popChange(5:6,3:4) = eye(2); % Ageing of type 1
popChange(7:8,7:8) = eye(2); % Ageing of type 2
popChange(9:16,1:8) = -eye(8); % Deaths (single)
popChange(17,3:4) = -1; % Union forming
popChange(18,[3,8]) = -1; % Union forming
popChange(19,[4,7]) = -1; % Union forming
popChange(20,7:8) = -1; % Union forming
popChange(17:20,9:12) = eye(4); % Union forming
popChange(21:24,9:12) = -eye(4); % Union divorce
popChange(21,3:4) = 1; % Union forming
popChange(22,[3,8]) = 1; % Union forming
popChange(23,[4,7]) = 1; % Union forming
popChange(24,7:8) = 1; % Union forming
popChange(25:2:32,9:12) = -eye(4); % Death within union
popChange(26:2:32,9:12) = -eye(4); % Death within union
popChange(25,4) = 1; % Death within union
popChange(26,3) = 1; % Death within union
popChange(27,8) = 1; % Death within union
popChange(28,3) = 1; % Death within union
popChange(29,4) = 1; % Death within union
popChange(30,7) = 1; % Death within union
popChange(31,8) = 1; % Death within union
popChange(32,7) = 1; % Death within union

% Preallocate rates vector
rates = zeros(32,1);

% Simulation parameters
count = 1;
t = 0;

% Outputs
if(fullOutput>0)
    pops = zeros(maxSteps,2);
    T = zeros(maxSteps,1);
    pops(1,:) = [N1,N2];
    T(1) = t;
end

% Main loop
while(count<maxSteps)
    
    % Set rates
    rates = setRates(rates,b,d,q,s1,s2,delta,tau,uMJ,uFJ,uMA,uFA,MJ1,FJ1,MA1,FA1,MJ2,FJ2,MA2,FA2,UM1F1,UM1F2,UM2F1,UM2F2);
    
    % Find next event
    sumRates = sum(rates);
    r = rand*sumRates;
    ratesCumSum = rates(1);
    event = 1;
    while(ratesCumSum<r)
        event = event+1;
        ratesCumSum = ratesCumSum + rates(event);
    end
    
    % Update populations based on event
    MJ1 = MJ1 + popChange(event,1);
    FJ1 = FJ1 + popChange(event,2);
    MA1 = MA1 + popChange(event,3);
    FA1 = FA1 + popChange(event,4);
    MJ2 = MJ2 + popChange(event,5);
    FJ2 = FJ2 + popChange(event,6);
    MA2 = MA2 + popChange(event,7);
    FA2 = FA2 + popChange(event,8);
    UM1F1 = UM1F1 + popChange(event,9);
    UM1F2 = UM1F2 + popChange(event,10);
    UM2F1 = UM2F1 + popChange(event,11);
    UM2F2 = UM2F2 + popChange(event,12);
    
    % If above culling threshold, then cull
    N1 = MJ1 + FJ1 + MA1 + FA1 + 2*UM1F1 + UM1F2 + UM2F1;
    N2 = MJ2 + FJ2 + MA2 + FA2 + 2*UM2F2 + UM1F2 + UM2F1;
    if((N1+N2)>cullingThreshold)
        MJ1 = round(MJ1*survivingProportion);
        FJ1 = round(FJ1*survivingProportion);
        MA1 = round(MA1*survivingProportion);
        FA1 = round(FA1*survivingProportion);
        MJ2 = round(MJ2*survivingProportion);
        FJ2 = round(FJ2*survivingProportion);
        MA2 = round(MA2*survivingProportion);
        FA2 = round(FA2*survivingProportion);
        UM1F1 = round(UM1F1*survivingProportion);
        UM1F2 = round(UM1F2*survivingProportion);
        UM2F1 = round(UM2F1*survivingProportion);
        UM2F2 = round(UM2F2*survivingProportion);
        N1 = MJ1 + FJ1 + MA1 + FA1 + 2*UM1F1 + UM1F2 + UM2F1;
        N2 = MJ2 + FJ2 + MA2 + FA2 + 2*UM2F2 + UM1F2 + UM2F1;
    end
    
    % Update outputs
    count = count+1;
    if(fullOutput>0)
        pops(count,:) = [N1,N2];
        if(sumRates==0)
            T(count) = T(count-1);
            break
        end
        T(count) = T(count-1) - log(rand)/sumRates;
    end
    
    % Exit if one is extinct
    if(N1==0 || N2==0 || sumRates==0)
        break;
    end
end

% Update outputs
if(fullOutput>0)
    pops = pops(1:count,:);
    T = T(1:count);
else
    T = NaN;
    pops = [N1,N2];
end

function rates = setRates(rates,b,d,q,s1,s2,delta,tau,uMJ,uFJ,uMA,uFA,MJ1,FJ1,MA1,FA1,MJ2,FJ2,MA2,FA2,UM1F1,UM1F2,UM2F1,UM2F2)

NU = max(1E-30,MA1 + FA1 + MA2 + FA2);
N = MJ1 + FJ1 + MJ2 + FJ2 + NU + 2*(UM1F1 + UM1F2 + UM2F1 + UM2F2);

% Births
dd = max(0,b*(1-q*N));
rates(1) = dd*(s1*UM1F1 + (s1+s2)*(UM1F2+UM2F1)/4); % birth of MJ1
rates(2) = dd*((1-s1)*UM1F1 + (1-(s1+s2)/2)*(UM1F2+UM2F1)/2); % birth of FJ1
rates(3) = dd*(s2*UM2F2 + (s1+s2)*(UM1F2+UM2F1)/4); % birth of MJ2
rates(4) = dd*((1-s2)*UM2F2 + (1-(s1+s2)/2)*(UM1F2+UM2F1)/2); % birth of FJ2

% Ageing
rates(5) = tau*MJ1;
rates(6) = tau*FJ1;
rates(7) = tau*MJ2;
rates(8) = tau*FJ2;

% Deaths (single)
rates(9) = uMJ*MJ1;
rates(10) = uFJ*FJ1;
rates(11) = uMA*MA1;
rates(12) = uFA*FA1;
rates(13) = uMJ*MJ2;
rates(14) = uFJ*FJ2;
rates(15) = uMA*MA2;
rates(16) = uFA*FA2;

% Union forming
rates(17) = delta*MA1*FA1/NU;
rates(18) = delta*MA1*FA2/NU;
rates(19) = delta*MA2*FA1/NU;
rates(20) = delta*MA2*FA2/NU;

% Union divorce
rates(21) = d*UM1F1;
rates(22) = d*UM1F2;
rates(23) = d*UM2F1;
rates(24) = d*UM2F2;

% Death within union
rates(25) = uMA*UM1F1;
rates(26) = uFA*UM1F1;
rates(27) = uMA*UM1F2;
rates(28) = uFA*UM1F2;
rates(29) = uMA*UM2F1;
rates(30) = uFA*UM2F1;
rates(31) = uMA*UM2F2;
rates(32) = uFA*UM2F2;
