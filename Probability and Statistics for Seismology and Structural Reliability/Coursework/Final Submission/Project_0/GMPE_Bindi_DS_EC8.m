function [sa,sigma_Inter,sigma_Intra,sigma_S2S,sigma_total]=GMPE_Bindi_DS_EC8(rupType,period,magnitude,distance,vs30,faultType,table1)
% D. Bindi, M. Massa, L. Luzi, et al (2014).
% Pan-European ground-motion prediction equations for the average horizontal...
%  component of PGA, PGV, and 5%-damped PSA at spectral periods up to 3.0s...
%  using the RESORCE dataset.
% Bull Earthquake Eng.
% results are in g (Sa and PGA) and in cm/s (PGV)

%% INPUT
% rupType:  1=Table 1(Rjb)  2=Table 3(Rhypo)
% period: -1=PGV, 0=PGA, T=0.02-3s
% magnitude        : Magnitude value
% distance         : Rjb Rhypo
% vs30             : S-waves velocity in the first 30 m
% faultType        : 0=undetermined 1=normal 2=reverse 3=strike-slip

%% input examples
% distance = [0.01:0.01:100]';
% rupType = 1;
% period = 2;
% magnitude = 5*ones(size(distance));
% vs30 = 800*ones(size(distance));
% faultType = 3*ones(size(distance));

%% parameters
if rupType == 1
%     load('table_bindi_DS_EC8.mat', 'table1')
    Values = table1; clear table1
elseif rupType == 2
%     load('table_bindi_DS_EC8.mat', 'table3')
    Values = table3; clear table3
end
Names  = {'T'; 'e1'; 'c1'; 'c2'; 'h'; 'c3'; 'b1'; 'b2'; 'b3'; 'classA'; 'classB'; 'classC'; 'classD'; 'sofN'; 'sofR'; 'sofS'; 'sofU'; 'sigma_Inter'; 'sigma_Intra'; 'sigma_S2S'; 'sigma_total';};
paras = table(Values','RowNames',Names);

%% parameter interpolation between periods
parasArray  = table2array(paras);
if period > 0
    index0      = find(parasArray(1,:)==0)+1;
    ParaTarPeriod  = interp1(log(parasArray(1,index0:end))',parasArray(2:end,index0:end)',log(period));
    clear Names Values index0 parasArray
elseif period == 0
    ParaTarPeriod = parasArray(2:end,2);
elseif period == -1
    ParaTarPeriod = parasArray(2:end,1);
end
%% Definition of the GMPE parameters
e1 = ParaTarPeriod(1);
c1 = ParaTarPeriod(2);
c2 = ParaTarPeriod(3);
h = ParaTarPeriod(4); 
c3 = ParaTarPeriod(5); 
b1 = ParaTarPeriod(6); 
b2 = ParaTarPeriod(7); 
b3 = ParaTarPeriod(8);
classA = ParaTarPeriod(9);
classB = ParaTarPeriod(10);
classC = ParaTarPeriod(11);
classD = ParaTarPeriod(12);
sofN = ParaTarPeriod(13); 
sofR = ParaTarPeriod(14); 
sofS = ParaTarPeriod(15);
sofU = ParaTarPeriod(16); 
sigma_Inter = ParaTarPeriod(17); 
sigma_Intra = ParaTarPeriod(18); 
sigma_S2S = ParaTarPeriod(19); 
sigma_total = ParaTarPeriod(20);

Mref = 5.5;
Mh = 6.75;
Rref = 1;
Vref = 800;
%% equation 2
Fd = (c1 + c2 * (magnitude - Mref)) .* log10(sqrt(distance.^2+h^2) / Rref) - c3 * (sqrt(distance.^2+h^2) - Rref);

%% equation 3
Fm = b1 * (magnitude - Mh) + b2 * (magnitude - Mh).^2;
Fm(magnitude>Mh) = b3 * (magnitude(magnitude>Mh) - Mh);

%% Fs and Fsof in equation 1
Fs = zeros(size(vs30));
Fs(vs30>=800) = classA;
Fs(vs30<800 & vs30>=360) = classB;
Fs(vs30<360 & vs30>=180) = classC;
Fs(vs30<180) = classD;

Fsof = zeros(size(faultType));
Fsof(faultType==0) = sofU;
Fsof(faultType==1) = sofN;
Fsof(faultType==2) = sofR;
Fsof(faultType==3) = sofS;

%% equation 1
sa = 10.^(e1 + Fd + Fm + Fs + Fsof);

if period ~= -1 
    sa = sa/980.665;
end

end