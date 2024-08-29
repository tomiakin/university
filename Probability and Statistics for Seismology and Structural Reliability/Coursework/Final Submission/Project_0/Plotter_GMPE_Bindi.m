%% Input
magnitude = 6;
Tn = 0;
R = 0:0.01:300;
Vs30 = 400;
FaultType = 'normal';
load('table_bindi_DS_EC8.mat', 'table1')

%% Processing

if strcmp(FaultType,'strike-slip')==1
    faultType=3;
elseif strcmp(FaultType,'normal')==1
    faultType=1;
elseif strcmp(FaultType,'reverse')==1
    faultType=2;
end

[gmp_central_estimate2,~,~,~,sigma_total]=GMPE_Bindi_DS_EC8(1,Tn,magnitude,R,Vs30,faultType,table1);

%% Plot
loglog(R,gmp_central_estimate2,'k'); hold on
loglog(R,10.^(log10(gmp_central_estimate2)+sigma_total),'k--'); hold on
loglog(R,10.^(log10(gmp_central_estimate2)-sigma_total),'k--'); hold on

xlabel('Distance [km]')
ylabel('IM [g]')
