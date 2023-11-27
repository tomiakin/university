
% Akkar and Bommer (2010) - Mw = 5.0-7.6 (constant sigma model)
% log10(PGV/PGA/PSA) = b1 + b2*M + b3*M^2 + (b4+b5*M)*log10(sqrt(Rjb^2+b6^2)) + b7*Ss + b8*Sa + b9*Fn + b10*Fr
% Ss = 1: soft soil site index (Vs30 < 360 m/s)
% Sa = 1: firm soil site index (360 < Vs30 < 750 m/s)
% Fn = 1: normal faulting index
% Fr = 1: reverse faulting index
% coef = [Tn b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 sigma-intra sigma-inter sigma-total]
% PGA: Tn = 0 and PGV: Tn = -1
clear
clc

%Tn   = [0 0.1 0.2 0.3 0.4 0.5 0.75 1.0 1.5 2.0 3.0];

Tn = 0;

m = 5.5;
R = 0:0.01:200;
Fault_Type = 'normal';

Vs30 = 1000;

load GMPEcoef_AB10
for ii = 1:length(Tn);
    coef1(ii,:) = coef_AB10(find(Tn(ii)==coef_AB10(:,1)),:);
end

if strcmp(Fault_Type,'strike-slip')==1
    FMech = [0 0];
elseif strcmp(Fault_Type,'normal')==1
    FMech = [1 0];
elseif strcmp(Fault_Type,'reverse')==1
    FMech = [0 1];
end

soil1 = [0 0]; % [Ss,Sa]
if Vs30 <= 360;
    soil1(1) = 1;
elseif Vs30 > 360 && Vs30 <= 750;
    soil1(2) = 1;
end


for i=1:length(R)
    gmp_central_estimate(i) = 10.^(coef1(:,2)       + ...
        coef1(:,3) * m   +  ...
        coef1(:,4) * m^2 +  ...
        (coef1(:,5)+coef1(:,6)*m).*log10(sqrt((R(i))^2+coef1(:,7).^2)) + ...
        coef1(:,8) * soil1(1) + ...
        coef1(:,9) * soil1(2) + ...
        coef1(:,10)*FMech(1)  + ...
        coef1(:,11)*FMech(2))/981; % (g)
end

figure(1)
subplot(221)
plot(R,gmp_central_estimate,'k','linewidth',3)
hold on
plot(R,10.^(log10(gmp_central_estimate)+coef1(:,14)),':r','linewidth',3)
hold on
plot(R,10.^(log10(gmp_central_estimate)-coef1(:,14)),':r','linewidth',3)
xlabel('Joyner-Boore Distance (km)')
ylabel('S_a (g)')
axis square
xlim([0 200])
hold on
title('Linear space')
grid on
subplot(222)
loglog(R,gmp_central_estimate,'k','linewidth',3)
hold on
loglog(R,10.^(log10(gmp_central_estimate)+coef1(:,14)),':r','linewidth',3)
hold on
loglog(R,10.^(log10(gmp_central_estimate)-coef1(:,14)),':r','linewidth',3)
xlabel('Joyner-Boore Distance (km)')
ylabel('S_a (g)')
axis square
xlim([0 200])
title('Loglog space')
grid on
%%
%%
%%
clear gmp_central_estimate
m = 5.5;
R = 30;
Fault_Type = 'normal';

Vs30 = 1000;

load GMPEcoef_AB10
for ii = 1:length(Tn);
    coef1(ii,:) = coef_AB10(find(Tn(ii)==coef_AB10(:,1)),:);
end

if strcmp(Fault_Type,'strike-slip')==1
    FMech = [0 0];
elseif strcmp(Fault_Type,'normal')==1
    FMech = [1 0];
elseif strcmp(Fault_Type,'reverse')==1
    FMech = [0 1];
end

soil1 = [0 0]; % [Ss,Sa]
if Vs30 <= 360;
    soil1(1) = 1;
elseif Vs30 > 360 && Vs30 <= 750;
    soil1(2) = 1;
end


for i=1:length(R)
    gmp_central_estimate(i) = 10.^(coef1(:,2)       + ...
        coef1(:,3) * m   +  ...
        coef1(:,4) * m^2 +  ...
        (coef1(:,5)+coef1(:,6)*m).*log10(sqrt((R(i))^2+coef1(:,7).^2)) + ...
        coef1(:,8) * soil1(1) + ...
        coef1(:,9) * soil1(2) + ...
        coef1(:,10)*FMech(1)  + ...
        coef1(:,11)*FMech(2))/981; % (g)
end

figure(1)
subplot(221)
plot([R R],[0.000001 5],'k','linewidth',2)
imx=0.000001:0.001:5;
hold on
plot(R+normpdf(log10(imx),log10(gmp_central_estimate),coef1(:,14))./imx./max(normpdf(log10(imx),log10(gmp_central_estimate),coef1(:,14))),imx,'b','linewidth',2)

xlabel('Joyner-Boore Distance (km)')
ylabel('S_a (g)')
axis square
xlim([1 200])
ylim([0 0.5])
hold on
title('Linear space')
grid on
set(gca,'FontSize',16)

subplot(222)
loglog([R R],[0.000001 5],'k','linewidth',2)
imx=0.000001:0.001:5;
hold on
loglog(R+normpdf(log10(imx),log10(gmp_central_estimate),coef1(:,14))./imx./max(normpdf(log10(imx),log10(gmp_central_estimate),coef1(:,14))),imx,'b','linewidth',2)

xlabel('Joyner-Boore Distance (km)')
ylabel('S_a (g)')
axis square
xlim([1 200])
ylim([0.001 1])
title('Loglog space')
set(gca,'FontSize',16)
grid on
%%
figure(1)
subplot(2,2,3:4)
imx=0.000001:0.001:5;
hold on
plot(imx,1-normcdf(log10(imx),log10(gmp_central_estimate),coef1(:,14)),'b','linewidth',2)

xlabel('S_a (g)')
ylabel('P(S_a>s_a|R,M_W)')
axis square

xlim([0 0.5])
grid on
box on
set(gca,'FontSize',16)