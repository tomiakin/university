%% Monte Carlo sampling of the Gutenberg-Richter inverse model
clear
clc

Mmin    = 3.9;
Mmax    = 5.3;
Delta_M = 0.1;

a       = 3.91;
b       = 1.04;

% Annual occurrence rate of earthquakes with M >= Mmin
Lambda = 10^(a-b*Mmin);

Tmax = 10000; % Duration of the simulated catalog (years)

time = 0;
count = 0;
f = waitbar(0,'Please wait...');
while time <= Tmax
    
    % First simulation (time of the next earthquake)
    % This is the inverse of the Exponential distribution
    iatime = expinv(rand(1,1),1/Lambda);
    time = time + iatime;
    
    count = count + 1;
    waitbar(time/Tmax,f,[num2str(round(time/Tmax*10*100)/10),' %'])
    
    % Second simulation (magnitude value)
    m = - log10(10^(-b*Mmin)-rand(1,1)*(10^(-b*Mmin)-10^(-b*Mmax)))/b; % CDF - Mmin and Mmax
    
    syncat(count,:) = [time m];
end
close(f)
%% Let's now see how the simulated earthquake catalog looks like
Delta_M = 0.1;

years = syncat(:,1);
magnitude = syncat(:,2);

ThreshYear = [min(years) max(years)];

Minterval0 = 3.5:Delta_M:Mmax;
Minterval  = Mmin:Delta_M:Mmax;

index = find(years>ThreshYear(1,1));
magnitude = magnitude(index,:);

N0 = zeros(length(Minterval0)-1,1);
for i = 1:length(Minterval0)-1
    N0(i) = length(find(magnitude >= Minterval0(i)))/(ThreshYear(2)-ThreshYear(1)+1);
end

N = zeros(length(Minterval)-1,1);
for i = 1:length(Minterval)-1
    N(i) = length(find(magnitude >= Minterval(i)))/(ThreshYear(2)-ThreshYear(1)+1);
end

N0(N0==0)=NaN;
N(N==0)=NaN;

% tmp = polyfit(Minterval(1:end-1)+Delta_M/2,log10(N'),1);
% GRFit(1,1)=tmp(2);
% GRFit(1,2)=-tmp(1);

% Calculate the b-value (maximum likelihood)
MminCat=min(magnitude(magnitude>Mmin));
MmeanCat=mean(magnitude(magnitude>Mmin));
GRFit(1,2) = (1/(MmeanCat-(MminCat-Delta_M/2)))*log10(exp(1));
% Calculate the a-value
GRFit(1,1) = log10(N(1)) + GRFit(1,2) * MminCat;


semilogy(Minterval0(1:end-1)+Delta_M/2,N0,'ks',[Mmin Mmax],10.^(GRFit(1,1)-GRFit(1,2)*[Mmin Mmax]),'r','linewidth',3,'markersize',10); hold on;
semilogy(Minterval(1:end-1)+Delta_M/2,N,'gs',[Mmin Mmax],10.^(GRFit(1,1)-GRFit(1,2)*[Mmin Mmax]),'r','linewidth',3,'markersize',10); hold on;
hold on
plot([Mmin Mmin],[0.01 8],'b--','linewidth',2)
hold on
plot([Mmax Mmax],[0.01 8],'b--','linewidth',2)

axis square; xlabel('Magnitude (Mw)'); ylabel('\lambda = Occurrence Rate'); 
axis([3 7 0.001 100]); grid on;
title ('Gutenberg-Richter Distribution')
txt11 = (['Characteristic values: ','Mw_{min}= ',num2str(Mmin)]);
text(4,18,txt11,'FontWeight','bold','FontSize',14)
txt12 = (['a=',num2str(round(GRFit(1,1)*100)/100),' & b=',num2str(round(GRFit(1,2)*100)/100)]);
text(4,12,txt12,'FontWeight','bold','FontSize',14)
text(4, 0.005,'Log_{10}(\lambda) = a - b Mw','FontWeight','bold','FontSize',14)
text(4, 0.002,['\lambda(Mw_{min}) = ',num2str(round(10^(GRFit(1,1)-GRFit(1,2)*Mmin)*100)/100)],'FontWeight','bold','FontSize',14)
set(gca,'FontSize',16,'Layer','Top','XTick',...
    [2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7],'YMinorTick','on','YScale','log')