%% ZONE 5 AKA ZONE 2

clear
clc
close all

load SUB_CATALOG

seismico_source = 2;
magnitude = SUB_CATALOG(seismico_source).CATALOG.Mw;
year      = SUB_CATALOG(seismico_source).CATALOG.Year;

Mmin = 4.5; % Min that will cause damage to structure, change in future?
Mmax = max(magnitude);
Delta_M = 0.2;
ThreshYear = [1800 max(year)];

%%
Minterval0 = 3.5:Delta_M:Mmax;
Minterval  = Mmin:Delta_M:Mmax;

index = find(year>ThreshYear(1,1));
magnitude = magnitude(index,:);

Mmintoshow = min(magnitude)
Mmaxtoshow = max(magnitude)

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