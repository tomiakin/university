clear
clc
close all
load SUB_CATALOG
seismico_source = 2; % change based on zone
magnitude = SUB_CATALOG(seismico_source).CATALOG.Mw;

%% Define the Magnitude increment
deltaM = 0.1;

%% Distribution of magnitude
Mmin    = 4.5; % same for both zones
Mmax    = max(magnitude);
Delta_M = 0.1;
b       = 0.98; % changes based on the 'b' value from GR_zone

mm = linspace(Mmin,Mmax,100);

G  = (1-10.^(-b*(mm-Mmin)))./(1-10.^(-b*(Mmax-Mmin)));

Mi  = Mmin + deltaM/2 : deltaM : Mmax - deltaM/2;
PMi = zeros(size(Mi));

for i=1:length(Mi)
    PMi(i) = interp1(mm',G',Mi(i)+deltaM/2) - interp1(mm',G',Mi(i)-deltaM/2);
end

figure
subplot(121)

plot(mm,G,'k','linewidth',2)
xlabel('M_W - Magnitude')
ylabel('CDF')
axis square
grid on
set(gca,'FontSize',16)

subplot(122)
bar(Mi,PMi,'facecolor','r')
xlabel('M_W - Magnitude')
ylabel('p(M=m|M\geqM_{min})')
axis square
grid on
set(gca,'FontSize',16)
