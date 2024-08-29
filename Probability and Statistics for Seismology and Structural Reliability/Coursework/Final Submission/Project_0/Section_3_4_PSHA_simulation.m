%% Initialisation of the workspace
clear % the "clear" command clear the workspace (imporant to clean the memory)
clc   % the "clc" clear the window command from previous command (only visual)
close all

%% POINT OF INTEREST and value of the VS30 of the soil
X = [35.36	25.12]; %[Latitude Longitude]
% X = [38 14.5];

Vs30 = 300; % m/s shear wave velocity of the soil under the site of interest
%% vector of IM to consider for the hazard curve
imx=0:0.01:4; %g

%% Vibration periods (T) to be considered for the spectral accelerations
Tnn   = [0 0.1 0.2 0.3 0.4 0.5 0.75 1.0 1.5 2.0 3.0]; % seconds

%% Return Period values to be consider
Tr = [30 50 475 2475]; % corresponding to 81%, 63%, 10% and 2% in 50 years

%% Upload the file of the seismic zones.
load Final_zones
SZ=Final_zones;

%% Below you need to defien the valus of Mmin, Mmax, a and b you 
%  calculated from the G&R
a    = [3.33 4.33]; % here you need to repleace your values
b    = [0.78 0.98]; % here you need to repleace your values

Mwmin = [4.5 4.5];      % here you need to repleace your values
Mwmax = [7.5 6.4];      % here you need to repleace your values

Lambda_min = 10.^(a-b.*Mwmin);

FaultType = {'normal', 'normal'}; % it can be 'normal', 'reverse', 'strike-slip'. Here you need to repleace your types ONE TYPE FOR EACH ZONE

%% loading the coefficients for the GMPE
load GMPEcoef_AB10
coef1 = zeros(length(Tnn),length(coef_AB10(1,:)));
for ii = 1:length(Tnn)
    coef1(ii,:) = coef_AB10(Tnn(ii)==coef_AB10(:,1),:);
end

%% Number of years for the simulation
Num = 10000; % years

%% Simulation-based procedure
plotter = 1;

time = 0;
counttmp = 0;

if plotter==1
figure(1)
plot([SZ.X],[SZ.Y])
end

wb = waitbar(0,'Please wait...');

while time < Num
    
    % Simulation of the inter-arrival time
%     iatime = expinv(rand(1,3),1./Lambda_min);
    iatime = -1./Lambda_min .*log(1-rand(1,length(SZ)));

    % Identification of the seismic zone where the earthquake occurs
    [~,sz] = min(iatime);
    
    % Calculation of the absolute time of the event
    time   = time + iatime(sz);
    
    % To show a waiting bar
    waitbar(time/Num,wb,[num2str(round(time/Num*10*100)/10),' %'])
    
    % Counter for the earthquakes
    counttmp = counttmp + 1;
    
    % GR relationship - Simulation of the magnitude level
    m = - log10(10^(-b(sz)*Mwmin(sz))-rand(1,1)*(10^(-b(sz)*Mwmin(sz))-10^(-b(sz)*Mwmax(sz))))/b(sz); % CDF - Mmin and Mmax
    
    % Earthquake location - Simulation of the epicentre
    in_event = 0;
    while in_event == 0
        lon = rand(1,1)*(max(SZ(sz).X) - min(SZ(sz).X)) + min(SZ(sz).X);
        lat = rand(1,1)*(max(SZ(sz).Y) - min(SZ(sz).Y)) + min(SZ(sz).Y);
        in_event = inpolygon(lon,lat,SZ(sz).X,SZ(sz).Y);
    end
    
    if plotter ==1
    hold on
    plot(lon,lat,'r.')
    end
    
    repi = deg2km(distance(lat,lon,X(1),X(2))); % Epicentral distance
    Ri  = distance_converter(repi);             % Conversion to RJB
    
    % GMPE
    eps = randn(length(Tnn),1); % simulation of the error for the GMPE for all the Spectral accelerations
    
    %% Fault mechanism term
    if strcmp(FaultType{sz},'strike-slip')==1
        FMech = [0 0];
    elseif strcmp(FaultType{sz},'normal')==1
        FMech = [1 0];
    elseif strcmp(FaultType{sz},'reverse')==1
        FMech = [0 1];
    end
    %% Soil Term
    soil1 = [0 0]; % [Ss,Sa]
    if Vs30 <= 360
        soil1(1) = 1;
    elseif Vs30 > 360 && Vs30 <= 750
        soil1(2) = 1;
    end
    %%
    
    a1=0.458;
    a2=-0.0549;
    a3=1.046;
    a4=-0.0361;
    a5=-1.297;
    a6=-0.138;
    a7=0.105;
    delta = (1-exp(-(a1+a2*m)*Ri^(a3+a4*m)))*exp(a5+a6*m+a7*m^2);
    r = max([Ri - delta, 0]);
    
    gmp_central_estimate = 10.^(coef1(:,2)       + ...
        coef1(:,3) * m   +  ...
        coef1(:,4) * m^2 +  ...
        (coef1(:,5)+coef1(:,6)*m).*log10(sqrt(r^2+coef1(:,7).^2)) + ...
        coef1(:,8) * soil1(1) + ...
        coef1(:,9) * soil1(2) + ...
        coef1(:,10)*FMech(1)  + ...
        coef1(:,11)*FMech(2))/981; % (g)
    
    % Simulation of the expected intensity measure from the lognormal
    % distribution of the GMPE
    gmp_estimate = log10(gmp_central_estimate)+ coef1(:,14).*eps;
    
    % Store results
    catalog1(counttmp,:) = [ii time m lat lon Ri sz]; 
    catalog2(counttmp,:) = 10.^(gmp_estimate');
    catalog3(counttmp,:) = eps';
    
end
close(wb)
[hazard_values, ~] = sort(catalog2);
%% PSHA results
Final_time = Num;
lambda = flip((1:counttmp)./Final_time);

figure
loglog(hazard_values,lambda,'linewidth',2);
ylim([10^-4 2])
xlim([0.01 2])
xlabel('S_a - Spectral acceleration [g]'); ylabel('\lambda - mean annual rate'); 
grid on; axis square
set(gca,'FontSize',16,'Layer','Top')
hold on
plot(repmat([0.0001 3],length(Tr),1)',[1./Tr' 1./Tr']',':k','linewidth',2)
legend([repmat('T = ',length(Tnn'),1),num2str(Tnn')],'Location','northeastoutside');

%%
for i=1:length(Tr)
[~,index] = min(abs(lambda-1./Tr(i)));
Response_Spectra(i,:)=hazard_values(index,:);
end

figure
subplot(121)
plot(Tnn,Response_Spectra,'x-','linewidth',2)
xlabel('T (s)'); ylabel('S_a (g)'); 
grid on; axis square
set(gca,'FontSize',16,'Layer','Top')
legend(num2str(Tr'))

subplot(122)
loglog(Tnn,Response_Spectra,'x-','linewidth',2)
xlabel('T (s)'); ylabel('S_a (g)'); 
grid on; axis square
set(gca,'FontSize',16,'Layer','Top')
legend(num2str(Tr'))

%% Deaggregation
clear catalog10 catalog20
load slipcolor
clear PMF mag_bin_center dist_bin
% Disaggregation binning
dmag = 0.2; %0.125; 
ddist = 20; %10; 
Rmax = 300; %km

T = 3;
Sa = 0.01;

[~, ind_sort] = sort(catalog2(:,Tnn==T));
catalog10 = catalog1(ind_sort,:);
catalog20 = catalog2(ind_sort,:);

mag_bin = min(Mwmin):dmag:max(Mwmax);
for i = 1:length(mag_bin)-1
    mag_bin_center(1,i) = (mag_bin(1,i) + mag_bin(1,i+1))/2;
end

dist_bin = 0:ddist:Rmax;
for i = 1:length(dist_bin)-1
    dist_bin_center(1,i) = (dist_bin(1,i) + dist_bin(1,i+1))/2;
end

Bin_deagg_num = zeros(length(mag_bin)-1,length(dist_bin)-1);
Bin_deagg_den = zeros(length(mag_bin)-1,length(dist_bin)-1);

for ii = 1:length(mag_bin)-1
    for jj = 1:length(dist_bin)-1
%         Bin_deagg_den(ii,jj) = length(find(catalog1(:,3) >= mag_bin(ii) & catalog1(:,3) < mag_bin(ii+1) & catalog1(:,6) >= dist_bin(jj) & catalog1(:,6) < dist_bin(jj+1)));
        Bin_deagg_num(ii,jj) = length(find(catalog20(:,Tnn==T)>Sa & catalog10(:,3) >= mag_bin(ii) & catalog10(:,3) < mag_bin(ii+1) & catalog10(:,6) >= dist_bin(jj) & catalog10(:,6) < dist_bin(jj+1)));       
    end
end

% Define the probability mass function.
PMF = Bin_deagg_num./sum(sum(Bin_deagg_num));

figure
for lm = 1:length(mag_bin)-1
    for ld = 1:length(dist_bin)-1
        Xp = [dist_bin_center(ld)-ddist/2 dist_bin_center(ld)+ddist/2 dist_bin_center(ld)+ddist/2 dist_bin_center(ld)-ddist/2];
        Yp = [mag_bin_center(lm)-dmag/2 mag_bin_center(lm)-dmag/2 mag_bin_center(lm)+dmag/2 mag_bin_center(lm)+dmag/2];
        patch(Xp,Yp,PMF(lm,ld))
    end
end
colormap(slipcolor)
colorbar
axis tight
axis square
xlabel('Distance [km]')
ylabel('Magnitude')

figure
load slipcolor
b=bar3(mag_bin_center,PMF,'g');
xlabel('Distance [km]')
ylabel('Magnitude')
zlabel('PMF')
diff=dist_bin_center(2)-dist_bin_center(1);
for k = 1:length(b)
  xData = b(k).XData;
  zData = b(k).ZData;
  set(b(k), 'XData', (xData-k).*diff+dist_bin_center(k), ...
            'CData', zData, 'FaceColor', 'interp');
%    set(b(k), 'XData', (xData-k).*diff+dist_bin_center(k), ...
%             'CData', zData, 'FaceColor', 'interp');      
end
colormap(slipcolor)
colorbar
view([150.9 44.4]);
axis tight
% set('XTick',dist_bin_center,'XTickLabel',num2cell(dist_bin_center))
set(gca,'XDir',...
   'reverse','YDir','reverse','XTick',dist_bin_center,'XTickLabel',num2cell(dist_bin_center))
