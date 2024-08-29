%% Initialisation of the workspace
clear % the "clear" command clear the workspace (imporant to clean the memory)
clc   % the "clc" clear the window command from previous command (only visual)
close all

%% POINT OF INTEREST and value of the VS30 of the soil
% Point 1 - Inside one of the seismic zone
X = [35.36	25.12]; %[Latitude Longitude]

% To acquire the Vs30: https://earthquake.usgs.gov/data/vs30/
Vs30 = 300; % m/s shear wave velocity of the soil under the site of interest RANGE WAS 300-360
%% vector of IM to consider for the hazard curve
imx=0:0.01:4; %g

%% Vibration periods (T) to be considered for the spectral accelerations
% Tnn   = [0 0.1 0.2 0.3 0.4 0.5 0.75 1.0 1.5 2.0 3.0]; % seconds 
Tnn   = [0];

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

FaultType = {'normal', 'normal'}; % it can be 'normal', 'reverse', 'strike-slip'. Here you need to repleace your types ONE TYPE FOR EACH ZONE

%% Discretisation of the mangnitude and of the geometry
deltaM   = [0.1 0.1];     % Magnitude discretization (the smaller the slower)
DeltaX   = [25 25];        % Number of division for the distance (the larger the slower)

%% loading the coefficients for the GMPE for all vibration periods
load GMPEcoef_AB10
coef1 = zeros(length(Tnn),length(coef_AB10(1,:)));
for ii = 1:length(Tnn)
    coef1(ii,:) = coef_AB10(Tnn(ii)==coef_AB10(:,1),:);
end

%% Initialisation of the output vector
Sa_out = zeros(length(Tnn),length(Tr));

%% Double for-loop for each seismic zone 
Dea(length(Tnn)).Deaggregation(length(a)).D = [];

for ti = 1:length(Tnn)
    Tn = Tnn(ti);
    disp(['*** Period = ',num2str(Tn),' s ***'])
    %%
    HAZARD_CURVES_SZ        = zeros(length(X),length(imx));
    HAZARD_CURVES_SZ_weight = zeros(length(X),length(imx));
    
    Lambda_min = zeros(length(X),1);
    
    for i=1:length(a)
        disp(['*** SEISMIC ZONE = ',num2str(i),' ***'])
        %% Discretisation of G&R
        Mmax = Mwmax(i);
        Mmin = Mwmin(i);
        A    = a(i);
        B    = b(i);
        mm = linspace(Mmin,Mmax,100);
        G  = (1-10.^(-B*(mm-Mmin)))./(1-10.^(-B*(Mmax-Mmin)));
        Mi  = Mmin + deltaM(i)/2 : deltaM(i) : Mmax - deltaM(i)/2;
        PMi = zeros(size(Mi));
        for j=1:length(Mi)
            PMi(j) = interp1(mm',G',Mi(j)+deltaM(i)/2) - interp1(mm',G',Mi(j)-deltaM(i)/2); % pdf of the magnitude
        end
        clear mm G j
        
        %% Discretisation of the seismic area and calculation of distances
        L = max(SZ(i).X)-min(SZ(i).X);
        H = max(SZ(i).Y)-min(SZ(i).Y);
        
        Delta_x = DeltaX(i);
        Delta_y = round(H/L*Delta_x);
        
        xg = min(SZ(i).X):(max(SZ(i).X)-min(SZ(i).X))/Delta_x:max(SZ(i).X);
        yg = min(SZ(i).Y):(max(SZ(i).Y)-min(SZ(i).Y))/Delta_y:max(SZ(i).Y);
        
        XG = zeros(length(xg)*length(yg),1);
        YG = zeros(length(xg)*length(yg),1);
        
        counter = 0;
        
        clear ii jj
        for ii=1:length(xg)
            for jj=1:length(yg)
                [in,on]=inpolygon(xg(ii),yg(jj),SZ(i).X,SZ(i).Y);
                if in>0 || on>0
                    counter = counter + 1;
                    XG(counter)=xg(ii);
                    YG(counter)=yg(jj);
                end
            end
        end
        XG(counter+1:end,:)=[];
        YG(counter+1:end,:)=[];
        
        % Distance on the Earth surface in degrees
        arclen = distance(YG,XG,X(1),X(2));
        % Conversion from degrees to km and from Epicentral Distance to Rjb
        Ri  = distance_converter(deg2km(arclen));
        
        % PMF of the distances (equal probability)
        PRi = 1/length(Ri); % pdf of the distance 
        
        clear L H Delta_x Delta_y xg yg ii jj in on counter arclen
        
        %% Fault mechanism term
        Fault_Type = FaultType(i);
        if strcmp(Fault_Type,'strike-slip')==1
            FMech = [0 0];
        elseif strcmp(Fault_Type,'normal')==1
            FMech = [1 0];
        elseif strcmp(Fault_Type,'reverse')==1
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
        counter = 0;
        
        Dgmpe=zeros(length(Ri),length(Mi),length(imx)); % initialisation deaggregation
        
        for ir=1:length(Ri)
            for jm=1:length(Mi)
                counter = counter+1;
                
                m = Mi(jm);
                
                % Function to prepare the distance for the GMPE
                a1=0.458;
                a2=-0.0549;
                a3=1.046;
                a4=-0.0361;
                a5=-1.297;
                a6=-0.138;
                a7=0.105;
                delta = (1-exp(-(a1+a2*m)*Ri(ir)^(a3+a4*m)))*exp(a5+a6*m+a7*m^2);
                r = max([Ri(ir) - delta, 0]);
                
                % Application of the GMPE
                gmp_central_estimate = 10.^(coef1(Tnn==Tn,2)       + ...
                    coef1(Tnn==Tn,3) * m   +  ...
                    coef1(Tnn==Tn,4) * m^2 +  ...
                    (coef1(Tnn==Tn,5)+coef1(Tnn==Tn,6)*m).*log10(sqrt(r^2+coef1(Tnn==Tn,7).^2)) + ...
                    coef1(Tnn==Tn,8) * soil1(1) + ...
                    coef1(Tnn==Tn,9) * soil1(2) + ...
                    coef1(Tnn==Tn,10)*FMech(1)  + ...
                    coef1(Tnn==Tn,11)*FMech(2))/981; % (g)
                
                if counter==1
                    Pgmpe=(1-normcdf(log10(imx),log10(gmp_central_estimate),coef1(Tnn==Tn,14))).*PRi.*PMi(jm);
                else
                    Pgmpe=Pgmpe+(1-normcdf(log10(imx),log10(gmp_central_estimate),coef1(Tnn==Tn,14))).*PRi.*PMi(jm);
                end
                
                Dgmpe(ir,jm,:)=normpdf(log(imx)./log(10),log(gmp_central_estimate)./log(10),coef1(Tnn==Tn,14))./imx.*PRi.*PMi(jm);
            end
        end
        
        Lambda_min(i) = round(10^(A-B*Mmin)*100)/100; % *100 and /100 is to round to the first 2 significative digits
        HAZARD_CURVES_SZ(i,:) = Pgmpe;
        
        if i==1
            HAZARD = HAZARD_CURVES_SZ(i,:)*Lambda_min(i);
        else
            HAZARD = HAZARD + HAZARD_CURVES_SZ(i,:)*Lambda_min(i);
        end
        Dea(ti).Deaggregation(i).R=Ri;
        Dea(ti).Deaggregation(i).M=Mi;
        Dea(ti).Deaggregation(i).D=Dgmpe*Lambda_min(i);
    end
    %%
    loglog(imx,HAZARD,'linewidth',2)
    hold all
    Sa_out(ti,:) = interp1(HAZARD,imx,1./Tr);
    clear HAZARD HAZARD_CURVES_SZ
    ylim([10^-4 2])
    xlim([0.01 2])
    xlabel('S_a - Spectral acceleration [g]'); ylabel('\lambda - mean annual rate'); 
    grid on; axis square
    set(gca,'FontSize',16,'Layer','Top')
end
hold on
plot(repmat([0.0001 3],length(Tr),1)',[1./Tr' 1./Tr']',':k','linewidth',2)

legend([repmat('T = ',length(Tnn'),1),num2str(Tnn')],'Location','northeastoutside');

%% Plot of the Uniform hazard spectra
figure
subplot(121)
plot(Tnn,Sa_out,'x-','linewidth',2)
xlabel('T (s)'); ylabel('S_a (g)'); 
grid on; axis square
set(gca,'FontSize',16,'Layer','Top')
legend(num2str(Tr'))

subplot(122)
loglog(Tnn,Sa_out,'x-','linewidth',2)
xlabel('T (s)'); ylabel('S_a (g)'); 
grid on; axis square
set(gca,'FontSize',16,'Layer','Top')
legend(num2str(Tr'))

%% Computation and plot of the deaggregation
figure
load slipcolor
% Definition of the value for which we want to perform the deaggregation
Tdea = 3; % seconds
index1=find(Tnn==Tdea);
Sadea = 0.01; % g
[~,index2]=min(abs(imx-Sadea));

% Disaggregation binning
dmag = 0.2; %0.125; 
ddist = 20; %10; 
Rmax = 300; %km

mag_bin = min(Mwmin):dmag:max(Mwmax);
mag_bin_center = zeros(1,length(mag_bin)-1);
for i = 1:length(mag_bin)-1
    mag_bin_center(1,i) = (mag_bin(1,i) + mag_bin(1,i+1))/2;
end

dist_bin = 0:ddist:Rmax;
dist_bin_center = zeros(1,length(dist_bin)-1);
for i = 1:length(dist_bin)-1
    dist_bin_center(1,i) = (dist_bin(1,i) + dist_bin(1,i+1))/2;
end

Bin_deagg = zeros(length(mag_bin)-1,length(dist_bin)-1);

for ii = 1:length(mag_bin)-1
    for jj = 1:length(dist_bin)-1
        
        for i=1:length(a)
            RR=Dea(index1).Deaggregation(i).R;
            MM=Dea(index1).Deaggregation(i).M;
            DD=Dea(index1).Deaggregation(i).D;
            if i==1
                Bin_deagg(ii,jj) = sum(sum(DD(RR >= dist_bin(jj) & RR < dist_bin(jj+1), MM >= mag_bin(ii) & MM < mag_bin(ii+1),index2))); %equal probability case
            else
                Bin_deagg(ii,jj) = Bin_deagg(ii,jj) + sum(sum(DD(RR >= dist_bin(jj) & RR < dist_bin(jj+1), MM >= mag_bin(ii) & MM < mag_bin(ii+1),index2))); %equal probability case
            end
        end
        
    end
end

for lm = 1:length(mag_bin)-1
    for ld = 1:length(dist_bin)-1
        pdf = Bin_deagg(lm,ld)./sum(sum(Bin_deagg));
        Xp = [dist_bin_center(ld)-ddist/2 dist_bin_center(ld)+ddist/2 dist_bin_center(ld)+ddist/2 dist_bin_center(ld)-ddist/2];
        Yp = [mag_bin_center(lm)-dmag/2 mag_bin_center(lm)-dmag/2 mag_bin_center(lm)+dmag/2 mag_bin_center(lm)+dmag/2];
        patch(Xp,Yp,pdf)
    end
end
colormap(slipcolor)
colorbar
axis tight
axis square
xlabel('Distance [km]')
ylabel('Magnitude')
%%
figure
load slipcolor
b=bar3(mag_bin_center,Bin_deagg./sum(sum(Bin_deagg)),'g');
xlabel('Distance [km]')
ylabel('Magnitude')
zlabel('PMF')
diff=dist_bin_center(2)-dist_bin_center(1);
for k = 1:length(b)
  xData = b(k).XData;
  zData = b(k).ZData;
%   set(b(k), 'XData', (xData-k).*diff+dist_bin_center(k), ...
%             'CData', zData, 'FaceColor', 'interp');
   set(b(k), 'XData', (xData-k).*diff+dist_bin_center(k), ...
            'CData', zData, 'FaceColor', 'interp');      
end
colormap(slipcolor)
colorbar
view([150.9 44.4]);
axis tight
set(gca,'XDir',...
   'reverse','YDir','reverse','XTick',dist_bin_center,'XTickLabel',num2cell(dist_bin_center))
