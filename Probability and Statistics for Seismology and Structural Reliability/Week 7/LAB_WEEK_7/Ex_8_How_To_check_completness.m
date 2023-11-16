%% This script check the completness for each zone in the
%  sub catalog

Mmin = 4;

load SUB_CATALOG

for j=1:length(SUB_CATALOG) % HENCE WHY THERE ARE 3 PLOTS
    
    magnitude    = SUB_CATALOG(j).CATALOG.Mw; 
    year         = SUB_CATALOG(j).CATALOG.Year;
    minimum_year = min(year);
    maximum_year = max(year);
    
    T=minimum_year:5:maximum_year;
    
    n = zeros(length(T),1);
    for i=2:length(T)
        n(i)=length(find(magnitude >= Mmin & year>T(i-1) & year<T(i)));
    end
    
    figure
    title(['Zone ',num2str(j)])
    plot(T,cumsum(n),'k','linewidth',2)
    grid on
    axis square; xlabel('Time (Years)'); ylabel('Cumulative number');
    set(gca,'FontSize',16)
end