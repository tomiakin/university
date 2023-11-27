%% This script check the completness for each zone in the for different MAGNITUDES
%  sub catalog
load SUB_CATALOG

Mmin_vector = [1 2 3 4 5 6];

for k=1:length(Mmin_vector)
    
    Mmin = Mmin_vector(k);
    
    for j=1:length(SUB_CATALOG)
        
        magnitude    = SUB_CATALOG(j).CATALOG.Mw;
        year         = SUB_CATALOG(j).CATALOG.Year;
        minimum_year = min(year);
        maximum_year = max(year);
        
        T=minimum_year:5:maximum_year;
        
        n = zeros(length(T),1);
        for i=2:length(T)
            n(i)=length(find(magnitude >= Mmin & year>T(i-1) & year<T(i)));
        end
        
        figure(j)
        title(['Zone ',num2str(j)])
        hold on
        loglog(T,cumsum(n),'linewidth',2)
        grid on
        axis square; xlabel('Time (Years)'); ylabel('Cumulative number');
        set(gca,'FontSize',16)
        if k==length(Mmin_vector)
            legend({num2str(Mmin_vector')},'location','northwest')
            box on
        end
    end
    
end