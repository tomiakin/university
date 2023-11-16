function Mw = magnitude_converter(value,M_Type)
% https://en.wikipedia.org/wiki/Seismic_magnitude_scales
% https://earthquake.usgs.gov/data/comcat/data-eventterms.php#magType
% https://www.usgs.gov/natural-hazards/earthquake-hazards/science/magnitude-types?qt-science_center_objects=0#qt-science_center_objects

%%
if strcmp(M_Type,'Mw') || strcmp(M_Type,'MW') || strcmp(M_Type,'Mww') || strcmp(M_Type,'Mwc') || strcmp(M_Type,'mwc') || strcmp(M_Type,'Mwb') || strcmp(M_Type,'mwb') || strcmp(M_Type,'Mwr') || strcmp(M_Type,'mwr') || strcmp(M_Type,'Mw(mB)') || strcmp(M_Type,'MwMwp') || strcmp(M_Type,'Mwp') || strcmp(M_Type,'mw') || strcmp(M_Type,'mww') || strcmp(M_Type,'Mjma')
    Mw = value;
elseif strcmp(M_Type,'Md') || strcmp(M_Type,'md') || strcmp(M_Type,'MD') || strcmp(M_Type,'Mc')
    Mw = 0.93 * value + 0.35;
elseif strcmp(M_Type,'ML') || strcmp(M_Type,'Ml') || strcmp(M_Type,'M1') || strcmp(M_Type,'ml') || strcmp(M_Type,'mL') || strcmp(M_Type,'MLV') || strcmp(M_Type,'MLv')
    Mw1 = value + 0.17;
    Ms  = 1.27*(value - 1) - 0.016*value^2;
    if Ms > 6.2
        Mw2 = 0.99 * Ms + 0.08;
    else
        Mw2 = 0.67 * Ms + 2.07;
    end
    Mw = mean([Mw1, Mw2]);
elseif strcmp(M_Type,'Ms') || strcmp(M_Type,'MS') || strcmp(M_Type,'MSZ') || strcmp(M_Type,'Ms20') || strcmp(M_Type,'Mh') || strcmp(M_Type,'uk') || strcmp(M_Type,'UK') || strcmp(M_Type,'m') || strcmp(M_Type,'M') || strcmp(M_Type,'') || strcmp(M_Type,'MG') || strcmp(M_Type,'ms') || strcmp(M_Type,'ms1mx') || strcmp(M_Type,'msmle') || strcmp(M_Type,'mpv') || strcmp(M_Type,'Mm') || strcmp(M_Type,'Msl') || strcmp(M_Type,'Ms7') || strcmp(M_Type,'Ms_20') || strcmp(M_Type,'Msz') || strcmp(M_Type,'Ms1') || strcmp(M_Type,'mmm') || strcmp(M_Type,'Mt') || strcmp(M_Type,'mt') || strcmp(M_Type,'Mws')
    if value > 6.2
        Mw = 0.99 * value + 0.08;
    else
        Mw = 0.67 * value + 2.07;
    end
elseif strcmp(M_Type,'Me') || strcmp(M_Type,'ME')
    E = 10^((value+2.9)*3/2);
    Ms = (log10(E)-11.8)/1.5;
    if Ms > 6.2
        Mw = 0.99 * Ms + 0.08;
    else
        Mw = 0.67 * Ms + 2.07;
    end
elseif strcmp(M_Type,'Mi') || strcmp(M_Type,'Mwp')
    Mw = 0.89 * value + 0.57;
elseif strcmp(M_Type,'Mb') || strcmp(M_Type,'mB') || strcmp(M_Type,'MB') || strcmp(M_Type,'mb') || strcmp(M_Type,'Mlg') || strcmp(M_Type,'mb_Lg') || strcmp(M_Type,'mb_lg') || strcmp(M_Type,'fa') || strcmp(M_Type,'Mfa') || strcmp(M_Type,'mbLg') || strcmp(M_Type,'mblg') || strcmp(M_Type,'MN') || strcmp(M_Type,'mb1') || strcmp(M_Type,'mb1mx') || strcmp(M_Type,'mbmle') || strcmp(M_Type,'mbtmp') || strcmp(M_Type,'mB_BB') || strcmp(M_Type,'Mn') || strcmp(M_Type,'Mbn') 
    Mw1 = 0.85 * value + 1.03;
    Mw2 = 1.21 * value - 0.76;
    Mw = mean([Mw1, Mw2]);
end
end