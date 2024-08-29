function RJB = distance_converter(R_epicentral)
    RJB = -3.5525+0.8845.*R_epicentral;
    RJB(RJB<0)=0;
end