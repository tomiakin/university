function pdf = lognormal_pdf(x,eta,beta,base)

% base can be either exp(1) or 10

if base == exp(1)
    pdf = 1./(x.*beta.*sqrt(2*pi)).*exp(-(log(x)-log(eta)).^2./(2.*beta^2));
    %  pdf_check = normpdf(log(x),log(eta),beta)./x;
elseif base ==10
    pdf = 1./(x.*beta.*sqrt(2*pi)).*exp(-(log10(x)-log10(eta)).^2./(2.*beta^2));
%     pdf_check = normpdf(log(x)./log(10),log(eta)./log(10),beta)./x;
end

end