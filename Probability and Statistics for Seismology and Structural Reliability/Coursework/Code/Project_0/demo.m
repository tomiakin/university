clear
clc
close all
% mu_log = log(mu^2/(sqrt(sigma^2 + mu^2))); 
% sigma_log = sqrt( log( (sigma^2)/(mu^2) + 1 ) );
%% Your data
% NORMAL s2
normal_like_s2 = sampleLike(chi_sq_norm_s2, df_s2);
% LOGNORMAL s2
lognormal_like_s2 = sampleLike(chi_sq_logn_s2, df_s2);
disp([' For s2 the sample likelihood for the normal is ' num2str(normal_like_s2) ', lognormal is ' num2str(lognormal_like_s2)] );



