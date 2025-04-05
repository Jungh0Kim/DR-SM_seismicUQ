% Jungho Kim
% Compute lognormal parameters based on Mu, Std

function [lambda, zeta] = Logn_para(Mu, Std)

lambda = log(Mu^2 / sqrt( Std^2 + Mu^2));
zeta = sqrt(log( Std^2 / Mu^2 + 1 ));

% [Mu_check, Var_check] = lognstat(lambda,zeta);    % check

end % function end


