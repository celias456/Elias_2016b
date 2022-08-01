function [ree] = solve_ree(lambda,rho_f,delta,alpha) 
% Calculates REE solutions
% Input:
%   lambda: Interest rate semi-elasticity of Money Demand
%   alpha: Fundamentals equation constant
%   delta: Fundamentals equation parameter on time
%   rho_f: Fundamentals equation parameter on lagged fundamentals term
%   sigma_epsilon: Standard deviation of white noise term in Fundamentals equation

% Output:
%   ree: 3x1 matrix of ree solutions
%     ree(1) = REE constant (agent 1)
%     ree(2) = REE coefficient on T (time) for agent 1
%     ree(3) = REE coefficient on f_t (fundamentals) for agent 1
% -------------------------------------------------------------------------


% Find REE solutions
c = 1/(1+lambda-lambda*rho_f); %coefficient on f_t
b = lambda*c*delta; %coefficient on time
a = lambda*(b+c*(alpha+delta)); %constant

ree = [a;b;c];


end   