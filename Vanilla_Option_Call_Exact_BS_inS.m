function [V S Delta Gamma] = Vanilla_Option_Call_Exact_BS_inS(sigma,dtau)
% function [V S ] = Vanilla_Option_Call_Exact_BS_inS(sigma)

% Parameters
r = 0.05;
% sigma = 0.2;
Ko = 100;
T = 1;

t = 0;
tau = T-t;
dx = 0.5;
Sint = 100;
Smin = 0;
Smax = exp(2)*100;
S = Smin:dx:Smax;
nx = length(S);
u  = zeros(nx,1);
% for i = 1:nx
%   d1 = 1/(sigma*sqrt(tau))*(x(i) + 0.5*sigma^2*tau + 0.5*sigma^2*tau);
%   d2 = 1/(sigma*sqrt(tau))*(x(i) + 0.5*sigma^2*tau - 0.5*sigma^2*tau);
%   N_1_der = 1/sqrt(2*pi)*exp(-1/2.*d1);
%   Gamma(i,1) = N_1_der/(x(i)*sigma*sqrt(tau));
%   N1 = normcdf(d1);
%   NN1(i,1) = normcdf(d1);
%   N2 = normcdf(d2);
%   u(i,1) = Ko*exp(x(i)+0.5*sigma^2*tau)*N1 - Ko*N2;
% end
% V = exp(-r*tau)*u;
% S = Ko*exp(x-(r-0.5*sigma^2)*tau);
% Delta = NN1;

% Corresponding exact solution for the transform x = ln(S/Sint); \tau = T-t;
for i = 1:nx
  d1 = 1/(sigma*sqrt(tau))*(log(S(i)/Ko) + (r + sigma^2/2)*tau);

  d2 = d1 - sigma*sqrt(tau);

    N_1_der = 1/sqrt(2*pi)*exp(-1/2*d1^2);
   Gamma(i,1) = N_1_der/(S(i)*sigma*sqrt(tau));
  N1 = normcdf(d1);
   NN1(i,1) = normcdf(d1);
  N2 = normcdf(d2);
  u(i,1) = S(i)*N1 - Ko*exp(-r*tau)*N2;
end
V = u;
Delta = NN1;
% figure;
% plot(S,V); axis([0 300 0 200]); hold on
% %
% % Pay-off
% Vpo = max(S-Ko,0);
% %figure;
% plot(S,Vpo);
return