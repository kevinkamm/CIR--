function [x,y,dfCIR1] = sim_CIR1(params,T,dW1,dW2)
%SIM_CIR1 simulates the CIR- model for given parameters using Euler
% scheme with \sqrt{\max(r_{i-1},0)}
%    Input:
%       params (8 x 1 array): $[\phi_1^x,...,\phi_1^y,...\phi_3^y,x_t0,y_t0]$
%       t0 (double): initial time
%       T (double): terminal time
%    Output:
%       x (N x M array): the simulated paths of the CIR process x in 
%                        r=x-y
%       y (N x M array): the simulated paths of the CIR process y in 
%                        r=x-y
%       dfCIR1 (N x M array): the CIR- model discount factor 
%                             $D(0,t_i)=\exp(\int_0^{t_i}{r_s ds})$ 
%                             using trapezoidal rule

% number of time steps and simulations, resp.
N = size(dW1,1)+1; % plus 1 because increments have one time step less
M = size(dW1,2);

% parameters for $x_t$
phi1x = params(1);
phi2x = params(2);
phi3x = params(3);
x0    = params(7);

% parameters for $y_t$
phi1y = params(4);
phi2y = params(5);
phi3y = params(6);
y0    = params(8);

% retrieve coeffiecents for $x_t$
kx     = 2*phi2x - phi1x;
sigmaxSquared = 2*(phi1x*phi2x - phi2x^2);
thetax = phi3x*sigmaxSquared/(2*kx);
sigmax = sqrt(sigmaxSquared);

% retrieve coeffiecents for $y_t$
ky     = 2*phi2y - phi1y;
sigmaySquared = 2*(phi2y^2 - phi1y*phi2y);
thetay = phi3y*sigmaySquared/(2*ky);
sigmay = sqrt(sigmaySquared);

% time step
dt=T/(N-1);

% solution arrays
x = zeros(N,M);
y = zeros(N,M);

% initial values
x(1,:) = x0;
y(1,:) = y0;

for i = 2:1:N
    x(i,:) = x(i-1,:) +...
             kx.*(thetax - x(i-1,:)).*dt +...
             (sigmax.*sqrt(max(x(i-1,:),0))).*dW1(i-1,:);
    y(i,:) = y(i-1,:) +...
             ky.*(thetay - y(i-1,:)).*dt +...
             (sigmay.*sqrt(max(y(i-1,:),0))).*dW2(i-1,:);  
end
r=x-y;
dfCIR1=exp(-cumtrapz(r,1).*dt);
% dfCIR1=exp(-cumsum(r,1).*dt);