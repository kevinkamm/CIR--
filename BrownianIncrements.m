function [dW1,dW2] = BrownianIncrements(T,N,M)
%%BROWNIANINCREMENTS computes two independent increments of Brownian
% motions with M paths and N-1 time steps.
%   Input:
%       T (1 x 1 double): terminal time
%       N (1 x 1 int): number of time increments
%       M (1 x 1 int): number of simulations
%   Output:
%       dW1 (N-1 x M array): Brownian increment for CIR process x
%       dW2 (N-1 x M array): Brownian increment for CIR process y

% time step
dt=T/(N-1);

% 2 Brownian increments
dW = sqrt(dt).*randn(N-1,2*M);
dW1 = dW(:,1:M);
dW2 = dW(:,M+1:2*M);

end