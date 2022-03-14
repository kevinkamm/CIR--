function Psi=determShift(params,P0TMarket,t,T)
%DETERMSHIFT calculates the deterministic shift for CIR-
%    Input:
%       params (8x1 array): params= $[\phi_1^x,...,\phi_1^y,...\phi_3^y,x_0,y_0]$
%       P0TMarket (gridded interpolant): linear interpolation of market
%                                        discount factor
%       t (p x 1 array): contains the current times
%       T (1 x q array): contains the maturities
%    Output:
%       Psi (p x q array): deterministic shift
%
% See also P0T_Market, P0T_CIR1.

if sum(t,'all')==0
    Psi = P0TMarket(T)./P0T_CIR1(params,T);
else
    Psi = P0TMarket(T)./P0TMarket(t) .*...
          P0T_CIR1(params,t)./P0T_CIR1(params,T);
end
end