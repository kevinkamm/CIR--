function P0TMarket=P0T_Market(marketTimes,marketDF,method)
%%P0T_Market creates function of the gridded interpolation of the market 
% discount factors evaluated for future evaluation.
%   Input:
%       marketTimes (p x 1 array): maturities of market discount factors
%       marketDF (p x 1 array): values of market discount factors
%       method (string): 'linear' (default), 'spline'
%   Output:
%       P0TMarket (function):
%   Usage:
%       P0T_Market(marketTimes,marketDF): uses method='linear'
%       P0T_Market(marketTimes,marketDF,method): uses given method
%
% See also griddedInterpolant.

if nargin<3
    method='linear';
end
    P0TMarket=griddedInterpolant(marketTimes,marketDF,method);
end