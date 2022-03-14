function strikeSwap = atmStrikes(maturitySwap,tenorSwap,P0TMarket)
%%ATMSTRIKES computes the at-the-money strikes for given market zero curve.
% We assume annual settlements.
%    Input:
%       maturitySwap (p x 1 array): maturities of swaption
%       tenorSwap (1 x q array): tenors of swaption
%       P0TMarket (griddedInterpolant): interpolated market zero coupon 
%                                       prices
%    Output:
%       strikeSwap (p x q array): atm strikes of swaption 

% we assume annual payments 
paymentTimes=1:1:tenorSwap(end);
TN=maturitySwap+paymentTimes;
Teval=[maturitySwap,TN];

zcMarket=P0TMarket(Teval);

% foward swap rate = atm strikes
strikeSwap=(zcMarket(:,1)-zcMarket(:,2:end))./(cumsum(zcMarket(:,2:end),2));
strikeSwap=strikeSwap(:,tenorSwap);
end