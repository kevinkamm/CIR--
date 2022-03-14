function price=swaption_matrix(params,x,y,modelTimes,dfCIR1,...
                               strikeSwap,maturitySwap,tenorSwap,...
                               P0TMarket,...
                               swapType,...
                               varargin)
%SWAPTION_MATRIX calculates the expectation of the simulated swaption
% prices
%    Input:
%       params (8x1 array): params= $[\phi_1^x,...,\phi_1^y,...\phi_3^y,x_0,y_0]$
%       x (N x M array): contains the simulated paths of the CIR process x
%       y (N x M array): contains the simulated paths of the CIR process y
%       modelTimes (N x 1 array): contains the timeline of the model
%       dfCIR1 (N x M array): contains the CIR- discount factors
%       strikeSwap (p x q array): strikes(= forward swap rates) of the 
%                                 swaption contract
%       maturitySwap (p x 1 array): contains the maturities of the contract
%       tenorSwap (1 x q array): contains the time to length of time till
%                                expiry
%       P0TMarket (gridded interpolant): linear interpolation of market
%                                        discount factor
%       swapType (int): equal to 1 for payer swaption and -1 for receiver
%    Output:
%       price (p x q array): expectation of the swaption prices
%
% See also P0T_Market.

% swaptionType='payer';
% for kk=1:2:length(varargin)
%     switch varargin{kk}
%         case 'swaptionType'
%             swaptionType=varargin{kk+1};
%     end
% end
% we assume annual payments 
paymentTimes=1:1:tenorSwap(end);

% for setting correct tau: T=maturity+tenor
T=maturitySwap+paymentTimes;
% calculate Zero coupon prices
matInd=zeros(length(maturitySwap),1);
for i=1:1:length(matInd)
    matInd(i) = find(modelTimes<=maturitySwap(i),1,'last');
end
% find tenor indices in paymentTimes
tenInd=tenorSwap; % true because annual payments

% zero coupon price for CIR--: maturity till expiry
zcPriceCIR1 = PtT_CIR1(params,maturitySwap,T,matInd,x,y);
determShift_tT = determShift(params,P0TMarket,maturitySwap,T);
zcPriceCIR2 = determShift_tT .* zcPriceCIR1;

% discount factors for CIR--: today till maturity
determShift_0T0 = determShift(params,P0TMarket,0,maturitySwap);
dfCIR2 = determShift_0T0 .* reshape(dfCIR1(matInd,:),[],1,size(x,2));

% accrual factor or present value of a basis point
SnN=cumsum(zcPriceCIR2,2);

% forward swap rate
swapRate=(1-zcPriceCIR2)./SnN;

K=strikeSwap;
% compute discounted swaption price
price=mean(dfCIR2.*max(0,swapType.*(swapRate(:,tenInd,:)-K)).*SnN(:,tenInd,:),3);
% switch swaptionType
%     case 'receiver'
%     % receiver swaption
%     price=mean(dfCIR2.*max(0,K-swapRate(:,tenInd,:)).*SnN(:,tenInd,:),3);
%     case 'payer'
%     % payer swaption
%     price=mean(dfCIR2.*max(0,swapRate(:,tenInd,:)-K).*SnN(:,tenInd,:),3);
%     otherwise 
%         error('Unknown swaption type: %s given, needs to be receiver or payer',swaptionType)
% end

end