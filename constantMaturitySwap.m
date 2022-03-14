function [price,varargout]=constantMaturitySwap(params,x,y,modelTimes,dfCIR1,...
                               K,Tn,TN,c,...
                               P0TMarket,...
                               swapType,...
                               varargin)
%CONSTANTMATURITYSWAP calculates today's price of a CMS with first exercise
% T_n, maturity T_N and index of the underlying swap rates c.
% If K is not empty, price will be computed with given fixed rate,
% otherwise the CMS rate will be computed and the price should be near zero
% prices.
%    Input:
%       params (8x1 array): params= $[\phi_1^x,...,\phi_1^y,...\phi_3^y,x_0,y_0]$
%       x (N x M array): contains the simulated paths of the CIR process x
%       y (N x M array): contains the simulated paths of the CIR process y
%       modelTimes (N x 1 array): contains the timeline of the model
%       dfCIR1 (N x M array): contains the CIR- discount factors
%       K (double or empty): strikes(= forward swap rates) of the 
%                                 swaption contract
%       Tn (double): contains the maturities of the contract
%       TN (double): contains the time till expiry
%       P0TMarket (gridded interpolant): linear interpolation of market
%                                        discount factor
%       swapType (int): equal to 1 for payer swaption and -1 for receiver
%    Output:
%       price (p x q array): expectation of the swaption prices
%
% See also P0T_Market.

% TODO vectorize in Tn,TN,c

%% Initialize
% settlement dates
T=Tn:1:TN;
% compute discount factor with determ shift
dfCIR2=determShift(params,P0TMarket,0,modelTimes)'.*dfCIR1;

    function [SnN,RnN]=swapRate(Tn,TN)
        %SWAPRATE computes accrual factor and forward swap rate assuming
        % annual payments.
        %   Input:
        %       Tn (double): t=Tn initial time
        %       TN (double): last payment date
        %   Output:
        %       SnN (1 x M array): simulated accrual factors
        %       RnN (1 x M array): simulated swap rates
        
        indTn=find(modelTimes<=Tn,1,'last');
        Tpayment=Tn+1:1:TN;
        zcPriceCIR1 = PtT_CIR1(params,Tn,Tpayment,indTn,x,y);
        determShift_tT = determShift(params,P0TMarket,Tn,Tpayment);
        zcPriceCIR2 = determShift_tT .* zcPriceCIR1;
        
        % accrual factor or present value of a basis point
        SnN=sum(zcPriceCIR2,2);

        % forward swap rate
        RnN=squeeze((1-zcPriceCIR2(:,end,:))./SnN);
        SnN=squeeze(SnN);
    end
%%
if isempty(K)
    K=0;
    for iT = 1:1:length(T)-1
        [~,RnN]=swapRate(T(iT),T(iT)+c);
        indT=find(modelTimes<=T(iT),1,'last');
        K=K+dfCIR2(indT,:).*RnN';
    end 
    K=mean(K,2)./sum(P0TMarket(T(1:end-1)));
end
price=0;
for iT = 1:1:length(T)-1
    [~,RnN]=swapRate(T(iT),T(iT)+c);
    indT=find(modelTimes<=T(iT),1,'last');
    price=price+dfCIR2(indT,:).*swapType.*(RnN'-K);
end

price=mean(price,2);
if nargout>1
    varargout{1}=K;
end

end