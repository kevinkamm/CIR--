function price=bermudanSwaption(params,x,y,modelTimes,dfCIR1,...
                               K,Tn,TN,...
                               P0TMarket,...
                               swapType,...
                               varargin)
%BERMUDANSWAPTION calculates today's price of a "T_N nc T_n" Bermudan
% swaption with fixed rate K and annual exercise dates.
% prices
%    Input:
%       params (8x1 array): params= $[\phi_1^x,...,\phi_1^y,...\phi_3^y,x_0,y_0]$
%       x (N x M array): contains the simulated paths of the CIR process x
%       y (N x M array): contains the simulated paths of the CIR process y
%       modelTimes (N x 1 array): contains the timeline of the model
%       dfCIR1 (N x M array): contains the CIR- discount factors
%       K (double): strikes(= forward swap rates) of the 
%                                 swaption contract
%       Tn (double): contains the maturities of the contract
%       TN (double): contains the time to expiry
%       P0TMarket (gridded interpolant): linear interpolation of market
%                                        discount factor
%       swapType (int): equal to 1 for payer swaption and -1 for receiver
%    Output:
%       price (p x q array): expectation of the swaption prices
%
% See also P0T_Market.

%% Initialize
% early exercise dates
TE=Tn:1:TN-1;
% compute discount factor with determ shift
dfCIR2=determShift(params,P0TMarket,0,modelTimes)'.*dfCIR1;

%% backward induction via LSMC
Nb=floor(log(size(x,2))/2);% number of accurate basis functions is O(log M)
T=[0,TE];
for i=length(T):-1:2
    tInd0=find(modelTimes<=T(i-1),1,'last');
    tInd1=find(modelTimes<=T(i),1,'last');
    [SnN,RnN]=swapRate(T(i),TN);
    payoff=max(swapType.*(RnN-K),0).*SnN;
    if i==length(T)
        % discount final payoff for next iteration
        % SnN equal to P(T_{N-1},T_N)
        % RnN equal to (1-P(T_{N-1},T_N))/P(T_{N-1},T_N)
        V=payoff;
    else
        % a decision is only required while in-the-money
        ITMind=find(payoff>0);
        
        if isempty(ITMind)
            % no update
            % should most likely never occur
            continue;
        end
        % continuation value via LSMC
        B=basis(RnN(ITMind),Nb);
        b=V(ITMind);
        lambda=B\b;
        c=B*lambda;
        % compare continuation to current payoff
        exerciseInd=find(payoff(ITMind)>c);
        ITMexerciseInd=ITMind(exerciseInd);
        V(ITMind(exerciseInd))=payoff(ITMexerciseInd);
%         V(ITMind)=max(payoff(ITMind),c);
    end
    %discount value for next iteration
    V=(dfCIR2(tInd1,:)./dfCIR2(tInd0,:))'.*V;
end
price = mean(V,1);
%% swap rate
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
%% basis functions: polynomials
    function B=basis(x,Nb)
       %BASIS evaluates the first Nb basis functions at x.
       % We choose the polynomials as basis functions.
       %    Input:
       %        x (ndarray): points of evaluation
       %        Nb (int): number of polynomials
       x=reshape(x,[1,size(x)]);
       B=ones([Nb,size(x)]);
       B(2:end,:)=repmat(x,[Nb-1,ones(1,sum(ndims(x)>1))]);
       B=cumprod(B,1);
       B=permute(B,[2:ndims(B),1]);
       B=squeeze(B);
    end
end