function SOV=GramCharlier_0T(params,...
                             P0TMarket,...
                             maturitySwap,...
                             tenorSwap,...
                             strikeSwap,...
                             S,...
                             order,...
                             swapType,...
                             varargin)
%%GRAMCHARLIER_0T calculates swaption prices for given parameters in CIR--
% model.
%    Input:
%       params (8x1 array): params= $[\phi_1^x,...,\phi_1^y,...\phi_3^y,x_t0,y_t0]$
%       strikeSwap (p x q array): strikes of swaption 
%       maturitySwap (p x 1 array): maturities of swaption
%       tenorSwap (1 x q array): tenors of swaption
%       P0TMarket (griddedInterpolant): interpolated market zero coupon 
%                                       prices
%       marketPrice (p x q array): market prices of swaption
%       S (s1 x s2 cell array): contains subset sums for Gram-Charlier
%       order (1 x 1 int): contains the order of the Gram-Charlier approx
%       varargin (cell array): contains optional name-value pairs
%       swapType (int): equal to 1 for payer swaption and -1 for receiver
%    Output:
%       SOV (order-2 x p x q array): Swaption value for all orders up to
%                                    order
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

%% Compute swap moments
sm = swapMoments(params,maturitySwap,T,tenorSwap,P0TMarket,S,strikeSwap,order,swapType);
%% Compute cumulants
Ptemp = P0TMarket(maturitySwap)'.^([1:order]');
Cn=cumulants(sm).*Ptemp;
%% Compute coefficients
qn=expansionCoefficients(Cn);
%% Compute expansion
if ~(all(isreal(sqrt(Cn(2,:))),'all')&& all(isreal(sm),'all'))
%     warning('Complex Riccati')
    SOV=10^6.*ones([order-2,size(strikeSwap)]);
else
% C1=reshape(Cn(1,:,:),[size(Cn,2),size(Cn,3)]);
% C2=reshape(Cn(2,:,:),[size(Cn,2),size(Cn,3)]);
% tempC=C1./sqrt(C2)
C1=Cn(1,:,:);
C2=Cn(2,:,:);
tempC=C1./sqrt(C2);
he=hermite(reshape(tempC,length(maturitySwap),length(tenorSwap)),1:1:order-2);
SOV=C1.*normcdf(tempC)+...
    sqrt(C2).*normpdf(tempC).*...
    (1+cumsum((-1).^[3:order]'.*qn(3:order,:,:).*he,1));
% SOV=reshape(SOV,order-2,length(maturitySwap),length(tenorSwap));
end
end
function sm=swapMoments(params,T0,T,tenorSwap,P0TMarket,S,strikes,order,swapType)
% parameters for $x_t$
phi1x = params(1);
phi2x = params(2);
phi3x = params(3);
x0   = params(7);

% parameters for $y_t$
phi1y = params(4);
phi2y = params(5);
phi3y = params(6);
y0   = params(8);

tau=cat(2,zeros(length(T0),1),T-T0);%tenor+mat - mat

[Ax,Bx] = RiccatiCIR(phi1x,phi2x,phi3x,tau);
[Ay,By] = RiccatiCIR(phi1y,phi2y,phi3y,tau);
P0TMCIR1temp = (P0TMarket(T)./P0T_CIR1(params,T));
P0T0Mtemp = P0TMarket(T0)./P0T_CIR1(params,T0);
P0T0CIRM1temp = (P0T_CIR1(params,T0)./P0TMarket(T0));
sm=zeros([order,size(strikes)]);
for iO=1:1:order
    for iTen=1:1:size(strikes,2)
        Scurr=S{iTen,iO}; %difficult to vectorize because different sizes
%         switch swaptionType
%             case 'receiver'
%                 % receiver swaption
%                 a=[-P0T0Mtemp,strikes(:,iTen).*ones(1,tenorSwap(iTen))];
%                 a(:,end)=a(:,end)+1;
%                 a(:,2:end)=a(:,2:end).*P0TMCIR1temp(:,1:tenorSwap(iTen));
%             case 'payer'
%                 % payer swaption
%                 a=[P0T0Mtemp,-strikes(:,iTen).*ones(1,tenorSwap(iTen))];
%                 a(:,end)=a(:,end)-1;
%                 a(:,2:end)=a(:,2:end).*P0TMCIR1temp(:,1:tenorSwap(iTen));
%             otherwise 
%                     error('Unknown swaption type: %s given, needs to be receiver or payer',swaptionType)
%         end
        a=[swapType.*P0T0Mtemp,-swapType.*strikes(:,iTen).*ones(1,tenorSwap(iTen))];
        a(:,end)=a(:,end)-swapType.*1;
        a(:,2:end)=a(:,2:end).*P0TMCIR1temp(:,1:tenorSwap(iTen));

        k=Scurr';
        kT=reshape(k,1,size(k,1),size(k,2)); %mat,ak exponent,permutations
        count=(factorial(iO)./prod(factorial(kT),2));
        ak=count.*prod(a.^kT,2);
        ax=prod(Ax(:,1:tenorSwap(iTen)+1).^kT(1,:,:),2);
        bx=sum(Bx(:,1:tenorSwap(iTen)+1).*kT(1,:,:),2);
        ay=prod(Ay(:,1:tenorSwap(iTen)+1).^kT(1,:,:),2);
        by=sum(By(:,1:tenorSwap(iTen)+1).*kT(1,:,:),2);
        tau=T0;
        % bond moments
        [Mx,Nx] = RiccatiCIRTerminal(phi1x,phi2x,phi3x,ax,bx,tau);
        [My,Ny] = RiccatiCIRTerminal(phi1y,phi2y,phi3y,ay,by,tau);
        sm(iO,:,iTen) = P0T0CIRM1temp.^iO.*...
                           sum(ak.*Mx.*exp(-Nx.*x0).*My.*exp(Ny.*y0),3)...
                           ./P0T_CIR1(params,T0);
    end
end
end
function c=cumulants(moments)
%%Cumulants computes the cumulants for given moments
%   Input:
%       moments (nd array): sucessive moments starting at first moment
%   Output:
%       c (size(moments) array): cumulants
    c=zeros(size(moments));
    for i=1:1:size(moments,1)
        switch i
            case 1
                c(i,:,:)=moments(1,:,:);
            case 2
                c(i,:,:)=moments(2,:,:)-moments(1,:,:).^2;
            case 3
                c(i,:,:)=moments(3,:,:)-3.*moments(1,:,:).*moments(2,:,:)+2.*moments(1,:,:).^3;
            case 4
                c(i,:,:)=moments(4,:,:)-4.*moments(1,:,:).*moments(3,:,:)-3.*moments(2,:,:).^2+...
                    12.*moments(1,:,:).^2.*moments(2,:,:)-6.*moments(1,:,:).^4;
            case 5
                c(i,:,:)=moments(5,:,:)-5.*moments(1,:,:).*moments(4,:,:)-...
                    10.*moments(2,:,:).*moments(3,:,:)+20.*moments(1,:,:).^2.*moments(3,:,:)+...
                    30.*moments(1,:,:).*moments(2,:,:).^2-60.*moments(1,:,:).^3.*moments(2,:,:)+...
                    24.*moments(1,:,:).^5;
            case 6 
                c(i,:,:)=moments(6,:,:)-6.*moments(1,:,:).*moments(5,:,:)-...
                        15.*moments(2,:,:).*moments(4,:,:)+...
                        30.*moments(1,:,:).^2.*moments(4,:,:)-...
                        10.*moments(3,:,:).^2+...
                        120.*moments(1,:,:).*moments(2,:,:).*moments(3,:,:)-...
                        120.*moments(1,:,:).^3.*moments(3,:,:)+...
                        30.*moments(2,:,:).^3-...
                        270.*moments(1,:,:).^2.*moments(2,:,:).^2+...
                        360.*moments(1,:,:).^4.*moments(2,:,:)-...
                        120.*moments(1,:,:).^6;
            case 7
                c(i,:,:)=moments(7,:,:)-7.*moments(1,:,:).*moments(6,:,:)-...
                        21.*moments(2,:,:).*moments(5,:,:)-...
                        35.*moments(3,:,:).*moments(4,:,:)+...
                        140.*moments(1,:,:).*moments(3,:,:).^2-...
                        630.*moments(1,:,:).*moments(2,:,:).^3+...
                        210.*moments(1,:,:).*moments(2,:,:).*moments(4,:,:)-...
                        1260.*moments(1,:,:).^2.*moments(2,:,:).*moments(3,:,:)+...
                        42.*moments(1,:,:).^2.*moments(5,:,:)+...
                        2520.*moments(1,:,:).^3.*moments(2,:,:).^2-...
                        210.*moments(1,:,:).^3.*moments(4,:,:)+...
                        210.*moments(2,:,:).^2.*moments(3,:,:)+...
                        840.*moments(1,:,:).^4.*moments(3,:,:)-...
                        2520.*moments(1,:,:).^5.*moments(2,:,:)+...
                        720.*moments(1,:,:).^7;
            otherwise
                error('Cumulant undefined')
        end
    end
end
function qn=expansionCoefficients(Cn)
%%ExpansionCoefficients computes the q's for the Gram-Charlier expansion
%   Input:
%       Cn (nd array): cumulants with sucessive order in first dimension
%   Output:
%       qn (size(Cn) array): expansion coefficients
%   Usage:
%       qn=expansionCoefficients(Cn)
qn=zeros(size(Cn));
for i=1:1:size(qn,1)
    switch i
        case 0
            qn(1,:,:)=1;
        case {1,2}
            
        case 3
            qn(3,:,:)=Cn(3,:,:)./(factorial(3).* Cn(2,:,:).^(3/2));
        case 4
            qn(4,:,:)=Cn(4,:,:)./(factorial(4).* Cn(2,:,:).^(2));
        case 5
            qn(5,:,:)=Cn(5,:,:)./(factorial(5).* Cn(2,:,:).^(5/2));
        case 6
            qn(6,:,:)=(Cn(6,:,:)+10.*Cn(3,:,:).^2)./(factorial(6).* Cn(2,:,:).^(3));
        case 7
            qn(7,:,:)=(Cn(7,:,:)+35.*Cn(3,:,:).*Cn(4,:,:))./(factorial(7).* Cn(2,:,:).^(7/2));
        otherwise
            error('Unknown order')
    end
end
end
function he=hermite(x,orders)
%%Hermite computes the probabilists hermite polynomials at x for given
% order.
%   Input:
%       x (nd array): evaluation points
%       orders (dx1 or 1xd array): orders to of hermite polynomials
%   Output:
%       he (d x size(x) array): evaluated hermite polynomials
xPowers=repmat(reshape(x,[1 size(x)]),[max(orders),ones(1,ndims(x))]);
xPowers=cumprod(xPowers,1);
function he0=hermite0()
    he0=ones(size(x));
end
function he1=hermite1()
    he1=xPowers(1,:,:);
end
function he2=hermite2()
    he2=xPowers(2,:,:)-1;
end
function he3=hermite3()
    he3=xPowers(3,:,:)-3.*xPowers(1,:,:);
end
function he4=hermite4()
    he4=xPowers(4,:,:)-6.*xPowers(2,:,:)+3;
end
function he5=hermite5()
    he5=xPowers(5,:,:)-10.*xPowers(3,:,:)+15.*xPowers(1,:,:);
end
function he6=hermite6()
    he6=xPowers(6,:,:)-15.*xPowers(4,:,:)+45.*xPowers(2,:,:)-15;
end
function he7=hermite7()
    he7=xPowers(7,:,:)-21.*xPowers(5,:,:)+105.*xPowers(3,:,:)-105.*xPowers(1,:,:);
end
he=zeros([length(orders),size(x)]);
for i=1:1:length(orders)
    switch orders(i)
        case 0
            he(i,:,:)=hermite0();
        case 1
            he(i,:,:)=hermite1();
        case 2
            he(i,:,:)=hermite2();
        case 3
            he(i,:,:)=hermite3();
        case 4
            he(i,:,:)=hermite4();
        case 5
            he(i,:,:)=hermite5();
        case 6
            he(i,:,:)=hermite6();
        case 7
            he(i,:,:)=hermite7();
        otherwise
            error('Unknown Hermite Polynomial')
    end
end
end