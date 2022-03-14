function zcPrice = PtT_CIR1(params,t,T,tind,x,y)
%%PTT_CIR1 computes zero-coupon price of CIR- model for given parameters
%    Input:
%       params (8 x 1 array): $[\phi_1^x,...,\phi_1^y,...\phi_3^y,x_t0,y_t0]$
%       t (p x 1): contains the current times
%       T (1 x n array or pxn): contains the maturities
%       tind (p x 1 int array): indices such that model-times matches t
%       x (N x M array): contains the simulated paths of the CIR process x
%       y (N x M array): contains the simulated paths of the CIR process y
%       modelTimes (N x 1 array): contains the timeline of the model
%    Output:
%       zcPrice (p x n x M array): contains the future discount factors DtT 
%                                  for given initial time t and maturities 
%                                  T. To get the future zero coupon price 
%                                  use mean(zcPrice,3)

% tind=zeros(size(t));
% for i=1:1:length(tind)
%     tind(i) = find(modelTimes<=t(i),1,'last');
% end

xt=x(tind,:);
xt=reshape(xt,size(xt,1),1,size(xt,2));
yt=y(tind,:);
yt=reshape(yt,size(yt,1),1,size(yt,2));
tau=T-t;

% parameters for $x_t$
phi1x = params(1);
phi2x = params(2);
phi3x = params(3);

% parameters for $y_t$
phi1y = params(4);
phi2y = params(5);
phi3y = params(6);

[Ax,Bx] = RiccatiCIR(phi1x,phi2x,phi3x,tau);
[Ay,By] = RiccatiCIR(phi1y,phi2y,phi3y,tau);
zcPrice = Ax.*exp(-Bx.*xt).*Ay.*exp(By.*yt);
end