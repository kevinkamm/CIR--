function f = fminGC(params,strikeSwap,maturitySwap,tenorSwap,...
                  P0TMarket,marketPrice,S,...
                  order,swapType,varargin)
%FMINGC provides the objective function for CIR- model calibration
%    Input:
%       marketPrice (p x q array): market swaption matrix
%       modelPrice (p x q array): model swaption matrix
%       swapType (int): equal to 1 for payer swaption and -1 for receiver
%    Output:
%       f (double): contains the value of the least squares problem
%
compMode='subMatrix';
calMode=order;
coords={};
for kk=1:2:length(varargin)
    switch varargin{kk}
        case 'compMode'
            compMode=varargin{kk+1};
        case 'coordinates'
            coords=varargin{kk+1};
        case 'calMode'
            calMode=varargin{kk+1};
    end
end
calMode=calMode-2; % because first two orders don't exist
switch compMode
    case 'subMatrix'
        modelPrice=GramCharlier_0T(params,...
                     P0TMarket,...
                     maturitySwap,...
                     tenorSwap,...
                     strikeSwap,...
                     S,...
                     order,...
                     swapType);
         mPrice=marketPrice;
    case 'coordinates'
        matInd=coords{1};
        tenInd=coords{2};
        modelPrice=zeros(order-2,length(matInd),1);
        mPrice=zeros(length(matInd),1);
        for k=1:1:length(matInd)
            i=matInd(k);
            j=tenInd(k);
            modelPrice(:,k)=GramCharlier_0T(params,...
                     P0TMarket,...
                     maturitySwap(i),...
                     tenorSwap(j),...
                     strikeSwap(i,j),...
                     S(j,:),...
                     order,...
                     swapType);
            mPrice(k)=marketPrice(i,j);
        end
    otherwise
        error('Wrong computational mode, %s given but subMatrix or coordinates required',compMode)
end
mPrice=reshape(mPrice,[1,size(mPrice)]);
% f = sum(abs(mPrice-modelPrice).^2,'all');
f = sum((modelPrice(calMode,:,:)./mPrice -1).^2,'all');
end
