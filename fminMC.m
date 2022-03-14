function f = fminMC(params,T,dW1,dW2,strikeSwap,maturitySwap,tenorSwap,...
                  P0TMarket,modelTimes,marketPrice,swapType,varargin)
%FMINGC provides the objective function for CIR- model calibration
%    Input:
%       marketPrice (p x q array): market swaption matrix
%       modelPrice (p x q array): model swaption matrix
%    Output:
%       f (double): contains the value of the least squares problem
%
compMode='subMatrix';
coords={};
for kk=1:2:length(varargin)
    switch varargin{kk}
        case 'compMode'
            compMode=varargin{kk+1};
        case 'coordinates'
            coords=varargin{kk+1};
    end
end

[x,y,dfCIR1] = sim_CIR1(params,T,dW1,dW2);
switch compMode
    case 'subMatrix'
        modelPrice=swaption_matrix(params,...
                        x,y,modelTimes,dfCIR1,...
                        strikeSwap,maturitySwap,...
                        tenorSwap,...
                        P0TMarket,...
                        swapType);
         mPrice=marketPrice;
    case 'coordinates'
        matInd=coords{1};
        tenInd=coords{2};
        modelPrice=zeros(order-3,length(matInd),1);
        mPrice=zeros(length(matInd),1);
        for k=1:1:length(matInd)
            i=matInd(k);
            j=tenInd(k);
            modelPrice(:,k)=swaption_matrix(params,...
                                x,y,modelTimes,dfCIR1,...
                                strikeSwap(i,j),maturitySwap(i),...
                                tenorSwap(j),...
                                P0TMarket,...
                                swapType);
            mPrice(k)=marketPrice(i,j);
        end
    otherwise
        error('Wrong computational mode, %s given but subMatrix or coordinates required',compMode)
end
                   
% f = sum(abs(marketPrice-mPrice).^2,'all');
f = sum((modelPrice./mPrice -1).^2,'all');
end
