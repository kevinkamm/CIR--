function [params, ctimes, errors, methods]=calibrationGC(strikeSwap,...
                                                maturitySwap,...
                                                tenorSwap,...
                                                P0TMarket,...
                                                marketPrice,...
                                                S,...
                                                order,...
                                                swapType,...
                                                varargin)
%%CALIBRATIONGC calibrates the CIR-- model to the selected swaption prices
%    Input:
%       strikeSwap (p x q array): strikes of swaption 
%       maturitySwap (p x 1 array): maturities of swaption
%       tenorSwap (1 x q array): tenors of swaption
%       P0TMarket (griddedInterpolant): interpolated market zero coupon 
%                                       prices
%       marketPrice (p x q array): market prices of swaption
%       S (s1 x s2 cell array): contains subset sums for Gram-Charlier
%       order (1 x 1 int): contains the order of the Gram-Charlier approx
%       swapType (int): equal to 1 for payer swaption and -1 for receiver
%       varargin (cell array): contains optional name-value pairs
%           'swaptionType': 'payer' (default), 'receiver'
%           'compMode': 'subMatrix' (default), 'coordinates'
%           'coordinates': cell array
%    Output:
%       params (d x 1 cell array): calibrated parameters
%       ctimes (d x 1 cell array): computational times
%       errors (d x 1 cell array): calibration error, value of fminGC
%       methods (d x 1 cell array): calibration methods or which initial
%                                   datum
%
% See also fmincon, ga, GramCharlier_0T, linconstr, subsetSum, fminGC

compMode='subMatrix';
calMode=order;
coords={};
saveParam=false;
loadParam=false;
paramDir='';
saveDir='Results/Matlab/';
for kk=1:2:length(varargin)
    switch varargin{kk}
        case 'compMode'
            compMode=varargin{kk+1};
        case 'coordinates'
            coords=varargin{kk+1};
        case 'calMode'
            calMode=varargin{kk+1};
        case 'saveParam'
            saveParam=varargin{kk+1};
        case 'loadParam'
            loadParam=varargin{kk+1};
        case 'paramDir'
            paramDir=varargin{kk+1};
    end
end
if saveParam
    mkDir([saveDir,paramDir]);
end
% linear constraint matrix
[A,b] = linconstr();

% optimization options
tolCon=1e-12;
optsGa = optimoptions('ga',...
                    'ConstraintTolerance',tolCon,...
                    'UseParallel',true);
opts = optimoptions('fmincon',...
                    'MaxFunctionEvaluations',10000,...
                    'MaxIterations',1000,...
                    'TolCon',tolCon,...
                    'UseParallel',true);

% params= $[\phi_1^x,...,\phi_1^y,...\phi_3^y,x_0,y_0]$
% lower and upper bounds
lb=1e-4*ones(1,8);
% Feller conditions
 lb(3)=1;lb(6)=1;
% for convenience
ub=1*ones(1,8);
ub(3)=2;ub(6)=2; 
% ub(3)=4;ub(6)=4; 
ub(7)= 1e-1;
ub(8) = 1e-1;

% Initialize outputs
params={};
ctimes={};
errors={};
methods={};
% Global optimization
disp('Start global optimization')
if loadParam && exist([saveDir,paramDir,'/','ga.mat'])
    tempGa=load([saveDir,paramDir,'/','ga.mat']);
    paramGa=tempGa.paramGa;
    errGa=tempGa.errGa;
    ctimeGa=tempGa.ctimeGa;
else
    ticGa=tic;
    [paramGa,errGa] = ga(@(params)fminGC(params,strikeSwap,maturitySwap,tenorSwap,...
                      P0TMarket,marketPrice,S,...
                      order,swapType,...
                      'compMode',compMode,...
                      'coordinates',coords,...
                      'calMode',calMode),...
                 8,...
                 A,b,...
                 [],[],...
                 lb,ub,...
                 [],optsGa);
    ctimeGa=toc(ticGa);
end
params{end+1}=paramGa;
ctimes{end+1}=ctimeGa;
errors{end+1}=errGa;
methods{end+1}='ga';
fprintf('Finished global optimization with err=%1.3e after %3.3f seconds\n',errGa,ctimeGa)
if saveParam
    save([saveDir,paramDir,'/','ga.mat'],'paramGa','errGa','ctimeGa');
end

% Initial datum ga
disp('Start local optimization with ga input')
if loadParam && exist([saveDir,paramDir,'/','gaFmin.mat'])
    tempGaFmin=load([saveDir,paramDir,'/','gaFmin.mat']);
    paramGaFmin=tempGaFmin.paramGaFmin;
    errGaFmin=tempGaFmin.errGaFmin;
    ctimeGaFmin=tempGaFmin.ctimeGaFmin;
else
ticGaFmin=tic;
[paramGaFmin,errGaFmin] = fmincon(@(params)fminGC(params,strikeSwap,maturitySwap,tenorSwap,...
                  P0TMarket,marketPrice,S,...
                  order,swapType,...
                  'compMode',compMode,...
                  'coordinates',coords,...
                  'calMode',calMode),...
             paramGa,...
             A,b,...
             [],[],...
             lb,ub,...
             [],...
             opts);
ctimeGaFmin=toc(ticGaFmin)+ctimeGa;
end
params{end+1}=paramGaFmin;
ctimes{end+1}=ctimeGaFmin;
errors{end+1}=errGaFmin;
methods{end+1}='ga + fmincon';
fprintf('Finished local optimization with err=%1.3e after %3.3f seconds\n',errGaFmin,ctimeGaFmin)
if saveParam
    save([saveDir,paramDir,'/','gaFmin.mat'],'paramGaFmin','errGaFmin','ctimeGaFmin');
end
% Initial datum (ub+lb)/2
% disp('Start local optimization with (ub+lb)/2')
% ticFmin1=tic;
% [paramFmin1,errFmin1] = fmincon(@(params)fminGC(params,strikeSwap,maturitySwap,tenorSwap,...
%                   P0TMarket,marketPrice,S,...
%                   order,swapType,...
%                   'compMode',compMode,...
%                   'coordinates',coords,...
%                   'calMode',calMode),...
%              (lb+ub)./2,...
%              A,b,...
%              [],[],...
%              lb,ub,...
%              [],...
%              opts);
% ctimeFmin1=toc(ticFmin1);
% fprintf('Finished local optimization with err=%1.3e after %3.3f seconds\n',errFmin1,ctimeFmin1)

% Initial datum (ub+lb)/10
disp('Start local optimization with admissible (ub+lb)./10')
if loadParam && exist([saveDir,paramDir,'/','fmin2.mat'])
    tempFmin2=load([saveDir,paramDir,'/','fmin2.mat']);
    paramFmin2=tempFmin2.paramFmin2;
    errFmin2=tempFmin2.errFmin2;
    ctimeFmin2=tempFmin2.ctimeFmin2;
else
    ticFmin2=tic;
    tempParam=(lb+ub)./10;
    tempParam(2)=tempParam(2).*0.95;
    tempParam(4)=tempParam(4).*0.95;
    [paramFmin2,errFmin2] = fmincon(@(params)fminGC(params,strikeSwap,maturitySwap,tenorSwap,...
                      P0TMarket,marketPrice,S,...
                      order,swapType,...
                      'compMode',compMode,...
                      'coordinates',coords,...
                      'calMode',calMode),...
                 tempParam,...
                 A,b,...
                 [],[],...
                 lb,ub,...
                 [],...
                 opts);
    ctimeFmin2=toc(ticFmin2);
end
params{end+1}=paramFmin2;
ctimes{end+1}=ctimeFmin2;
errors{end+1}=errFmin2;
methods{end+1}='fmincon 2';
fprintf('Finished local optimization with err=%1.3e after %3.3f seconds\n',errFmin2,ctimeFmin2)
if saveParam
    save([saveDir,paramDir,'/','fmin2.mat'],'paramFmin2','errFmin2','ctimeFmin2');
end
% Initial datum (ub+lb)/5
disp('Start local optimization with admissible (ub+lb)./20')
if loadParam && exist([saveDir,paramDir,'/','fmin3.mat'])
    tempFmin3=load([saveDir,paramDir,'/','fmin3.mat']);
    paramFmin3=tempFmin3.paramFmin3;
    errFmin3=tempFmin3.errFmin3;
    ctimeFmin3=tempFmin3.ctimeFmin3;
else
    ticFmin3=tic;
    tempParam=(lb+ub)./20;
    tempParam(2)=tempParam(2).*0.95;
    tempParam(4)=tempParam(4).*0.95;
    [paramFmin3,errFmin3] = fmincon(@(params)fminGC(params,strikeSwap,maturitySwap,tenorSwap,...
                      P0TMarket,marketPrice,S,...
                      order,swapType,...
                      'compMode',compMode,...
                      'coordinates',coords,...
                      'calMode',calMode),...
                 tempParam,...
                 A,b,...
                 [],[],...
                 lb,ub,...
                 [],...
                 opts);
    ctimeFmin3=toc(ticFmin3);
end
params{end+1}=paramFmin3;
ctimes{end+1}=ctimeFmin3;
errors{end+1}=errFmin3;
methods{end+1}='fmincon 3';
fprintf('Finished local optimization with err=%1.3e after %3.3f seconds\n',errFmin3,ctimeFmin3)
if saveParam
    save([saveDir,paramDir,'/','fmin3.mat'],'paramFmin3','errFmin3','ctimeFmin3');
end

% Outputs
% params={paramGa,paramGaFmin,paramFmin2};
% ctimes={ctimeGa,ctimeGa+ctimeGaFmin,ctimeFmin2};
% errors={errGa,errGaFmin,errFmin2};
% methods={'ga','ga + fmincon','fmincon 2'};
end
function mkDir(dir)
    if exist(dir)==0
        mkdir(dir);
    end
end
function delDir(dir)
    if exist(dir)==7
        rmdir(dir,'s');
    end
end