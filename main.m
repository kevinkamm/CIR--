clear all;close all;fclose('all');
pool = gcp('nocreate');
if isempty(pool)
    % pool = parpool('threads'); % multithreading
    pool = parpool('local'); % multiprocessing
end
rng(0);
%% Data sets
% Choose the dataset and swaption type
dataIndex=2;
swaptionType='Receiver';
% swaptionType='Payer';
switch swaptionType
    case 'Receiver'
        swapType=-1;
    case 'Payer'
        swapType=1;
end
dataNames={'TermStructure_EUR_20191230','TermStructure_EUR_20201130'};
dataNamesSwap={'Swaption_EUR_20191231','Swaption_EUR_20201211'};
dataNamesBermudanSwap={['bermudanSwaption',swaptionType,'_EUR_20191230'],...
                       ['bermudanSwaption',swaptionType,'_EUR_20201130']};
dataNamesCMS={'cmsRate_EUR_20191230','cmsRate_EUR_20201130'};
%% Load the term structure
dataDir='TermStructure';
dataName=dataNames{dataIndex};
dataFile = [dataDir,'\',dataName,'.xls'];
sheet = 'Term_Str_Valori';
disp('Loading Termstructure data')

dataset = xlsread(dataFile,sheet);

% zero rates converted from percent to decimal
marketZeroRates = dataset(:,1)/100;
% market discount factor in decimals
marketDF = dataset(:,2);
% times of market discount factors and rates in years
marketTimes = dataset(:,3);

%% Load the swaption data
disp('Loading Swaption data')
dataDirSwap='Swaption';
dataNameSwap=dataNamesSwap{dataIndex};
dataFileSwap = [dataDirSwap,'\',dataNameSwap,'.xls'];
sheetSwap = 'Principale';

strikeSwap = xlsread(dataFileSwap,sheetSwap,'R14:V20');
maturitySwap = xlsread(dataFileSwap,sheetSwap,'Q14:Q20');
tenorSwap = xlsread(dataFileSwap,sheetSwap,'R13:V13');
marketSwapPrice = xlsread(dataFileSwap,sheetSwap,'J14:N20');
volSwap=xlsread(dataFileSwap,sheetSwap,'J28:P36');
volMatSwap=xlsread(dataFileSwap,sheetSwap,'I28:I36');
volTenorSwap=xlsread(dataFileSwap,sheetSwap,'J27:P27');
%% Load the bermudan swaption data
dataDirBermudanSwap='BermudanSwaption';
dataNameBermudanSwap=dataNamesBermudanSwap{dataIndex};
dataFileBermudanSwap = [dataDirBermudanSwap,'\',dataNameBermudanSwap,'.xlsx'];
strikeBermudanSwap = xlsread(dataFileBermudanSwap,'Atm Strike','B2:E6')./100;
maturityBermudanSwap = xlsread(dataFileBermudanSwap,'NPV','A2:A6');
tenorBermudanSwap = xlsread(dataFileBermudanSwap,'NPV','B1:E1');
hw1BermudanSwapPrice = xlsread(dataFileBermudanSwap,'NPV','B2:E6')./100;
%% Load CMS data
dataDirCMS='CMS';
dataNameCMS=dataNamesCMS{dataIndex};
dataFileCMS = [dataDirCMS,'\',dataNameCMS,'.xlsx'];
effectiveCMS = xlsread(dataFileCMS,'Sheet1','A2:A10');
tenorCMS = xlsread(dataFileCMS,'Sheet1','B2:B10');
indexCMS = xlsread(dataFileCMS,'Sheet1','C2:C10');
marketRateCMS = xlsread(dataFileCMS,'Sheet1','D2:D10')./100;
disp('Finished loading data')
%% Market Data interpolation
P0TMarket=P0T_Market(marketTimes,marketDF,'spline');
strikeSwap = atmStrikes(maturitySwap,tenorSwap,P0TMarket);
%% Select swaption prices for calibration by maturity and tenor
compMode='subMatrix';
% compMode='coordinates';
switch compMode
    case 'subMatrix'
        coords={};%dont change
        p=3:1:length(maturitySwap)-1;
        q=5;

%         p=3:1:length(maturitySwap)-1;
%         q=4;         

%         p=3:1:length(maturitySwap)-1;
%         q=3;

%         p=3:1:length(maturitySwap)-1;
%         q=2;
        
%         p=3:1:length(maturitySwap)-1;
%         q=1;

        maturitySwapTemp=maturitySwap(p);
        tenorSwapTemp=tenorSwap(q);
        strikeSwapTemp=strikeSwap(p,q);
        marketPrice=marketSwapPrice(p,q)
        
    case 'coordinates'
        maturitySwapTemp=maturitySwap;
        tenorSwapTemp=tenorSwap;
        strikeSwapTemp=strikeSwap;
        marketPrice=marketSwapPrice;
        
        % co terminal: large maturity with small tenor
        matInd=length(maturitySwap)-2:-1:2;
        tenInd=1:1:length(tenorSwap)-1;
        
%         matInd=3:1:length(maturitySwap);
%         tenInd=1:1:length(tenorSwap);
%         matInd=2:4;
%         tenInd=2:4;
%         matInd=1:4;
%         tenInd=1:4;
        coords={matInd,tenInd};
        matTenInd=sub2ind(size(marketSwapPrice),matInd,tenInd);
        marketSwapPrice(matTenInd)
    otherwise
        error('Wrong computational mode, %s given but subMatrix or coordinates required',compMode)
end

%% Simulate Brownian motions
% Choose the parameters for simulation
disp('Simulate Brownian motions for MC')
t0 = 0;
T = marketTimes(end);
% T=max(maturitySwapTemp,[],'all')+max(tenorSwapTemp,[],'all');
dt=1/256;
% dt=1/512;
N = ceil((T-t0)/dt)+ 1;
modelTimes = linspace(t0,T,N);
M = 1e4; 
% M=1e5;
[dW1,dW2] = BrownianIncrements(T,N,M);

%% Calibration with Gram Charlier
% Choose order: 3 up to 7 possible
order = 7;
% Choose which GC orders are taken into account in fmin
% calMode = order;
calMode = [3,5,7];

fileDates=extract(dataName,digitsPattern);
if strcmp(compMode,'coordinates')
    temp=num2str(cell2mat(coords'));
    fileName=sprintf('CIR%s_%s_M%s_T%s_GC%s_coords',...
                 fileDates{1},...
                 swaptionType,...
                 replace(temp(1,:),' ',''),...
                 replace(temp(2,:),' ',''),...
                 replace(num2str(calMode),' ',''));
else
    fileName=sprintf('CIR%s_%s_M%s_T%s_GC%s',...
                 fileDates{1},...
                 swaptionType,...
                 replace(num2str(maturitySwapTemp'),' ',''),...
                 replace(num2str(tenorSwapTemp),' ',''),...
                 replace(num2str(calMode),' ',''));
end
saveParam=true;
loadParam=true;

% compute subset sums
S = cell(length(tenorSwapTemp),order);
for iTen = 1:1:length(tenorSwapTemp)
    for iOrd=1:1:order
        S{iTen,iOrd} = subsetSum(tenorSwapTemp(iTen)+1,iOrd);
    end
end
fprintf('Calibration with Gram-Charlier orders: %s\n',num2str(calMode))
[paramsGC, ctimesGC, errorsGCFmin, methodsGC]=calibrationGC(strikeSwapTemp,...
                  maturitySwapTemp,tenorSwapTemp,...
                  P0TMarket,marketPrice,S,order,...
                  swapType,...
                  'compMode',compMode,...
                  'coordinates',coords,...
                  'calMode',calMode,...
                  'paramDir',fileName,...
                  'saveParam',saveParam,...
                  'loadParam',loadParam);
%% Compute results with GC parameters
% Evaluate Gram Charlier prices
fprintf('Compute Gram Charlier prices\n')
sovGC=cell(length(paramsGC),order-2);
for iP=1:1:length(paramsGC)
%     fprintf('Gram Charlier prices with params %d\n',iP)
    temp=GramCharlier_0T(paramsGC{iP},...
                     P0TMarket,...
                     maturitySwapTemp,...
                     tenorSwapTemp,...
                     strikeSwapTemp,...
                     S,...
                     order,...
                     swapType);
    for iO=3:1:order
        if strcmp(compMode,'coordinates')
            sovGC{iP,iO-2}=temp(iO-2,matTenInd);
        else
            sovGC{iP,iO-2}=reshape(temp(iO-2,:),...
                length(maturitySwapTemp),length(tenorSwapTemp));
        end
    end
end
sovGCFull=cell(length(paramsGC),order-2);
Sfull = cell(length(tenorSwap),order);
for iTen = 1:1:length(tenorSwap)
    for iOrd=1:1:order
        Sfull{iTen,iOrd} = subsetSum(tenorSwap(iTen)+1,iOrd);
    end
end

for iP=1:1:length(paramsGC)
%     fprintf('Gram Charlier prices with params %d\n',iP)
    temp=GramCharlier_0T(paramsGC{iP},...
                             P0TMarket,...
                             maturitySwap,...
                             tenorSwap,...
                             strikeSwap,...
                             Sfull,...
                             order,...
                             swapType);
    for iO=3:1:order
        sovGCFull{iP,iO-2}=reshape(temp(iO-2,:),...
                length(maturitySwap),length(tenorSwap));
    end
end
% Evaluate Monte Carlo prices
fprintf('Compute Monte Carlo prices\n')
sovMC=cell(length(paramsGC),1);
for iP=1:1:length(paramsGC)
%     fprintf('Monte Carlo prices with params %d\n',iP)
    [x,y,dfCIR1] = sim_CIR1(paramsGC{iP},T,dW1,dW2);
    temp=swaption_matrix(paramsGC{iP},x,y,modelTimes,dfCIR1,...
                       strikeSwapTemp,maturitySwapTemp,tenorSwapTemp,...
                       P0TMarket,swapType);
    if strcmp(compMode,'coordinates')
        sovMC{iP}=temp(matTenInd);
    else
        sovMC{iP}=temp;
    end
end
sovMCFull=cell(length(paramsGC),1);
for iP=1:1:length(paramsGC)
%     fprintf('Monte Carlo prices with params %d\n',iP)
    [x,y,dfCIR1] = sim_CIR1(paramsGC{iP},T,dW1,dW2);
    temp=swaption_matrix(paramsGC{iP},x,y,modelTimes,dfCIR1,...
                       strikeSwap,maturitySwap,tenorSwap,...
                       P0TMarket,swapType);
    sovMCFull{iP}=temp;
end
%% Errors with GC calibration
% we use an L1-error divided by number of entries
if strcmp(compMode,'coordinates')
    currMarketPrice=marketPrice(matTenInd);
else
    currMarketPrice=marketPrice;
end
% Errors GC vs MC
disp('Compute errors: GC vs MC')
errGCMC=cell(size(sovGC));
for iP=1:1:size(sovGC,1)
    for iO=1:1:size(sovGC,2)
        errGCMC{iP,iO}=mean(abs(sovGC{iP,iO}-sovMC{iP}),'all');
        fprintf('Mean-L1 Error GC%d vs MC with params %d: %3.3e\n',iO+2,iP,errGCMC{iP,iO})
    end
end
% Errors GC vs Market
disp('Compute errors: GC vs Market')
errGCM=cell(size(sovGC));
for iP=1:1:size(sovGC,1)
    for iO=1:1:size(sovGC,2)
        errGCM{iP,iO}=mean(abs(sovGC{iP,iO}-currMarketPrice),'all');
        fprintf('Mean-L1 Error GC%d vs Market with params %d: %3.3e\n',iO+2,iP,errGCM{iP,iO})
    end
end
% Errors MC vs Market
disp('Compute errors: MC vs Market')
errMCM=cell(length(paramsGC),1);
for iP=1:1:size(sovGC,1)
        errMCM{iP,1}=mean(abs(sovMC{iP,1}-currMarketPrice),'all');
        fprintf('Mean-L1 Error MC vs Market with params %d: %3.3e\n',iP,errMCM{iP,1})
end
%% Calibration with Monte Carlo
paramsMC={};
ctimesMC={};
errorsMCFmin={};
methodsMC={};
% fprintf('Calibration with Monte Carlo\n')
% [paramsMC, ctimesMC, errorsMCFmin, methodsMC]=calibrationMC(T,dW1,dW2,...
%                                 strikeSwapTemp,maturitySwapTemp,tenorSwapTemp,...
%                                 P0TMarket,modelTimes,marketPrice,...
%                                 swapType,...
%                                 'compMode',compMode,...
%                                 'coordinates',coords);

%% Bermudan swaption
disp('Compute Bermudan Swaption Prices')
priceBSwaption=zeros(length(paramsGC),length(maturityBermudanSwap),length(tenorBermudanSwap));
for iP = 1:1:length(paramsGC)
    for iM = 1:1:size(priceBSwaption,2)
        for iT = 1:1:size(priceBSwaption,3)
            priceBSwaption(iP,iM,iT)=...
                bermudanSwaption(paramsGC{iP},x,y,modelTimes,dfCIR1,...
                                   strikeBermudanSwap(iM,iT),...
                                   maturityBermudanSwap(iM),...
                                   maturityBermudanSwap(iM)+tenorBermudanSwap(iT),...
                                   P0TMarket,...
                                   swapType);
        end
    end
end
errBSwaption=abs(priceBSwaption-reshape(hw1BermudanSwapPrice,[1,size(hw1BermudanSwapPrice)]));
%%
% hw1BermudanSwapPrice
% squeeze(priceBSwaption(3,:,:))
squeeze(errBSwaption(3,:,:))
mean(errBSwaption(3,:,:),'all')
squeeze(errBSwaption(4,:,:))
mean(errBSwaption(4,:,:),'all')
%% CMS rates
disp('Compute CMS par rates')
CMSrate=zeros(length(paramsGC),length(marketRateCMS));
CMSprice=zeros(size(CMSrate));
for iP=1:1:length(paramsGC)
    for i=1:1:length(CMSrate)
    [CMSprice(iP,i),CMSrate(iP,i)]=constantMaturitySwap(paramsGC{iP},x,y,modelTimes,dfCIR1,...
                                   [],effectiveCMS(i),effectiveCMS(i)+tenorCMS(i),indexCMS(i),...
                                   P0TMarket,...
                                   swapType);
    end
end
errCMSrate=abs(CMSrate-marketRateCMS');
errCMSrate(3,:)
errCMSrate(4,:)
%% Create output
disp('Saving figures and creating output')
output(fileName,...
        paramsGC, ctimesGC, errorsGCFmin, methodsGC,calMode,...
        sovGC,sovMC,...
        errGCMC,errGCM,errMCM,...
        marketTimes,marketDF,marketZeroRates,...
        maturitySwapTemp,tenorSwapTemp,marketPrice,strikeSwapTemp,...
        maturitySwap,tenorSwap,marketSwapPrice,strikeSwap,...
        maturityBermudanSwap,tenorBermudanSwap,strikeBermudanSwap,...
        hw1BermudanSwapPrice,priceBSwaption,errBSwaption,...
        effectiveCMS,tenorCMS,indexCMS,marketRateCMS,CMSrate,errCMSrate,...
        volSwap,volMatSwap,volTenorSwap)
disp('Done')