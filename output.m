function output(fileName,...
                paramsGC, ctimesGC, errorsGCFmin,methodsGC,calMode,...
                sovGC,sovMC,...
                errGCMC,errGCM,errMCM,...
                marketTimes,marketDF,marketZeroRates,...
                maturitySwapTemp,tenorSwapTemp,marketPrice,strikeSwapTemp,...
                maturitySwap,tenorSwap,marketSwapPrice,strikeSwap,...
                maturityBermudanSwap,tenorBermudanSwap,strikeBermudanSwap,...
                hw1BermudanSwapPrice,priceBSwaption,errBSwaption,...
                effectiveCMS,tenorCMS,indexCMS,marketRateCMS,CMSrate,errCMSrate,...
                volSwap,volMatSwap,volTenorSwap)
fclose('all');

numDict={'One','Two','Three','Four','Five','Six','Seven','Eight','Nine'};

temp=split(fileName,'_');
yearStr=extract(temp{1},digitsPattern);
yearStr=yearStr{1};
year=str2num(yearStr(1:4));
month=str2num(yearStr(5:6));
day=str2num(yearStr(7:end));
switch yearStr
    case '20191230'
        identifier='A';
    case '20201130'
        identifier='B';
    otherwise
        identifier='Z';
end
            
picType='eps';
saveParam='epsc';

root=[pwd, '\' ,'Results'];
pdfRoot=[root,'\','Pdf'];
tempPath=[pdfRoot,'\','temp'];
copyPath=[pdfRoot,'\',fileName];
templatePath=[tempPath,'\','template', '.','tex'];
outputFilePath=[copyPath,'\','template','.','pdf';...
            copyPath,'\','template','.','tex'];
copyFilePath=[copyPath,'\',fileName,'.','pdf';...
              copyPath,'\',fileName,'.','tex'];
inputPath=[tempPath,'\','input','.','tex'];
          
mkDir(pdfRoot);
% delDir(tempPath);
mkDir(tempPath);
cleanDir(tempPath,{'template.tex'});
delDir(copyPath)
mkDir(copyPath)
delFile(inputPath);

inputFile=fopen(inputPath,'a');
% Head
fprintf(inputFile,...
        '\\section{CIR-{}- model}\n');
%% Calibration Parameters
fprintf(inputFile,...
        '\\subsection{Calibration Parameters}\n');
fprintf(inputFile,...
    '\\paragraph*{Calibration Parameters for GC using orders %s}\\hfill\\\\\n',num2str(calMode));
[latexFilePath,latexCommand]=saveCalParamGC('CIRParam');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:1:size(latexCommand,1)
    fprintf(inputFile,...
            '\t%s\\hfill\\\\\n',latexCommand(iLC,:));
end
% % Euler Parameter
% fprintf(inputFile,...
%         '\\subsection{Euler Parameter}\n');
% [latexFilePath,latexCommand]=saveEulerParam('Euler');
% fprintf(inputFile,...
%         '\t\\input{%s}\n',changeSlash(latexFilePath));
% fprintf(inputFile,...
%         '\t%s\n',latexCommand);
%% Computational times
fprintf(inputFile,...
        '\\subsection{Computational Times}\n');
fprintf(inputFile,...
        '\\paragraph*{Calibration Times for GC}\\hfill\\\\\n\\noindent\n');
[latexFilePath,latexCommand]=saveCalTimesGC('CompTimes');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
fprintf(inputFile,...
        '\t%s\n\\hfill\\\\\n',latexCommand);
%% Errors
fprintf(inputFile,...
        '\\subsection{Errors}\n');
fprintf(inputFile,...
        '\\paragraph*{Errors with parameters from GC calibration}\\hfill\\\\\n\\noindent\n');
[latexFilePath,latexCommand]=saveErrGCParam('Errors');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
fprintf(inputFile,...
        '\t%s\n\\hfill\\\\\n',latexCommand);
%% Swaption Results
fprintf(inputFile,...
        '\\subsection{Swaption surface}\n');
fprintf(inputFile,...
        '\\subsubsection{Only GC-calibrated values}\n');
% Market
fprintf(inputFile,...
        '\\paragraph*{Market}\\hfill\\\\\n\\noindent\n');
[latexFilePath,latexCommand]=saveSwaptionMarket('Swaptions');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:2:size(latexCommand,2)
    fprintf(inputFile,...
        '\\begin{table}\n');
    fprintf(inputFile,...
        '\\centering\n');
    fprintf(inputFile,...
            '\t%s\\hfill\\\\\n',latexCommand{iLC});
    fprintf(inputFile,...
        '\\caption{%s}\n',latexCommand{iLC+1});
    fprintf(inputFile,...
        '\\end{table}\n');
end
fprintf(inputFile,...
        '\\FloatBarrier\n');
% GC
fprintf(inputFile,...
        '\\paragraph*{Gram-Charlier}\\hfill\\\\\n\\noindent\n');
[latexFilePath,latexCommand]=saveSwaptionGC('Swaptions');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:2:size(latexCommand,2)
    fprintf(inputFile,...
        '\\begin{table}');
    fprintf(inputFile,...
        '\\centering');
    fprintf(inputFile,...
            '\t%s\\hfill\\\\\n',latexCommand{iLC});
    fprintf(inputFile,...
        '\\caption{%s}',latexCommand{iLC+1});
    fprintf(inputFile,...
        '\\end{table}');
end
fprintf(inputFile,...
        '\\FloatBarrier\n');
% MC
fprintf(inputFile,...
        '\\paragraph*{Monte-Carlo}\\hfill\\\\\n\\noindent\n');
[latexFilePath,latexCommand]=saveSwaptionMC('Swaptions');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:2:size(latexCommand,2)
    fprintf(inputFile,...
        '\\begin{table}\n');
    fprintf(inputFile,...
        '\\centering\n');
    fprintf(inputFile,...
            '\t%s\\hfill\\\\\n',latexCommand{iLC});
    fprintf(inputFile,...
        '\\caption{%s}\n',latexCommand{iLC+1});
    fprintf(inputFile,...
        '\\end{table}\n');
end
fprintf(inputFile,...
        '\\FloatBarrier\n');
%% Bermudan Swaption Results
fprintf(inputFile,...
        '\\section{Bermudan swaption prices}\n');
% CIR-- prices
fprintf(inputFile,...
        '\\paragraph*{CIR-{}-}\\hfill\\\\\n\\noindent\n');
[latexFilePath,latexCommand]=saveBermudanSwaptionCIR('BSwaptions');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:2:size(latexCommand,2)
    fprintf(inputFile,...
        '\\begin{table}\n');
    fprintf(inputFile,...
        '\\centering\n');
    fprintf(inputFile,...
            '\t%s\\hfill\\\\\n',latexCommand{iLC});
    fprintf(inputFile,...
        '\\caption{%s}\n',latexCommand{iLC+1});
    fprintf(inputFile,...
        '\\end{table}\n');
end
fprintf(inputFile,...
        '\\FloatBarrier\n');
% HW1
fprintf(inputFile,...
        '\\paragraph*{HW1}\\hfill\\\\\n\\noindent\n');
[latexFilePath,latexCommand]=saveBermudanSwaptionHW1('BSwaptions');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:2:size(latexCommand,2)
    fprintf(inputFile,...
        '\\begin{table}\n');
    fprintf(inputFile,...
        '\\centering\n');
    fprintf(inputFile,...
            '\t%s\\hfill\\\\\n',latexCommand{iLC});
    fprintf(inputFile,...
        '\\caption{%s}\n',latexCommand{iLC+1});
    fprintf(inputFile,...
        '\\end{table}\n');
end
fprintf(inputFile,...
        '\\FloatBarrier\n');
% Errors
fprintf(inputFile,...
        '\\paragraph*{Errors CIR-{}- vs HW1}\\hfill\\\\\n\\noindent\n');
[latexFilePath,latexCommand]=saveBermudanSwaptionErrors('BSwaptions');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:2:size(latexCommand,2)
    fprintf(inputFile,...
        '\\begin{table}\n');
    fprintf(inputFile,...
        '\\centering\n');
    fprintf(inputFile,...
            '\t%s\\hfill\\\\\n',latexCommand{iLC});
    fprintf(inputFile,...
        '\\caption{%s}\n',latexCommand{iLC+1});
    fprintf(inputFile,...
        '\\end{table}\n');
end
fprintf(inputFile,...
        '\\FloatBarrier\n');
% Strikes
fprintf(inputFile,...
        '\\paragraph*{Strikes}\\hfill\\\\\n\\noindent\n');
[latexFilePath,latexCommand]=saveBermudanSwaptionStrikes('BSwaptions');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:2:size(latexCommand,2)
    fprintf(inputFile,...
        '\\begin{table}\n');
    fprintf(inputFile,...
        '\\centering\n');
    fprintf(inputFile,...
            '\t%s\\hfill\\\\\n',latexCommand{iLC});
    fprintf(inputFile,...
        '\\caption{%s}\n',latexCommand{iLC+1});
    fprintf(inputFile,...
        '\\end{table}\n');
end
fprintf(inputFile,...
        '\\FloatBarrier\n');
%% CMS
fprintf(inputFile,...
        '\\section{Bermudan swaption prices}\n');
[latexFilePath,latexCommand]=saveCMSCIR('CMS');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:2:size(latexCommand,2)
    fprintf(inputFile,...
        '\\begin{table}\n');
    fprintf(inputFile,...
        '\\centering\n');
    fprintf(inputFile,...
            '\t%s\\hfill\\\\\n',latexCommand{iLC});
    fprintf(inputFile,...
        '\\caption{%s}\n',latexCommand{iLC+1});
    fprintf(inputFile,...
        '\\end{table}\n');
end
fprintf(inputFile,...
        '\\FloatBarrier\n');
%% Appendix
fprintf(inputFile,...
        '\\appendix\n');
% Market Data
fprintf(inputFile,...
        '\\section{Market Data}\n');
% Zero Coupon Curve
fprintf(inputFile,...
        '\\subsection{Zero Coupon Curve}\n');
[latexFilePath,latexCommand]=saveMarketData('MarketData');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
fprintf(inputFile,...
        '\\begin{table}\n');
fprintf(inputFile,...
        '\t%s\n',latexCommand);
fprintf(inputFile,...
        '\\centering\n');
fprintf(inputFile,...
        '\\caption{Market data containing the zero rate curve and zero coupon curve at %d/%d/%d.}\n',day,month,year);
fprintf(inputFile,...
        '\\label{tab:market_data%s}\n',identifier);
fprintf(inputFile,...
        '\\end{table}\n');
% Swaption
fprintf(inputFile,...
        '\\subsection{Swaption Data}\n');
[latexFilePath,latexCommand]=saveStrikeSwaption('MarketData');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
fprintf(inputFile,...
        '\\begin{table}\n');
fprintf(inputFile,...
        '\t%s\n',latexCommand);
fprintf(inputFile,...
        '\\centering\n');
fprintf(inputFile,...
        '\\caption{Market data containing the swaption strikes at %d/%d/%d.}\n',day,month,year);
fprintf(inputFile,...
        '\\label{tab:strike_swaption%s}\n',identifier);
fprintf(inputFile,...
        '\\end{table}\n');
[latexFilePath,latexCommand]=saveMarketSwaption('MarketData');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
fprintf(inputFile,...
        '\\begin{table}\n');
fprintf(inputFile,...
        '\t%s\n',latexCommand);
fprintf(inputFile,...
        '\\centering\n');
fprintf(inputFile,...
        '\\caption{Market data containing the swaption prices at %d/%d/%d.}\n',day,month,year);
fprintf(inputFile,...
        '\\label{tab:market_swaption%s}\n',identifier);
fprintf(inputFile,...
        '\\end{table}\n');
[latexFilePath,latexCommand]=saveVolSwaption('MarketData');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
fprintf(inputFile,...
        '\\begin{table}\n');
fprintf(inputFile,...
        '\t%s\n',latexCommand);
fprintf(inputFile,...
        '\\centering\n');
fprintf(inputFile,...
        '\\caption{Market data containing the swaption volatilities in bps at %d/%d/%d.}\n',day,month,year);
fprintf(inputFile,...
        '\\label{tab:market_swaption%s}\n',identifier);
fprintf(inputFile,...
        '\\end{table}\n');
fclose(inputFile);
%% Compile Latex
currFolder=cd(tempPath);
str1=sprintf('pdflatex %s',templatePath);
[returncode, ~] =system(str1);
% system(str1)
cd(currFolder);

copyfile(tempPath,copyPath);
% Renaming files
for iFile=1:1:size(copyFilePath,1)
    movefile(outputFilePath(iFile,:),copyFilePath(iFile,:))
end
% Delete auxiliary latex files
delete([copyPath,'\','template*.*']);
%% Functions
    function [latexFilePath,latexCommands]=saveCMSCIR(saveAt)
        latexCommands={};
        if strcmp(saveAt, '')
            latexFilePath='cmsCIR.tex';
        else
            latexFilePath=[saveAt,'\','cmsCIR.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
        file = fopen([tempPath,'\',latexFilePath],'a');
        for iP=1:1:size(paramsGC,2)
            latexCommand=['\cmsCIR',identifier,numDict{iP}];
            latexCommands{end+1}=latexCommand;
            fprintf(file,...
                    '\\newcommand{%s}{\n',latexCommand);
            fprintf(file,...
                    '\\renewcommand{\\arraystretch}{1}\n');
            fprintf(file,...
                    '\\begin{tabular}{|*{%d}{c}|}\n',6);
            fprintf(file,...
                    '\\hline\n');
            fprintf(file,...
                mat2Table('',...
                          {'Effective Date','Tenor','Index','Market CMS Rate','Model CMS Rate','Abs Error'},...
                          {},...
                          [effectiveCMS,tenorCMS,indexCMS,marketRateCMS,CMSrate(iP,:)',errCMSrate(iP,:)'],...
                          {'','%3.3g'},...
                          {'',''},...
                          'headHook','\\hline\n'));
            fprintf(file,...
                    '\\\\\\hline\n');
            fprintf(file,...
                    '\\end{tabular}\n');
            fprintf(file,...
                    '}\n');
            latexCommand=['\cmsCIRCaption',identifier,numDict{iP}];
            latexCommands{end+1}=latexCommand;
            fprintf(file,...
                    '\\newcommand{%s}{\n',latexCommand);
            fprintf(file,...
                    'CMS rates with parameters %s',methodsGC{iP});
            fprintf(file,...
                    '}\n');
        end
        fclose(file);
    end
    function [latexFilePath,latexCommands]=saveBermudanSwaptionErrors(saveAt)
        latexCommands={};
        if strcmp(saveAt, '')
            latexFilePath='bSwaptionErr.tex';
        else
            latexFilePath=[saveAt,'\','bSwaptionErr.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
        file = fopen([tempPath,'\',latexFilePath],'a');
        for iP=1:1:size(paramsGC,2)
            latexCommand=['\bSwaptionErr',identifier,numDict{iP}];
            latexCommands{end+1}=latexCommand;
            fprintf(file,...
                    '\\newcommand{%s}{\n',latexCommand);
            fprintf(file,...
                    '\\renewcommand{\\arraystretch}{1}\n');
            fprintf(file,...
                    '\\begin{tabular}{|c|*{%d}{c}|}\n',length(tenorBermudanSwap));
            fprintf(file,...
                    '\\hline\n');
            fprintf(file,...
                mat2Table('\\diagbox{Maturity}{Tenor}',...
                          tenorBermudanSwap,...
                          maturityBermudanSwap,...
                          squeeze(errBSwaption(iP,:,:)),...
                          {'%d','%1.3e'},...
                          {'',''},...
                          'headHook','\\hline\n'));
            fprintf(file,...
                    '\\\\\\hline\n');
            fprintf(file,...
                    '\\end{tabular}\n');
            fprintf(file,...
                    '}\n');
            latexCommand=['\bSwaptionErrCaption',identifier,numDict{iP}];
            latexCommands{end+1}=latexCommand;
            fprintf(file,...
                    '\\newcommand{%s}{\n',latexCommand);
            fprintf(file,...
                    'Absolute errors of Bermudan swaption with parameters %s and HW1. Mean error is equal to %3.3g',methodsGC{iP},mean(errBSwaption(iP,:,:),'all'));
            fprintf(file,...
                    '}\n');
        end
        fclose(file);
    end
    function [latexFilePath,latexCommands]=saveBermudanSwaptionStrikes(saveAt)
        latexCommands={};
        if strcmp(saveAt, '')
            latexFilePath='bSwaptionStrikes.tex';
        else
            latexFilePath=[saveAt,'\','bSwaptionStrikes.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
        file = fopen([tempPath,'\',latexFilePath],'a');
        latexCommand=['\bSwaptionStrikes',identifier];
        latexCommands{end+1}=latexCommand;
        fprintf(file,...
                '\\newcommand{%s}{\n',latexCommand);
        fprintf(file,...
                '\\renewcommand{\\arraystretch}{1}\n');
        fprintf(file,...
                '\\begin{tabular}{|c|*{%d}{c}|}\n',length(tenorBermudanSwap));
        fprintf(file,...
                '\\hline\n');
        fprintf(file,...
            mat2Table('\\diagbox{Maturity}{Tenor}',...
                      tenorBermudanSwap,...
                      maturityBermudanSwap,...
                      strikeBermudanSwap.*100,...
                      {'%d','%3.3g'},...
                      {'','\\,\\%%'},...
                      'headHook','\\hline\n'));
        fprintf(file,...
                '\\\\\\hline\n');
        fprintf(file,...
                '\\end{tabular}\n');
        fprintf(file,...
                '}\n');
        latexCommand=['\bSwaptionStrikesCaption',identifier];
        latexCommands{end+1}=latexCommand;
        fprintf(file,...
                '\\newcommand{%s}{\n',latexCommand);
        fprintf(file,...
                'Bermudan swaption strikes');
        fprintf(file,...
                '}\n');
        fclose(file);
    end
    function [latexFilePath,latexCommands]=saveBermudanSwaptionHW1(saveAt)
        latexCommands={};
        if strcmp(saveAt, '')
            latexFilePath='bSwaptionHW.tex';
        else
            latexFilePath=[saveAt,'\','bSwaptionHW.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
        file = fopen([tempPath,'\',latexFilePath],'a');
        latexCommand=['\bSwaptionHW',identifier];
        latexCommands{end+1}=latexCommand;
        fprintf(file,...
                '\\newcommand{%s}{\n',latexCommand);
        fprintf(file,...
                '\\renewcommand{\\arraystretch}{1}\n');
        fprintf(file,...
                '\\begin{tabular}{|c|*{%d}{c}|}\n',length(tenorBermudanSwap));
        fprintf(file,...
                '\\hline\n');
        fprintf(file,...
            mat2Table('\\diagbox{Maturity}{Tenor}',...
                      tenorBermudanSwap,...
                      maturityBermudanSwap,...
                      hw1BermudanSwapPrice.*100,...
                      {'%d','%3.3g'},...
                      {'','\\,\\%%'},...
                      'headHook','\\hline\n'));
        fprintf(file,...
                '\\\\\\hline\n');
        fprintf(file,...
                '\\end{tabular}\n');
        fprintf(file,...
                '}\n');
        latexCommand=['\bSwaptionHWCaption',identifier];
        latexCommands{end+1}=latexCommand;
        fprintf(file,...
                '\\newcommand{%s}{\n',latexCommand);
        fprintf(file,...
                'Bermudan swaption with HW1 model');
        fprintf(file,...
                '}\n');
        fclose(file);
    end
    function [latexFilePath,latexCommands]=saveBermudanSwaptionCIR(saveAt)
        latexCommands={};
        if strcmp(saveAt, '')
            latexFilePath='bSwaption.tex';
        else
            latexFilePath=[saveAt,'\','bSwaption.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
        file = fopen([tempPath,'\',latexFilePath],'a');
        for iP=1:1:size(paramsGC,2)
            latexCommand=['\bSwaption',identifier,numDict{iP}];
            latexCommands{end+1}=latexCommand;
            fprintf(file,...
                    '\\newcommand{%s}{\n',latexCommand);
            fprintf(file,...
                    '\\renewcommand{\\arraystretch}{1}\n');
            fprintf(file,...
                    '\\begin{tabular}{|c|*{%d}{c}|}\n',length(tenorBermudanSwap));
            fprintf(file,...
                    '\\hline\n');
            fprintf(file,...
                mat2Table('\\diagbox{Maturity}{Tenor}',...
                          tenorBermudanSwap,...
                          maturityBermudanSwap,...
                          squeeze(priceBSwaption(iP,:,:)).*100,...
                          {'%d','%3.3g'},...
                          {'','\\,\\%%'},...
                          'headHook','\\hline\n'));
            fprintf(file,...
                    '\\\\\\hline\n');
            fprintf(file,...
                    '\\end{tabular}\n');
            fprintf(file,...
                    '}\n');
            latexCommand=['\bSwaptionCaption',identifier,numDict{iP}];
            latexCommands{end+1}=latexCommand;
            fprintf(file,...
                    '\\newcommand{%s}{\n',latexCommand);
            fprintf(file,...
                    'Bermudan swaption with parameters %s',methodsGC{iP});
            fprintf(file,...
                    '}\n');
        end
        fclose(file);
    end
    function [latexFilePath,latexCommands]=saveSwaptionMarket(saveAt)
        latexCommands={};
        if strcmp(saveAt, '')
            latexFilePath='swaptionMarket.tex';
        else
            latexFilePath=[saveAt,'\','swaptionMarket.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
        file = fopen([tempPath,'\',latexFilePath],'a');
        latexCommand=['\swaptionMarket',identifier];
        latexCommands{end+1}=latexCommand;
        fprintf(file,...
                '\\newcommand{%s}{\n',latexCommand);
        fprintf(file,...
                '\\renewcommand{\\arraystretch}{1}\n');
        fprintf(file,...
                '\\begin{tabular}{|c|*{%d}{c}|}\n',length(tenorSwapTemp));
        fprintf(file,...
                '\\hline\n');
        fprintf(file,...
            mat2Table('\\diagbox{Maturity}{Tenor}',...
                      tenorSwapTemp,...
                      maturitySwapTemp,...
                      marketPrice,...
                      {'%d','%1.3e'},...
                      {'',''},...
                      'headHook','\\hline\n'));
        fprintf(file,...
                '\\\\\\hline\n');
        fprintf(file,...
                '\\end{tabular}\n');
        fprintf(file,...
                '}\n');
        latexCommand=['\swaptionMarketcaption',identifier];
        latexCommands{end+1}=latexCommand;
        fprintf(file,...
                '\\newcommand{%s}{\n',latexCommand);
        fprintf(file,...
                'Market Swaption surface');
        fprintf(file,...
                '}\n');
        fclose(file);
    end
    function [latexFilePath,latexCommands]=saveSwaptionMC(saveAt)
        latexCommands={};
        if strcmp(saveAt, '')
            latexFilePath='swaptionMC.tex';
        else
            latexFilePath=[saveAt,'\','swaptionMC.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
        file = fopen([tempPath,'\',latexFilePath],'a');
        for iP=1:1:size(paramsGC,2)
            latexCommand=['\swaptionMC',identifier,numDict{iP}];
            latexCommands{end+1}=latexCommand;
            fprintf(file,...
                    '\\newcommand{%s}{\n',latexCommand);
            fprintf(file,...
                    '\\renewcommand{\\arraystretch}{1}\n');
            fprintf(file,...
                    '\\begin{tabular}{|c|*{%d}{c}|}\n',length(tenorSwapTemp));
            fprintf(file,...
                    '\\hline\n');
            fprintf(file,...
                mat2Table('\\diagbox{Maturity}{Tenor}',...
                          tenorSwapTemp,...
                          maturitySwapTemp,...
                          sovMC{iP,1},...
                          {'%d','%1.3e'},...
                          {'',''},...
                          'headHook','\\hline\n'));
            fprintf(file,...
                    '\\\\\\hline\n');
            fprintf(file,...
                    '\\end{tabular}\n');
            fprintf(file,...
                    '}\n');
            latexCommand=['\swaptionMCcaption',identifier,numDict{iP}];
            latexCommands{end+1}=latexCommand;
            fprintf(file,...
                    '\\newcommand{%s}{\n',latexCommand);
            fprintf(file,...
                    'Swaption surface with parameters %s and MC',methodsGC{iP});
            fprintf(file,...
                    '}\n');
        end
        fclose(file);
    end
    function [latexFilePath,latexCommands]=saveSwaptionGC(saveAt)
        latexCommands={};
        if strcmp(saveAt, '')
            latexFilePath='swaptionGC.tex';
        else
            latexFilePath=[saveAt,'\','swaptionGC.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
        file = fopen([tempPath,'\',latexFilePath],'a');
        for iP=1:1:size(paramsGC,2)
            for iO=1:1:size(sovGC,2)
                latexCommand=['\swaptionGC',identifier,numDict{iP},numDict{iO+2}];
                latexCommands{end+1}=latexCommand;
                fprintf(file,...
                        '\\newcommand{%s}{\n',latexCommand);
                fprintf(file,...
                        '\\renewcommand{\\arraystretch}{1}\n');
                fprintf(file,...
                        '\\begin{tabular}{|c|*{%d}{c}|}\n',length(tenorSwapTemp));
                fprintf(file,...
                        '\\hline\n');
                fprintf(file,...
                    mat2Table('\\diagbox{Maturity}{Tenor}',...
                              tenorSwapTemp,...
                              maturitySwapTemp,...
                              sovGC{iP,iO},...
                              {'%d','%1.3e'},...
                              {'',''},...
                              'headHook','\\hline\n'));
                fprintf(file,...
                        '\\\\\\hline\n');
                fprintf(file,...
                        '\\end{tabular}\n');
                fprintf(file,...
                        '}\n');
                latexCommand=['\swaptionGCcaption',identifier,numDict{iP},numDict{iO+2}];
                latexCommands{end+1}=latexCommand;
                fprintf(file,...
                        '\\newcommand{%s}{\n',latexCommand);
                fprintf(file,...
                        'Swaption surface with parameters %s and GC%d',methodsGC{iP},iO+2);
                fprintf(file,...
                        '}\n');
            end
        end
        fclose(file);
    end
    function [latexFilePath,latexCommand]=saveErrGCParam(saveAt)
        latexCommand=['\errGCParam',identifier];
        if strcmp(saveAt, '')
            latexFilePath='errGCParam.tex';
        else
            latexFilePath=[saveAt,'\','errGCParam.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
        file = fopen([tempPath,'\',latexFilePath],'a');
        fprintf(file,...
                '\\newcommand{%s}{\n',latexCommand);
        fprintf(file,...
                '\\renewcommand{\\arraystretch}{1}\n');
        fprintf(file,...
                '\\begin{landscape}\n');
        fprintf(file,...
                '\\begin{tabular}{|c|*{%d}{c}|}\n',12);
        fprintf(file,...
                '\\hline\n');
        fprintf(file,...
            mat2Table('\\diagbox{Methods}{Errors}',...
                      {'Fmin','MC-M','GC3-MC','GC4-MC',...
                                          'GC5-MC','GC6-MC',...
                                          'GC7-MC','GC3-M',...
                                          'GC4-M','GC5-M',...
                                          'GC6-M','GC7-M'},...
                      methodsGC',...
                      cat(2,cell2mat(errorsGCFmin'),...
                            cell2mat(errMCM),....
                            cell2mat(errGCMC),...
                            cell2mat(errGCM)),...
                      {'%d','%1.2e'},...
                      {'',''},...
                      'headHook','\\hline\n'));
        fprintf(file,...
                '\\\\\\hline\n');
        fprintf(file,...
                '\\end{tabular}\n');
        fprintf(file,...
                '\\end{landscape}\n');
        fprintf(file,...
                '}\n');
        fclose(file);
    end
    function [latexFilePath,latexCommand]=saveCalParamGC(saveAt)
        latexCommand=['\calParamGC',identifier];
        if strcmp(saveAt, '')
            latexFilePath='calParamGC.tex';
        else
            latexFilePath=[saveAt,'\','calParamGC.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
        file = fopen([tempPath,'\',latexFilePath],'a');
        fprintf(file,...
                '\\newcommand{%s}{\n',latexCommand);
        fprintf(file,...
                '\\renewcommand{\\arraystretch}{1}\n');
        fprintf(file,...
                '\\begin{tabular}{|c|*{%d}{c}|}\n',8);
        fprintf(file,...
                '\\hline\n');
        head={'$\\phi_x^1$','$\\phi_x^2$','$\\phi_x^3$',...
              '$\\phi_y^1$','$\\phi_y^2$','$\\phi_y^3$',...
              '$x_0$','$y_0$'};
        fprintf(file,...
            mat2Table('\\diagbox{Methods}{Parameters}',...
                      head,...
                      methodsGC',...
                      cell2mat(paramsGC'),...
                      {'%d','%3.3g'},...
                      {'',''},...
                      'headHook','\\hline\n'));
        fprintf(file,...
                '\\\\\\hline\n');
        fprintf(file,...
                '\\end{tabular}\n');
        fprintf(file,...
                '}\n');
        fclose(file);
    end
    function [latexFilePath,latexCommand]=saveCalTimesGC(saveAt)
        latexCommand=['\calTimesGC',identifier];
        if strcmp(saveAt, '')
            latexFilePath='calTimesGC.tex';
        else
            latexFilePath=[saveAt,'\','calTimesGC.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
        file = fopen([tempPath,'\',latexFilePath],'a');
        fprintf(file,...
                '\\newcommand{%s}{\n',latexCommand);
        fprintf(file,...
                '\\renewcommand{\\arraystretch}{1}\n');
        fprintf(file,...
                '\\begin{tabular}{|c|*{%d}{c}|}\n',1);
        fprintf(file,...
                '\\hline\n');
        fprintf(file,...
            mat2Table('\\diagbox{Methods}{Comp Times}',...
                      {['GC with orders ',num2str(calMode)]},...
                      methodsGC',...
                      cell2mat(ctimesGC)',...
                      {'%d','%3.3g'},...
                      {'','\\,s'},...
                      'headHook','\\hline\n'));
        fprintf(file,...
                '\\\\\\hline\n');
        fprintf(file,...
                '\\end{tabular}\n');
        fprintf(file,...
                '}\n');
        fclose(file);
    end
    function [latexFilePath,latexCommand]=saveVolSwaption(saveAt)
        latexCommand=['\volSwaption',identifier];
        if strcmp(saveAt, '')
            latexFilePath='swaptionVol.tex';
        else
            latexFilePath=[saveAt,'\','swaptionVol.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
        file = fopen([tempPath,'\',latexFilePath],'a'); 
        fprintf(file,...
                '\\newcommand{%s}{\n',latexCommand);
        fprintf(file,...
                '\\renewcommand{\\arraystretch}{1}\n');
        fprintf(file,...
                '\\begin{tabular}{|c|*{%d}{c}|}\n',length(volTenorSwap));
        fprintf(file,...
                '\\hline\n');
        fprintf(file,...
            mat2Table('\\diagbox{Maturity}{Tenor}',...
                      volTenorSwap,...
                      volMatSwap,...
                      volSwap,...
                      {'%d','%3.3g'},...
                      {'',''},...
                      'headHook','\\hline\n'));
        fprintf(file,...
                '\\\\\\hline\n');
        fprintf(file,...
                '\\end{tabular}\n');
        fprintf(file,...
                '}\n');
        fclose(file);
    end
    function [latexFilePath,latexCommand]=saveMarketSwaption(saveAt)
        latexCommand=['\marketSwaption',identifier];
        if strcmp(saveAt, '')
            latexFilePath='swaptionMarket.tex';
        else
            latexFilePath=[saveAt,'\','swaptionMarket.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
        file = fopen([tempPath,'\',latexFilePath],'a'); 
        fprintf(file,...
                '\\newcommand{%s}{\n',latexCommand);
        fprintf(file,...
                '\\renewcommand{\\arraystretch}{1}\n');
        fprintf(file,...
                '\\begin{tabular}{|c|*{%d}{c}|}\n',length(tenorSwap));
        fprintf(file,...
                '\\hline\n');
        fprintf(file,...
            mat2Table('\\diagbox{Maturity}{Tenor}',...
                      tenorSwap,...
                      maturitySwap,...
                      marketSwapPrice,...
                      {'%d','%3.3e'},...
                      {'',''},...
                      'headHook','\\hline\n'));
        fprintf(file,...
                '\\\\\\hline\n');
        fprintf(file,...
                '\\end{tabular}\n');
        fprintf(file,...
                '}\n');
        fclose(file);
    end
    function [latexFilePath,latexCommand]=saveStrikeSwaption(saveAt)
        latexCommand=['\strikeSwaption',identifier];
        if strcmp(saveAt, '')
            latexFilePath='swaptionStrike.tex';
        else
            latexFilePath=[saveAt,'\','swaptionStrike.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
        file = fopen([tempPath,'\',latexFilePath],'a');
        fprintf(file,...
                '\\newcommand{%s}{\n',latexCommand);
        fprintf(file,...
                '\\renewcommand{\\arraystretch}{1}\n');
        fprintf(file,...
                '\\begin{tabular}{|c|*{%d}{c}|}\n',length(tenorSwap));
        fprintf(file,...
                '\\hline\n');
        fprintf(file,...
                mat2Table('\\diagbox{Maturity}{Tenor}',...
                          tenorSwap,...
                          maturitySwap,...
                          100*strikeSwap,...
                          {'%d','%3.3e'},...
                          {'','\\,\\%%'},...
                          'headHook','\\hline\n'));
        fprintf(file,...
                '\\\\\\hline\n');
        fprintf(file,...
                '\\end{tabular}\n');
        fprintf(file,...
                '}\n');
        fclose(file);
    end    
    function [latexFilePath,latexCommand]=saveMarketData(saveAt)
        latexCommand=['\marketData',identifier];
        if strcmp(saveAt, '')
            latexFilePath='marketData.tex';
        else
            latexFilePath=[saveAt,'\','marketData.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
        file = fopen([tempPath,'\',latexFilePath],'a');
        fprintf(file,...
                '\\newcommand{%s}{\n',latexCommand);
        fprintf(file,...
                '\\renewcommand{\\arraystretch}{1}\n');
        fprintf(file,...
                '\\begin{tabular}{|*{%d}{c}|}\n',3);
        fprintf(file,'\\hline\n');
        fprintf(file,...
                mat2Table('',...
                          {'Maturity (in years)','Zero rate (in \\%%)','Zero-coupon price'},...
                          marketTimes,...
                          [100*marketZeroRates,marketDF],...
                          {'%g','%3.9e'},...
                          {'',''},...
                          'headHook','\\hline\n'));
        fprintf(file,'\\\\\\hline\n');
        fprintf(file,...
                '\\end{tabular}\n');
        fprintf(file,...
                '}\n');
        fclose(file);
    end
    function latexFilePath=saveFigures(figures,saveAt,figName)
        if strcmp(figName,'')
            figName='fig';
        end
        if strcmp(saveAt, '')
            latexFilePath='figure.tex';
            latexFolderPath=tempPath;
        else
            latexFilePath=[saveAt,'\','figure.tex'];
            latexFolderPath=[tempPath,'\',saveAt];
            mkDir([tempPath,'\',saveAt]);
        end
        file = fopen([tempPath,'\',latexFilePath],'a');
        for i=1:1:length(figures)
            picName=sprintf('%s_%d',figName,i);
            picPath = [latexFolderPath,'\',picName,'.',picType];
            if strcmp(saveAt, '')
                relPicPath=picName;
            else
                relPicPath=[saveAt,'\',picName];
            end
            saveas(figures(i),picPath,saveParam);
            fprintf(file,...
                '\\begin{landscape}\n');
            fprintf(file,...
                '\\includegraphics[width=.95\\columnwidth]{%s}\n',...
                changeSlash(relPicPath));
            fprintf(file,...
                '\\end{landscape}\n');
        end
        fclose(file);
    end


end
function str=changeSlash(str)
    for i=1:1:length(str)
        if strcmp(str(i),'\')
            str(i)='/';
        end
    end
end
function delFile(file)
    if exist(file)
        delete(file);
    end
end
function delDir(dir)
    if exist(dir)==7
        rmdir(dir,'s');
    end
end
function cleanDir(mdir,except)
    except{end+1}='.';
    except{end+1}='..';
    for d = dir(mdir).'
      if ~any(strcmp(d.name,except))
          if d.isdir
              rmdir([d.folder,'\',d.name],'s');
          else
              delete([d.folder,'\',d.name]);
          end
      end
    end
end
function mkDir(dir)
    if exist(dir)==0
        mkdir(dir);
    end
end
function latexStr=mat2Table(corner,head,index,body,formatSpec,unit,varargin)
%     percent='\\,\\%%';
    headHook='';
    bodyHook='';
    for k=1:1:length(varargin)
        switch varargin{k}
            case 'headHook'
                headHook=varargin{k+1};
            case 'bodyHook'
                bodyHook=varargin{k+1};
        end
    end

    latexStr='';
    function x=functor(x,formatSpec,unit)
        if ~ischar(x)
            if mod(abs(x),1)==0 
                x=num2str(x,'%d');
            else
                x=num2str(x,formatSpec);
            end
        end
        x=[x,unit];
    end
    if ~isempty(corner)
        cornerStr=functor(corner,formatSpec,'');
        latexStr=[latexStr,cornerStr,' & '];
    end
    if ~isempty(head)
        if iscell(head)
            headCell=cellfun(@(x)functor(x,formatSpec,''),head,...
                'UniformOutput',false);
        else
            headCell=arrayfun(@(x)functor(x,formatSpec,''),head,...
                'UniformOutput',false);
        end
        headStr=join(headCell,' & ');
        latexStr=[latexStr,headStr{1},'\\\\\n'];
    end
    latexStr=[latexStr,headHook];
    bodyCell=arrayfun(@(x)functor(x,formatSpec{2},unit{2}),body,...
        'UniformOutput',false);
    if ~isempty(index)
        if iscell(index)
            indexCell=cellfun(@(x)functor(x,formatSpec{1},unit{1}),index,...
                'UniformOutput',false);
        else
            indexCell=arrayfun(@(x)functor(x,formatSpec{1},unit{1}),index,...
                'UniformOutput',false);
        end
        bodyCell=cat(2,indexCell,bodyCell);
    end
    bodyStr=join(join(bodyCell,' & ',2),'\\\\\n',1);
    latexStr=[latexStr,bodyStr{1}];
    latexStr=[latexStr,bodyHook];
end