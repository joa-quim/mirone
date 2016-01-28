%#mex
function calc_omori(events_time,events_mag,h_mir_fig);
%   This was plot_llkstest ZMAP function. This part of ZMAP is better written
%   I assembled in this file all functions needed to do the job and it is here
%   also that the input parameters are requested via the inputdlg function
%   Mironified by J. Luis at 6-07-05

% WARNING: WHEN COMPILING NEEDS TO INCLUDE THE CALC_OMORI_AUX.M FILE AND IT WILL
%		   ALSO BE NECESSARY TO TEST THAT THE COMPILED VERSION (DLL) WORKS SINCE
%		   PREVIOUS VERSION WAS COMPILED WITH BOTH FILES MERGED.

% function plot_llkstest(a,time,timef,bootloops,maepi);
% Plots Ncum observed vs. Ncum modeled for specified time windows
% with choosing the model for the learning period and performs a KS-Test
%
% Input variables:      (in the original plot_llkstest function)
% a         : earthquake catalog
% time      : learning period fo fit Omori parameters
% timef     : forecast period
% bootloops : Number of bootstraps
% maepi     : mainshock
%
% J.Woessner
% last update: 20.07.04

% $Id: calc_omori.m 4270 2014-01-11 22:20:41Z j $

% Get input parameters
prompt = {'Enter length of learning period (days)','Enter number of bootstraps:',...
        'Enter model number (1:pck, 2:pckk, 3:ppckk, 4:ppcckk:'};
tit    = 'Parameters';
def    = {'50','50','1'};
answer = inputdlg(prompt,tit,[1 50; 1 50; 1 50],def);
if isempty(answer);    return;     end
time = str2double(answer{1});
bootloops = str2double(answer{2});
nMod = str2double(answer{3});
timef = 1;      % They use this as default

%[m_main main] = max(a(:,6));
% date_matlab = datenum(floor(a(:,3)),a(:,4),a(:,5),a(:,8),a(:,9),zeros(size(a,1),1)); 
% date_main = datenum(floor(maepi(3)),maepi(4),maepi(5),maepi(8),maepi(9),0); 
% time_aftershock = date_matlab-date_main;

[maepi,ind_magMax] = max(events_mag);
time_aftershock = (events_time - events_time(ind_magMax)) * 365;    % Wrong for leap-years

% Aftershock catalog
vSel1 = (time_aftershock(:) > 0);
tas = time_aftershock(vSel1); 
%eqcatalogue = a(vSel1,:);
eqcatalogue = [events_time(vSel1) events_mag(vSel1)];   % I have join them otherwise is too complicated

% Estimation of Omori parameters from learning period
l = (tas <= time);
time_as = tas(l);
if (isempty(time_as))
    warndlg('There are no events inside the learning period. Quiting.','Warning')
    return
end
% Times up to the forecast time
lf = (tas <= time+timef);
time_asf = tas(lf);
time_asf = sort(time_asf);          % ACHO QUE ISTO NAO E PRECISO

% Select biggest aftershock earliest in time, but more than 1 day after
% mainshock and in learning period
mAfLearnCat = eqcatalogue(l,:);     clear eqcatalogue;
ft_c = 1/365;               % Time not considered to find biggest aftershock
%vSel = ((mAfLearnCat(:,3) > maepi(:,3)+ft_c) & (mAfLearnCat(:,3) <= maepi(:,3)+time/365));
vSel = ( (mAfLearnCat(:,1) > events_time(ind_magMax)+ft_c) & (mAfLearnCat(:,1) <= events_time(ind_magMax)+time/365) );
mCat = mAfLearnCat(vSel,:);
vSel = (mCat(:,2) == max(mCat(:,2)));
vBigAf = mCat(vSel,:);
if (length(mCat) > 1)
    vSel = ( fix(vBigAf(:,1)) == min(fix(vBigAf(:,1))) );
    vBigAf = vBigAf(vSel,:);
end

%date_biga = datenum(floor(vBigAf(3)),vBigAf(4),vBigAf(5),vBigAf(8),vBigAf(9),0); 
%fT1 = date_biga - date_main; % Time of big aftershock
if (isempty(vBigAf))
    warndlg('No big after-shock found (too short swarm?). Quiting','Warning')
    return
end
fT1 = (vBigAf(1) - events_time(ind_magMax)) * 365;      % Time of big aftershock

% Calculate uncertainty and mean values of p,c,and k
[mMedModF, mStdL, loopout] = brutebootloglike_a2(time_as, time_asf, bootloops,fT1,nMod);
pmed1 = mMedModF(1,1);      pmed2 = mMedModF(1,3);
cmed1 = mMedModF(1,5);      cmed2 = mMedModF(1,7);
kmed1 = mMedModF(1,9);      kmed2 = mMedModF(1,11);

% Compute model according to model choice
[pval1, pval2, cval1, cval2, kval1, kval2, fAIC, fL] = bruteforceloglike_a2(time_as,fT1,nMod);

% Start plotting
h_fig = figure('Numbertitle','off','Name','Aftershock modelling fit','Visible','off','Color',get(0,'factoryUicontrolBackgroundColor'));
uimenu('Label','Save file','callback',{@calSave,h_fig,h_mir_fig});

% Plot the forecast ... 
cumnrf = (1:length(time_as))'; 
cumnr_modelf = [];
if (nMod == 1)
    for i=1:length(time_as)
        if pval1 ~= 1
            cm = kval1/(pval1-1)*(cval1^(1-pval1)-(time_as(i)+cval1)^(1-pval1));
        else
            cm = kval1*log(time_as(i)/cval1+1);
        end
        cumnr_modelf = [cumnr_modelf; cm];
    end
else
    for i=1:length(time_as)
        if (time_as(i) <= fT1)
            if pval1 ~= 1
                cm = kval1/(pval1-1)*(cval1^(1-pval1)-(time_as(i)+cval1)^(1-pval1));
            else
                cm = kval1*log(time_as(i)/cval1+1);
            end
            cumnr_modelf = [cumnr_modelf; cm];
        else
            if (pval1 ~= 1 & pval2 ~= 1)
                cm = kval1/(pval1-1)*(cval1^(1-pval1)-(time_as(i)+cval1)^(1-pval1))+ kval2/(pval2-1)*(cval2^(1-pval2)-(time_as(i)-fT1+cval2)^(1-pval2));
            else
                cm = kval1*log(time_as(i)/cval1+1) + kval2*log((time_as(i)-fT1)/cval2+1);
            end
            cumnr_modelf = [cumnr_modelf; cm];
        end
    end
end; % End of if on nMod
time_as=sort(time_as);
cumnr_modelf=sort(cumnr_modelf);

pf1 =  plot(time_as,cumnr_modelf,'g-.','Linewidth',2,'Tag','Model');
hold on
pf2 =  plot(time_as,cumnrf, 'b-','Linewidth',2,'Tag','Orig');
paf = plot(fT1, 0,'h','MarkerFaceColor',[1 1 0],'MarkerSize',12,'MarkerEdgeColor',[0 0 0] );

% Calculate KSTEST2 as a measure of the goodness of fit
[H,P,KSSTAT] = calc_omori_aux('kstest2', cumnr_modelf,cumnrf);

% Calculate RMS
i=(1:1:length(time_as))';
fRMS = (sum((i-cumnr_modelf).^2)/length(i))^0.5;

% Round values for output
pval1 = round(100*pval1)/100;       pval2 = round(100*pval2)/100;
cval1 = round(1000*cval1)/1000;     cval2 = round(1000*cval2)/1000;
kval1 = round(10*kval1)/10;         kval2 = round(10*kval2)/10;
pmed1 = round(100*pmed1)/100;       mStdL(1,1) = round(100*mStdL(1,1))/100;
pmed2 = round(100*pmed2)/100;       mStdL(1,2) = round(100*mStdL(1,2))/100;
cmed1 = round(1000*cmed1)/1000;     mStdL(1,3) = round(1000*mStdL(1,3))/1000;
cmed2 = round(1000*cmed2)/1000;     mStdL(1,4) = round(1000*mStdL(1,4))/1000;
kmed1 = round(10*kmed1)/10;         mStdL(1,5) = round(100*mStdL(1,5))/100;
kmed2 = round(10*kmed2)/10;         mStdL(1,6)= round(100*mStdL(1,6))/100;
fRMS = round(100*fRMS)/100;

% Get Y limits for positioning texts
yy = get(gca,'ylim');

if (nMod == 1)
    string1=['p = ' num2str(pval1) '; c = ' num2str(cval1) '; k = ' num2str(kval1) ];
    string3=['pm = ' num2str(pmed1) '+-' num2str(mStdL(1,1)) '; cm = ' num2str(cmed1) '+-' num2str(mStdL(1,3)) '; km = ' num2str(kmed1) '+-' num2str(mStdL(1,5))];
    text(max(time_asf)*0.05,yy(2)*0.9,string1,'FontSize',10);
    text(max(time_asf)*0.05,yy(2)*0.8,string3,'FontSize',10);
elseif (nMod == 2)
    string1=['p = ' num2str(pval1) '; c = ' num2str(cval1) '; k1 = ' num2str(kval1) '; k2 = ' num2str(kval2) ];
    string3=['pm = ' num2str(pmed1) '+-' num2str(mStdL(1,1)) '; cm = ' num2str(cmed1) '+-' num2str(mStdL(1,3)) '; km1 = ' num2str(kmed1) '+-' num2str(mStdL(1,5)) '; km2 = ' num2str(kmed2) '+-' num2str(mStdL(1,6))];
    text(max(time_asf)*0.05,yy(2)*0.9,string1,'FontSize',10);
    text(max(time_asf)*0.05,yy(2)*0.8,string3,'FontSize',10);
elseif (nMod == 3)
    string1=['p1 = ' num2str(pval1) '; c = ' num2str(cval1) '; k1 = ' num2str(kval1) ];
    string2=['p2 = ' num2str(pval2) '; k2 = ' num2str(kval2) ];
    string3=['pm1 = ' num2str(pmed1) '+-' num2str(mStdL(1,1)) '; cm = ' num2str(cmed1) '+-' num2str(mStdL(1,3)) '; km1 = ' num2str(kmed1) '+-' num2str(mStdL(1,5))];
    string4=['pm2 = ' num2str(pmed2) '+-' num2str(mStdL(1,2)) '; km2 = ' num2str(kmed2) '+-' num2str(mStdL(1,6))];
    text(max(time_asf)*0.05,yy(2)*0.9,string1,'FontSize',10);
    text(max(time_asf)*0.05,yy(2)*0.85,string2,'FontSize',10);
    text(max(time_asf)*0.05,yy(2)*0.8,string3,'FontSize',10);
    text(max(time_asf)*0.05,yy(2)*0.75,string4,'FontSize',10);
else
    string1=['p1 = ' num2str(pval1) '; c1 = ' num2str(cval1) '; k1 = ' num2str(kval1) ];
    string2=['p2 = ' num2str(pval2) '; c2 = ' num2str(cval2) '; k2 = ' num2str(kval2) ];
    string3=['pm1 = ' num2str(pmed1) '+-' num2str(mStdL(1,1)) '; cm1 = ' num2str(cmed1) '+-' num2str(mStdL(1,3)) '; km1 = ' num2str(kmed1) '+-' num2str(mStdL(1,5))];
    string4=['pm2 = ' num2str(pmed2) '+-' num2str(mStdL(1,2)) '; cm2 = ' num2str(cmed2) '+-' num2str(mStdL(1,4)) '; km2 = ' num2str(kmed2) '+-' num2str(mStdL(1,6))];
    text(max(time_asf)*0.05,yy(2)*0.9,string1,'FontSize',10);
    text(max(time_asf)*0.05,yy(2)*0.85,string2,'FontSize',10);
    text(max(time_asf)*0.05,yy(2)*0.8,string3,'FontSize',10);
    text(max(time_asf)*0.05,yy(2)*0.75,string4,'FontSize',10);
end
string=['H = ' num2str(H) ' P = ' num2str(P) ' KS-Statistic) = ' num2str(KSSTAT)];
text(max(time_asf)*0.05,yy(2)*0.1,string,'FontSize',10);
text(max(time_asf)*0.05,yy(2)*0.05,['AIC = ' num2str(fAIC)],'FontSize',10);
text(max(time_asf)*0.05,yy(2)*0.15,['RMS = ' num2str(fRMS)],'FontSize',10);
% Legend
legend([pf2 pf1 paf],'Data',['Model ' num2str(nMod)],'Sec. AF',0);


% Calculate fits of different models
mRes = [];
% Modified Omori law (pck)
nMod = 1; [pval1, pval2, cval1, cval2, kval1, kval2, fAIC, fL] = bruteforceloglike_a2(time_as,fT1,nMod);
mRes = [mRes; nMod, pval1, pval2, cval1, cval2, kval1, kval2, fAIC, fL];
% MOL with secondary aftershock (pckk)
nMod = 2; [pval1, pval2, cval1, cval2, kval1, kval2, fAIC, fL] = bruteforceloglike_a2(time_as,fT1,nMod);
mRes = [mRes; nMod, pval1, pval2, cval1, cval2, kval1, kval2, fAIC, fL];
% MOL with secondary aftershock (ppckk)
nMod = 3; [pval1, pval2, cval1, cval2, kval1, kval2, fAIC, fL] = bruteforceloglike_a2(time_as,fT1,nMod);
mRes = [mRes; nMod, pval1, pval2, cval1, cval2, kval1, kval2, fAIC, fL];
% MOL with secondary aftershock (ppcckk)
nMod = 4; [pval1, pval2, cval1, cval2, kval1, kval2, fAIC, fL] = bruteforceloglike_a2(time_as,fT1,nMod);
mRes = [mRes; nMod, pval1, pval2, cval1, cval2, kval1, kval2, fAIC, fL];

% Select best fitting model by AIC
vSel = (mRes(:,8) == min(mRes(:,8)));
mRes = mRes(vSel,:);
if (length(mRes(:,1)) > 1)
    vSel = (mRes(:,1) == min(mRes(:,1)));
    mRes = mRes(vSel,:);
end
% Model to use for bootstrapping as of lowest AIC to observed data
nMod = mRes(1,1);

text(max(time_asf)*0.05,yy(2)*0.2,['Info: Best model is ' num2str(nMod)],'FontSize',10);
% Figure settings
set(gca,'Fontsize',12,'Linewidth',2)
xlabel('Time [Days after mainshock]','Fontsize',12,'Fontweight','bold')
set(h_fig,'Visible','on')

% -------------------------------------------------------------------------------
function calSave(obj,eventdata,h_fig,h_mir_fig)
% Save the curves in file. [x_dat y_dat y_mod]
	handles_mir = guidata(h_mir_fig);		% Get the Mirone handles structure

	[FileName,PathName] = put_or_get_file(handles_mir, ...
		{'*.dat;*.DAT', 'Gutt Rich file (*.dat,*.DAT)'; '*.*', 'All Files (*.*)'},'Select File name','put','.dat');
	if isequal(FileName,0),		return,		end
	fname = [PathName FileName];

	% Fish the data
	h_mod = findobj(h_fig,'Tag','Model');		h_dat = findobj(h_fig,'Tag','Orig');
	x_mod = get(h_mod,'XData');					y_mod = get(h_mod,'YData');
	x_dat = get(h_dat,'XData');					y_dat = get(h_dat,'YData');

	fid = fopen(fname, 'w');
	if (fid < 0)    errordlg(['Can''t open file:  ' fname],'Error');    return;     end
	fprintf(fid,'# Time (days)\tCum number\tCum number (model)\n');
	fprintf(fid,'%.7f\t%d\t%.2f\n',[x_dat; y_dat; y_mod]);
	fclose(fid);

% -----------------------------------------------------------------------------
function [mMedModF, mStdL, loopout] = brutebootloglike_a2(time_as, time_asf, bootloops,fT1, nMod);
% Bootstrap analysis of Omori parameters calculated by bruteforce.m
% (p1,p2,c1,c2,k1,k2)-pair is mean of the bootstrap values by determining the mean cumulative number modeled a end of the learning period
% Standard deviations are calculated as the 2nd moment, not to rely fully on normal distributions
%
% Input parameters:
%   time_as     Delay times [days] of learning period
%   time_asf    Delay times [days] until end of forecast period
%   bootloops   Number of bootstraps
%   fT1         Time of biggest aftershock in learning period
%   nMod        Model to fit data, three models including a secondary aftershock sequence.
%               Different models have varying amount of free parameters
%               before (p1,c1,k1) and after (p2,c2,k2) the aftershock occurence
%               1: modified Omori law (MOL): 3 free parameters
%                  p1=p2,c1=c2,k1=k2
%               2: MOL with one secondary aftershock sequence:4 free parameters
%                  p1=p2,c1=c2,k1~=k2
%               3: MOL with one secondary aftershock sequence:5 free parameters
%                  p1~=p2,c1=c2,k1~=k2
%               4: MOL with one secondary aftershock sequence:6 free parameters
%                  p1~=p2,c1~=c2,k1~=k2
%
% Output parameters:
%  mMedModF :  Result matrix including the values for the mean forecast at end of forecast period
%  mStdL    :  Uncertainties of fit to the data in learning period
%  loopout     contains all results
%
% Samuel Neukomm / S. Wiemer / J. Woessner
% last update: 05.08.03

time_as = sort(time_as);
n = length(time_as);
loopout = []; 
% Initialize random seed
rand('seed',sum(100*clock));
hWaitbar1 = aguentabar(0,'title','Bootstrapping...');
set(hWaitbar1,'Numbertitle','off','Name','Bootstap Omori parameters')
for (j = 1:bootloops)
    randnr = ceil(rand(n,1)*n);
    i = (1:n)';
    newtas(i,:) = time_as(randnr(i),:); % bootstrap sample
    newtas = sort(newtas);
    [pv1, pv2, cv1, cv2, kv1, kv2, fAIC, fL] = bruteforceloglike_a2(newtas, fT1, nMod);
    loopout = [loopout; pv1, pv2, cv1, cv2, kv1, kv2, fAIC, fL];
    aguentabar(j/bootloops)
end

% New version: Choose mean (p,c,k)-variables by modelling the cumulative number at end of 
% the learning period

% 2nd moment i.e. Standard deviations
[pstd1] = calc_StdDev(loopout(:,1));    [pstd2] = calc_StdDev(loopout(:,2));
[cstd1] = calc_StdDev(loopout(:,3));    [cstd2] = calc_StdDev(loopout(:,4));
[kstd1] = calc_StdDev(loopout(:,5));    [kstd2] = calc_StdDev(loopout(:,6));

% Uncertainties of fit
mStdL = [pstd1 pstd2 cstd1 cstd2 kstd1 kstd2];

% Compute best fitting pair of variates
loopout = [ loopout , loopout(:,1)*0];
for (j = 1:length(loopout(:,1)))
    cumnr = (1:length(time_asf))'; 
    cumnr_model = [];     
    pval1 = loopout(j,1);    pval2 = loopout(j,2);
    cval1 = loopout(j,3);    cval2 = loopout(j,4);
    kval1 = loopout(j,5);    kval2 = loopout(j,6);
    if (nMod == 1)
        for (i=1:length(time_asf))
            if (pval1 ~= 1)
                cm = kval1/(pval1-1)*(cval1^(1-pval1)-(time_asf(i)+cval1)^(1-pval1));
            else
                cm = kval1*log(time_asf(i)/cval1+1);
            end
            cumnr_model = [cumnr_model; cm];
        end
        loopout(j,9) = max(cumnr_model);
    else
        for (i=1:length(time_asf))
            if (time_asf(i) <= fT1)
                if (pval1 ~= 1)
                    cm = kval1/(pval1-1)*(cval1^(1-pval1)-(time_asf(i)+cval1)^(1-pval1));
                else
                    cm = kval1*log(time_asf(i)/cval1+1);
                end
                cumnr_model = [cumnr_model; cm];
            else
                if (pval1 ~= 1 & pval2 ~= 1)
                    cm = kval1/(pval1-1)*(cval1^(1-pval1)-(time_asf(i)+cval1)^(1-pval1))+ kval2/(pval2-1)*(cval2^(1-pval2)-(time_asf(i)-fT1+cval2)^(1-pval2));
                else
                    cm = kval1*log(time_asf(i)/cval1+1) + kval2*log((time_asf(i)-fT1)/cval2+1);
                end
                cumnr_model = [cumnr_model; cm];
            end
        end
        loopout(j,9) = max(cumnr_model);
    end; % End of if on nMod
end

[Y in] = sort(loopout(:,9));
loops = loopout(in,:);

% Mean values
vMean = abs(loops(:,9)-mean(loops(:,9)));
nMean = (find(vMean == min(vMean)));

if (length(nMean(:,1)) > 1)     nMean = nMean(1,1);  end
pMean1 = loops(nMean,1);        pMean2 = loops(nMean,2);
cMean1 = loops(nMean,3);        cMean2 = loops(nMean,4);
kMean1 = loops(nMean,5);        kMean2 = loops(nMean,6);

mMedModF = [pMean1, pstd1, pMean2, pstd2, cMean1, cstd1, cMean2, cstd2, kMean1, kstd1, kMean2, kstd2];

% -----------------------------------------------
function [fStdDev] = calc_StdDev(vDistribution)
% Computes the standard deviation of a non-parameterized distribution 
%
% Input:   vDistribution -> Vector containing the distribution
% Output:  fStdDev       -> Standard deviation of the given distribution
%
% Danijel Schorlemmer - March 10, 2003

nNumElements = length(vDistribution);       % Get number of elements of distribution

% Compute standard deviation as the second moment of the given distribution
vDist = vDistribution.^2 .* (1/nNumElements);
fStdDev = sqrt(nansum(vDist)-(nanmean(vDistribution))^2);

% Check result for imaginary part
% Explanation: For a large vector vDistribution with all the same numbers
% the argument of the sqrt-coomand may be not zero although it should. 
% Thus the sqrt outputs a complex number, although it should be zero! 
% This is fixed now by the following if command.
if (~isreal(fStdDev) && ~isnan(fStdDev))    fStdDev = 0;     end

% --------------------------------------------------------------------------------------
function [pv1, pv2, cv1, cv2, kv1, kv2, fAIC, fL] = bruteforceloglike_a2(tas, fT1, nMod)
% Calculates by a constrained grid search the parameters of the modified Omori formula 
% using the log likelihood function by Ogata (1983) and calculates the best model by the
% corrected AIC (Burnham & Anderson(2002)
%
% Input parameters:
%   tas     Delay time of entire sequence [days]
%   fT1     Time of large aftershock after mainshock in days
%   nMod    Model number for choosing amount of free parameters
%
% Output parameters:
%   pv1         p value before large aftershock
%   pv2         p value after large aftershock
%   cv1         c value before large aftershock
%   cv2         c value after large aftershock
%   kv1         k1 value before large aftershock
%   kv2         k2 value after large aftershock
%   fAIC        Akaike Information Criterion value
%   fL          Maximum likelihood estunation
%
% J. Woessner / S. Neukomm
% last update: 20.04.04

options = struct('ActiveConstrTol', [], 'DerivativeCheck', [], 'Diagnostics', [], 'DiffMaxChange', [], ...
    'DiffMinChange', [], 'Display', 'none', 'GoalsExactAchieve', [],'GradConstr', [], 'GradObj', [], ...
    'Hessian', [], 'HessMult', [], 'HessPattern', [], 'HessUpdate', [], 'Jacobian', [], 'JacobMult', [], ...
    'JacobPattern', [], 'LargeScale', [], 'LevenbergMarquardt', [], 'LineSearchType', [], ...
    'MaxFunEvals', 400, 'MaxIter', 500, 'MaxPCGIter', [], 'MaxSQPIter', 500, 'MeritFunction', [], ...
    'MinAbsMax', [], 'NonlEqnAlgorithm', [], 'OutputFcn', [], 'Preconditioner', [], 'PrecondBandWidth', [], ...
    'ShowStatusWindow', [], 'Simplex', [], 'TolCon', [], 'TolFun', 1e-04, 'TolPCG', [], ...
    'TolX', [], 'TypicalX', []);

% Starting values
fPmin = 0.2;    fPmax = 2.7;    fCmin = 0.01;   fCmax = 5;
fKmin1 = 10;    fKmax1 = 5000;  fKmin2 = 10;    fKmax2 = 4000;

if (nMod == 1)      % Modified Omori law.  3 free parameters: p, c , k
    fPar = 3;
    vStartValues = [1.1 0.5 50];
    %[vValues, fL] = fmincon_m('bruteloglike', vStartValues, [], [], [], [],...
    [vValues, fL] = calc_omori_aux('fmincon_m',@bruteloglike, vStartValues, [], [], [], [],...
        [fPmin fCmin fKmin1 ], [fPmax fCmax fKmax1 ], [], options, tas);
    pv1 = vValues(1);    pv2 = vValues(1);    cv1 = vValues(2);
    cv2 = vValues(2);    kv1 = vValues(3);    kv2 = vValues(3);
elseif (nMod == 2)      % 4 free parameters: p, c , k1, k2
    fPar = 4;
    vStartValues = [1.1 0.5 50 50];
    %[vValues, fL] = fmincon_m('bruteloglike_pck2', vStartValues, [], [], [], [],...
    [vValues, fL] = calc_omori_aux('fmincon_m', @bruteloglike_pck2, vStartValues, [], [], [], [],...
        [fPmin fCmin fKmin1 fKmin2], [fPmax fCmax fKmax1 fKmax2], [], options, tas, fT1);
    pv1 = vValues(1);    pv2 = vValues(1);    cv1 = vValues(2);
    cv2 = vValues(2);    kv1 = vValues(3);    kv2 = vValues(4);
elseif (nMod == 3)      % 5 free parameters: p1,p2,c,k1,k2
    fPar = 5;
    vStartValues = [1.1 1.1 0.5 0.5 50 50];
    %[vValues, fL] = fmincon_m('bruteloglike_p2ck2', vStartValues, [], [], [], [],...
    [vValues, fL] = calc_omori_aux('fmincon_m', @bruteloglike_p2ck2, vStartValues, [], [], [], [],...
        [fPmin fPmin fCmin fCmin fKmin1 fKmin2], [fPmax fPmax fCmax fCmax fKmax1 fKmax2], [], options, tas, fT1);
    pv1 = vValues(1);    pv2 = vValues(2);    cv1 = vValues(3);
    cv2 = vValues(3);    kv1 = vValues(5);    kv2 = vValues(6);
else                    % 6 free parameters: p1,p2,c1, c2,k1,k2
    fPar = 6;
    vStartValues = [1.1 1.1 0.5 0.5 50 50];
    %[vValues, fL] = fmincon_m('bruteloglike_p2c2k2', vStartValues, [], [], [], [],...
    [vValues, fL] = calc_omori_aux('fmincon_m', @bruteloglike_p2c2k2, vStartValues, [], [], [], [],...
        [fPmin fPmin fCmin fCmin fKmin1 fKmin2], [fPmax fPmax fCmax fCmax fKmax1 fKmax2], [], options, tas, fT1);
    pv1 = vValues(1);    pv2 = vValues(2);    cv1 = vValues(3);
    cv2 = vValues(4);    kv1 = vValues(5);    kv2 = vValues(6);
end  
% corrected Akaike Information Criterion
[fk,nX]=size(tas);
fAIC = -2*(-fL)+2*fPar+2*fPar*(fPar+1)/(fk-fPar-1);

%----------------------------------------------------------------------------------
function fL = bruteloglike(vValues,time_as)
% This function calculates the log likelihood 
% function for the modeled aftershock sequence and the maximum likelihood estimate for k, c and p
% Reference: Ogata, Estimation of the parameters in the modified Omori formula 
% for aftershock sequences by  the maximum likelihood procedure, J. Phys. Earth, 1983
% (Formula 6)
%
% J. Woessner
% last update: 29.07.03

p = vValues(1);             c = vValues(2);     k = vValues(3);
fTstart = min(time_as);     fTend = max(time_as);   % Setting start end end time

if (p ~= 1)
    fAcp = ((fTend+c).^(1-p)-(fTstart+c).^(1-p))./(1-p);
    %cumnr_model = k/(p-1)*(c^(1-p)-(c+time_as(i)).^(1-p)); % integrated form of MOL
else
    fAcp = log(fTend+c)-log(fTstart+c);
    %cumnr_model = k*log(time_as(i)/c+1); % integrated form of MOL
end
% rms = (sum((i-cumnr_model).^2)/length(i))^0.5; % RMS between observed data and MOL
% Log likelihood 
nNumEvents = length(time_as);
fL = -(nNumEvents*log(k)-p*sum(log(time_as+c))-k*fAcp);

%----------------------------------------------------------------------------------
function fL = bruteloglike_pck2(vValues,tas,fT1)
% Function to calculate the log likelihood function of an Omori law including one
% secondary aftershock at time fT1. Assume p and c constant for the entire sequence,
% but different k's before and after fT1
%
% Incoming variables:
% vValues : Starting values for p,c,k1,k2
% time_as : Vector of aftershock times from mainshock time
% fT1     : time after mainshock of the large aftershock
% 
% Outgoing:
% fL  : Log likelihood function value
%
% J. Woessner
% last update: 05.08.03

p = vValues(1);     c = vValues(2);
k1 = vValues(3);    k2 = vValues(4);

% Calculate different time lengths
vSel = (tas > fT1);

% Select the two time periods
vTperiod1 = tas(~vSel,:);
vTperiod2 = tas(vSel,:);

% Setting start end end time
fTstart = min(tas);
fTend = max(tas);

%% Calculate the likelihood function for the mainshock sequence up to the large aftershock time ft1
if p ~= 1
    fAcp = ((fT1+c).^(1-p)-(fTstart+c).^(1-p))./(1-p);
    %cumnr_model = k/(p-1)*(c^(1-p)-(c+time_as(i)).^(1-p)); % integrated form of MOL
else
    fAcp = log(fT1+c)-log(fTstart+c);
    %cumnr_model = k*log(time_as(i)/c+1); % integrated form of MOL
end
% rms = (sum((i-cumnr_model).^2)/length(i))^0.5; % RMS between observed data and MOL
% Log likelihood first period 
nEvents = length(vTperiod1);
fL_per1 = nEvents*log(k1)-p*sum(log(vTperiod1+c))-k1*fAcp;


%% Calculate the likelihood function for the sequence after the large aftershock time fT1
% Some shortcuts
fT2 = min(vTperiod2); % Staring time of events after large aftershock
fTerm1 = sum(log(k1*(vTperiod2+c).^(-p)+k2*(vTperiod2-fT1+c).^(-p)));
fpsup = 1-p;
if (p~=1)
    fTerm2a = k1/fpsup*((fTend+c).^fpsup-(fT1+c).^fpsup);
    fTerm2b = k2/fpsup*((fTend-fT1+c).^fpsup-c.^fpsup);
    fTerm2 = fTerm2a + fTerm2b;
else
    fTerm2 = k1*(log(fTend+c)-log(fT1+c))+k2*(log(fTend-fT1+c)-log(c));
end
% Log likelihood second period 
fL_per2 = fTerm1-fTerm2;

% Add upp likelihoods
fL = -(fL_per1+fL_per2);

%----------------------------------------------------------------------------------
function fL = bruteloglike_p2ck2(vValues,tas,fT1)
% Function to calculate the log likelihood function of an Omori law including one
% secondary aftershock at time fT1. c constant , p and k different before and after fT1
%
% Incoming variables:
% vValues : Starting values for p,c,k1,k2
% time_as : Vector of aftershock times from mainshock time
% fT1     : time after mainshock of the large aftershock
% 
% Outgoing:
% fL  : Log likelihood function value
%
% J. Woessner
% last update: 05.08.03

p1 = vValues(1);    p2 = vValues(2);
c1 = vValues(3);    c2 = vValues(3);
k1 = vValues(5);    k2 = vValues(6);

% Calculate different time lengths
vSel = (tas > fT1);
% Select the two time periods
vTperiod1 = tas(~vSel,:);
vTperiod2 = tas(vSel,:);

% Setting start end end time
fTstart = min(tas);
fTend = max(tas);

%% Calculate the likelihood function for the mainshock sequence up to the large aftershock time ft1
if (p1 ~= 1)
    fAcp = ((fT1+c1).^(1-p1)-(fTstart+c1).^(1-p1))./(1-p1);
else
    fAcp = log(fT1+c1)-log(fTstart+c1);
end
% Log likelihood first period 
nEvents = length(vTperiod1);
fL_per1 = nEvents*log(k1)-p1*sum(log(vTperiod1+c1))-k1*fAcp;


%% Calculate the likelihood function for the sequence after the large aftershock time fT1
% Some shortcuts
fpsup1 = 1-p1;
fpsup2 = 1-p2;

fTerm1 = sum(log(k1*(vTperiod2+c1).^(-p1)+k2*(vTperiod2-fT1+c2).^(-p2)));

if (p1~=1 && p2~=2)
    fTerm2a = k1/fpsup1*((fTend+c1).^fpsup1-(fT1+c1).^fpsup1);
    fTerm2b = k2/fpsup2*((fTend-fT1+c2).^fpsup2-c2.^fpsup2);
    fTerm2 = fTerm2a + fTerm2b;
elseif (p1==1 && p2==1)
    fTerm2 = k1*(log(fTend+c1)-log(fT1+c1))+k2*(log(fTend-fT1+c2)-log(c2));
elseif (p1~=1 && p2==1)
    fTerm2a = k1/fpsup1*((fTend+c1).^fpsup1-(fT1+c1).^fpsup1);
    fTerm2b = k2*(log(fTend-fT1+c2)-log(c2));
    fTerm2 = fTerm2a + fTerm2b;
else % (p1==1 & p2~=1)
    fTerm2a = k1*(log(fTend+c1)-log(fT1+c1));
    fTerm2b = k2/fpsup2*((fTend-fT1+c2).^fpsup2-c2.^fpsup2);
    fTerm2 = fTerm2a + fTerm2b;
end;
% Log likelihood second period 
fL_per2 = fTerm1-fTerm2;

% Add upp likelihoods
fL = -(fL_per1+fL_per2);

%----------------------------------------------------------------------------------
function fL = bruteloglike_p2c2k2(vValues,tas,fT1)
% Function to calculate the log likelihood function of an Omori law including one
% secondary aftershock at time fT1. p, c and k different before and after fT1
%
% Incoming variables:
% vValues : Starting values for p,c,k1,k2
% time_as : Vector of aftershock times from mainshock time
% fT1     : time after mainshock of the large aftershock
% 
% Outgoing:
% fL  : Log likelihood function value
%
% J. Woessner
% last update: 05.08.03

p1 = vValues(1);    p2 = vValues(2);
c1 = vValues(3);    c2 = vValues(4);
k1 = vValues(5);    k2 = vValues(6);

% Calculate different time lengths
vSel = (tas > fT1);
% Select the two time periods
vTperiod1 = tas(~vSel,:);
vTperiod2 = tas(vSel,:);

% Setting start end end time
fTstart = min(tas);
fTend = max(tas);

%% Calculate the likelihood function for the mainshock sequence up to the large aftershock time ft1
if (p1 ~= 1)
    fAcp = ((fT1+c1).^(1-p1)-(fTstart+c1).^(1-p1))./(1-p1);
else
    fAcp = log(fT1+c1)-log(fTstart+c1);
end
% Log likelihood first period 
nEvents = length(vTperiod1);
fL_per1 = nEvents*log(k1)-p1*sum(log(vTperiod1+c1))-k1*fAcp;


%% Calculate the likelihood function for the sequence after the large aftershock time fT1
% Some shortcuts
fpsup1 = 1-p1;      fpsup2 = 1-p2;

fTerm1 = sum(log(k1*(vTperiod2+c1).^(-p1)+k2*(vTperiod2-fT1+c2).^(-p2)));
if (p1~=1 && p2~=2)
    fTerm2a = k1/fpsup1*((fTend+c1).^fpsup1-(fT1+c1).^fpsup1);
    fTerm2b = k2/fpsup2*((fTend-fT1+c2).^fpsup2-c2.^fpsup2);
    fTerm2 = fTerm2a + fTerm2b;
elseif (p1==1 && p2==1)
    fTerm2 = k1*(log(fTend+c1)-log(fT1+c1))+k2*(log(fTend-fT1+c2)-log(c2));
elseif (p1~=1 && p2==1)
    fTerm2a = k1/fpsup1*((fTend+c1).^fpsup1-(fT1+c1).^fpsup1);
    fTerm2b = k2*(log(fTend-fT1+c2)-log(c2));
    fTerm2 = fTerm2a + fTerm2b;
else % (p1==1 & p2~=1)
    fTerm2a = k1*(log(fTend+c1)-log(fT1+c1));
    fTerm2b = k2/fpsup2*((fTend-fT1+c2).^fpsup2-c2.^fpsup2);
    fTerm2 = fTerm2a + fTerm2b;
end;
% Log likelihood second period 
fL_per2 = fTerm1-fTerm2;

% Add upp likelihoods
fL = -(fL_per1+fL_per2);

% --------------------------------------------------------------------
function y = nanmean(x)
% Average or mean ignoring NaNs.

if (isempty(x))
	y = NaN;    return
end

% Replace NaNs with zeros.
nans = isnan(x);
i = find(nans);
x(i) = zeros(size(i));

% count terms in sum over first non-singleton dimension
dim = find(size(x) > 1);
if (isempty(dim)),	dim = 1;
else				dim = dim(1);
end
count = sum(~nans,dim);

% Protect against a column of all NaNs
i = find(count == 0);
count(i) = 1;
y = sum(x,dim)./count;
y(i) = NaN;

% -----------------------------------------------
function y = nansum(x)
% Sum ignoring NaNs.

% Replace NaNs with zeros.
nans = isnan(x);
i = find(nans);
x(i) = zeros(size(i));

% Protect against an entire column of NaNs
y = sum(x);     i = find(all(nans));
y(i) = i + NaN;

