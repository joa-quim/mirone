%#mex
function calc_omori(events_time,events_mag,h_mir_fig);
%   This was plot_llkstest ZMAP function. This part of ZMAP is better written
%   I assembled in this file all functions needed to do the job and it is here
%   also that the input parameters are requested via the inputdlg function
%   Mironified by J. Luis at 6-07-05

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
[H,P,KSSTAT] = kstest2(cumnr_modelf,cumnrf);

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
hWaitbar1 = waitbar(0,'Bootstrapping...');
set(hWaitbar1,'Numbertitle','off','Name','Bootstap Omori parameters')
for (j = 1:bootloops)
    randnr = ceil(rand(n,1)*n);
    i = (1:n)';
    newtas(i,:) = time_as(randnr(i),:); % bootstrap sample
    newtas = sort(newtas);
    [pv1, pv2, cv1, cv2, kv1, kv2, fAIC, fL] = bruteforceloglike_a2(newtas, fT1, nMod);
    loopout = [loopout; pv1, pv2, cv1, cv2, kv1, kv2, fAIC, fL];
    waitbar(j/bootloops)
end
close(hWaitbar1)

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
function [fStdDev] = calc_StdDev(vDistribution);
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
if (~isreal(fStdDev) & ~isnan(fStdDev))    fStdDev = 0;     end

% --------------------------------------------------------------------------------------
function [H,pValue,KSstatistic] = kstest2(x1 , x2 , alpha , tail)
%KSTEST2 Two-sample Kolmogorov-Smirnov goodness-of-fit hypothesis test.
%   H = KSTEST2(X1,X2,ALPHA,TAIL) performs a Kolmogorov-Smirnov (K-S) test 
%   to determine if independent random samples, X1 and X2, are drawn from 
%   the same underlying continuous population. ALPHA and TAIL are optional
%   scalar inputs: ALPHA is the desired significance level (default = 0.05); 
%   TAIL indicates the type of test (default = 0). H indicates the result of
%   the hypothesis test:
%      H = 0 => Do not reject the null hypothesis at significance level ALPHA.
%      H = 1 => Reject the null hypothesis at significance level ALPHA.
% 
%   Let F1(x) and F2(x) be the empirical distribution functions from sample 
%   vectors X1 and X2, respectively. The 2-sample K-S test hypotheses and 
%   test statistic are:
%
%   Null Hypothesis: F1(x) = F2(x) for all x
%      For TAIL =  0 (2-sided test), alternative: F1(x) not equal to F2(x).
%      For TAIL =  1 (1-sided test), alternative: F1(x) > F2(x).
%      For TAIL = -1 (1-sided test), alternative: F1(x) < F2(x).
%
%   For TAIL = 0, 1, and -1, the test statistics are T = max|F1(x) - F2(x)|,
%   T = max[F1(x) - F2(x)], and T = max[F2(x) - F1(x)], respectively.
%
%   The decision to reject the null hypothesis occurs when the significance 
%   level, ALPHA, equals or exceeds the P-value.
%
%   X1 and X2 are row or column vectors of lengths N1 and N2, respectively, 
%   and represent random samples from some underlying distribution(s). 
%   Missing observations, indicated by NaN's (Not-a-Number), are ignored.
%
%   [H,P] = KSTEST2(...) also returns the asymptotic P-value P.
%
%   [H,P,KSSTAT] = KSTEST2(...) also returns the K-S test statistic KSSTAT
%   defined above for the test type indicated by TAIL.
%
%   The asymptotic P-value becomes very accurate for large sample sizes, and
%   is believed to be reasonably accurate for sample sizes N1 and N2 such 
%   that (N1*N2)/(N1 + N2) >= 4.

% Author(s): R.A. Baker, 08/14/98
% Copyright 1993-2002 The MathWorks, Inc. 
% $Revision: 1.5 $   $ Date: 1998/01/30 13:45:34 $

%
% References:
%   (1) Massey, F.J., "The Kolmogorov-Smirnov Test for Goodness of Fit",
%         Journal of the American Statistical Association, 46 (March 1956), 68-77.
%   (2) Miller, L.H., "Table of Percentage Points of Kolmogorov Statistics",
%         Journal of the American Statistical Association, (March 1951), 111-121.
%   (3) Conover, W.J., "Practical Nonparametric Statistics", 
%         John Wiley & Sons, Inc., 1980.
%   (4) Press, W.H., et. al., "Numerical Recipes in C", 
%         Cambridge University Press, 1992.
 
if (nargin < 2)    error(' At least 2 inputs are required.');   end

% Ensure each sample is a VECTOR.
[rows1 , columns1]  =  size(x1);
[rows2 , columns2]  =  size(x2);

if (rows1 ~= 1) & (columns1 ~= 1) 
    error(' Sample ''X1'' must be a vector.');
end

if (rows2 ~= 1) & (columns2 ~= 1) 
    error(' Sample ''X2'' must be a vector.');
end

% Remove missing observations indicated by NaN's, and 
% ensure that valid observations remain.
x1  =  x1(~isnan(x1));
x2  =  x2(~isnan(x2));
x1  =  x1(:);
x2  =  x2(:);

if (isempty(x1))   error(' Sample vector ''X1'' is composed of all NaN''s.');   end
if (isempty(x2))   error(' Sample vector ''X2'' is composed of all NaN''s.');   end

% Ensure the significance level, ALPHA, is a scalar 
% between 0 and 1 and set default if necessary.
if (nargin >= 3) & ~isempty(alpha)
   if prod(size(alpha)) > 1
      error(' Significance level ''Alpha'' must be a scalar.');
   end
   if (alpha <= 0 | alpha >= 1)
      error(' Significance level ''Alpha'' must be between 0 and 1.'); 
   end
else
   alpha  =  0.05;
end

% Ensure the type-of-test indicator, TAIL, is a scalar integer from 
% the allowable set {-1 , 0 , 1}, and set default if necessary.
if (nargin >= 4) & ~isempty(tail)
   if prod(size(tail)) > 1
      error(' Type-of-test indicator ''Tail'' must be a scalar.');
   end
   if (tail ~= -1) & (tail ~= 0) & (tail ~= 1)
      error(' Type-of-test indicator ''Tail'' must be -1, 0, or 1.');
   end
else
   tail  =  0;
end

% Calculate F1(x) and F2(x), the empirical (i.e., sample) CDFs.
binEdges   = [-inf ; sort([x1;x2]) ; inf];

binCounts1 = histc (x1, binEdges);
binCounts2 = histc (x2, binEdges);

sumCounts1 = cumsum(binCounts1)./sum(binCounts1);
sumCounts2 = cumsum(binCounts2)./sum(binCounts2);

sampleCDF1 = sumCounts1(1:end-1);
sampleCDF2 = sumCounts2(1:end-1);

% Compute the test statistic of interest.
switch tail
   case  0      %  2-sided test: T = max|F1(x) - F2(x)|.
      deltaCDF = abs(sampleCDF1 - sampleCDF2);
   case -1      %  1-sided test: T = max[F2(x) - F1(x)].
      deltaCDF = sampleCDF2 - sampleCDF1;
   case  1      %  1-sided test: T = max[F1(x) - F2(x)].
      deltaCDF = sampleCDF1 - sampleCDF2;
end

KSstatistic = max(deltaCDF);

% Compute the asymptotic P-value approximation and accept or
% reject the null hypothesis on the basis of the P-value.
n1     = length(x1);
n2     = length(x2);
n      = n1 * n2 /(n1 + n2);
lambda = max((sqrt(n) + 0.12 + 0.11/sqrt(n)) * KSstatistic , 0);

if (tail ~= 0)        % 1-sided test.
    pValue = exp(-2 * lambda * lambda);
else                % 2-sided test (default).
    %  Use the asymptotic Q-function to approximate the 2-sided P-value.
   j = [1:101]';
   pValue  =  2 * sum((-1).^(j-1).*exp(-2*lambda*lambda*j.^2));
   if (pValue < 0)  pValue = 0; end
   if (pValue > 1)  pValue = 1; end
end

H = (alpha >= pValue);

% --------------------------------------------------------------------------------------
function [pv1, pv2, cv1, cv2, kv1, kv2, fAIC, fL] = bruteforceloglike_a2(tas, fT1, nMod);
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
    [vValues, fL] = fmincon_m(@bruteloglike, vStartValues, [], [], [], [],...
        [fPmin fCmin fKmin1 ], [fPmax fCmax fKmax1 ], [], options, tas);
    pv1 = vValues(1);    pv2 = vValues(1);    cv1 = vValues(2);
    cv2 = vValues(2);    kv1 = vValues(3);    kv2 = vValues(3);
elseif (nMod == 2)      % 4 free parameters: p, c , k1, k2
    fPar = 4;
    vStartValues = [1.1 0.5 50 50];
    %[vValues, fL] = fmincon_m('bruteloglike_pck2', vStartValues, [], [], [], [],...
    [vValues, fL] = fmincon_m(@bruteloglike_pck2, vStartValues, [], [], [], [],...
        [fPmin fCmin fKmin1 fKmin2], [fPmax fCmax fKmax1 fKmax2], [], options, tas, fT1);
    pv1 = vValues(1);    pv2 = vValues(1);    cv1 = vValues(2);
    cv2 = vValues(2);    kv1 = vValues(3);    kv2 = vValues(4);
elseif (nMod == 3)      % 5 free parameters: p1,p2,c,k1,k2
    fPar = 5;
    vStartValues = [1.1 1.1 0.5 0.5 50 50];
    %[vValues, fL] = fmincon_m('bruteloglike_p2ck2', vStartValues, [], [], [], [],...
    [vValues, fL] = fmincon_m(@bruteloglike_p2ck2, vStartValues, [], [], [], [],...
        [fPmin fPmin fCmin fCmin fKmin1 fKmin2], [fPmax fPmax fCmax fCmax fKmax1 fKmax2], [], options, tas, fT1);
    pv1 = vValues(1);    pv2 = vValues(2);    cv1 = vValues(3);
    cv2 = vValues(3);    kv1 = vValues(5);    kv2 = vValues(6);
else                    % 6 free parameters: p1,p2,c1, c2,k1,k2
    fPar = 6;
    vStartValues = [1.1 1.1 0.5 0.5 50 50];
    %[vValues, fL] = fmincon_m('bruteloglike_p2c2k2', vStartValues, [], [], [], [],...
    [vValues, fL] = fmincon_m(@bruteloglike_p2c2k2, vStartValues, [], [], [], [],...
        [fPmin fPmin fCmin fCmin fKmin1 fKmin2], [fPmax fPmax fCmax fCmax fKmax1 fKmax2], [], options, tas, fT1);
    pv1 = vValues(1);    pv2 = vValues(2);    cv1 = vValues(3);
    cv2 = vValues(4);    kv1 = vValues(5);    kv2 = vValues(6);
end  
% corrected Akaike Information Criterion
[fk,nX]=size(tas);
fAIC = -2*(-fL)+2*fPar+2*fPar*(fPar+1)/(fk-fPar-1);

%----------------------------------------------------------------------------------
function fL = bruteloglike(vValues,time_as);
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
function fL = bruteloglike_pck2(vValues,tas,fT1);
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
function fL = bruteloglike_p2ck2(vValues,tas,fT1);
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

if (p1~=1 & p2~=2)
    fTerm2a = k1/fpsup1*((fTend+c1).^fpsup1-(fT1+c1).^fpsup1);
    fTerm2b = k2/fpsup2*((fTend-fT1+c2).^fpsup2-c2.^fpsup2);
    fTerm2 = fTerm2a + fTerm2b;
elseif (p1==1 & p2==1)
    fTerm2 = k1*(log(fTend+c1)-log(fT1+c1))+k2*(log(fTend-fT1+c2)-log(c2));
elseif (p1~=1 & p2==1)
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
function fL = bruteloglike_p2c2k2(vValues,tas,fT1);
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
if (p1~=1 & p2~=2)
    fTerm2a = k1/fpsup1*((fTend+c1).^fpsup1-(fT1+c1).^fpsup1);
    fTerm2b = k2/fpsup2*((fTend-fT1+c2).^fpsup2-c2.^fpsup2);
    fTerm2 = fTerm2a + fTerm2b;
elseif (p1==1 & p2==1)
    fTerm2 = k1*(log(fTend+c1)-log(fT1+c1))+k2*(log(fTend-fT1+c2)-log(c2));
elseif (p1~=1 & p2==1)
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

% -------------------------------------------------------------------------------------------------
function [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon_m(FUN,X,A,B,Aeq,Beq,LB,UB,NONLCON,options,varargin)
%FMINCON Finds a constrained minimum of a function of several variables.
%   FMINCON solves problems of the form:
%       min F(X)  subject to:  A*X  <= B, Aeq*X  = Beq (linear constraints)
%        X                       C(X) <= 0, Ceq(X) = 0   (nonlinear constraints)
%                                LB <= X <= UB            
%
%   Examples
%     FUN can be specified using @:
%        X = fmincon(@humps,...)
%     In this case, F = humps(X) returns the scalar function value F of the HUMPS function
%     evaluated at X.
%

%   Copyright 1990-2003 The MathWorks, Inc. 

defaultopt = struct('Display','final','LargeScale','on', ...
   'TolX',1e-6,'TolFun',1e-6,'TolCon',1e-6,'DerivativeCheck','off',...
   'Diagnostics','off',...
   'GradObj','off','GradConstr','off',...
   'HessMult',[],...% HessMult [] by default
   'Hessian','off','HessPattern','sparse(ones(numberOfVariables))',...
   'MaxFunEvals','100*numberOfVariables',...
   'MaxSQPIter',Inf,...
   'DiffMaxChange',1e-1,'DiffMinChange',1e-8,...
   'PrecondBandWidth',0,'TypicalX','ones(numberOfVariables,1)',...
   'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
   'TolPCG',0.1,'MaxIter',400,'OutputFcn',[]);

large = 'large-scale';
medium = 'medium-scale';

if nargin < 10, options=[];
   if nargin < 9, NONLCON=[];
      if nargin < 8, UB = [];
         if nargin < 7, LB = [];
            if nargin < 6, Beq=[];
               if nargin < 5, Aeq =[];
               end, end, end, end, end, end

lenVarIn = length(varargin);
XOUT=X(:);
numberOfVariables=length(XOUT);

verbosity = 0;

% Set to column vectors
B = B(:);
Beq = Beq(:);
l = LB(:);
u = UB(:);
lFinite = l(~isinf(l));
uFinite = u(~isinf(u));

meritFunctionType = 0;

gradflag = 0;
hessflag = 0;
gradconstflag = 0;

% Convert to inline function as needed
if ~isempty(FUN)  % will detect empty string, empty matrix, empty cell array
   [funfcn, msg] = optimfcnchk1(FUN,'fmincon',length(varargin),gradflag,hessflag);
end

confcn{1} = '';      % ESTA
[rowAeq,colAeq]=size(Aeq);
hessflag = 0;
OUTPUT.algorithm = medium;       % ESTA

lenvlb=length(l);
lenvub=length(u);

% Ensure starting point lies within bounds
i=1:lenvlb;
lindex = XOUT(i)<l(i);
if any(lindex)
    XOUT(lindex)=l(lindex)+1e-4; 
end
i=1:lenvub;
uindex = XOUT(i)>u(i);
if any(uindex)
    XOUT(uindex)=u(uindex);
end
X(:) = XOUT;

% Evaluate function
GRAD=zeros(numberOfVariables,1);
HESS = [];

try
  f = feval(funfcn{3},X,varargin{:});
catch
  errmsg = sprintf('%s\n%s\n\n%s',...
     'FMINCON cannot continue because user supplied objective function', ...
     ' failed with the following error:', lasterr);
  error(errmsg);
end

% Evaluate constraints
c=[]; ceq =[];
cGRAD = zeros(numberOfVariables,length(c));
ceqGRAD = zeros(numberOfVariables,length(ceq));

% call algorithm
[X,FVAL,lambda,EXITFLAG,OUTPUT,GRAD,HESSIAN]=...
  nlconst1(funfcn,X,l,u,full(A),B,full(Aeq),Beq,confcn,options,defaultopt, ...
  verbosity,gradflag,gradconstflag,hessflag,meritFunctionType,...
  f,GRAD,HESS,c,ceq,cGRAD,ceqGRAD,varargin{:});
LAMBDA=lambda;

% ----------------------------------------------------------------------------
function [x,FVAL,lambda_out,EXITFLAG,OUTPUT,GRADIENT,HESS]= ...
    nlconst1(funfcn,x,lb,ub,Ain,Bin,Aeq,Beq,confcn,OPTIONS,defaultopt,...
    verbosity,gradflag,gradconstflag,hessflag,meritFunctionType,...
    fval,gval,Hval,ncineqval,nceqval,gncval,gnceqval,varargin)
%NLCONST Helper function to find the constrained minimum of a function 
%   of several variables. Called by FMINCON, FGOALATTAIN FSEMINF and FMINIMAX.
%
%   Copyright 1990-2003 The MathWorks, Inc.
%   $Revision: 1.41 $  $Date: 2003/02/10 20:48:36 $

%   Calls OPTEVAL.
%
%   meritFunctionType==5 for fseminf
%                    ==1 for fminimax & fgoalattain (can use 0, but 1 is default)
%                    ==0 for fmincon

% Initialize some parameters
FVAL=[]; lambda_out=[]; OUTPUT=[]; lambdaNLP = []; GRADIENT = []; 

% Handle the output
outputfcn = defaultopt.OutputFcn;

iter = 0;
XOUT = x(:);
% numberOfVariables must be the name of this variable
numberOfVariables = length(XOUT);
SD = ones(numberOfVariables,1); 
Nlconst = 'nlconst';
bestf = Inf; 
if isempty(confcn{1})
    constflag = 0;
else
    constflag = 1;
end
stepsize = 1;
HESS=eye(numberOfVariables,numberOfVariables); % initial Hessian approximation.
done = false; 
EXITFLAG = 1;

% Get options
tolX = optimgetfast(OPTIONS,'TolX',defaultopt);
tolFun = optimgetfast(OPTIONS,'TolFun',defaultopt);
tolCon = optimgetfast(OPTIONS,'TolCon',defaultopt);
DiffMinChange = optimgetfast(OPTIONS,'DiffMinChange',defaultopt);
DiffMaxChange = optimgetfast(OPTIONS,'DiffMaxChange',defaultopt);
if DiffMinChange >= DiffMaxChange
    errmsg =sprintf('%s %0.5g %s %0.5g%s\n%s\n',...
    'DiffMinChange options parameter is',DiffMinChange,'and DiffMaxChange is',...
    DiffMaxChange,'.','DiffMinChange must be strictly less than DiffMaxChange.');
    error(errmsg)
end
DerivativeCheck = strcmp(optimgetfast(OPTIONS,'DerivativeCheck',defaultopt),'on');
typicalx = optimgetfast(OPTIONS,'TypicalX',defaultopt) ;
if ischar(typicalx)
   if isequal(lower(typicalx),'ones(numberofvariables,1)')
      typicalx = ones(numberOfVariables,1);
   else
      error('Option ''TypicalX'' must be a numeric value if not the default.')
   end
end
maxFunEvals = optimgetfast(OPTIONS,'MaxFunEvals',defaultopt);
maxIter = optimgetfast(OPTIONS,'MaxIter',defaultopt);
% In case the defaults were gathered from calling: optimset('fmincon'):
if ischar(maxFunEvals)
    if isequal(lower(maxFunEvals),'100*numberofvariables')
        maxFunEvals = 100*numberOfVariables;
    else
        error('Option ''MaxFunEvals'' must be an integer value if not the default.')
    end
end

% For SQP problem create default sqpopt to pass to qpsub. 
% (Values not needed until we fix default MaxSQPIter in R14.)
sqpopt = []; 

% Handle bounds as linear constraints
arglb = ~isinf(lb);
lenlb=length(lb); % maybe less than numberOfVariables due to old code
argub = ~isinf(ub);
lenub=length(ub);
boundmatrix = eye(max(lenub,lenlb),numberOfVariables);
if nnz(arglb) > 0     
    lbmatrix = -boundmatrix(arglb,1:numberOfVariables);% select non-Inf bounds 
    lbrhs = -lb(arglb);
else
    lbmatrix = []; lbrhs = [];
end
if nnz(argub) > 0
    ubmatrix = boundmatrix(argub,1:numberOfVariables);
    ubrhs=ub(argub);
else
    ubmatrix = []; ubrhs=[];
end 

% Update constraint matrix and right hand side vector with bound constraints.
A = [lbmatrix;ubmatrix;Ain];
B = [lbrhs;ubrhs;Bin];
if isempty(A)
    A = zeros(0,numberOfVariables); B=zeros(0,1);
end
if isempty(Aeq)
    Aeq = zeros(0,numberOfVariables); Beq=zeros(0,1);
end

% Used for semi-infinite optimization:
s = nan; POINT =[]; NEWLAMBDA =[]; LAMBDA = []; NPOINT =[]; FLAG = 2;
OLDLAMBDA = [];

x(:) = XOUT;  % Set x to have user expected size
% Compute the objective function and constraints
f = fval;
nceq = nceqval; ncineq = ncineqval;  % nonlinear constraints only
nc = [nceq; ncineq];
c = [ Aeq*XOUT-Beq; nceq; A*XOUT-B; ncineq];

% Get information on the number and type of constraints.
non_eq = length(nceq);
non_ineq = length(ncineq);
[lin_eq,Aeqcol] = size(Aeq);
[lin_ineq,Acol] = size(A);  % includes upper and lower bounds
eq = non_eq + lin_eq;
ineq = non_ineq + lin_ineq;
ncstr = ineq + eq;
% Boolean inequalitiesExist = true if and only if there exist either
% finite bounds or linear inequalities or nonlinear inequalities. 
% Used only for printing indices of active inequalities at the solution
inequalitiesExist = any(arglb) || any(argub) || size(Ain,1) > 0 || non_ineq > 0;

% Compute the initial constraint violation.
ga=[abs(c( (1:eq)' )) ; c( (eq+1:ncstr)' ) ];
if (~isempty(c))    mg=max(ga);
else                mg = 0;
end

if (isempty(f))
    error('FUN must return a non-empty objective function.')
end

OLDX=XOUT;
OLDC=c; OLDNC=nc;
OLDgf=zeros(numberOfVariables,1);
gf=zeros(numberOfVariables,1);
OLDAN=zeros(ncstr,numberOfVariables);
LAMBDA=zeros(ncstr,1);
lambdaNLP = zeros(ncstr,1);
numFunEvals=1;
numGradEvals=1;
% Just for not having anoying warnings of false errors from the compiler
OLDF = [];  MATX = [];  ACTIND = [];    infeasIllPosedMaxSQPIter = [];

how = ''; 
optimError = []; % In case we have convergence in 0th iteration, this needs a value.
%---------------------------------Main Loop-----------------------------
while ~done 
   %----------------GRADIENTS----------------
   
   if ~gradconstflag || ~gradflag || DerivativeCheck
      % Finite Difference gradients (even if just checking analytical)
      POINT = NPOINT; 
      len_nc = length(nc);
      ncstr =  lin_eq + lin_ineq + len_nc;     
      FLAG = 0; % For semi-infinite

      % Compute finite difference gradients
      %
      if (DerivativeCheck || ~gradconstflag)            % ESTA
         [gf,gnc,NEWLAMBDA,OLDLAMBDA,s]=finiteDifferences(XOUT,x,funfcn,confcn,lb,ub,f,nc,...
                               [],DiffMinChange,DiffMaxChange,typicalx,[],'all',...
                               LAMBDA,NEWLAMBDA,OLDLAMBDA,POINT,FLAG,s,varargin{:});
         gnc = gnc'; % nlconst requires the transpose of the Jacobian
      else
         % no derivative check and user-supplied constraint gradients.             
         % Estimate objective gradient only.
         gf=finiteDifferences(XOUT,x,funfcn,confcn,lb,ub,f,[],...   
                         DiffMinChange,DiffMaxChange,typicalx,[],...
                         [],[],[],[],[],[],varargin{:});
      end
      
      FLAG = 1; % For semi-infinite
      numFunEvals = numFunEvals + numberOfVariables;
   end  
   
   % Now add in Aeq, and A
   if ~isempty(gnc)
      gc = [Aeq', gnc(:,1:non_eq), A', gnc(:,non_eq+1:non_ineq+non_eq)];
   elseif ~isempty(Aeq) || ~isempty(A)
      gc = [Aeq',A'];
   else
      gc = zeros(numberOfVariables,0);
   end
   AN=gc';
   
   % Iteration 0 is handled separately below
   if (iter > 0) % All but 0th iteration ----------------------------------------
       % Compute the first order KKT conditions.
       if strcmp(funfcn{2},'fseminf'), lambdaNLP = NEWLAMBDA; end
       normgradLag = norm(gf + AN'*lambdaNLP,inf);
       normcomp = norm(lambdaNLP(eq+1:ncstr).*c(eq+1:ncstr),inf);
       if isfinite(normgradLag) && isfinite(normcomp)
           optimError = max(normgradLag, normcomp);
       else
           optimError = inf;
       end
       feasError  = mg;
       optimScal = 1; feasScal = 1; 
              
       %-------------TEST CONVERGENCE---------------
       
       if (optimError < tolFun*optimScal && feasError < tolCon*feasScal)
           EXITFLAG = 1;
           done = true;

           if inequalitiesExist
              % Report active inequalities
              [activeLb,activeUb,activeIneqLin,activeIneqNonlin] = ...
                  activeInequalities(c,tolCon,arglb,argub,lin_eq,non_eq,size(Ain));           
           end   
       elseif (max(abs(SD)) < 2*tolX || abs(gf'*SD) < 2*tolFun ) && (mg < tolCon || infeasIllPosedMaxSQPIter) 
           % The algorithm can make no more progress.  If feasible, compute 
           % the new up-to-date Lagrange multipliers (with new gradients) 
           % and recompute the KKT error.  Then output appropriate termination
           % message.
           if (mg < tolCon)
               lambdaNLP(:,1) = 0;
               [Q,R] = qr(AN(ACTIND,:)');
               ws = warning('off');
               lambdaNLP(ACTIND) = -R\Q'*gf;
               warning(ws);
               lambdaNLP(eq+1:ncstr) = max(0,lambdaNLP(eq+1:ncstr));
               if strcmp(funfcn{2},'fseminf'), lambdaNLP = NEWLAMBDA; end
               normgradLag = norm(gf + AN'*lambdaNLP,inf);
               normcomp = norm(lambdaNLP(eq+1:ncstr).*c(eq+1:ncstr),inf);
               if isfinite(normgradLag) && isfinite(normcomp)
                   optimError = max(normgradLag, normcomp);
               else
                   optimError = inf;
               end
               optimScal = 1;
               if optimError < tolFun*optimScal
                   EXITFLAG = 1;
               elseif max(abs(SD)) < 2*tolX
                   EXITFLAG = 1;
               else 
                   EXITFLAG = 1;
               end 
               
               if inequalitiesExist
                  % Report active inequalities
                  [activeLb,activeUb,activeIneqLin,activeIneqNonlin] = ...
                      activeInequalities(c,tolCon,arglb,argub,lin_eq,non_eq,size(Ain));  
               end
           else                         % if mg < tolCon
               EXITFLAG = -1;
           end                          % if mg < tolCon
           done = true;
       else % continue
           % NEED=[LAMBDA>0] | G>0
           if numFunEvals > maxFunEvals
               XOUT = MATX;
               f = OLDF;
               EXITFLAG = 0;
               done = true;
           end
           if iter > maxIter
               XOUT = MATX;
               f = OLDF;
               EXITFLAG = 0;
               done = true;
           end
       end 
   else % ------------------------0th Iteration----------------------------------
       
   end % if iter > 0
   
   % Continue if termination criteria do not hold or it is the 0th iteration-------------------------------------------
   if ~done 
      how=''; 
      iter = iter + 1;

      %-------------SEARCH DIRECTION---------------
      % For equality constraints make gradient face in 
      % opposite direction to function gradient.
      for i=1:eq 
         schg=AN(i,:)*gf;
         if schg>0
            AN(i,:)=-AN(i,:);
            c(i)=-c(i);
         end
      end
   
      if numGradEvals>1  % Check for first call    
         if meritFunctionType~=5,   
            NEWLAMBDA=LAMBDA; 
         end
         [ma,na] = size(AN);
         GNEW=gf+AN'*NEWLAMBDA;
         GOLD=OLDgf+OLDAN'*LAMBDA;
         YL=GNEW-GOLD;
         sdiff=XOUT-OLDX;

         % Make sure Hessian is positive definite in update.
         if YL'*sdiff<stepsize^2*1e-3
            while YL'*sdiff<-1e-5
               [YMAX,YIND]=min(YL.*sdiff);
               YL(YIND)=YL(YIND)/2;
            end
            if YL'*sdiff < (eps*norm(HESS,'fro'));
               how=' Hessian modified twice';
               FACTOR=AN'*c - OLDAN'*OLDC;
               FACTOR=FACTOR.*(sdiff.*FACTOR>0).*(YL.*sdiff<=eps);
               WT=1e-2;
               if max(abs(FACTOR))==0; FACTOR=1e-5*sign(sdiff); end
               while YL'*sdiff < (eps*norm(HESS,'fro')) && WT < 1/eps
                  YL=YL+WT*FACTOR;
                  WT=WT*2;
               end
            else
               how=' Hessian modified';
            end
         end
         
         %----------Perform BFGS Update If YL'S Is Positive---------
         if YL'*sdiff>eps
             HESS=HESS +(YL*YL')/(YL'*sdiff)-((HESS*sdiff)*(sdiff'*HESS'))/(sdiff'*HESS*sdiff);
             % BFGS Update using Cholesky factorization  of Gill, Murray and Wright.
             % In practice this was less robust than above method and slower. 
             %   R=chol(HESS); 
             %   s2=R*S; y=R'\YL; 
             %   W=eye(numberOfVariables,numberOfVariables)-(s2'*s2)\(s2*s2') + (y'*s2)\(y*y');
             %   HESS=R'*W*R;
         else
            how=' Hessian not updated';
         end
      else % First call
         OLDLAMBDA=repmat(eps+gf'*gf,ncstr,1)./(sum(AN'.*AN')'+eps);
         ACTIND = 1:eq;     
      end % if numGradEvals>1
      numGradEvals=numGradEvals+1;
   
      LOLD=LAMBDA;
      OLDAN=AN;
      OLDgf=gf;
      OLDC=c;
      OLDF=f;
      OLDX=XOUT;
      XN=zeros(numberOfVariables,1);
      GT =c;
      HESS = (HESS + HESS')*0.5;
   
      [SD,lambda,exitflagqp,outputqp,howqp,ACTIND] ...
         = qpsub(HESS,gf,AN,-GT,[],[],XN,eq,-1,Nlconst,size(AN,1),numberOfVariables,OPTIONS,sqpopt,ACTIND);
    
      lambdaNLP(:,1) = 0;
      lambdaNLP(ACTIND) = lambda(ACTIND);
      lambda((1:eq)') = abs(lambda( (1:eq)' ));
      ga=[abs(c( (1:eq)' )) ; c( (eq+1:ncstr)' ) ];
      if (~isempty(c))  mg = max(ga);
      else              mg = 0;
      end

      if strncmp(howqp,'ok',2); 
          howqp =''; 
      end
      if ~isempty(how) && ~isempty(howqp) 
          how = [how,'; '];
      end

      LAMBDA=lambda((1:ncstr)');
      OLDLAMBDA=max([LAMBDA';0.5*(LAMBDA+OLDLAMBDA)'])' ;

      %---------------LINESEARCH--------------------
      MATX=XOUT;
      MATL = f+sum(OLDLAMBDA.*(ga>0).*ga) + 1e-30;

      infeasIllPosedMaxSQPIter = strcmp(howqp,'infeasible') || strcmp(howqp,'ill posed') || strcmp(howqp,'MaxSQPIter');
     % This merit function looks for improvement in either the constraint
     % or the objective function unless the sub-problem is infeasible in which
     % case only a reduction in the maximum constraint is tolerated.
     % This less "stringent" merit function has produced faster convergence in
     % a large number of problems.
     if mg > 0
        MATL2 = mg;
     elseif f >=0 
        MATL2 = -1/(f+1);
     else 
        MATL2 = 0;
     end
     if ~infeasIllPosedMaxSQPIter && f < 0
        MATL2 = MATL2 + f - 1;
     end
      if mg < eps && f < bestf
         bestf = f;
         bestx = XOUT;
         bestHess = HESS;
         bestgrad = gf;
         bestlambda = lambda;
         bestmg = mg;
         bestOptimError = optimError;
      end
      MERIT = MATL + 1;
      MERIT2 = MATL2 + 1; 
      stepsize=2;
      while  ((MERIT2 > MATL2) && (MERIT > MATL) && numFunEvals < maxFunEvals)
         stepsize=stepsize/2;
         if stepsize < 1e-4,  
            stepsize = -stepsize;          
         end
         XOUT = MATX + stepsize*SD;
         x(:)=XOUT; 
      
         f = feval(funfcn{3},x,varargin{:});
         if constflag
            [nctmp,nceqtmp] = feval(confcn{3},x,varargin{:});
            nctmp = nctmp(:); nceqtmp = nceqtmp(:);
         else
            nctmp = []; nceqtmp=[];
         end
            
         nc = [nceqtmp(:); nctmp(:)];
         c = [Aeq*XOUT-Beq; nceqtmp(:); A*XOUT-B; nctmp(:)];  

         numFunEvals = numFunEvals + 1;
         ga=[abs(c( (1:eq)' )) ; c( (eq+1:length(c))' )];
         if ~isempty(c)
            mg=max(ga);
         else
            mg = 0;
         end

         MERIT = f+sum(OLDLAMBDA.*(ga>0).*ga);
        if mg > 0
           MERIT2 = mg;
        elseif f >=0 
           MERIT2 = -1/(f+1);
        else 
           MERIT2 = 0;
        end
        if ~infeasIllPosedMaxSQPIter && f < 0
           MERIT2 = MERIT2 + f - 1;
        end
                                                                                                                                                                                                                            end  % line search loop
      %------------Finished Line Search-------------
   
      mf=abs(stepsize);
      LAMBDA=mf*LAMBDA+(1-mf)*LOLD;

      x(:) = XOUT;
      numFunEvals=numFunEvals+1;
   
      % Evaluate constraint gradients
      switch confcn{1}
      case 'fun'
         gnceq=[]; gncineq=[];
      case ''
         nctmp=[]; nceqtmp =[];
         gncineq = zeros(numberOfVariables,length(nctmp));
         gnceq = zeros(numberOfVariables,length(nceqtmp));
      otherwise
         error('Undefined calltype in FMINCON');
      end
      gnc_user = [gnceq, gncineq];
      gc = [Aeq', gnceq, A', gncineq];
   
   end % if ~done   
end % while ~done


% Update 
numConstrEvals = numGradEvals;

% Gradient is in the variable gf
GRADIENT = gf;

% If a better solution was found earlier, use it:
if (f > bestf )
   XOUT = bestx;
   f = bestf;
   HESS = bestHess;
   GRADIENT = bestgrad;
   lambda = bestlambda;
   mg = bestmg;
   gf = bestgrad;
   optimError = bestOptimError;
end

FVAL = f;
x(:) = XOUT;

OUTPUT.iterations = iter;
OUTPUT.funcCount = numFunEvals;
OUTPUT.stepsize = stepsize;
OUTPUT.algorithm = 'medium-scale: SQP, Quasi-Newton, line-search';
OUTPUT.firstorderopt = optimError;
OUTPUT.cgiterations = [];

[lin_ineq,Acol] = size(Ain);  % excludes upper and lower

lambda_out.lower=zeros(lenlb,1);
lambda_out.upper=zeros(lenub,1);

lambda_out.eqlin = lambdaNLP(1:lin_eq);
ii = lin_eq ;
lambda_out.eqnonlin = lambdaNLP(ii+1: ii+ non_eq);
ii = ii+non_eq;
lambda_out.lower(arglb) = lambdaNLP(ii+1 :ii+nnz(arglb));
ii = ii + nnz(arglb) ;
lambda_out.upper(argub) = lambdaNLP(ii+1 :ii+nnz(argub));
ii = ii + nnz(argub);
lambda_out.ineqlin = lambdaNLP(ii+1: ii + lin_ineq);
ii = ii + lin_ineq ;
lambda_out.ineqnonlin = lambdaNLP(ii+1 : end);

%--------------------------------------------------------------------------
function [activeLb,activeUb,activeIneqLin,activeIneqNonlin] = ...
    activeInequalities(c,tol,arglb,argub,linEq,nonlinEq,linIneq)
% ACTIVEINEQUALITIES returns the indices of the active inequalities
% and bounds.
% INPUT:
% c                 vector of constraints and bounds (see nlconst main code)
% tol               tolerance to determine when an inequality is active
% arglb, argub      boolean vectors indicating finite bounds (see nlconst
%                   main code)
% linEq             number of linear equalities
% nonlinEq          number of nonlinear equalities
% linIneq           number of linear inequalities
%
% OUTPUT
% activeLB          indices of active lower bounds
% activeUb          indices of active upper bounds  
% activeIneqLin     indices of active linear inequalities
% activeIneqNonlin  indices of active nonlinear inequalities
%

% We check wether a constraint is active or not using '< tol'
% instead of '<= tol' to be onsistent with nlconst main code, 
% where feasibility is checked using '<'.
finiteLb = nnz(arglb);      % number of finite lower bounds
finiteUb = nnz(argub);      % number of finite upper bounds

indexFiniteLb = find(arglb); % indices of variables with LB
indexFiniteUb = find(argub); % indices of variables with UB

% lower bounds
i = linEq + nonlinEq; % skip equalities

% Boolean vector that indicates which among the finite bounds is active
activeFiniteLb = abs(c(i + 1 : i + finiteLb)) < tol;

% indices of the finite bounds that are active
activeLb = indexFiniteLb(activeFiniteLb);

% upper bounds
i = i + finiteLb;

% Boolean vector that indicates which among the finite bounds is active
activeFiniteUb = abs(c(i + 1 : i + finiteUb)) < tol;

% indices of the finite bounds that are active
activeUb = indexFiniteUb(activeFiniteUb);

% linear inequalities
i = i + finiteUb;
activeIneqLin = find(abs(c(i + 1 : i + linIneq)) < tol); 
% nonlinear inequalities
i = i + linIneq;
activeIneqNonlin = find(abs(c(i + 1 : end)) < tol);   

% --------------------------------------------------------------------------------
function [gradf,cJac,NEWLAMBDA,OLDLAMBDA,s] = finiteDifferences(xCurrent,...
                  xOriginalShape,funfcn,confcn,lb,ub,fCurrent,cCurrent,...
                  YDATA,DiffMinChange,DiffMaxChange,typicalx,finDiffType,...
                  variables,LAMBDA,NEWLAMBDA,OLDLAMBDA,POINT,FLAG,s,varargin)
%
%  [gradf,cJac,NEWLAMBDA,OLDLAMBDA,s] = FINITEDIFFERENCES(xCurrent,...
%                  xOriginalShape,funfcn,confcn,lb,ub,fCurrent,cCurrent,...
%                  YDATA,DiffMinChange,DiffMaxChange,typicalx,finDiffType,...
%                  variables,LAMBDA,NEWLAMBDA,OLDLAMBDA,POINT,FLAG,s,varargin)
%
% computes the finite-difference gradients of the objective and
% constraint functions.
%
%  gradf = FINITEDIFFERENCES(xCurrent,xOriginalShape,funfcn,[],lb,ub,fCurrent,...
%                  [],YDATA,DiffMinChange,DiffMaxChange,typicalx,finDiffType,...
%                  variables,[],[],[],[],[],[],varargin)
%
% computes the finite-difference gradients of the objective function.
%
%
% INPUT:
% xCurrent              Point where gradient is desired
% xOriginalShape        Shape of the vector of variables supplied by the user
%                       (The value of xOriginalShape is NOT used)
% funfcn, confcn        Cell arrays containing info about objective and 
%                       constraints
% lb, ub                Lower and upper bounds
% fCurrent, cCurrent    Values at xCurrent of the function and the constraints 
%                       to be differentiated. Note that fCurrent can be a scalar 
%                       or a vector. 
% DiffMinChange, 
% DiffMaxChange         Minimum and maximum values of perturbation of xCurrent 
% finDiffType           Type of finite difference desired (only forward 
%                       differences implemented so far)
% variables             Variables w.r.t which we want to differentiate. Possible 
%
% LAMBDA,NEWLAMBDA,
% OLDLAMBDA,POINT,
% FLAG,s                Parameters for semi-infinite constraints
% varargin              Problem-dependent parameters passed to the objective and 
%                       constraint functions
%
% OUTPUT:
% gradf                 If fCurrent is a scalar, gradf is the finite-difference 
%                       gradient of the objective; if fCurrent is a vector,
%                       gradf is the finite-difference Jacobian  
% cJac                  Finite-difference Jacobian of the constraints
% NEWLAMBDA,
% OLDLAMBDA,s           Parameters for semi-infinite constraints
%
%   Copyright 1990-2003 The MathWorks, Inc.

[fNumberOfRows,fNumberOfColumns] = size(fCurrent);

% Make sure that fCurrent is either a scalar or a column vector
% (to ensure that the given fCurrent and the computed fplus 
% will have the same shape)
if fNumberOfColumns ~= 1
   fCurrent(:) = fCurrent;
else
   functionIsScalar = (fNumberOfRows == 1);
end

numberOfVariables = length(xCurrent); 

% nonEmptyLowerBounds = true if lb is not empty, false if it's empty;
% analogoulsy for nonEmptyUpperBound
nonEmptyLowerBounds = ~isempty(lb);
nonEmptyUpperBounds = ~isempty(ub);

% Make sure xCurrent and typicalx are column vectors so that the 
% operation max(abs(xCurrent),abs(typicalx)) won't error
xCurrent = xCurrent(:); typicalx = typicalx(:);
% Value of stepsize suggested in Trust Region Methods, Conn-Gould-Toint, section 8.4.3
CHG = sqrt(eps)*sign(xCurrent).*max(abs(xCurrent),abs(typicalx));
%
% Make sure step size lies within DiffminChange and DiffMaxChange
%
CHG = sign(CHG+eps).*min(max(abs(CHG),DiffMinChange),DiffMaxChange);
len_cCurrent = length(cCurrent);
cJac = zeros(len_cCurrent,numberOfVariables);  % For semi-infinite

if nargout < 3
   NEWLAMBDA=[]; OLDLAMBDA=[]; s=[];
   if (nargout == 1)      cJac=[];   end
end

% allVariables = true/false if finite-differencing wrt to one/all variables
allVariables = false;
if ischar(variables)
   if strcmp(variables,'all')
      variables = 1:numberOfVariables;
      allVariables = true;
   else
      error('Unknown value of input ''variables''.')
   end
end

% Preallocate gradf for speed 
if functionIsScalar
   gradf = zeros(numberOfVariables,1);
elseif allVariables % vector-function and gradf estimates full Jacobian 
   gradf = zeros(fNumberOfRows,numberOfVariables); 
else % vector-function and gradf estimates one column of Jacobian
   gradf = zeros(fNumberOfRows,1);
end
   
for gcnt=variables
   if (gcnt == numberOfVariables)      FLAG = -1;    end
   temp = xCurrent(gcnt);
   xCurrent(gcnt)= temp + CHG(gcnt);
         
   if (nonEmptyLowerBounds && isfinite(lb(gcnt))) || (nonEmptyUpperBounds && isfinite(ub(gcnt)))
      % Enforce bounds while finite-differencing.
      % Need lb(gcnt) ~= ub(gcnt), and lb(gcnt) <= temp <= ub(gcnt) to enforce bounds.
      % (If the last qpsub problem was 'infeasible', the bounds could be currently violated.)
      if (lb(gcnt) ~= ub(gcnt)) && (temp >= lb(gcnt)) && (temp <= ub(gcnt)) 
          if  ((xCurrent(gcnt) > ub(gcnt)) || (xCurrent(gcnt) < lb(gcnt))) % outside bound ?
              CHG(gcnt) = -CHG(gcnt);
              xCurrent(gcnt)= temp + CHG(gcnt);
              if (xCurrent(gcnt) > ub(gcnt)) || (xCurrent(gcnt) < lb(gcnt)) % outside other bound ?
                  [newchg,indsign] = max([temp-lb(gcnt), ub(gcnt)-temp]);  % largest distance to bound
                  if newchg >= DiffMinChange
                      CHG(gcnt) = ((-1)^indsign)*newchg;  % make sure sign is correct
                      xCurrent(gcnt)= temp + CHG(gcnt);
                  else
                      errmsg = sprintf('%s %d %s\n%s\n%s %0.5g%s\n',...
                          'Distance between lower and upper bounds, in dimension',gcnt,', is too small to compute', ...
                          'finite-difference approximation of derivative. Increase distance between these', ...
                          'bounds to be at least',2*DiffMinChange,'.');
                      error(errmsg)
                  end          
              end
          end
      end
   end % of 'if isfinite(lb(gcnt)) || isfinite(ub(gcnt))'
   
   xOriginalShape(:) = xCurrent;
   fplus = feval(funfcn{3},xOriginalShape,varargin{:});
   % YDATA: Only used by lsqcurvefit, which has no nonlinear constraints
   % (the only type of constraints we do finite differences on: bounds 
   % and linear constraints do not require finite differences) and thus 
   % no needed after evaluation of constraints
   if ~isempty(YDATA)    
      fplus = fplus - YDATA;
   end
   % Make sure it's in column form
   fplus = fplus(:);

   if functionIsScalar
      gradf(gcnt,1) =  (fplus-fCurrent)/CHG(gcnt);
   elseif allVariables % vector-function and gradf estimates full Jacobian 
      gradf(:,gcnt) = (fplus-fCurrent)/CHG(gcnt);
   else % vector-function and gradf estimates only one column of Jacobian
      gradf = (fplus-fCurrent)/CHG(gcnt);
   end

   if ~isempty(cJac) % Constraint gradient required
      % This is necessary in case confcn is empty, then in the comparison 
      % below confcn{2} would be out of range
      if ~isempty(confcn{1}) 
         [ctmp,ceqtmp] = feval(confcn{3},xOriginalShape,varargin{:});
         cplus = [ceqtmp(:); ctmp(:)];
      end
      % Next line used for problems with varying number of constraints
      if ~isempty(cplus)
         cJac(:,gcnt) = (cplus - cCurrent)/CHG(gcnt); 
      end           
   end
    xCurrent(gcnt) = temp;
end % for 

% --------------------------------------------------------------------------------
function [X,lambda,exitflag,output,how,ACTIND] = ...
    qpsub(H,f,A,B,lb,ub,X,neqcstr,verbosity,caller,ncstr, ...
    numberOfVariables,options,defaultopt,ACTIND)
%QP Quadratic programming subproblem. Handles qp and constrained
%   linear least-squares as well as subproblems generated from NLCONST.
%
%   X=QP(H,f,A,b) solves the quadratic programming problem:
%
%            min 0.5*x'Hx + f'x   subject to:  Ax <= b 
%             x    
%

%   Copyright 1990-2003 The MathWorks, Inc. 

% Define constant strings
NewtonStep = 'Newton';
NegCurv = 'negative curvature chol';   
ZeroStep = 'zero step';
SteepDescent = 'steepest descent';
Conls = 'lsqlin';
Lp = 'linprog';
Qp = 'quadprog';
Qpsub = 'qpsub';
how = 'ok'; 

exitflag = 1;
output = [];
iterations = 0;
if nargin < 15, ACTIND = [];
    if nargin < 13, options = []; 
    end
end

lb=lb(:); ub = ub(:);

if isempty(neqcstr), neqcstr = 0; end

LLS = 0;
if strcmp(caller, Conls)
    LLS = 1;
    [rowH,colH]=size(H);
    numberOfVariables = colH;
end
if strcmp(caller, Qpsub)
    normalize = -1;
else
    normalize = 1;
end

simplex_iter = 0;
if  norm(H,'inf')==0 || isempty(H), is_qp=0; else, is_qp=1; end

if LLS==1
    is_qp=0;
end

normf = 1;
if normalize > 0
    % Check for lp
    if (~is_qp && ~LLS)
        normf = norm(f);
        if (normf > 0)  f = f./normf;   end
    end
end

% Handle bounds as linear constraints
arglb = ~eq(lb,-inf);
lenlb=length(lb); % maybe less than numberOfVariables due to old code
if nnz(arglb) > 0     
    lbmatrix = -eye(lenlb,numberOfVariables);
    A=[A; lbmatrix(arglb,1:numberOfVariables)]; % select non-Inf bounds
    B=[B;-lb(arglb)];
end

argub = ~eq(ub,inf);
lenub=length(ub);
if (nnz(argub) > 0)
    ubmatrix = eye(lenub,numberOfVariables);
    A=[A; ubmatrix(argub,1:numberOfVariables)];
    B=[B; ub(argub)];
end 

% Bounds are treated as constraints: Reset ncstr accordingly
ncstr=ncstr + nnz(arglb) + nnz(argub);

% Figure out max iteration count
% For linprog/quadprog/lsqlin/qpsub problems, use 'MaxIter' for this.
% For nlconst (fmincon, etc) problems, use 'MaxSQPIter' for this.

if (isempty(options) || isempty(options.MaxSQPIter))
  maxiter = 10*max(numberOfVariables,ncstr-neqcstr);
else 
  maxiter = options.MaxSQPIter;
end 

% Used for determining threshold for whether a direction will violate a constraint.
normA = ones(ncstr,1);
if (normalize > 0)
    for i=1:ncstr
        n = norm(A(i,:));
        if (n ~= 0)
            A(i,:) = A(i,:)/n;
            B(i) = B(i)/n;
            normA(i,1) = n;
        end
    end
else 
    normA = ones(ncstr,1);
end
errnorm = 0.01*sqrt(eps); 

tolDep = 100*numberOfVariables*eps;      
lambda = zeros(ncstr,1);
eqix = 1:neqcstr;

% Modifications for warm-start.
% Do some error checking on the incoming working set indices.
ACTCNT = length(ACTIND);
if isempty(ACTIND)
    ACTIND = eqix;
elseif neqcstr > 0
    i = max(find(ACTIND<=neqcstr));
    if isempty(i) || i > neqcstr % safeguard which should not occur
        ACTIND = eqix;
    elseif i < neqcstr
        % A redundant equality constraint was removed on the last
        % SQP iteration.  We will go ahead and reinsert it here.
        numremoved = neqcstr - i;
        ACTIND(neqcstr+1:ACTCNT+numremoved) = ACTIND(i+1:ACTCNT);
        ACTIND(1:neqcstr) = eqix;
    end
end
aix = zeros(ncstr,1);
aix(ACTIND) = 1;
ACTCNT = length(ACTIND);
ACTSET = A(ACTIND,:);

% Check that the constraints in the initial working set are not
% dependent and find an initial point which satisfies the initial working set.
indepInd = 1:ncstr;
remove = [];
if ACTCNT > 0 && normalize ~= -1
    % call constraint solver
    [Q,R,A,B,X,Z,how,ACTSET,ACTIND,ACTCNT,aix,eqix,neqcstr,ncstr, ...
            remove,exitflag]= ...
        eqnsolv(A,B,eqix,neqcstr,ncstr,numberOfVariables,LLS,H,X,f, ...
        normf,normA,verbosity,aix,ACTSET,ACTIND,ACTCNT,how,exitflag); 
    
    if ~isempty(remove)
        indepInd(remove)=[];
        normA = normA(indepInd);
    end
    
    if strcmp(how,'infeasible')
        % Equalities are inconsistent, so X and lambda have no valid values
        % Return original X and zeros for lambda.
        ACTIND = indepInd(ACTIND);
        output.iterations = iterations;
        exitflag = -1;
        return
    end
    
    err = 0;
    if neqcstr >= numberOfVariables
        err = max(abs(A(eqix,:)*X-B(eqix)));
        if (err > 1e-8)  % Equalities not met
            how='infeasible';
            exitflag = -1;
            % Equalities are inconsistent, X and lambda have no valid values
            % Return original X and zeros for lambda.
            ACTIND = indepInd(ACTIND);
            output.iterations = iterations;
            return
        else % Check inequalities
            if (max(A*X-B) > 1e-8)
                how = 'infeasible';
                % was exitflag = 8; 
                exitflag = -1;
            end
        end
        if is_qp
            actlambda = -R\(Q'*(H*X+f));
        elseif LLS
            actlambda = -R\(Q'*(H'*(H*X-f)));
        else
            actlambda = -R\(Q'*f);
        end
        lambda(indepInd(ACTIND)) = normf * (actlambda ./normA(ACTIND));
        ACTIND = indepInd(ACTIND);
        output.iterations = iterations;
        return
    end
    
    % Check whether in Phase 1 of feasibility point finding. 
    if (verbosity == -2)
        cstr = A*X-B; 
        mc=max(cstr(neqcstr+1:ncstr));
        if (mc > 0)
            X(numberOfVariables) = mc + 1;
        end
    end
else 
    if ACTCNT == 0 % initial working set is empty 
        Q = eye(numberOfVariables,numberOfVariables);
        R = [];        Z = 1;
    else           % in Phase I and working set not empty
        [Q,R] = qr(ACTSET');
        Z = Q(:,ACTCNT+1:numberOfVariables);
    end   
end

% Find Initial Feasible Solution 
cstr = A*X-B;
if ncstr > neqcstr
    mc = max(cstr(neqcstr+1:ncstr));
else
    mc = 0;
end
if mc > eps
    quiet = -2;
    optionsPhase1 = []; % Use default options in phase 1
    ACTIND2 = 1:neqcstr;
    A2=[[A;zeros(1,numberOfVariables)],[zeros(neqcstr,1);-ones(ncstr+1-neqcstr,1)]];
    [XS,lambdaS,exitflagS,outputS,howS,ACTIND2] = ...
        qpsub([],[zeros(numberOfVariables,1);1],A2,[B;1e-5], ...
        [],[],[X;mc+1],neqcstr,quiet,Qpsub,size(A2,1),numberOfVariables+1, ...
        optionsPhase1,defaultopt,ACTIND2);
    slack = XS(numberOfVariables+1);
    X=XS(1:numberOfVariables);
    cstr=A*X-B;
    if slack > eps 
        if slack > 1e-8 
            how='infeasible';
            exitflag = -1;
        else
            how = 'overly constrained';
            exitflag = -1;
        end
        lambda(indepInd) = normf * (lambdaS((1:ncstr)')./normA);
        ACTIND = 1:neqcstr;
        ACTIND = indepInd(ACTIND);
        output.iterations = iterations;
        return
    else
        % Initialize active set info based on solution of Phase I.
        %      ACTIND = ACTIND2(find(ACTIND2<=ncstr));
        ACTIND = 1:neqcstr;
        ACTSET = A(ACTIND,:);
        ACTCNT = length(ACTIND);
        aix = zeros(ncstr,1);
        aix(ACTIND) = 1;
        if ACTCNT == 0
            Q = zeros(numberOfVariables,numberOfVariables);
            R = [];
            Z = 1;
        else
            [Q,R] = qr(ACTSET');
            Z = Q(:,ACTCNT+1:numberOfVariables);
        end
    end
end

if ACTCNT >= numberOfVariables - 1  
    simplex_iter = 1; 
end
[m,n]=size(ACTSET);

if (is_qp)
    gf=H*X+f;
    [SD, dirType] = compdir(Z,H,gf,numberOfVariables,f);
elseif (LLS)
    HXf=H*X-f;
    gf=H'*(HXf);
    HZ= H*Z;
    [mm,nn]=size(HZ);
    if mm >= nn
        [QHZ, RHZ] =  qr(HZ,0);
        Pd = QHZ'*HXf;
        % Now need to check which is dependent
        if min(size(RHZ))==1 % Make sure RHZ isn't a vector
            depInd = find( abs(RHZ(1,1)) < tolDep);
        else
            depInd = find( abs(diag(RHZ)) < tolDep );
        end  
    end
    if mm >= nn && isempty(depInd) % Newton step
        SD = - Z*(RHZ(1:nn, 1:nn) \ Pd(1:nn,:));
        dirType = NewtonStep;
    else % steepest descent direction
        SD = -Z*(Z'*gf);
        dirType = SteepDescent;
    end
else % lp
    gf = f;
    SD=-Z*Z'*gf;
    dirType = SteepDescent; 
    if norm(SD) < 1e-10 && neqcstr
        % This happens when equality constraint is perpendicular to objective function f.x.
        actlambda = -R\(Q'*(gf));
        lambda(indepInd(ACTIND)) = normf * (actlambda ./ normA(ACTIND));
        ACTIND = indepInd(ACTIND);
        output.iterations = iterations;
        return;
    end
end

oldind = 0; 

%--------------Main Routine-------------------

while iterations < maxiter
    iterations = iterations + 1;
    if isinf(verbosity)
      curr_out = sprintf('Iter: %5.0f, Active: %5.0f, step: %s, proc: %s',iterations,ACTCNT,dirType,how);
        disp(curr_out); 
    end
    
    % Find distance we can move in search direction SD before a 
    % constraint is violated.
    % Gradient with respect to search direction.
    GSD=A*SD;
    
    % Note: we consider only constraints whose gradients are greater
    % than some threshold. If we considered all gradients greater than 
    % zero then it might be possible to add a constraint which would lead to
    % a singular (rank deficient) working set. The gradient (GSD) of such
    % a constraint in the direction of search would be very close to zero.
    indf = find((GSD > errnorm * norm(SD))  &  ~aix);
    
    if isempty(indf) % No constraints to hit
        STEPMIN=1e16;
        dist=[]; ind2=[]; ind=[];
    else % Find distance to the nearest constraint
        dist = abs(cstr(indf)./GSD(indf));
        [STEPMIN,ind2] =  min(dist);
        ind2 = find(dist == STEPMIN);
        % Bland's rule for anti-cycling: if there is more than one 
        % blocking constraint then add the one with the smallest index.
        ind=indf(min(ind2));
    end
    
    %----------------Update X---------------------
    
    % Assume we do not delete a constraint
    delete_constr = 0;   
    
    if ~isempty(indf) && isfinite(STEPMIN) % Hit a constraint
        if strcmp(dirType, NewtonStep)
            % Newton step and hit a constraint: LLS or is_qp
            if STEPMIN > 1  % Overstepped minimum; reset STEPMIN
                STEPMIN = 1;
                delete_constr = 1;
            end
            X = X+STEPMIN*SD;
        else
            % Not a Newton step and hit a constraint: is_qp or LLS or maybe lp
            X = X+STEPMIN*SD;  
        end              
    else %  isempty(indf) | ~isfinite(STEPMIN)
        % did not hit a constraint
        if strcmp(dirType, NewtonStep)
            % Newton step and no constraint hit: LLS or maybe is_qp
            STEPMIN = 1;   % Exact distance to the solution. Now delete constr.
            X = X + SD;
            delete_constr = 1;
        else % Not a Newton step: is_qp or lp or LLS
            
            if (~is_qp && ~LLS) || strcmp(dirType, NegCurv) % LP or neg def (implies is_qp)
                % neg def -- unbounded
                if norm(SD) > errnorm
                    if normalize < 0
                        STEPMIN=abs((X(numberOfVariables)+1e-5)/(SD(numberOfVariables)+eps));
                    else 
                        STEPMIN = 1e16;
                    end
                    X=X+STEPMIN*SD;
                    how='unbounded'; 
                    % was exitflag = 5; 
                    exitflag = -1;
                else % norm(SD) <= errnorm
                    how = 'ill posed';
                    % was exitflag = 6; 
                    exitflag = -1;
                end
                ACTIND = indepInd(ACTIND);
                output.iterations = iterations;
                return
            else % singular: solve compatible system for a solution: is_qp or LLS
                if is_qp
                    projH = Z'*H*Z; 
                    Zgf = Z'*gf;
                    projSD = pinv(projH)*(-Zgf);
                else % LLS
                    projH = HZ'*HZ; 
                    Zgf = Z'*gf;
                    projSD = pinv(projH)*(-Zgf);
                end
                
                % Check if compatible
                if norm(projH*projSD+Zgf) > 10*eps*(norm(projH) + norm(Zgf))
                    % system is incompatible --> it's a "chute": use SD from compdir
                    % unbounded in SD direction
                    if norm(SD) > errnorm
                        if normalize < 0
                            STEPMIN=abs((X(numberOfVariables)+1e-5)/(SD(numberOfVariables)+eps));
                        else 
                            STEPMIN = 1e16;
                        end
                        X=X+STEPMIN*SD;
                        how='unbounded'; 
                        % was exitflag = 5;
                        exitflag = -1;
                    else % norm(SD) <= errnorm
                        how = 'ill posed';
                        %was exitflag = 6;
                        exitflag = -1;
                    end
                    ACTIND = indepInd(ACTIND);
                    output.iterations = iterations;
                    return
                else % Convex -- move to the minimum (compatible system)
                    SD = Z*projSD;
                    if gf'*SD > 0
                        SD = -SD;
                    end
                    dirType = 'singular';
                    % First check if constraint is violated.
                    GSD=A*SD;
                    indf = find((GSD > errnorm * norm(SD))  &  ~aix);
                    if isempty(indf) % No constraints to hit
                        STEPMIN=1;
                        delete_constr = 1;
                        dist=[]; ind2=[]; ind=[];
                    else % Find distance to the nearest constraint
                        dist = abs(cstr(indf)./GSD(indf));
                        [STEPMIN,ind2] =  min(dist);
                        ind2 = find(dist == STEPMIN);
                        % Bland's rule for anti-cycling: if there is more than one 
                        % blocking constraint then add the one with the smallest index.
                        ind=indf(min(ind2));
                    end
                    if STEPMIN > 1  % Overstepped minimum; reset STEPMIN
                        STEPMIN = 1;
                        delete_constr = 1;
                    end
                    X = X + STEPMIN*SD; 
                end
            end % if ~is_qp | smallRealEig < -eps
        end % if strcmp(dirType, NewtonStep)
    end % if ~isempty(indf)& isfinite(STEPMIN) % Hit a constraint
    
    %----Check if reached minimum in current subspace-----
    
    if delete_constr
        % Note: only reach here if a minimum in the current subspace found
        %       LP's do not enter here.
        if ACTCNT>0
            if is_qp
                rlambda = -R\(Q'*(H*X+f));
            elseif LLS
                rlambda = -R\(Q'*(H'*(H*X-f)));
                % else: lp does not reach this point
            end
            actlambda = rlambda;
            actlambda(eqix) = abs(rlambda(eqix));
            indlam = find(actlambda < 0);
            if (~length(indlam)) 
                lambda(indepInd(ACTIND)) = normf * (rlambda./normA(ACTIND));
                ACTIND = indepInd(ACTIND);
                output.iterations = iterations;
                return
            end
            % Remove constraint
            lind = find(ACTIND == min(ACTIND(indlam)));
            lind = lind(1);
            ACTSET(lind,:) = [];
            aix(ACTIND(lind)) = 0;
            [Q,R]=qrdelete(Q,R,lind);
            ACTIND(lind) = [];
            ACTCNT = length(ACTIND);
            simplex_iter = 0;
            ind = 0;
        else % ACTCNT == 0
            output.iterations = iterations;
            return
        end
        delete_constr = 0;
    end
    
    % If we are in the Phase-1 procedure check if the slack variable
    % is zero indicating we have found a feasible starting point.
    if normalize < 0
        if X(numberOfVariables,1) < eps
            ACTIND = indepInd(ACTIND);
            output.iterations = iterations;
            return;
        end
    end   
    
    % Calculate gradient w.r.t objective at this point
    if is_qp
        gf=H*X+f;
    elseif LLS % LLS
        gf=H'*(H*X-f);
        % else gf=f still true.
    end
    
    % Update constraints
    cstr = A*X-B;
    cstr(eqix) = abs(cstr(eqix));
    if max(cstr) > 1e5 * errnorm
        if max(cstr) > norm(X) * errnorm 
            how='unreliable'; 
            exitflag = -1;
        end
    end
    
    %----Add blocking constraint to working set----
    if ind % Hit a constraint
        aix(ind)=1;
        CIND = length(ACTIND) + 1;
        ACTSET(CIND,:)=A(ind,:);
        ACTIND(CIND)=ind;
        [m,n]=size(ACTSET);
        [Q,R] = qrinsert(Q,R,CIND,A(ind,:)');
        ACTCNT = length(ACTIND);
    end
    if ~simplex_iter
        [m,n]=size(ACTSET);
        Z = Q(:,m+1:n);
        if ACTCNT == numberOfVariables - 1, simplex_iter = 1; end
        oldind = 0; 
    else
        
        %---If Simplex Alg. choose leaving constraint---
        rlambda = -R\(Q'*gf);
        if isinf(rlambda(1)) && rlambda(1) < 0 
            fprintf('         Working set is singular; results may still be reliable.\n');
            [m,n] = size(ACTSET);
            rlambda = -(ACTSET + sqrt(eps)*randn(m,n))'\gf;
        end
        actlambda = rlambda;
        actlambda(eqix)=abs(actlambda(eqix));
        indlam = find(actlambda<0);
        if length(indlam)
            if STEPMIN > errnorm
                % If there is no chance of cycling then pick the constraint 
                % which causes the biggest reduction in the cost function. 
                % i.e the constraint with the most negative Lagrangian 
                % multiplier. Since the constraints are normalized this may 
                % result in less iterations.
                [minl,lind] = min(actlambda);
            else
                % Bland's rule for anti-cycling: if there is more than one 
                % negative Lagrangian multiplier then delete the constraint
                % with the smallest index in the active set.
                lind = find(ACTIND == min(ACTIND(indlam)));
            end
            lind = lind(1);
            ACTSET(lind,:) = [];
            aix(ACTIND(lind)) = 0;
            [Q,R]=qrdelete(Q,R,lind);
            Z = Q(:,numberOfVariables);
            oldind = ACTIND(lind);
            ACTIND(lind) = [];
            ACTCNT = length(ACTIND);
        else
            lambda(indepInd(ACTIND))= normf * (rlambda./normA(ACTIND));
            ACTIND = indepInd(ACTIND);
            output.iterations = iterations;
            return
        end
    end %if ACTCNT<numberOfVariables
    
    %----------Compute Search Direction-------------      
    if (is_qp)
        Zgf = Z'*gf; 
        if ~isempty(Zgf) && (norm(Zgf) < 1e-15)
            SD = zeros(numberOfVariables,1); 
            dirType = ZeroStep;
        else
            [SD, dirType] = compdir(Z,H,gf,numberOfVariables,f);
        end
    elseif (LLS)
        Zgf = Z'*gf;
        HZ = H*Z;
        if (norm(Zgf) < 1e-15)
            SD = zeros(numberOfVariables,1);
            dirType = ZeroStep;
        else
            HXf=H*X-f;
            gf=H'*(HXf);
            [mm,nn]=size(HZ);
            if mm >= nn
                [QHZ, RHZ] =  qr(HZ,0);
                Pd = QHZ'*HXf;
                % SD = - Z*(RHZ(1:nn, 1:nn) \ Pd(1:nn,:));
                % Now need to check which is dependent
                if min(size(RHZ))==1 % Make sure RHZ isn't a vector
                    depInd = find( abs(RHZ(1,1)) < tolDep);
                else
                    depInd = find( abs(diag(RHZ)) < tolDep );
                end  
            end
            if mm >= nn && isempty(depInd) % Newton step
                SD = - Z*(RHZ(1:nn, 1:nn) \ Pd(1:nn,:));
                dirType = NewtonStep;
            else % steepest descent direction
                SD = -Z*(Z'*gf);
                dirType = SteepDescent;
            end
        end
    else % LP
        if ~simplex_iter
            SD = -Z*(Z'*gf);
            gradsd = norm(SD);
        else
            gradsd = Z'*gf;
            if  (gradsd > 0)    SD = -Z;
            else                SD = Z;
            end
        end
        if abs(gradsd) < 1e-10 % Search direction null
            % Check whether any constraints can be deleted from active set.
            % rlambda = -ACTSET'\gf;
            if ~oldind
                rlambda = -R\(Q'*gf);
                ACTINDtmp = ACTIND; Qtmp = Q; Rtmp = R;
            else
                % Reinsert just deleted constraint.
                ACTINDtmp = ACTIND;
                ACTINDtmp(lind+1:ACTCNT+1) = ACTIND(lind:ACTCNT);
                ACTINDtmp(lind) = oldind;
                [Qtmp,Rtmp] = qrinsert(Q,R,lind,A(oldind,:)');
            end
            actlambda = rlambda;
            actlambda(1:neqcstr) = abs(actlambda(1:neqcstr));
            indlam = find(actlambda < errnorm);
            lambda(indepInd(ACTINDtmp)) = normf * (rlambda./normA(ACTINDtmp));
            if ~length(indlam)
                ACTIND = indepInd(ACTIND);
                output.iterations = iterations;
                return
            end
            cindmax = length(indlam);
            cindcnt = 0;
            m = length(ACTINDtmp);
            while (abs(gradsd) < 1e-10) && (cindcnt < cindmax)
                cindcnt = cindcnt + 1;
                lind = indlam(cindcnt);
                [Q,R]=qrdelete(Qtmp,Rtmp,lind);
                Z = Q(:,m:numberOfVariables);
                if m ~= numberOfVariables
                    SD = -Z*Z'*gf;
                    gradsd = norm(SD);
                else
                    gradsd = Z'*gf;
                    if  (gradsd > 0)    SD = -Z;
                    else                SD = Z;
                    end
                end
            end
            if abs(gradsd) < 1e-10  % Search direction still null
                ACTIND = indepInd(ACTIND);
                output.iterations = iterations;
                return;
            else
                ACTIND = ACTINDtmp;
                ACTIND(lind) = [];
                aix = zeros(ncstr,1);
                aix(ACTIND) = 1;
                ACTCNT = length(ACTIND);
                ACTSET = A(ACTIND,:);
            end
            lambda = zeros(ncstr,1);
        end
    end % if is_qp
end % while 

if iterations >= maxiter
    exitflag = 0;
    how = 'MaxSQPIter';
end

output.iterations = iterations;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Q,R,A,B,X,Z,how,ACTSET,ACTIND,ACTCNT,aix,eqix,neqcstr,ncstr,remove,exitflag]= ...
    eqnsolv(A,B,eqix,neqcstr,ncstr,numberOfVariables,LLS,H,X,f,normf, ...
    normA,verbosity,aix,ACTSET,ACTIND,ACTCNT,how,exitflag)
% EQNSOLV Helper function for QPSUB.
%    Checks whether the working set is linearly independent and
%    finds a feasible point with respect to the working set constraints.
%    If the equalities are dependent but not consistent, warning
%    messages are given. If the equalities are dependent but consistent, 
%    the redundant constraints are removed and the corresponding variables 
%    adjusted.

% set tolerances
tolDep = 100*numberOfVariables*eps;      
tolCons = 1e-10;

Z=[]; remove =[];

% First see if the equality constraints form a consistent system.
[Qa,Ra,Ea]=qr(A(eqix,:));

% Form vector of dependent indices.
if min(size(Ra))==1 % Make sure Ra isn't a vector
    depInd = find( abs(Ra(1,1)) < tolDep);
else
    depInd = find( abs(diag(Ra)) < tolDep );
end
if neqcstr > numberOfVariables
    depInd = [depInd; ((numberOfVariables+1):neqcstr)'];
end      

if ~isempty(depInd)    % equality constraints are dependent
    how='dependent';
    exitflag = 1;
    bdepInd =  abs(Qa(:,depInd)'*B(eqix)) >= tolDep ;
    
    if any( bdepInd ) % Not consistent
        how='infeasible';   
        exitflag = -1;
    else % the equality constraints are consistent
        % Delete the redundant constraints
        % By QR factoring the transpose, we see which columns of A' (rows of A) move to the end
        [Qat,Rat,Eat]=qr(A(eqix,:)');        
        [i,j] = find(Eat); % Eat permutes the columns of A' (rows of A)
        remove = i(depInd);
        numDepend = nnz(remove);
        A(eqix(remove),:)=[];
        B(eqix(remove))=[];
        neqcstr = neqcstr - numDepend;
        ncstr = ncstr - numDepend;
        eqix = 1:neqcstr;
        aix(remove) = [];
        ACTIND(1:numDepend) = [];
        ACTIND = ACTIND - numDepend;      
        ACTSET = A(ACTIND,:);
        ACTCNT = ACTCNT - numDepend;
    end % consistency check
end % dependency check

% Now that we have done all we can to make the equality constraints
% consistent and independent we will check the inequality constraints
% in the working set.  First we want to make sure that the number of 
% constraints in the working set is only greater than or equal to the
% number of variables if the number of (non-redundant) equality 
% constraints is greater than or equal to the number of variables.
if ACTCNT >= numberOfVariables
    ACTCNT = max(neqcstr, numberOfVariables-1);
    ACTIND = ACTIND(1:ACTCNT);
    ACTSET = A(ACTIND,:);
    aix = zeros(ncstr,1);
    aix(ACTIND) = 1;
end

% Now check to see that all the constraints in the working set are linearly independent.
if ACTCNT > neqcstr
    [Qat,Rat,Eat]=qr(ACTSET');
    
    % Form vector of dependent indices.
    if min(size(Rat))==1            % Make sure Rat isn't a vector
        depInd = find( abs(Rat(1,1)) < tolDep);
    else
        depInd = find( abs(diag(Rat)) < tolDep );
    end
    
    if ~isempty(depInd)
        [i,j] = find(Eat);          % Eat permutes the columns of A' (rows of A)
        remove2 = i(depInd);
        removeEq   = remove2(find(remove2 <= neqcstr));
        removeIneq = remove2(find(remove2 > neqcstr));
        
        if (~isempty(removeEq))     % Just take equalities as initial working set.
            ACTIND = 1:neqcstr; 
        else                        % Remove dependent inequality constraints.
            ACTIND(removeIneq) = [];
        end
        aix = zeros(ncstr,1);
        aix(ACTIND) = 1;
        ACTSET = A(ACTIND,:);
        ACTCNT = length(ACTIND);
    end  
end

[Q,R]=qr(ACTSET');
Z = Q(:,ACTCNT+1:numberOfVariables);

if ~strcmp(how,'infeasible') && ACTCNT > 0
    % Find point closest to the given initial X which satisfies working set constraints.
    minnormstep = Q(:,1:ACTCNT) * ...
        ((R(1:ACTCNT,1:ACTCNT)') \ (B(ACTIND) - ACTSET*X));
    X = X + minnormstep; 
    % Sometimes the "basic" solution satisfies Aeq*x= Beq 
    % and A*X < B better than the minnorm solution. Choose the one
    % that the minimizes the max constraint violation.
    err = A*X - B;
    err(eqix) = abs(err(eqix));
    if any(err > eps)
        Xbasic = ACTSET\B(ACTIND);
        errbasic = A*Xbasic - B;
        errbasic(eqix) = abs(errbasic(eqix));
        if (max(errbasic) < max(err))   X = Xbasic;     end
    end
end

% --------------------------------------------------------------------------------
function [SD, dirType] = compdir(Z,H,gf,nvars,f);
% COMPDIR Computes a search direction in a subspace defined by Z.  [SD,
% dirType] = compdir(Z,H,gf,nvars,f) returns a search direction for the
% subproblem 0.5*Z'*H*Z + Z'*gf. Helper function for NLCONST. SD is Newton
% direction if possible. SD is a direction of negative curvature if the
% Cholesky factorization of Z'*H*Z fails. If the negative curvature
% direction isn't negative "enough", SD is the steepest descent direction.
% For singular Z'*H*Z, SD is the steepest descent direction even if small,
% or even zero, magnitude.

%   Copyright 1990-2003 The MathWorks, Inc.

% Define constant strings
Newton = 'Newton';                     % Z'*H*Z positive definite
NegCurv = 'negative curvature chol';   % Z'*H*Z indefinite
SteepDescent = 'steepest descent';     % Z'*H*Z (nearly) singular

dirType = [];
% Compute the projected Newton direction if possible
projH = Z'*H*Z;
[R, p] = chol(projH);
if ~p  % positive definite: use Newton direction
    SD = - Z*(R \ ( R'\(Z'*gf)));
    dirType = Newton;
else % not positive definite
    [L,sneg] = choltrap(projH);
    if ~isempty(sneg) & sneg'*projH*sneg < -sqrt(eps) % if negative enough
        SD = Z*sneg;
        dirType = NegCurv;
    else % Not positive definite, not negative definite "enough" so use steepest descent direction
        stpDesc = - Z*(Z'*gf);
        % ||SD|| may be (close to) zero, but qpsub handles that case
        SD = stpDesc;
        dirType = SteepDescent;
    end %   
end % ~p  (positive definite)

% Make sure it is a descent direction
if gf'*SD > 0
    SD = -SD;
end

%-----------------------------------------------
function [L,sneg] = choltrap(A)
% CHOLTRAP Compute Cholesky factor or direction of negative curvature.
%     [L, SNEG] = CHOLTRAP(A) computes the Cholesky factor L, such that
%     L*L'= A, if it exists, or returns a direction of negative curvature
%     SNEG for matrix A when A is not positive definite. If A is positive
%     definite, SNEG will be []. 
%
%     If A is singular, it is possible that SNEG will not be a direction of
%     negative curvature (but will be nonempty). In particular, if A is
%     positive semi-definite, SNEG will be nonempty but not a direction of
%     negative curvature. If A is indefinite but singular, SNEG may or may
%     not be a direction of negative curvature.

sneg = [];
n = size(A,1);
L = eye(n);
tol = 0;    % Dividing by sqrt of small number isn't a problem 
for k=1:n-1
    if A(k,k) <= tol
        elem = zeros(length(A),1); 
        elem(k,1) = 1;
        sneg = L'\elem;
        return;
    else
        L(k,k) = sqrt(A(k,k));
        s = k+1:n;
        L(s,k) = A(s,k)/L(k,k);
        A(k+1:n,k+1:n) =  A(k+1:n,k+1:n)  - tril(L(k+1:n,k)*L(k+1:n,k)');   
    end
end
if A(n,n) <= tol
    elem = zeros(length(A),1); 
    elem(n,1) = 1;
    sneg = L'\elem;
else
    L(n,n) = sqrt(A(n,n));
end

%------------------------------------------------------------------
function value = optimgetfast(options,name,defaultopt)
%OPTIMGETFAST Get OPTIM OPTIONS parameter with no error checking so fast.

if (~isempty(options))      value = options.(name);
else                        value = [];
end

if (isempty(value))    value = defaultopt.(name);   end

% -----------------------------------------------------------------------
function [allfcns,msg] = optimfcnchk1(funstr,caller,lenVarIn,gradflag,hessflag,constrflag)
% OPTIMFCNCHK Pre- and post-process function expression for FUNCHK.

%   Copyright 1990-2002 The MathWorks, Inc. 
%   $Revision: 1.7 $  $Date: 2002/03/12 20:36:19 $

% Initialize
if nargin < 6
    constrflag = 0;
    if nargin < 5
        hessflag = 0;
        if nargin < 4
            gradflag = 0;
        end,end,end

msg='';
nonlconmsg =  ['NONLCON must be a function.'];
allfcns = {};
funfcn = [];
gradfcn = [];
hessfcn = [];
calltype = 'fun';

% {fun}
if (isa(funstr, 'cell') & length(funstr)==1)    % take the cellarray apart: we know it is nonempty
    [funfcn, msg] = fcnchk1(funstr{1},lenVarIn);
    if (~isempty(msg))      error(msg);    end
    % {fun,[]}      
elseif isa(funstr, 'cell') & length(funstr)==2 & isempty(funstr{2})
    [funfcn, msg] = fcnchk1(funstr{1},lenVarIn);
    if (~isempty(msg))      error(msg);    end
    % {fun, grad}   
elseif isa(funstr, 'cell') & length(funstr)==2 % and ~isempty(funstr{2})
    
    [funfcn, msg] = fcnchk1(funstr{1},lenVarIn);
    if (~isempty(msg))      error(msg);    end
    [gradfcn, msg] = fcnchk1(funstr{2},lenVarIn);
    if (~isempty(msg))      error(msg);    end
    calltype = 'fun_then_grad';
    if (~gradflag)
        calltype = 'fun';
    end
    % {fun, [], []}   
elseif isa(funstr, 'cell') & length(funstr)==3 & ~isempty(funstr{1}) & isempty(funstr{2}) & isempty(funstr{3})
    [funfcn, msg] = fcnchk1(funstr{1},lenVarIn);
    if (~isempty(msg))      error(msg);    end
    % {fun, grad, hess}   
elseif isa(funstr, 'cell') & length(funstr)==3 & ~isempty(funstr{2}) & ~isempty(funstr{3})
    [funfcn, msg] = fcnchk1(funstr{1},lenVarIn);
    if (~isempty(msg))      error(msg);    end
    [gradfcn, msg] = fcnchk1(funstr{2},lenVarIn);
    if (~isempty(msg))      error(msg);    end
    [hessfcn, msg] = fcnchk1(funstr{3},lenVarIn);
    if (~isempty(msg))      error(msg);    end
    calltype = 'fun_then_grad_then_hess';
    if ~hessflag & ~gradflag
        calltype = 'fun';
    elseif hessflag & ~gradflag
        calltype = 'fun';
    elseif ~hessflag & gradflag
        calltype = 'fun_then_grad';
    end
    
    % {fun, grad, []}   
elseif isa(funstr, 'cell') & length(funstr)==3 & ~isempty(funstr{2}) & isempty(funstr{3})
    [funfcn, msg] = fcnchk1(funstr{1},lenVarIn);
    if (~isempty(msg))      error(msg);    end
    [gradfcn, msg] = fcnchk1(funstr{2},lenVarIn);
    if (~isempty(msg))      error(msg);    end
    calltype = 'fun_then_grad';
    if ~gradflag
        calltype = 'fun';
    end
    
elseif ~isa(funstr, 'cell')  %Not a cell; is a string expression, function name string or inline object
    [funfcn, msg] = fcnchk1(funstr,lenVarIn);
    if (~isempty(msg))      error(msg);    end
else
    errmsg = sprintf('%s\n%s', ...
        'FUN must be a function or an inline object;', ...
        ' or, FUN may be a cell array that contains these type of objects.');
    error(errmsg)
end

allfcns{1} = calltype;      allfcns{2} = caller;
allfcns{3} = funfcn;        allfcns{4} = gradfcn;
allfcns{5} = hessfcn;

% ----------------------------------------------------
function [f,msg] = fcnchk1(fun,varargin)
msg = '';
if isstr(fun)
    f = strtrim(fun);
elseif isa(fun,'function_handle') 
    f = fun; 
else
    f = '';
    msg = ['FUN must be a function, a valid string expression, ', ...
            sprintf('\n'),'or an inline function object.'];
end

%------------------------------------------
function s1 = strtrim(s)
%STRTRIM Trim spaces from string.
if isempty(s)
    s1 = s;
else            % remove leading and trailing blanks (including nulls)
    c = find(s ~= ' ' & s ~= 0);
    s1 = s(min(c):max(c));
end

% --------------------------------------------------------------------
function y = nanmean(x)
%NANMEAN Average or mean ignoring NaNs.
%   Copyright 1993-2003 The MathWorks, Inc. 

if (isempty(x))     y = NaN;    return;     end     % Check for empty input.

% Replace NaNs with zeros.
nans = isnan(x);
i = find(nans);
x(i) = zeros(size(i));

% count terms in sum over first non-singleton dimension
dim = find(size(x)>1);
if (isempty(dim))   dim = 1;
else                dim = dim(1);
end
count = sum(~nans,dim);

% Protect against a column of all NaNs
i = find(count==0);
count(i) = 1;
y = sum(x,dim)./count;
y(i) = NaN;

% -----------------------------------------------
function y = nansum(x)
%NANSUM Sum ignoring NaNs.
%   Copyright 1993-2002 The MathWorks, Inc. 

% Replace NaNs with zeros.
nans = isnan(x);
i = find(nans);
x(i) = zeros(size(i));

% Protect against an entire column of NaNs
y = sum(x);     i = find(all(nans));
y(i) = i + NaN;
