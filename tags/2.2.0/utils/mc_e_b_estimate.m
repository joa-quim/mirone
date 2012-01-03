function mc_e_b_estimate(selt,ud)
%   This was the ZMAP bdiff2 script.
%   The credit for what is donne here goes obviously to the ZMAP team.
%   However, I deserve also some credit for beeing able to decipher among
%   the ZMAP programing mess of scripts transmiting data trough global variables
%   and for code cleaning (converting all scripts to functions) .
%   Mironified by J. Luis at 6-07-05
%
% This routine estimates the b-value of a curve automatically
% The b-value curve is differenciated and the point
% of the magnitude of completeness is marked. The b-value will be calculated
% using this point and the point half way toward the high
% magnitude end of the b-value curve.
%
% Stefan Wiemer 1/95
% last update: J.Woessner, 27.08.04

% Changes record:
% 02.06.03: Added choice of EMR-method to calculate Mc
% 06.10.03: Added choice of MBS-method to calculate Mc (fixed at 5)
%           Added bootstrap choice
% 28.07.04: Many changes: Now able to do computatios with all functions
%           available in calc_Mc

fs10 = 10;

% Create the input interface
%
% when run from timeplot.m selt=in and it an input menu is created
% this initiates a call back, where selt =  ca, and we go directly
% to the calculations, skipping the input menu

if (strcmp(selt, 'in'))
    ud.nBstSample = 100;    ud.fMccorr = 0;     % Default values
    
    ud.nPar = figure('Name','Mc Input Parameter','NumberTitle','off', 'MenuBar','none','units','points', ...
        'Visible','on', 'Position',[200 200 250 130],'Color',get(0,'factoryUicontrolBackgroundColor'));

    % Get list of Mc computation possibilities
    [labelList2] = calc_Mc;
    ud.popup=uicontrol('Style','popup', 'Units','normalized', 'Position',[0.027 0.659 0.934 0.127], 'String',labelList2,...
        'BackGroundColor','w', 'value',1);

    % Editable fields
    ud.field1 = uicontrol('BackGroundColor','w','Style','edit','Units','normalized',...
        'Position',[.817 .491 .141 .121], 'String',num2str(ud.nBstSample),...
        'CallBack',{@edit_bootstr,ud});

    ud.field2 =uicontrol('BackGroundColor','w','Style','edit','Units','normalized', 'Position',[.817 .329 .141 .121],...
        'String',num2str(ud.fMccorr),...
        'CallBack',{@edit_mcCorr,ud});

    % Buttons
    ud.bBst_button = uicontrol('Style','checkbox', 'string','Uncertianty by boostrapping',...
        'Units','normalized','Position',[.027 .468 .453 .087]);

    uicontrol('Style','Pushbutton', 'Units','normalized','Position',[.763 .064 .198 .133],...
        'Callback','close;','String','Cancel');

    uicontrol('Style','Pushbutton','Units','normalized','Position',[.336 .064 .198 .133],'String','Continue', ...
        'CallBack',{@button_go,ud});

    %%%  Test fields
    uicontrol('Style','text','Units','normalized','Position',[.081 .829 .73 .114],'FontSize',12,...
        'FontWeight','bold','String','Maximum likelihood estimation');
    uicontrol('Style','text','Units','normalized','Position',[.577 .491 .228 .104],'FontSize',10,...
        'FontWeight','bold','String','Bootstraps');
    uicontrol('Style','text','Units','normalized','Position',[.52 .347 .288 .092],'FontSize',10,...
        'FontWeight','bold','String','Mc correction');

    set(ud.nPar,'visible','on','UserData',ud);
    
end

% selt = ca after input menu is run and parameters have been set
if (selt == 'ca')           % build figure for plot
    fBinning = 0.1;

    bfig = figure('Units','normalized','NumberTitle','off','Name','Frequency-magnitude distribution',...
        'visible','on', 'pos',[0.3 0.3 0.42 0.6],'Color','w');
    ud2.h_fig = bfig;   ud2.h_mir_fig = ud.h_mir_fig;
    ud2.events_mag = ud.events_mag;     % Copy them here because ud will be destroyed with the input params figure
    ud2.events_time = ud.events_time;

    options = uimenu('Label','ZTools');
    uimenu(options,'Label','Estimate recurrence time/probability','callback',{@plorem,ud2});
    %uimenu(options,'Label','Manual fit of b-value','callback','bfitnew(newcat)');
    uimenu(options,'Label','Save values to file','callback',{@calSave,ud2});

    maxmag = ceil(10*max(ud.events_mag))/10;
    mima = min(ud.events_mag);
    hStep = 0.1;
    if (mima > 100),    hStep = 1;      end     % For hydrophone SL mags
    if (mima > 0),      mima = 0;       end

    % bval contains the number of events in each bin
    % bvalsum is the cum. sum in each bin
    % bval2 is number events in each bin, in reverse order
    % bvalsum3 is reverse order cum. sum.
    % xt3 is the step in magnitude for the bins == .1

    [bval,xt2] = histo_m('hist',ud.events_mag,(mima:hStep:maxmag));
    bvalsum = cumsum(bval);             % N for M <=
    bval2 = bval(length(bval):-1:1);
    bvalsum3 = cumsum(bval(length(bval):-1:1));    % N for M >= (counted backwards)
    xt3 = (maxmag:-hStep:mima);
    
    ud2.xt3 = xt3;    ud2.bvalsum3 = bvalsum3;     % Save those for Saving purposes
    backg_ab = log10(bvalsum3);

    %figure(bfig);   %delete(gca);delete(gca); delete(gca); delete(gca)
    axHnd = axes('position',[0.15, 0.20, 0.8, 0.75]);    % plot Freq-Mag curves

    % plot the cum. sum in each bin

    pl = semilogy(xt3,bvalsum3,'sb','LineWidth',1,'MarkerSize',6,'MarkerFaceColor','w','MarkerEdgeColor','k');
    hold on
    pl1 = semilogy(xt3,bval2,'^b','LineWidth',1,'MarkerSize',4,'MarkerFaceColor',[0.7 0.7 .7],'MarkerEdgeColor','k');

    % CALCULATE the diff in cum sum from the previous bin

    xlabel('Magnitude','FontWeight','bold','FontSize',12)
    ylabel('Cumulative Number','FontWeight','bold','FontSize',12)
    set(axHnd,'visible','on','FontSize',12,'FontWeight','bold','LineWidth',1,...
        'TickDir','out','Ticklength',[0.02 0.02],'Box','on','Tag','cufi','color','w')
    set(bfig,'visible','on','pointer','watch');   pause(0.1)

    % Estimate the b value
    % calculates max likelihood b value(bvml) && WLS(bvls)

    % SET DEFAULTS TO BE ADDED INTERACTIVLY LATER
    Nmin = 10;

    if (length(ud.events_mag) >= Nmin)     % enough events??
        % Added to obtain goodness-of-fit to powerlaw value
        %mcperc_ca3(ud.events_mag);     % O QUE E QUE RETORNA? COMENTEI E NAO FEZ DIFERENCA

        fMc = calc_Mc(ud.events_mag, ud.inpr1, fBinning, ud.fMccorr);
        l = (ud.events_mag >= fMc-(fBinning/2));
        if (length(ud.events_mag(l)) >= Nmin)
            [fMeanMag, fBValue, fStd_B, fAValue] = calc_bmemag(ud.events_mag(l), fBinning);
        else
            fMc = NaN;  fBValue = NaN;  fStd_B = NaN;   fAValue = NaN;
        end

        fStd_A = NaN;   fStd_Mc = NaN;      % Set standard deviation of a-value to nan;

        % Bootstrap uncertainties
        if (ud.bBst_button == 1)               % Check Mc from original catalog
            l = (ud.events_mag >= fMc-(fBinning/2));
            if length(ud.events_mag(l)) >= Nmin
                [fMc, fStd_Mc, fBValue, fStd_B, fAValue, fStd_A, vMc, mBvalue] = ...
                    calc_McBboot(ud.events_mag, fBinning, ud.nBstSample, ud.inpr1, Nmin, ud.fMccorr);
            else
                fMc = NaN; fStd_Mc = NaN; fBValue = NaN; fStd_B = NaN; fAValue = NaN; fStd_A = NaN;
            end
        else        % Set standard deviation of a-value to nan;
            fStd_A = NaN;   fStd_Mc = NaN;
        end
    else
        fMc = NaN; fStd_Mc = NaN; fBValue = NaN; fStd_B = NaN; fAValue = NaN; fStdDevB = NaN; fStdDevMc = NaN;
    end

    % calculate limits of line to plot for b value line
    % For ZMAP
    magco = fMc;
    index_low = find(xt3 < magco+.05 & xt3 > magco-.05);

    mag_hi = xt3(1);
    mz = xt3 <= mag_hi & xt3 >= magco-.0001;
    mag_zone = xt3(mz);

    % PLOTS an 'x' in the point of Mc
    figure(bfig);       % Lets us be sure that it will polt in the due figure
    semilogy(xt3(index_low),bvalsum3(index_low)*1.5,'vb','LineWidth',1.0,'MarkerSize',7);
    text(xt3(index_low)+0.2,bvalsum3(index_low)*1.5,'Mc','FontWeight','bold','FontSize',fs10,'Color','b');

    % Set to correct method, maximum like or least squares
    if (ud.bBst_button == 0)
        sol_type = 'Maximum Likelihood Solution';
    else
        sol_type = 'Maximum Likelihood Estimate, Uncertainties by bootstrapping';
    end;
    bw = fBValue;   ud2.bw = bw;    %bvml;
    aw = fAValue;   ud2.aw = aw;    %avml;
    ew = fStd_B;                    %stanml;
    
    % create and draw a line corresponding to the b value
    p = [-1*bw aw];
    f = polyval(p,mag_zone);
    f = 10.^f;
    hold on
    semilogy(mag_zone,f,'r','LineWidth',1);           % plot linear fit to backg
    std_backg = ew;      % standard deviation of fit

    set(axHnd,'XLim',[min(ud.events_mag)-0.5  max(ud.events_mag)+0.5],'YLim',[0.9 length(ud.events_time+30)*2.5])

    tt1 = num2str(bw,3);
    tt2 = num2str(std_backg,1);
    tmc = num2str(magco,2);

    h2 = axes('position',[0 0 1 1]);
    set(h2,'visible','off');

    a0 = aw - log10(max(ud.events_time)-min(ud.events_time));

    if (ud.bBst_button == 0)
        text(.16, .06,['b-value = ',tt1,' +/- ',tt2,',  a value = ',num2str(aw,3) ', a value (annual) = ',...
                num2str(a0,3)],'FontSize',fs10,'FontWeight','normal');
        set(bfig,'PaperPosition',[0.5 0.5 4.0 5.5])
        text(.16, .09,sol_type,'FontSize',fs10 );
        text(.16, .03,['Magnitude of Completeness = ',tmc],'FontSize',fs10);
    else
        txt1=text(.16, .06,['b-value = ',num2str(round(100*fBValue)/100),' +/- ',num2str(round(100*fStd_B)/100),...
                ',  a value = ',num2str(aw,3) ',  a value (annual) = ', num2str(a0,3)],'FontSize',fs10);
        set(txt1,'FontWeight','normal')
        set(bfig,'PaperPosition',[0.5 0.5 4.0 5.5])
        text(.16, .09,sol_type,'FontSize',fs10 );
        text(.16, .03,['Magnitude of Completeness = ',tmc ' +/- ', num2str(round(100*fStd_Mc)/100)],'FontSize',fs10);
    end;

    set(bfig,'UserData',ud2,'pointer','arrow');
end

% -------------------------------------------------------------------------------
function edit_bootstr(obj,eventdata,ud)
ud = get(ud.nPar,'UserData');          % Get the updated version
ud.nBstSample = str2double(get(ud.field1,'String'));
set(ud.nPar,'UserData',ud)

% -------------------------------------------------------------------------------
function edit_mcCorr(obj,eventdata,ud)
ud = get(ud.nPar,'UserData');          % Get the updated version
ud.fMccorr = str2double(get(ud.field2,'String'));
set(ud.popup,'Value',1);
set(ud.nPar,'UserData',ud)

% -------------------------------------------------------------------------------
function button_go(obj,eventdata,ud)
ud = get(ud.nPar,'UserData');          % Get the updated version
ud.inpr1 = get(ud.popup,'Value');
ud.bBst_button = get(ud.bBst_button,'Value');
set(ud.nPar,'UserData',ud);
delete(ud.nPar);    pause(0.05)
mc_e_b_estimate('ca',ud)

% -------------------------------------------------------------------------------
function calSave(obj,eventdata,ud)
% Save G & R curve in file
	ud = get(ud.h_fig,'UserData');          % Get the updated version
	handles_mir = guidata(ud.h_mir_fig);       % Get the Mirone handles structure

	[FileName,PathName] = put_or_get_file(handles_mir, ...
		{'*.dat;*.DAT', 'Gutt Rich file (*.dat,*.DAT)'; '*.*', 'All Files (*.*)'},'Select File name','put','.dat');
	if isequal(FileName,0),		return,		end
	fname = [PathName FileName];


	fid = fopen(fname, 'w');
	if (fid < 0),    errordlg(['Can''t open file:  ' fname],'Error');    return;     end
	fprintf(fid,'%.1f\t%d\n',fliplr([ud.xt3; ud.bvalsum3]));    % Flip because Mags were in decreasing order
	fclose(fid);

% -------------------------------------------------------------------------------
function plorem(obj,eventdata,ud)
ud = get(ud.h_fig,'UserData');          % Get the updated version

tr2 = [];   tr2u = [];  tr2l = [];
si = [];    % Note: this is one more of those variables comming out from nowhere. I don't know what it should
            % contain, but the code works(?) without it. I'll let it here in case I will ever find out what it is.
t0b = min(ud.events_time);    teb = max(ud.events_time);
for (m = max(ud.events_mag)-1:0.1:max(ud.events_mag)+2)
    tr = (teb-t0b)  ./ (10^(ud.aw-ud.bw*m));
    %tru = (teb-t0b) ./ (10^(ud.aw-(ud.bw+si)*m));       % With a si =[] those 3 are equal ??
    %trl = (teb-t0b) ./ (10^(ud.aw-(ud.bw-si)*m));
    tr2 = [tr2 ; tr  m];                % This the only one that is actualy used. Big code mess.
    %tr2u = [tr2u ; tru  m];
    %tr2l = [tr2l ; trl  m];
end

h1 = figure('NumberTitle','off','color','w','Name','Recurrence Time');
pl = plot(tr2(:,2),tr2(:,1),'k','LineWidth',2);

set(gca,'box','on','TickDir','out','FontWeight','bold','Yscale','log','XGrid','on','YGrid','on',...
    'FontSize',12,'Linewidth',1.2,'Ticklength',[0.02 0.02],'color','w')
xlabel('Magnitude');        ylabel('Recurrence Time [yrs]')

h2 = figure('NumberTitle','off','color','w','Name','Annual Probability');
pl =  plot(tr2(:,2),1./tr2(:,1),'k','LineWidth',2);
set(gca,'box','on','TickDir','out','FontWeight','bold','Yscale','log','XGrid','on','YGrid','on',...
    'FontSize',12,'Linewidth',1,'Ticklength',[0.02 0.02],'color','w')
xlabel('Magnitude');        ylabel('Annual Probability')

% -------------------------------------------------------------------------------
function mcperc_ca3(events_mag)
% This is a comleteness determination test 

[bval,xt2] = histo_m('hist',events_mag,-2:0.1:6);
l = (max(find(bval == max(bval)))); 
magco0 = xt2(l);
dat = []; 

for (i = magco0-0.5:0.1:magco0+0.7)
   l = (events_mag >= i - 0.0499);      nu = length(events_mag(l));
   if (length(events_mag(l)) >= 25) 
      [mw bv2 stan2 av] = bmemag(events_mag(l));
      res2 = synthb_aut(events_mag,bv2,i,l);
      dat = [dat; i res2];
   else
      dat = [dat; i NaN];
   end
end

j =  min(find(dat(:,2) < 10 )); 
if (isempty(j) == 1),    Mc90 = NaN;
else                    Mc90 = dat(j,1); 
end
   
j =  min(find(dat(:,2) < 5 )); 
if (isempty(j) == 1),    Mc95 = NaN;
else                    Mc95 = dat(j,1); 
end

j =  min(find(dat(:,2) < 10 ));
if (isempty(j) == 1),       j = min(find(dat(:,2) < 15 ));      end
if (isempty(j) == 1),       j = min(find(dat(:,2) < 20 ));      end
if (isempty(j) == 1),       j = min(find(dat(:,2) < 25 ));      end
j2 = min(find(dat(:,2) == min(dat(:,2)) )); 

Mc = dat(j,1); 
magco = Mc; 
if (isempty(magco))
    magco = NaN;    prf = 100 -min(dat(:,2));
end

% -------------------------------------------------------------------------------
function res = synthb_aut(events_mag,bv2,i,l)
%This is synthetic219
%This program generates a synthetic catalog of given total number of events, b-value, minimum magnitude, 
%and magnitude increment.
% Yuzo Toya 2/1999

TN = length(events_mag);  %total number of events
B  = bv2 ;      %b-value
IM = i;         %starting magnitude (hypothetical Mc)
inc = 0.1;      %magnitude increment

% log10(N)=A-B*M
M=IM:inc:15;
N=10.^(log10(TN)-B*(M-IM));
aval=(log10(TN)-B*(0-IM));
N=round(N);

new = repmat(NaN,TN,1);

ct1  = min(find(N == 0)) - 1;
if isempty(ct1) == 1 ; ct1 = length(N); end

ctM=M(ct1);
count=0;
ct=0;
for (I=IM:inc:ctM)
   ct=ct+1;
   if I~=ctM;
      for sc=1:(N(ct)-N(ct+1));
         count=count+1;
         new(count)=I;
      end
   else   
      count=count+1;
      new(count)=I;
   end
end

PM = M(1:ct);
%PN = log10(N(1:ct));
N = N(1:ct); 
%le = length(events_mag(l));
[bval,xt2] = histo_m('hist',events_mag(l),PM);
b3 = fliplr(cumsum(fliplr(bval)));    % N for M >= (counted backwards)
res = sum(abs(b3 - N))/sum(b3)*100; 

% -------------------------------------------------------------------------------
function [mean_m1, b1, sig1, av2] =  bmemag(events_mag)
% function calculates the mean magnitute, the b value based 
% on the mean and the standart deviation 
% Stefan Wiemer 03/95

mima = min(events_mag);
if mima > 0;    mima = 0;   end

% calculate the mean magnitude, b(mean) and std 
%
n = length(events_mag);
mean_m1 = mean(events_mag);
b1 = (1 / (mean_m1-min(events_mag-0.05))) * log10(exp(1));      % VERIFICAR SE min(events_mag-0.05) == mima - .05
sig1 = (sum((events_mag-mean_m1).^2)) / (n*(n-1));
sig1 = sqrt(sig1);
sig1 = 2.30 * sig1 * b1^2;            % standard deviation
av2 = log10(n)+b1 * mima;

% -------------------------------------------------------------------------------
function [fMeanMag, fBValue, fStdDev, fAValue] =  calc_bmemag(events_mag, fBinning)
% Calculates the mean magnitute, the b-value based 
% on the maximum likelihood estimation, the a-value and the 
% standard deviation of the b-value
%
% Input parameters:
%   events_mag      Magnitude of events
%   fBinning        Binning of the earthquake magnitudes (default 0.1)
%
% Output parameters:
%   fMeanMag        Mean magnitude
%   fBValue         b-value
%   fStdDev         Standard deviation of b-value
%   fAValue         a-value
%
% Danijel Schorlemmer
% June 2, 2003

% Set the default value if not passed to the function
if (nargin == 1),   fBinning = 0.1;     end

% Calculate the minimum and mean magnitude, length of catalog
nLen = length(events_mag);
fMinMag = min(events_mag);
fMeanMag = mean(events_mag);
% Calculate the b-value (maximum likelihood)
fBValue = (1/(fMeanMag-(fMinMag-(fBinning/2)))) * log10(exp(1));   
% Calculate the standard deviation 
fStdDev = (sum((events_mag-fMeanMag).^2))/(nLen*(nLen-1));
fStdDev = 2.30 * sqrt(fStdDev) * fBValue^2;            
% Calculate the a-value
fAValue = log10(nLen) + fBValue * fMinMag;

%--------------------------------------------------------------------------------
function [fMc, fStd_Mc, fBvalue, fStd_B, fAvalue, fStd_A, vMc, mBvalue] = ...
    calc_McBboot(events_mag, fBinning, nSample, nMethod, nMinNum, fMcCorr)
% Calculate Mc, b-value and their uncertainties from bootstrapping
% Mc and b are actually the mean values of the empirical distribution
%
% Input parameters
% events_mag  : Magnitude of events
% fBinning    : Binning interval
% nSample     : Number of bootstraps to determine Mc
% nMethod     : Method to determine Mc (see calc_Mc)
% nMinNum     : Minimum number of events above Mc to calculate 
% fMcCorr     : Correction term for Mc
%
% Output parameters
% fMc     : Mean Mc of the bootstrap
% fStd_Mc : 2nd moment of Mc distribution
% fBvalue : Mean b-value of the bootstrap
% fStd_B  : 2nd moment of b-value distribution
% fAvalue : Mean a-value of the bootstrap
% fStd_A  : 2nd moment of a-value distribution
% vMc     : Vector of Mc-values
% mBvalue : Matrix of [fMeanMag fBvalue fStdDev fAvalue]
%
% J. Woessner 
% last update: 22.04.04

% Check input variables
if (nargin < 5)
    nMinNum = 50;    fMcCorr = 0;
elseif (nargin == 5)
    fMcCorr = 0;
end

% Initialize
vMc = [];   mBvalue = [];

% Reset randomizer
rand('state',sum(100*clock));
% Create bootstrap samples using bootstrap matlab toolbox
mMag_bstsamp = bootrsp(events_mag,nSample);
% Determine Mc uncertainty
for (nSamp=1:nSample)
    events_mag = mMag_bstsamp(:,nSamp);
    [fMc] = calc_Mc(events_mag, nMethod, fBinning, fMcCorr);
    vMc =  [vMc; fMc];
    % Select magnitude range and calculate b-value
    vSel = (events_mag >= fMc-fBinning/2);
    mCat = events_mag(vSel);
    % Check static for length of catalog
    nY = length(mCat);
    if (nY >= nMinNum)
        [fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
    else
        fMeanMag = NaN;        fBvalue = NaN;        fStdDev = NaN;        fAvalue = NaN;
    end;
    mBvalue = [mBvalue; fMeanMag fBvalue fStdDev fAvalue];
    
end;

% Calculate mean Mc and standard deviation
fStd_Mc = calc_StdDev(vMc);
fMc = nanmean(vMc);

% Calculate mean b-value and standard deviation
fStd_B = calc_StdDev(mBvalue(:,2));
fStd_A = calc_StdDev(mBvalue(:,4));
fBvalue = nanmean(mBvalue(:,2));
fAvalue = nanmean(mBvalue(:,4));

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
if (~isreal(fStdDev) && ~isnan(fStdDev)),    fStdDev = 0;     end

% --------------------------------------------------------------------
function y = nanmean(x)
% Average or mean ignoring NaNs.

	if (isempty(x)),     y = NaN;    return;     end     % Check for empty input.

	% Replace NaNs with zeros.
	nans = isnan(x);
	i = find(nans);
	x(i) = zeros(size(i));

	% count terms in sum over first non-singleton dimension
	dim = find(size(x)>1);
	if (isempty(dim)),  dim = 1;
	else				dim = dim(1);
	end
	count = sum(~nans,dim);

	% Protect against a column of all NaNs
	i = find(count==0);
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

% -----------------------------------------------
function out = bootrsp(in,B)
%   out=bootrsp(in,B)
%    
%   Bootstrap  resampling  procedure. 
%
%     Inputs:
%        in - input data 
%        B  - number of bootstrap resamples (default B=1)        
%     Outputs:
%       out - B bootstrap resamples of the input data  
%
%   For a vector input data of size [N,1], the  resampling 
%   procedure produces a matrix of size [N,B] with columns 
%   being resamples of the input vector.
%
%   For a matrix input data of size  [N,M], the resampling
%   procedure produces a 3D matrix of  size  [N,M,B]  with 
%   out(:,:,i), i = 1,...,B, being a resample of the input 
%   matrix.
%
%   Example:
%
%   out=bootrsp(randn(10,1),10);

%  Created by A. M. Zoubir and D. R. Iskander  May 1998
%
%  References:
%  Efron, B.and Tibshirani, R.  An Introduction to the Bootstrap.
%               Chapman and Hall, 1993.
%
%  Zoubir, A.M. Bootstrap: Theory and Applications. Proceedings 
%               of the SPIE 1993 Conference on Advanced  Signal 
%               Processing Algorithms, Architectures and Imple-
%               mentations. pp. 216-235, San Diego, July  1993.
%
%  Zoubir, A.M. and Boashash, B. The Bootstrap and Its Application
%               in Signal Processing. IEEE Signal Processing Magazine, 
%               Vol. 15, No. 1, pp. 55-76, 1998.

if (nargin == 0),    error('Provide input data'); end;
if (nargin == 1),    B = 1;  end;

s = size(in);     
if (length(s) > 2) 
    error('Input data can be a vector or a 2D matrix only'); 
end
if (min(s) == 1)  
    out = in(ceil(max(s)*rand(max(s),B)));    
else         
    out = in(ceil(s(1)*s(2)*rand(s(1),s(2),B))); 
end

% ----------------------------------------------------------------------------------
function [fMc] = calc_Mc(mCatalog, nMethod, fBinning, fMcCorrection)
% Calculates the magnitude of completeness for a given catalog
%
% Input parameters:
%   mCatalog       Earthquake catalog for determing the magnitude of completeness
%   nMethod        Method to determine the magnitude of completeness
%                  1: Maximum curvature
%                  2: Fixed Mc = minimum magnitude (Mmin)
%                  3: Mc90 (90% probability)
%                  4: Mc95 (95% probability)
%                  5: Best combination (Mc95 - Mc90 - maximum curvature)
%         ==>      6: Mc using EMR-method       REMOVED THIS OPTION BECAUSE IT DEPENDS ON THE OPTIMISATION TBX
%                  7: Mc due b using Shi & Bolt uncertainty
%                  8: Mc due b using bootstrap uncertainty
%                  9: Mc due b Cao-criterion
%   fBinning       Binning of catalog's magnitudes (default 0.1)
%   fMcCorrection  Correction term to be added to fMc (default 0)
%
% Output parameters:
%   fMc            Magnitude of completeness
%
% Special function
%   If called without any parameters, calc_Mc returns a string containing the names 
%   of all available Mc-determination routines
%
% Danijel Schorlemmer, Jochen Woessner: April 27, 2004

if (nargin == 0)
  fMc = ['1: Maximum curvature|' ...
         '2: Fixed Mc = minimum magnitude (Mmin)|' ...
         '3: Mc90 (90% probability)|' ...
         '4: Mc95 (95% probability)|' ...
         '5: Best combination (Mc95 - Mc90 - maximum curvature)|' ...
         '6: Mc due b using Shi & Bolt uncertainty|' ...
         '7: Mc due b using bootstrap uncertainty|' ...
         '8: Mc due b Cao-criterion'];
  return;
end;

% Check input variables
if (nargin < 3)
    fBinning = 0.1;     fMcCorrection = 0;
elseif (nargin == 4)
    fMcCorrection = 0;
end

% Init return variable
fMc = NaN;

if (nMethod == 1)               % Maximum curvature
    fMc = calc_McMaxCurvature(mCatalog); 
elseif (nMethod == 2  )         % Fixed Mc (Mc = Mmin)
    fMc = min(mCatalog); 
elseif (nMethod == 3)           % Automatic Mc90 
    [fDummy, fDummy, fMc] = calc_McBest(mCatalog, fBinning);  
elseif (nMethod == 4)           % Automatic Mc95 
    [fDummy, fMc, fDummy] = calc_McBest(mCatalog, fBinning);  
elseif (nMethod == 5  )         % Best combination (Mc95 - Mc90 - maximum curvature)
    [fMc, Mc95, Mc90] = calc_McBest(mCatalog, fBinning);  
    if (~isnan(Mc95)),          fMc = Mc95; 
    elseif (~isnan(Mc90)),     fMc = Mc90; 
    else                        fMc = calc_McMaxCurvature(mCatalog);
    end
elseif (nMethod == 6)           % Mc due b using Shi & Bolt uncertainty
    [fMc_shi, fBvalue_shi, fBStd_shi, fAvalue_shi, mBave] = calc_Mcdueb(mCatalog);
    fMc = fMc_shi;
elseif (nMethod == 7)           % Mc due b using bootstrap uncertainty
    nSample = 500;
    [fMc_bst, fBvalue_bst, fBStd_bst, fAvalue_bst, mBave] = calc_McduebBst(mCatalog, fBinning, 5, 50,nSample);
    fMc = fMc_bst;
else % nMethod == 8             % Mc due b Cao-criterion
    [fMc_cao, fBvalue_cao, fBStd_cao, fAvalue_cao] = calc_McduebCao(mCatalog);
    fMc = fMc_cao;
end

if (isempty(fMc)),   fMc = NaN;  end                     % Check fMc  
if (~isnan(fMc)),    fMc = fMc + fMcCorrection;  end     % Apply correction

%--------------------------------------------------------------------------------
function [fMc] = calc_McMaxCurvature(mCatalog)
% Determines the magnitude of completeness at the point of maximum
%   curvature of the frequency magnitude distribution
%
% Input parameter:
%   mCatalog        Earthquake catalog
%
% Output parameter:
%   fMc             Magnitude of completeness, nan if not computable
%
% Danijel Schorlemmer
% November 7, 2001

try         % Get maximum and minimum magnitudes of the catalog
    fMaxMagnitude = max(mCatalog);
    fMinMagnitude = min(mCatalog);
    if (fMinMagnitude > 0),    fMinMagnitude = 0;    end
  
    nNumberMagnitudes = round(fMaxMagnitude*10) + 1;       % Number of magnitudes units
  
    % Create a histogram over magnitudes
    [vHist, vMagBins] = hist(mCatalog, (fMinMagnitude:0.1:fMaxMagnitude));
  
    % Get the points with highest number of events -> maximum curvature  
    fMc = vMagBins(max(find(vHist == max(vHist))));
    if (isempty(fMc)),   fMc = nan;      end  
catch
    fMc = nan;
end

%--------------------------------------------------------------------------------
function [fMc, fMc95, fMc90] = calc_McBest(mCatalog, fBinning)
% First estimation of magnitude of completeness (maximum curvature) 
[vEvents, vMag] = hist(mCatalog, -2:0.1:6);
nSel = max(find(vEvents == max(vEvents))); 
fMcStart = vMag(nSel);

mData = [];     % Data container

% Magnitude increment
if (nargin == 1),    fBinning = 0.1;     end

for (nCnt = (fMcStart - 0.9):fBinning:(fMcStart + 1.5))
    vSel = (mCatalog > (nCnt - (fBinning/2)));
    nNumberEvents = length(mCatalog(vSel));
    if (nNumberEvents >= 25)
        [fDummy fBValue fDummy fDummy] =  bmemag(mCatalog(vSel));
        fStartMag = nCnt;           % Starting magnitude (hypothetical Mc)
   
        % log10(N)=A-B*M
        vMag = (fStartMag:fBinning:15); % Ending magnitude must be sufficiently high
        vNumber = 10.^(log10(nNumberEvents)-fBValue*(vMag - fStartMag));
        vNumber = round(vNumber);
    
        % Find the last bin with an event
        nLastEventBin = min(find(vNumber == 0)) - 1;
        if (isempty(nLastEventBin)),    nLastEventBin = length(vNumber);    end
    
        % Determine set of all magnitude bins with number of events > 0
        ct = round((vMag(nLastEventBin)-fStartMag)*(1/fBinning) + 1);
    
        PM=vMag(1:ct);
        vNumber = vNumber(1:ct); 
        [bval, vDummy] = hist(mCatalog(vSel),PM);
        b3 = fliplr(cumsum(fliplr(bval)));    % N for M >= (counted backwards)
        res2 = sum(abs(b3 - vNumber))/sum(b3)*100; 
        mData = [mData; nCnt res2];
    else     
        mData = [mData; nCnt nan];
    end
end

% Evaluation of results

% Is fMc90 available
nSel = min(find(mData(:,2) < 10)); 
if (isempty(nSel)),  fMc90 = nan;
else                fMc90 = mData(nSel,1); 
end

% Is fMc95 available
nSel = min(find(mData(:,2) < 5)); 
if (isempty(nSel)),      fMc95 = nan;
else                    fMc95 = mData(nSel,1); 
end

% ?????
j =  min(find(mData(:,2) < 10 ));
if (isempty(j)),    j = min(find(mData(:,2) < 15 )); end
if (isempty(j)),    j = min(find(mData(:,2) < 20 )); end
if (isempty(j)),    j = min(find(mData(:,2) < 25 )); end

fMc = mData(j,1);  
if (isempty(fMc)),   fMc = nan;  end

%--------------------------------------------------------------------------------
function [fMc, fBvalue, fBStd, fAvalue, mBave] = calc_Mcdueb(mCatalog, fBinning, nWindowSize, nMinNumberEvents)
% Calculate Mc using the function b-value vs. cut-off-magnitude
% Decision criterion for b and Mc: b_i-std_Shi(b_i) <= b_ave <= b_i+std_Shi(b_i)
%
% Relevant reference: Cao A., Gao, S.S., Temporal variation of seismic b-values 
% beneath northeastern Japan island arc, GRL, 29, 9, 2002
%
% Incoming variables:
% mCatalog         : EQ catalog
% fBinning         : Bin size
% nWindowSize      : Window size
% nMinNumberEvents : Minimum number of events
%
% Outgoing variables:
% fMc              : Magnitude of completeness
% fBStd            : Shi & Bolt deviation for b
% fBvalue          : b-value
% fAvalue          : a-value
%
% Author: J. Woessner
% last update: 04.06.03

% Check input
if nargin == 0, error('No catalog input'); end;
if nargin == 1, fBinning = 0.1; nWindowSize = 5; nMinNumberEvents = 50;     end;
if nargin == 2, nWindowSize = 5; nMinNumberEvents = 50; end
if nargin == 3, nMinNumberEvents = 50;  end;
if nargin > 4, disp('Too many arguments!'), return; end;

% Initialize
fMc = nan;      fBvalue = nan;      mBave = [];     mMcBA = [];

% Calculate b-with magnitude
[mBvalue] = calc_bwithmag(mCatalog, fBinning, nMinNumberEvents);
% Remove NANs
vSel = isnan(mBvalue(:,1));
mBvalue = mBvalue(~vSel,:);

% Use Shi & Bolt uncertainty to decide for Mc 
for nStep = 1:(length(mBvalue(:,1))-nWindowSize)
    fBave = mean(mBvalue(nStep:nStep+nWindowSize,1));
    mBave = [mBave; fBave mBvalue(nStep,1) mBvalue(nStep,2) mBvalue(nStep,3) mBvalue(nStep,4)];
    % Criterion: If fBave is in in between the error estimate of the b-value of the first cut-off magnitude
    % take it as guess
    if (fBave >= mBvalue(nStep,1)-mBvalue(nStep,2) && fBave <= mBvalue(nStep,1)+mBvalue(nStep,2))
        mMcBA = [mMcBA; fBave mBvalue(nStep,1) mBvalue(nStep,2) mBvalue(nStep,3) mBvalue(nStep,4)];
    end;
end;
% Create output
try
    fMc = mMcBA(1,5);    fBvalue = mMcBA(1,2);    fAvalue = mMcBA(1,4);    fBStd = mMcBA(1,3);
catch
    fMc = nan;    fBvalue = nan;    fAvalue = nan;    fBStd = nan;
end

%--------------------------------------------------------------------------------
function [fMc, fBvalue, fBStd, fAvalue, fSigmaLow, fSigmaHi, mBave, mBvalue] = ...
    calc_McduebBst(mCatalog, fBinning, nWindowSize, nMinNumberEvents, nSample)
% Calculate Mc using the function b-value vs. cut-off-magnitude: Bootstrap approach
% Decision criterion for b and Mc: b_i-std_Bst(b_i) <= b_ave <= b_i+std_Bst(b_i)

% Relevant reference: Cao A., Gao, S.S., Temporal variation of seismic b-values 
% beneath northeastern Japan island arc, GRL, 29, 9, 2002
%
% Incoming variables:
% mCatalog         : EQ catalog
% fBinning         : Bin size
% nWindowSize      : Window size
% nMinNumberEvents : Minimum number of events
% nSample          : Number of bootstrap samples
%
% Outgoing variables:
% fMc              : Magnitude of completeness
% fBvalue          : b-value
% fBStd            : 2nd moment of b-value-distribution (comparable to standard deviation) 
% fAvalue          : a-value
% fSigmaLow        : 16-percentile of b-value distribution
% fSigmaHi         : 84-percentile of b-value distribution
% mBave            : Result matrix for plotting (average values)
% mBvalue          : Result matrix for plotting
% Author: J. Woessner
% last update: 04.06.03

% Check input
if nargin == 0, error('No catalog input'); end;
if nargin == 1, fBinning = 0.1; nWindowSize = 5; nMinNumberEvents = 50; nSample = 100; 
    disp('Default Bin size: 0.1, Windowsize = 5, Minimum number of events: 50, Bootstrap samples = 100');end;
if nargin == 2, nWindowSize = 5; nMinNumberEvents = 50; nSample = 100; 
    disp('Default Windowsize = 5, Minimum number of events: 50, Bootstrap samples = 100');end;
if nargin == 3, nMinNumberEvents = 50; nSample = 100; disp('Default Bootstrap samples = 100');end;
if nargin == 4, nSample = 100; disp('Default Minimum number of events: 50, Bootstrap samples = 100');end;
if (nargin > 5),     disp('Too many arguments!'), return; end;

% Initialize
mBvalue = [];   mBvalue_bst = [];   mBave = []; mMcBA = [];

% Set fix values
fMinMag = min(mCatalog);
fMaxMag = max(mCatalog);

% Create bootstrap samples using bootstrap matlab toolbox
mMag_bstsamp = bootrsp(mCatalog,nSample);

% Calculate b-with magnitude
hWaitbar1 = aguentabar('title','Bootstaping ...');
for (fMag = fMinMag:fBinning:fMaxMag)
	for (nSamp = 1:nSample)
		mCatalog = mMag_bstsamp(:,nSamp);
		% Select magnitude range
		vSel = (mCatalog >= fMag-0.05);
		mCat = mCatalog(vSel);
		% Check for minimum number of events
		if (length(mCat) >= nMinNumberEvents)
			try
				[fMeanMag, fBValue, fStdDev, fAValue] =  calc_bmemag(mCat, fBinning);
				mBvalue_bst = [mBvalue_bst; fBValue fStdDev fAValue fMag];
			catch
				mBvalue_bst = [mBvalue_bst; nan nan nan fMag];
			end;
		else
			mBvalue_bst = [mBvalue_bst; nan nan nan fMag];
		end		% END of IF
	end			% END of FOR nSamp

	% Check for Nan and create output for [16 84]-percentile
	vSel = isnan(mBvalue_bst(:,1));
	mBvalue_bst_tmp = mBvalue_bst(~vSel,:);
	if (~isempty(mBvalue_bst_tmp(:,1)) && length(mBvalue_bst_tmp(:,1)) > 1)
		vSigma = prctile(mBvalue_bst_tmp(:,1),[16 84]);
	elseif (~isempty(mBvalue_bst_tmp(:,1)) && length(mBvalue_bst_tmp(:,1)) == 1)
		vSigma = prctile(mBvalue_bst_tmp(:,1),[16 84]);
		vSigma = vSigma';
	else
		vSigma = [nan nan];
	end;
	% Calculate 2nd moment
	if (~isempty(mBvalue_bst_tmp(:,1)))
		fStdBst = calc_StdDev(mBvalue_bst_tmp(:,1));
	else
		fStdBst = nan;
	end

	try             % mBvalue: b std_bolt(b) a Mc 16-perc 18-perc std(b_2nd moment)
		mBvalue = [mBvalue; nanmean(mBvalue_bst) vSigma fStdBst];
	catch
		mBvalue = [mBvalue; nan nan nan nan nan nan nan];
	end;
	mBvalue_bst =[];
	aguentabar(fMag/fMaxMag)
end     % END of FOR fMag
if (ishandle(hWaitbar1))	close(hWaitbar1),	end

% Use bootstrap percentiles to decide for Mc 
for (nStep = 1:(length(mBvalue(:,1))-nWindowSize))
    fBave = mean(mBvalue(nStep:nStep+nWindowSize,1));
    mBave = [mBave; fBave mBvalue(nStep,:)];
    % Criterion: If fBave is in in between the error estimate of the b-value of the first cut-off magnitude
    % take it as guess
    if (fBave >= mBvalue(nStep,5) && fBave <= mBvalue(nStep,6))
        mMcBA = [mMcBA; fBave mBvalue(nStep,:)];
    end;
end;

% Create output
try
    fMc = mMcBA(1,5);       fBvalue = mMcBA(1,2);       fAvalue = mMcBA(1,4);
    fBStd = mMcBA(1,8);     fSigmaLow = mMcBA(1,6);     fSigmaHi = mMcBA(1,7);
catch
    fMc = nan;    fBvalue = nan;    fAvalue = nan;    fBStd = nan;    fSigmaLow = nan;    fSigmaHi = nan;
end

%--------------------------------------------------------------------------------
function [fMc, fBvalue, fBStd, fAvalue] = calc_McduebCao(mCatalog, fBinning, nMinNumberEvents)
% Calculate Mc using the function b-value vs. cut-off-magnitude
% Decision criterion for Mc and b is: b_i - b_i-1 <= 0.03 as in reference
% 
% Reference: Cao A., Gao, S.S., Temporal variation of seismic b-values 
% beneath northeastern Japan island arc, GRL, 29, 9, 2002
%
% Incoming variables:
% mCatalog         : EQ catalog
% fBinning         : Bin size
% nMinNumberEvents : Minimum number of events
%
% Outgoing variables:
% fMc              : Magnitude of completeness
% fBStd            : Shi & Bolt standard deviation of b
% fBvalue          : b-value
% fAvalue          : a-value
%
% Author: J. Woessner
% last update: 04.06.03

% Check input
if nargin == 0, error('No catalog input'); end;
if nargin == 1, fBinning = 0.1; nMinNumberEvents = 50;  end;
if nargin == 2, nMinNumberEvents = 50;  end;
if nargin > 3, error('Too many arguments!'); end;

% Initialize
fMc = nan;  fBvalue = nan;  mMcBA = [];

% Calculate b-with magnitude
[mBvalue] = calc_bwithmag(mCatalog, fBinning, nMinNumberEvents);
% Remove NANs
vSel = isnan(mBvalue(:,1));
mBvalue = mBvalue(~vSel,:);

% Use Shi & Bolt uncertainty to decide for Mc 
for nStep = 2:(length(mBvalue(:,1)))
    % Criterion: If bi+1 - bi < 0.03, then use bi as b-value and cut-off magnitude as Mc
    if (mBvalue(nStep,1)- mBvalue(nStep-1,1)<= 0.03)
        mMcBA = [mMcBA; mBvalue(nStep,1) mBvalue(nStep,2) mBvalue(nStep,3) mBvalue(nStep,4)];
    end
end
% Create output
try
    fMc = mMcBA(1,4);    fBvalue = mMcBA(1,1);    fAvalue = mMcBA(1,3);    fBStd = mMcBA(1,2);
catch
    fMc = nan;    fBvalue = nan;    fAvalue = nan;    fBStd = nan;
end

%--------------------------------------------------------------------------------
function y = prctile(x,p)
% Percentiles of the sample in X.
%   Y = PRCTILE(X,P) returns a value that is greater than P percent
%   of the values in X. For example, if P = 50  Y is the median of X. 
%
%   P may be either a scalar or a vector.  If X is a matrix, the ith
%   row of Y is the P(i) percentile of each column of X.  If X is a
%   vector, Y has the same shape as P.

[prows pcols] = size(p);
if (prows ~= 1 && pcols ~= 1)
    error('P must be a scalar or a vector.');
end
if (any(p > 100) || any(p < 0))
    error('P must take values between 0 and 100');
end

if (~any(isnan(x)))
   y = prctilecol(x,p);
else                    % if there are NaNs, process each column
   if (size(x,1) == 1),  x = x';   end
   c = size(x,2);
   np = length(p);
   y = zeros(np,c);
   for (j=1:c)
      xx = x(:,j);
      xx = xx(~isnan(xx));
      if (isempty(xx)),  z = NaN;
      else                  z = prctilecol(xx,p);
      end
      y(:,j) = z(:);
   end
   
   % For a vector x, the orientation of y comes from p
   if (min(size(x)) == 1 && prows==1 && pcols > 1),     y = y';   end
end
      
function y = prctilecol(x,p)
xx = sort(x);
[m,n] = size(x);

if (m == 1 || n == 1)
    m = max(m,n);
	if (m == 1)
	   y = x*ones(length(p),1);
	   return;
	end
    q = 100*(0.5:m - 0.5)./m;
    xx = [min(x); xx(:); max(x)];
else
    q = 100*(0.5:m - 0.5)./m;
    xx = [min(x); xx; max(x)];
end
y = interp1([0 q 100],xx,p);

%--------------------------------------------------------------------------------
function [mBvalue] = calc_bwithmag(mCatalog, fBinning, nMinNumberevents)
% Calculate b-value depending on cut-off magnitude
%
% Incoming variables:
% mCatalog : Earthquake catalog
% fBinnig  : Binning interval
% nMinNumberevents : Minimum number of events
%
% Outgoing variables:
% mBvalue(:,1) : b-values ascending with magnitude
% mBvalue(:,2) : Standard deviation of b (Shi & Bolt, 1982) ascending with magnitude
% mBvalue(:,3) : a-values ascending with magnitude
% mBvalue(:,4) : Ascending magnitudes
%
% Author: J. Woessner
% last update: 04.06.03

% Check input
if nargin == 0, error('No catalog input'); end;
if nargin == 1, fBinning = 0.1;  end;
if nargin > 3, error('Too many arguments!');    end;

% Initialze
mBvalue = [];

% Set fix values
fMinMag = min(mCatalog);
fMaxMag = max(mCatalog);

for fMag=fMinMag:fBinning:fMaxMag
    % Select magnitude range
    vSel = (mCatalog >= fMag-0.05);
    mCat = mCatalog(vSel);
    
    % Check for minimum number of events
    if (length(mCat) >= nMinNumberevents)
        try
            [fMeanMag, fBValue, fStdDev, fAValue] =  calc_bmemag(mCat, fBinning);
            mBvalue = [mBvalue; fBValue fStdDev fAValue fMag];
        catch
            mBvalue = [mBvalue; nan nan nan fMag];
        end
    else
        mBvalue = [mBvalue; nan nan nan fMag];
    end
end
%--------------------------------------------------------------------------------
