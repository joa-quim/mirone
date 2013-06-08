function map = mkpj(n)
% MKPJ Returns a perceptually improved Jet colormap
%   MKPJ (N) returns an Nx3 colormap. 
%   usage: map=mkpj(n);
%
% JUSTIFICATION: rainbow, or spectrum color schemes are considered a poor
% choice for scientific data display by many in the scientific community
% (see for example reference 1 and 3) in that they introduce artifacts 
% that mislead the viewer. "The rainbow color map appears as if it’s separated
% into bands of almost constant hue, with sharp transitions between hues. 
% Viewers perceive these sharp transitions as sharp transitions in the data,
% even when this is not the casein how regularly spaced (interval) data are
% displayed (quoted from reference 1). The default Matlab Jet is no exception.
% Please see examples in my previous FEX submission "Perceptually improved
% colormaps".
% This submission is intended to share the results of my work to create a
% more perceptually balanced, Jet colormap. Please see output arguments section
% for more details.
% A full series of posts on improved colormaps will appear on my blog at
% http://mycarta.wordpress.com/
%
%   arguments: (input)
%   n - scalar specifying number of points in the colorbar. Maximum n=256
%       If n is not specified, defaults to 256
%
%       The new, more perceptually balanced, Jet color map. 
%       This explanation should be read in conjunction with the submitted figure 
%       "JET_idea_after_scaling_final_L_plt.png"
%       The idea to improve the default jet came to me after seeing the human
%       wavelenght discrimination curve from Gregory (1964) in the paper by 
%       Welland et al. (2006).
%       The top display in the figure is the standard Matlab Jet. Below it, 
%       second display, is the Luminance versus wavelenght profile. This is
%       essentially the same as plot you would find in the left plot of 
%       figure 3 in Rogowitz and Kalvin (2001). In the third display I reproduced
%       Gregory's curve. It's considering these two together that brought the
%       Eureka moment (though it took me a few experiments and some calculations
%       to line them up properly). The idea was: "can I use this curve to
%       dynamically stretch the rainbow where transitions are too sharp 
%       (basically around the cyan and yellow), compared to everywhere else?".
%       The answer is yes, it can be done. In practice this amounted to calculate
%       the "inverted" function in the fourth display. Assuming a distance of 1
%       between each of the 256 samples on the x axis in the original Jet colormap,
%       the function was used to resample it at non-integer distances (up
%       to 1.5 where the yellow is, unmodified at 1 where the cyan is, less
%       than 1 everywhere else) resulting in a greater number of samples in
%       the yellow area in the display compared to all the blue areas, with
%       the total number of samples staying at 256. The next step was to force
%       all these new samples back to a distance of 1, achieving a
%       continuously dynamic stretch. The resulting colormap is shown in
%       the fifth display and accompanied by its Luminance versus wavelenght
%       in the last display. This profile is noticeably more perceptually
%       balnced with gentler transitions and a compressive character. 
%       To me this is a very good result, even though it's not perfect. 
%       In an ideal world I would have paired one of the several human wavelenght 
%       discrimination curves found in Wyszecki and Stiles (2000) 
%       with in conjunction with the corresponding spectrum color functions
%       (many are referenced in the FEX submission "Spectral and XYZ Color
%       Functions"), but the problem I had was to identify a matching pair; 
%       a problem that in the end proved insurmountable.
%
%   arguments: (output)
%   map - the output colormap i
%  
%   Example: compare cape topography using the improved Jet vs. default Jet(clipped)
%     %  load cape;
%     %  imagesc(X); colormap(mkpj(128)); colorbar;
%     %  figure;
%     %  imagesc(X); colormap(mkpj(128)); colorbar;
%
%   Acknowledgements
%     For function architecture and code syntax I was inspired by:
%     Light Bartlein Color Maps 
%     www.mathworks.com/matlabcentral/fileexchange/17555
%     (and comments posted therein)
% 
%     A great way to learn more about improved colormaps and making colormaps:
%     MakeColorMap
%     www.mathworks.com/matlabcentral/fileexchange/17552
%     blogs.mathworks.com/videos/2007/11/15/practical-example-algorithm-development-for-making-colormaps/
%
%
%  References
%     1)  Borland, D. and Taylor, R. M. II (2007) - Rainbow Color Map (Still) 
%         Considered Harmful
%         IEEE Computer Graphics and Applications, Volume 27, Issue 2
%         Pdf paper included in submission
% 
%     2)  Kindlmann, G. Reinhard, E. and Creem, S. Face-based Luminance Matching
%         for Perceptual Colormap Generation
%         IEEE - Proceedings of the conference on Visualization '02
%         www.cs.utah.edu/~gk/papers/vis02/FaceLumin.pdf
% 
%     3)  Light, A. and Bartlein, P.J. (2004) - The end of the rainbow? 
%         Color schemes for improved data graphics.
%         EOS Transactions of the American Geophysical Union 85 (40)
%         Reprint of Article with Comments and Reply
%         http://geography.uoregon.edu/datagraphics/EOS/Light-and-Bartlein.pdf
% 
%     4)  Rogowitz, B.E. and  Kalvin, A.D. (2001) - The "Which Blair project":
%         a quick visual method for evaluating perceptual color maps. 
%         IEEE - Proceedings of the conference on Visualization ‘01
%         www.research.ibm.com/visualanalysis/papers/WhichBlair-Viz01Rogowitz_Kalvin._final.pdf
% 
%     5)  Rogowitz, B.E. and  Kalvin, A.D. - Why Should Engineers and Scientists
%         Be Worried About Color?
%         www.research.ibm.com/people/l/lloydt/color/color.HTM
% 
%     6)  Rogowitz, B.E. and  Kalvin, A.D. - How NOT to Lie with Visualization
%         www.research.ibm.com/dx/proceedings/pravda/truevis.htm
%
%     7)  Welland, M., Donnelly, N., and Menneer, T., (2006) - Are we properly
%         using our brains in seismic interpretation? - The Leading Edge; February
%         2006; v. 25; No. 2
%         http://tle.geoscienceworld.org/cgi/content/abstract/25/2/142
%
%     8)  Wyszecki, G. and Stiles W. S. (2000) - Color Science: Concepts and 
%         Methods, Quantitative Data and Formulae, 2nd Edition, John Wiley and Sons
%         http://ca.wiley.com/WileyCDA/WileyTitle/productCd-0471399183.html
%
%     9)  Gregory, R.L. (1966) Eye and Brain: The Psychology of Seeing
%         Fifth edition, 1997 - http://press.princeton.edu/titles/6016.html
% 
%
%  Author: Matteo Niccoli
%  http://mycarta.wordpress.com/
%  matteo@mycarta.ca
%  Release: 1.00
%  Release date: September 28 2011
%
%	As usual reworked for maximum simplicity/speed by Joaquim Luis

	if (nargin == 0)
		n = 256;
	elseif (nargin > 1),
		error('MKJP: wrong number of input args')
	end
	if (n > 256),		error('Maximum number of 256 points for colormap exceeded'),	end

JetI =	[
	0	0.0116	1
	0	0.030375	1
	0	0.0497	1
	0	0.069555	1
	0	0.089918	1
	0	0.11102	1
	0	0.13287	1
	0	0.15503	1
	0	0.17703	1
	0	0.19887	1
	0	0.2208	1
	0	0.24252	1
	0	0.26373	1
	0	0.28435	1
	0	0.30458	1
	0	0.32436	1
	0	0.3436	1
	0	0.36229	1
	0	0.3805	1
	0	0.39825	1
	0	0.41555	1
	0	0.43246	1
	0	0.44906	1
	0	0.46517	1
	0	0.48063	1
	0	0.49535	1
	0	0.50949	1
	0	0.52309	1
	0	0.53616	1
	0	0.54871	1
	0	0.56077	1
	0	0.57239	1
	0	0.58362	1
	0	0.59449	1
	0	0.60501	1
	0	0.61519	1
	0	0.62506	1
	0	0.63463	1
	0	0.64391	1
	0	0.65292	1
	0	0.6617	1
	0	0.67026	1
	0	0.6786	1
	0	0.68674	1
	0	0.6947	1
	0	0.70249	1
	0	0.71009	1
	0	0.71753	1
	0	0.72486	1
	0	0.73209	1
	0	0.7392	1
	0	0.74621	1
	0	0.75314	1
	0	0.76001	1
	0	0.7668	1
	0	0.77352	1
	0	0.78018	1
	0	0.7868	1
	0	0.79337	1
	0	0.79989	1
	0	0.80638	1
	0	0.81285	1
	0	0.8193	1
	0	0.82572	1
	0	0.83212	1
	0	0.83852	1
	0	0.84491	1
	0	0.85128	1
	0	0.85765	1
	0	0.86401	1
	0	0.87036	1
	0	0.87669	1
	0	0.88304	1
	0	0.88941	1
	0	0.89581	1
	0	0.90221	1
	0	0.90865	1
	0	0.91512	1
	0	0.9216	1
	0	0.92811	1
	0	0.93467	1
	0	0.94132	1
	0	0.94806	1
	0	0.95487	1
	0	0.96175	1
	0	0.9687	1
	3.129e-005	0.97685	0.99997
	0.00012415	0.98585	0.99988
	0.00027706	0.99346	0.99972
	0.0004885	0.99744	0.99951
	0.0033792	0.99848	0.99662
	0.010635	0.99928	0.98936
	0.019956	0.9998	0.98004
	0.029032	1	0.97097
	0.037192	1	0.96281
	0.045799	1	0.9542
	0.054762	1	0.94524
	0.063971	1	0.93603
	0.073375	1	0.92662
	0.083057	1	0.91694
	0.093071	1	0.90693
	0.10347	1	0.89653
	0.11427	1	0.88573
	0.12547	1	0.87453
	0.13715	1	0.86285
	0.14939	1	0.85061
	0.16228	1	0.83772
	0.17585	1	0.82415
	0.19009	1	0.80991
	0.20501	1	0.79499
	0.22072	1	0.77928
	0.23739	1	0.76261
	0.25479	1	0.74521
	0.27267	1	0.72733
	0.29097	1	0.70903
	0.30999	1	0.69001
	0.32946	1	0.67054
	0.34906	1	0.65094
	0.36861	1	0.63139
	0.38835	1	0.61165
	0.40818	1	0.59182
	0.42797	1	0.57203
	0.44764	1	0.55236
	0.46733	1	0.53267
	0.48694	1	0.51306
	0.5063	1	0.4937
	0.52534	1	0.47466
	0.54421	1	0.45579
	0.56282	1	0.43718
	0.581	1	0.419
	0.59864	1	0.40136
	0.6159	1	0.3841
	0.63275	1	0.36725
	0.6491	1	0.3509
	0.66491	1	0.33509
	0.68025	1	0.31975
	0.69515	1	0.30485
	0.70959	1	0.29041
	0.72358	1	0.27642
	0.73719	1	0.26281
	0.75039	1	0.24961
	0.76315	1	0.23685
	0.77545	1	0.22455
	0.78732	1	0.21268
	0.79881	1	0.20119
	0.80991	1	0.19009
	0.82062	1	0.17938
	0.83098	1	0.16902
	0.84099	1	0.15901
	0.8507	1	0.1493
	0.86012	1	0.13988
	0.86927	1	0.13073
	0.87816	1	0.12184
	0.88679	1	0.11321
	0.89518	1	0.10482
	0.90333	1	0.096675
	0.91125	1	0.088755
	0.91895	1	0.081048
	0.92646	1	0.073539
	0.93377	1	0.066226
	0.9409	1	0.059102
	0.94785	1	0.052151
	0.95464	1	0.045359
	0.96139	1	0.038605
	0.96805	1	0.031954
	0.9744	1	0.025599
	0.98027	1	0.019732
	0.98654	0.9997	0.01346
	0.99297	0.99892	0.0070311
	0.99798	0.99775	0.00202
	1	0.9963	0
	1	0.99303	0
	1	0.98732	0
	1	0.98067	0
	1	0.97455	0
	1	0.96933	0
	1	0.9642	0
	1	0.95914	0
	1	0.95417	0
	1	0.94928	0
	1	0.94447	0
	1	0.93973	0
	1	0.93506	0
	1	0.93046	0
	1	0.92594	0
	1	0.92147	0
	1	0.91704	0
	1	0.91264	0
	1	0.90828	0
	1	0.90396	0
	1	0.8997	0
	1	0.89552	0
	1	0.89139	0
	1	0.88731	0
	1	0.88326	0
	1	0.87921	0
	1	0.8752	0
	1	0.8712	0
	1	0.86722	0
	1	0.86326	0
	1	0.85932	0
	1	0.8554	0
	1	0.85149	0
	1	0.84759	0
	1	0.84371	0
	1	0.83983	0
	1	0.83596	0
	1	0.83209	0
	1	0.82822	0
	1	0.82435	0
	1	0.82047	0
	1	0.81654	0
	1	0.81259	0
	1	0.80862	0
	1	0.80461	0
	1	0.80058	0
	1	0.79652	0
	1	0.79244	0
	1	0.7883	0
	1	0.78409	0
	1	0.77983	0
	1	0.7755	0
	1	0.77109	0
	1	0.76657	0
	1	0.76196	0
	1	0.75724	0
	1	0.75241	0
	1	0.74742	0
	1	0.74232	0
	1	0.73706	0
	1	0.73162	0
	1	0.72596	0
	1	0.72009	0
	1	0.71399	0
	1	0.70763	0
	1	0.70099	0
	1	0.69406	0
	1	0.68682	0
	1	0.67923	0
	1	0.67125	0
	1	0.66283	0
	1	0.65397	0
	1	0.64465	0
	1	0.63488	0
	1	0.62462	0
	1	0.61382	0
	1	0.60241	0
	1	0.59032	0
	1	0.57748	0
	1	0.56384	0
	1	0.54932	0
	1	0.53384	0
	1	0.51719	0
	1	0.49929	0
	1	0.48026	0
	1	0.4602	0];

	% this is the kernel of MKPJ
	if (n == 256)		% Pre computed cmap. Nothing to do.
		map = JetI;
	else
		idx1 = linspace(1,n,size(JetI,1));
		map = interp1(idx1,JetI,1:n,'cubic');
	end
