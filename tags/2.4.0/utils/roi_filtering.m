function [Z,Z_rect,handles] = roi_filtering(handles, Z, head, Z_rect, r_c, mask, border)
% Do median filtering on either a rectangular or arbitrarely shaped ROI zone.

%	Copyright (c) 2004-2012 by J. Luis
%
% 	This program is part of Mirone and is free software; you can redistribute
% 	it and/or modify it under the terms of the GNU Lesser General Public
% 	License as published by the Free Software Foundation; either
% 	version 2.1 of the License, or any later version.
% 
% 	This program is distributed in the hope that it will be useful,
% 	but WITHOUT ANY WARRANTY; without even the implied warranty of
% 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% 	Lesser General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

	is_rect = 0;			done = 0;
	if (ischar(mask))		is_rect = 1;	end
	if (nargin == 6)		border = 'no';	end

	prompt = {'Enter number of filter rows' ,'Enter number of filter cols', 'Maximum allowed change (z units)'};
	def = {num2str(5) num2str(5) num2str(10)};
	resp  = inputdlg(prompt,'Median Filtering',[1 30; 1 30; 1 30],def);    pause(0.01)
	if (isempty(resp))		set(handles.figure1,'pointer','arrow'),		return,		end
	Z_rect = double(Z_rect);				% It has to be with medfilt2 of R13, but I'll hve to change it to R2006b

	Z_filt = img_fun('medfilt2',Z_rect,[str2double(resp{1}) str2double(resp{2})]);
	max_z_cut = str2double(resp{3});
	dife = (Z_rect - Z_filt);
	idx_dif_p = find(dife > max_z_cut);
	idx_dif_n = find(dife < -max_z_cut);
	Z_filt(idx_dif_p) = Z_rect(idx_dif_p) - max_z_cut;
	Z_filt(idx_dif_n) = Z_rect(idx_dif_n) + max_z_cut;
	clear idx_dif_p idx_dif_n dife;

	if (is_rect)
		if (strcmp(border,'no'))     % We don't need to do border re-interpolation
			hfc = fix(str2double(resp{1})/2);   % Half filter width. Don't use the filtered data on beyond
			hfr = fix(str2double(resp{2})/2);   % the transition zone. It may produce sharp edges.
			if (isa(Z,'single'))
				Z_rect = single(Z_filt);		clear Z_filt;
			elseif (isa(Z,'int16'))
				Z_rect = int16(Z_filt);			clear Z_filt;
			elseif (isa(Z,'uint16'))
				Z_rect = uint16(Z_filt);		clear Z_filt;
			else		% Should never come here
				Z_rect = single(Z_filt);		clear Z_filt;
			end
			Z(r_c(1)+hfr:r_c(2)-hfr,r_c(3)+hfc:r_c(4)-hfc) = Z_rect(hfr+1:end-hfr,hfc+1:end-hfc);
			done = 1;
		else            % User asked for border re-interpolation. We must create a mask grid
			mask = true(size(Z_rect));
		end
	end

	handles.Z_back = Z(r_c(1):r_c(2),r_c(3):r_c(4));    % For the undo op
	handles.r_c = r_c;

	if (done)		return,		end			% Rectangle without border re-interpolation

	Z_rect = smooth_roipoly_edge(head, handles.have_nans, Z, handles.Z_back, Z_filt, r_c, mask, 3);

	if (isa(Z,'single'))
		Z(r_c(1):r_c(2),r_c(3):r_c(4)) = single(Z_rect);
	elseif (isa(Z,'int16'))
		Z(r_c(1):r_c(2),r_c(3):r_c(4)) = int16(Z_rect);
	elseif (isa(Z,'uint16'))
		Z(r_c(1):r_c(2),r_c(3):r_c(4)) = uint16(Z_rect);
	else		% Should never come here
		Z(r_c(1):r_c(2),r_c(3):r_c(4)) = single(Z_rect);
	end
