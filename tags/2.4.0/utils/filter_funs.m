function img = filter_funs( handles, opt )
% image filters functions
% _______ A FAZER

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

	img = (get(handles.hImg,'CData'));
    
	switch opt
		case 'SUSAN'
			for (i=1:size(img,3)),		img(:,:,i) = susan(img(:,:,i),'-s','-3');    end
		case 'Median'
			if (ndims(img) == 2)
				img = img_fun('medfilt2',img,'indexed');
			else				% RGB image. Do't really know if the following is correct
				for (k = 1:3),			img(:,:,k) = img_fun('medfilt2',img(:,:,k));	end
			end
		case 'Adaptive'
			for (k = 1:size(img,3)),	img(:,:,k) = img_fun('adpmedian', img(:,:,k), 7);		end
		case 'STD'
			img = img_fun('stdfilt',img);
		case 'Min'
			for (k = 1:size(img,3)),	img(:,:,k) = ordfilt2(img(:,:,k), 1, ones(3), 'symmetric');		end
		case 'Max'
			for (k = 1:size(img,3)),	img(:,:,k) = ordfilt2(img(:,:,k), 9, ones(3), 'symmetric');		end
		case 'range'
			img = img_fun('rangefilt',img);
	end
	
	set(handles.hImg,'CData',img)
