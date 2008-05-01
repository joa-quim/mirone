function img = filter_funs( handles, opt )
% image filters functions
% _______ A FAZER

	img = (get(handles.hImg,'CData'));
    
switch opt
    case 'SUSAN'
        if (ndims(img) == 2)
            img = susan(img,'-s','-3');
        else
            for (i=1:3),    img(:,:,i) = susan(img(:,:,i),'-s','-3');    end
        end
    case 'Median'
        if (ndims(img) == 2)
            img = img_fun('medfilt2',img,'indexed');
        else				% RGB image. Do't really know if the following is correct
			for (k = 1:3)
            	img(:,:,k) = img_fun('medfilt2',img(:,:,k));
			end
        end
    case 'STD'
		img = img_fun('stdfilt',img);
    case 'range'
		img = img_fun('rangefilt',img);
end

set(handles.hImg,'CData',img)
