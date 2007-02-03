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
        else
            img = img_fun('medfilt2',img);
        end
end

set(handles.hImg,'CData',img)
