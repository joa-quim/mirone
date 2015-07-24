function [out1,out2,out3] = click_e_point(arg1,arg2,arg3)
%GINPUT Graphical input from mouse.
%   [X,Y] = GINPUT_POINTER(N) gets N points from the current axes and returns 
%   the X- and Y-coordinates in length N vectors X and Y.  The cursor
%   can be positioned using a mouse (or by using the Arrow Keys on some 
%   systems).  Data points are entered by pressing a mouse button
%   or any key on the keyboard except carriage return, which terminates
%   the input before N points are entered.
%
%   [X,Y] = GINPUT_POINTER gathers an unlimited number of points until the
%   return key is pressed.
% 
%   [X,Y,BUTTON] = GINPUT_POINTER(N) returns a third result, BUTTON, that 
%   contains a vector of integers specifying which mouse button was
%   used (1,2,3 from left) or ASCII numbers if a key on the keyboard
%   was used.
%
%   This version of ginput lets the user specify what pointer should be used while
%   in ginput mode.
%   [...] = GINPUT_POINTER(N,POINTER)     or  [...] = GINPUT_POINTER(POINTER)
%   where POINTER is either a valid pointer name or a 16x16 matrix deffining a custom pointer
%
%   [...] = GINPUT_POINTER(N,POINTER,PSHS)
%   where PSHS is a two element vector, sets the 'PointerShapeHotSpot' property

%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 5.32 $  $Date: 2002/04/14 13:20:49 $

%   In order to change as little as possible, I changed the original tests against nargin
%   to tests to nargs, where this new variable holds the values of the original nargin.
%   Coffeeright J. Luis 2003

PSHS = [];
if (nargin >= 2)
    if ~(ischar(arg2) || size(arg2,1) == 16)
        error('ginput error: second argument must be a string containing the pointer name');
    elseif size(arg2,1) == 16
        new_pointer = 'custom';
        nargs = 1;
        if (nargin == 3 && numel(arg3) == 2)    % Got a 'PointerShapeHotSpot'
            PSHS = arg3;
        end
    else
        new_pointer = arg2;
        nargs = 1;             % So that the rest of the code may work as if this change has not been made
    end
elseif (nargin == 1)
    if ischar(arg1)
        new_pointer = arg1;
        nargs = 0;             % So that the rest of the code may work as if this change has not been made
    elseif length(arg1) == 16
        new_pointer = 'custom';
        nargs = 0;
    else
        new_pointer = 'fullcrosshair';
        nargs = 1;
    end
else
    new_pointer = 'fullcrosshair';
    nargs = 0;
end

out1 = []; out2 = []; out3 = []; y = [];
if ~ispc
   tp = get(0,'TerminalProtocol');
else
   tp = 'micro';
end

if ~strcmp(tp,'none') && ~strcmp(tp,'x') && ~strcmp(tp,'micro')
%    the calls to trmginput give a lot anoying compiling warnings
%    if nargout == 1,
%       if nargs == 1,
%          out1 = trmginput(arg1);
%       else
%          out1 = trmginput;
%       end
%    elseif nargout == 2 | nargout == 0,
%       if nargs == 1,
%          [out1,out2] = trmginput(arg1);
%       else
%          [out1,out2] = trmginput;
%       end
%       if  nargout == 0
%          out1 = [ out1 out2 ];
%       end
%    elseif nargout == 3,
%       if nargs == 1,
%          [out1,out2,out3] = trmginput(arg1);
%       else
%          [out1,out2,out3] = trmginput;
%       end
%    end
else   
    fig = gcf;
    figure(fig);
    if nargs == 0
        how_many = -1;
        b = [];
    else
        how_many = arg1;
        b = [];
        if  (ischar(how_many) || size(how_many,1) ~= 1 || size(how_many,2) ~= 1 ...
                || ~(fix(how_many) == how_many) || how_many < 0)
            error('Requires a positive integer.')
        end
        if how_many == 0
            ptr_fig = 0;
            while (ptr_fig ~= fig)
                ptr_fig = get(0,'PointerWindow');
            end
            scrn_pt = get(0,'PointerLocation');
            loc = get(fig,'Position');
            pt = [scrn_pt(1) - loc(1), scrn_pt(2) - loc(2)];
            out1 = pt(1); y = pt(2);
        elseif how_many < 0
            error('Argument must be a positive integer.')
        end
    end
   
    % Remove figure button functions
    state = uisuspend_j(fig);
    if strcmp(new_pointer,'custom')              % Also changed here
        set(fig,'pointer',new_pointer,'PointerShapeCData',arg2)
        if (~isempty(PSHS)),    set(fig,'PointerShapeHotSpot',PSHS);    end
    else
        set(fig,'pointer',new_pointer);
    end
    fig_units = get(fig,'units');
    char = 0;
   
    while how_many ~= 0      % Use no-side effect WAITFORBUTTONPRESS
        waserr = 0;
        try
			keydown = wfbp(fig);
		catch
			waserr = 1;
        end
        if (waserr == 1)
            if(ishandle(fig))
                set(fig,'units',fig_units);
                uirestore_j(state, 'nochildren');
                error('Interrupted');
            else
                error('Interrupted by figure deletion');
            end
        end
      
        ptr_fig = get(0,'CurrentFigure');
        if (ptr_fig == fig)
            if keydown
                char = get(fig, 'CurrentCharacter');
                button = abs(get(fig, 'CurrentCharacter'));
                scrn_pt = get(0, 'PointerLocation');
                set(fig,'units','pixels')
                loc = get(fig, 'Position');
                pt = [scrn_pt(1) - loc(1), scrn_pt(2) - loc(2)];
                set(fig,'CurrentPoint',pt);
            else
                button = get(fig, 'SelectionType');
                if strcmp(button,'open')
                    button = b(length(b));
                elseif strcmp(button,'normal')
                    button = 1;
                elseif strcmp(button,'extend')
                    button = 2;
                elseif strcmp(button,'alt')
                    button = 3;
                else
                    error('Invalid mouse selection.')
                end
            end
            pt = get(gca, 'CurrentPoint');
         
            how_many = how_many - 1;
         
            if (char == 13) % & how_many ~= 0)
                % if the return key was pressed, char will == 13,
                % and that's our signal to break out of here whether
                % or not we have collected all the requested data points.  
                % If this was an early breakout, don't include
                % the <Return> key info in the return arrays.
                % We will no longer count it if it's the last input.
                break
            end
         
            out1 = [out1; pt(1,1)];
            y = [y; pt(1,2)];
            b = [b; button];
        end
    end
   
    uirestore_j(state, 'nochildren');
    set(fig,'units',fig_units);
   
    if nargout > 1
        out2 = y;
        if (nargout > 2),   out3 = b;   end
    else
        out1 = [out1 y];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = wfbp(fig)
%WFBP   Replacement for WAITFORBUTTONPRESS that has no side effects.

% Now wait for that buttonpress, and check for error conditions
waserr = 0;
try
	h=findall(fig,'type','uimenu','accel','C');   % Disabling ^C for edit menu so the only ^C is for
	set(h,'accel','');                            % interrupting the function.
	keydown = waitforbuttonpress;
	current_char = double(get(fig,'CurrentCharacter')); % Capturing the character.
	if ~isempty(current_char) && (keydown == 1)          % If the character was generated by the 
		if(current_char == 3)                       % current keypress AND is ^C, set 'waserr'to 1
			waserr = 1;                             % so that it errors out. 
		end
	end
    set(h,'accel','C');                                 % Set back the accelerator for edit menu.
catch
    waserr = 1;
end
drawnow;
if (waserr == 1)
    set(h,'accel','C');                                % Set back the accelerator if it errored out.
    error('Interrupted');
end

if nargout>0, key = keydown; end
