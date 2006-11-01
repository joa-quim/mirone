function uistate = uiclearmode_j(fig, varargin)
%UICLEARMODE Clears the current interactive figure mode;
%  UISTATE=UICLEARMODE(FIG) suspends the interactive properties of a 
%  figure window and returns the previous state in the structure
%  UISTATE.  This structure contains information about the figure's
%  WindowButton* functions and the cursor.  It also contains the 
%  ButtonDownFcn's for all children of the figure.
%
%  UISTATE=UICLEARMODE(FIG, FUNCTION [, ARGS]) suspends the
%  interactive properties of the figure FIG in two ways.
%  First, uiclearmode notifies the currently active mode, if
%  any are active, that the current mode should deinstall its
%  event handlers, such as its Figure WindowButtonDown
%  callbacks.  Next, UICLEARMODE installs the function
%  FUNCTION as the deinstaller for the new mode.  Finally, 
%  UICLEARMODE, like UISUSPEND, resets the WindowButton* 
%  functions and returns that information in a struct which
%  can be saved and passed back to UIRESTORE.
%
%  UISTATE=UICLEARMODE(FIG,'docontext',...) also suspends
%  uicontext menus
%
%  Example:
%      
%  The following function defines a new interactive mode that
%  cooperates with other modes such as plotedit and rotate3d.
%
%  That is, before myinteractivemode is activated, plot
%  editing or rotate3d is turned off.  If myinteractivemode
%  is active, then activating plot editing calls
%  myinteractive(fig,'off').  The calling syntax for
%  myinteractivemode is:
%
%     myinteractivemode(gcf,'on')   % display figure current
%                                   % point on mouse down
%   
%   function myinteractivemode(fig,newstate)
%   %MYINTERACTIVEMODE.M
%   persistent uistate;
%      switch newstate
%      case 'on'
%         disp('myinteractivemode: on');
%         uistate = uiclearmode(fig,'myinteractivemode',fig,'off');
%         set(fig,'UserData',uistate,...
%                 'WindowButtonDownFcn',...
%                 'get(gcbf,''CurrentPoint'')',...
%                 'Pointer','crosshair');
%      case 'off'
%         disp('myinteractivemode: off');
%         if ~isempty(uistate)
%            uirestore(uistate);
%            uistate = [];
%         end
%      end
%
%  See also UISUSPEND, UIRESTORE, SCRIBECLEARMODE.

%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 1.10.4.5 $ $Date: 2004/10/22 19:45:11 $

%  *** Undocumented Syntax ***
%  UISTATE=UICLEARMODE(FIG,'docontext_preserve',...) stores
%  uicontext menus in output state, but does not remove them 
%  from the current plot like 'docontext'
%  UISTATE=UICLEARMODE(FIG,'keepWatch',OLDPTR,...) treats watch
%  cursors specially by not reseting them and storing
%  OLDPTR in UISTATE instead.

%#function scribeclearmode

nargs = nargin;
docontext = 0;
docontext_preserve = 0;
keepWatch = false;
firstarg = 1;

% catch docontext flag
if nargs>1
    if isa(varargin{1},'char')
        if strcmp(varargin{1},'docontext')
            docontext=1;
        elseif strcmpi(varargin{1},'docontext_preserve')
            docontext_preserve=1;
        elseif strcmp(varargin{1},'keepWatch')
          keepWatch = true;
          firstarg = 3;
        end
    end
end

% fetch current uistate from figure using uisuspend 
% the second input false tells uisuspend to just return
% the struct.
old_uistate = uisuspend_j(fig, false);

if (docontext || docontext_preserve)
    scribeclearmode(fig, varargin{2:nargs-1});
else
    scribeclearmode(fig, varargin{firstarg:end});
end

% now let uisuspend handle the defaults
uistate = uisuspend_j(fig, true);

% treat watch cursors specially since they are telling users
% something time consuming is happening and we shouldn't reset
% to the default cursor until the action is done.
if keepWatch && strcmp(old_uistate.Pointer,'watch')
    set(fig, 'Pointer', old_uistate.Pointer);
    uistate.Pointer = varargin{2};
end

% make sure plot edit stays on if it was on
% plotedit(fig, 'setenabletools',old_uistate.ploteditEnable);       % COMMENTED

% extra gets and sets if doing the uicontextmenus
if docontext
    uistate.docontext=1;
    nulluicontext = [];
    set(uistate.Children, 'UIContextMenu', nulluicontext);
elseif docontext_preserve
    % Same thing as docontext, but don't null out menus
    uistate.docontext=1;
else  
    uistate.docontext=0;
end
