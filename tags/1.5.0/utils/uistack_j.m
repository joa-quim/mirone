function varargout = uistack_j(Handles,StackOpt,Step)
%UISTACK_J Restack objects children of AXES.
%   UISTACK_J(H) visually raises the visual stacking order of the objects, H. 
%   UISTACK_J(H,STACKOPT) where STACKOPT is 'up','down','top','bottom'
%   UISTACK_J(H,STACKOPT,STEP) where STEP is the distance to move 'up' and 'down'
%   applies the stacking option to the handles specified by H. All
%   handles, H, must have the same parent.

%   This is a hack of UISTACK that should only be used with handles children of
%   Image and Surface objects (e.g. NO UICONTROLS).
%   However, it corrects a bug present on the ML function that doesn't allow
%   restacking between different object types.

error(nargchk(0,1,nargout));
error(nargchk(1,3,nargin));
if ~all(ishandle(Handles)),
    error('Invalid Handles passed to UISTACK');
end
if (nargin == 1)
    StackOpt='up';
    Step=1;
end
if (nargin == 2),    Step=1;     end
if (Step < 0),       Step = 0;   end

Parent=get(Handles,{'Parent'});
Parent=[Parent{:}];
UParent=unique(Parent);
if length(UParent)>1,
    error('All handles passed to UISTACK must have the same Parent.');
end

% change stack order
NewOrder = getNewOrder(UParent,Handles, StackOpt, Step);
set(UParent,'Children',NewOrder);

if (nargout),    varargout{1}=allchild(UParent);    end

% --------------------------------------------------------------------
function NewOrder = getNewOrder(UParent,Handles, StackOpt, Step)

NOUSE = -1;

Children = allchild(UParent);
HandleLoc=find(ismember(Children,Handles));

switch StackOpt,
case 'up',
    NewOrder=[ones(Step,1).*NOUSE;Children];
    HandleLoc = HandleLoc + Step;
    for lp=1:length(Handles),
        Idx=HandleLoc(lp);
        NewOrder= [NewOrder(1:Idx-Step-1);NewOrder(Idx);NewOrder(Idx-Step:Idx-1);NewOrder(Idx+1:length(NewOrder))];
    end
    NewOrder(NewOrder == NOUSE) = [];
    
case 'down',
    NewOrder=[Children;ones(Step,1).*NOUSE];
    for lp=length(Handles):-1:1,
        Idx=HandleLoc(lp);
        NewOrder = [NewOrder(1:Idx-1);NewOrder(Idx+1:Idx+Step);NewOrder(Idx);NewOrder(Idx+Step+1:length(NewOrder))];
    end
    NewOrder(NewOrder == NOUSE) = [];
    
case 'top',
    Children(HandleLoc)=[];  
    NewOrder=[Handles;Children];
    
case 'bottom',
    Children(HandleLoc)=[];  
    NewOrder=[Children;Handles];    
    
otherwise,
    error('Invalid Stack option for UISTACK');      
    
end % switch

% If image, light, or surface object exists, put them at the bottom of the stack
hImg = findobj(UParent,'type','image');
hLight = findobj(UParent,'type','light');
hSurf = findobj(UParent,'type','surface');
hAx_bot = [hImg; hLight; hSurf];
if (~isempty(hAx_bot))
    NewOrder(NewOrder == hAx_bot) = [];
    NewOrder=[NewOrder; hAx_bot];
end
