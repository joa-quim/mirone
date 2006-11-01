% This is a adapted copy of a startup file from a program called Ospray.

% This is an example startup.m file for using M_GMT.  Use it directly,
% or if you already have a startup.m file, append the contents of this
% file to yours.  In any case, edit the stuff below as needed.


% Change the following command to specify the directories where you installed 
% the 'Mirone' directory.  Leave the '-end' part alone.
addpath('D:\mirone_devel\', '-end');


% You can make an M_GMT icon for your desktop this way:  Find matlab.exe
% on your hard drive.  Create a shortcut to it (by right-clicking on it)
% and put the shortcut on your desktop.  Rename the shortcut 'Mirone'.
% Right-click on the shortcut and choose 'Properties'.  In the box labeled
% 'Start in', type the directory name where this file, startup.m, is
% stored.  For example, if this file is C:\MyHome\startup.m, then change
% the 'Start in' box to say C:\MyHome .


% The following cd command is optional.  If you uncomment it (by removing 
% the '%' at the start of the line), it will make mirone start up in the 
% directory it specifies, so you don't have to hunt around to find your data
% files.  Edit the stuff between the quote marks to make the directory be 
% where you store your data:
%cd('C:\MyHome\MyDataDir\');


% If you have trouble with Matlab's 'splash box' (the box that says 
% 'Matlab/The language of technical computing/Version X.Y.Z ...') sticking on
% the screen when you first launch Matlab/m_gmt, then uncomment this line:
%pause(3)


% Launch mirone.
mirone
