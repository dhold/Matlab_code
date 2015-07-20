
%call to exit matlab
%if nargin == 0
setpref('StartupDirectory','LastWorkingDirectory',pwd)
%elseif ischar(starting_directory)
%setpref('StartupDirectory','LastWorkingDirectory',starting_directory) %lies to it
%end