function bool=istilde(n)
%ISTILDE  True if an output is supressed via the tilde operator.
%   ISTILDE(N) returns logical 1 (true) if the N-th output argument of the
%   calling function is ignored with the tilde operator, and logical 0 (false)
%   otherwise. The calling function is the function that calls the the function
%   in which ISTILDE is called.
%
%   ISTILDE() without an argument returns logical 1 (true) if any arguments of
%   the calling function are ignored with the tilde operator, and logical 0
%   (false) otherwise.
%
%   Example:    
%       function istildedemo
%           [A,B,C] = f();	% Prints false 
%           [A,~,C] = f();	% Prints true, large B not allocated
%           [~,B,~] = f();	% Prints false
%       end
%
%       function [A,B,C] = f()
%           isB = istilde(2)
%           if isB
%               B = [];     % Output must still be assigned
%           else
%               B = rand(100);
%           end
%           A = 1; C = 3;
%       end
%
%   Note:
%       Outputs that are ignored via the tilde operator must still be assigned,
%       as the example above illustrates. ISTILDE permits costly memory
%       allocations and computations to be avoided by detecting if a particular
%       output was supressed.
%
%       ISTILDE is not supported for functions called from Matlab's Command
%       Window and will return false in all cases in this situation.
%
%       ISTILDE is not a robust solution. This function opens and reads the
%       M-File of the calling function. If this file changes while a program is
%       running, e.g., a user makes an edit and saves, ISTILDE could give
%       erroneous results.
%
%   See also PUNCT

%   Andrew D. Horchler, adh9 @ case . edu, Created 5-4-13
%   Revision: 1.2, 5-6-13


% Check optional input
if nargin > 0
    if ~isnumeric(n) || ~isscalar(n) || ~isreal(n)
        error('istilde:InvalidDataType',...
              'Input must be a string or a real numeric scalar value.');
    end
    if ~isfinite(n) || n < 1 || n ~= floor(n)
        error('istilde:NonIntegerValue',...
              'Input must be an integer greater than or equal to one.');
    end
    if evalin('caller','nargout') < n
        bool = false;
        return;
    end
end

% Get full path of M-file of calling function
st = dbstack(2,'-completenames');
if isempty(st)
    % ISTILDE is not supported for functions called from the Command Window
    bool = false;
    return;
end

% Open file where call was made for reading (assume has not changed)
fid = fopen(st.file,'r');

% Set up cleanup function to close file in case something goes wrong
C = onCleanup(@()~isempty(fopen(fid))&&fclose(fid));

% Read just line of function call, parse out arguments between square brackets
str = textscan(fid,'%[[]%[^]]',1,'Headerlines',st.line-1);

% Close file
fclose(fid);

str = str{2}{1};
if isempty(str)
    % No output or one output not in square brackets, tilde not possible
    bool = false;
else
    % If a tilde is present among output arguments (may have spaces or '...')
    bool = ~isempty(regexp(str,...
        '(^|,)(\s*(\.\.\.)?\s*)*~(\s*(\.\.\.)?\s*)*(,|$)','once'));
    if bool && nargin > 0
        % Split up output arguments around commas
        str = regexp(str,'(\s*(\.\.\.)?\s*)*,(\s*(\.\.\.)?\s*)*','split');
        
        % Check if a tilde is present in N-th argument
        bool = strcmp(str(n),'~');
    end
end