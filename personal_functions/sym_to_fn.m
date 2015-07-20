function sym_to_fn(fn_name,sym_fn,fn_args)
% Converts a symbolic expression to a .m function which can take
% identically sized matrix / scalar inputs for the symbolic variables
% with only element operations happening between them, will also worth with
% things like cross(var1,var2) but not matrix products etc unless you 
% define a function "mat_prod" and put it in the symbolic expression. 
% 
%fn_name is name of the resulting functions **assumed to end in .m**
%sym_fn is the symbolic expression to be converted to a .m files
%fn_args should be a cell array with all the sym variables found in the
%function sym_fn
% Remember to delete these temporary .m files after you have used them and
% no longer require them.
% Author : David Holdaway
if nargin == 2
fn_args = symvar(sym_fn);
end
if exist(fn_name,'file')
delete(fn_name) %remove old function, don't call with fn_name a file you want for gods sake
end
temp = char(sym_fn);

%replace operations with element wise versions
L = strfind(temp,'*');
if ~isempty(L)
temp2 =temp(1:L(1)-1);
for k = 2:length(L)
temp2 = strcat(temp2,'.*',temp(L(k-1)+1:L(k)-1));
end
temp  = strcat(temp2,'.*',temp(L(end)+1:end));
end

L = strfind(temp,'/');
if ~isempty(L)
temp2 =temp(1:L(1)-1);
for k = 2:length(L)
temp2 = strcat(temp2,'./',temp(L(k-1)+1:L(k)-1));
end
temp  = strcat(temp2,'./',temp(L(end)+1:end));
end

L = strfind(temp,'^');
if ~isempty(L)
temp2 =temp(1:L(1)-1);
for k = 2:length(L)
temp2 = strcat(temp2,'.^',temp(L(k-1)+1:L(k)-1));
end
temp  = strcat(temp2,'.^',temp(L(end)+1:end));
end

fid = fopen(fn_name,'w');
fprintf(fid,strcat('function out = ',fn_name(1:end-2),'('));
%reason for end-2 is last two characters assumed to be .m
for k = 1:length(fn_args)-1
fprintf(fid,strcat(char(fn_args(k)),','));    
end
fprintf(fid,strcat(char(fn_args(end)),') \n'));  
fprintf(fid,'out ='); %print variable into a function to quicken subs
fprintf(fid,temp);
fprintf(fid,';');
fclose(fid);