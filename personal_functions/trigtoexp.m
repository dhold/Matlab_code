function out = trigtoexp(s)

out = char(sym(s));
for cses = 1:4
    
    switch cses
        case 1
L = strfind(out,'sin(');
        case 2
L = strfind(out,'cos(');
        case 3
L = strfind(out,'tan(');            
    end

for lp1 = 1:length(L)
    
    lft_cnt = 1; rgt_cnt=0;%number of brackets counted
    
    for lp2 = L(lp1)+4:length(out)
        if strcmp(out(lp2),'(')
            lft_cnt = lft_cnt+1;
        elseif strcmp(out(lp2),')')
            rgt_cnt = rgt_cnt+1;
            if rgt_cnt == lft_cnt
                break
            end
        end
    end
    exponnt = out(L(lp1)+4:(lp2-1));
    switch cses
        case 1
 new_str = strcat('(exp(i*',exponnt,') - exp(-i*',exponnt,'))/(2i)'); 
        case 2
new_str = strcat('(exp(i*',exponnt,') + exp(-i*',exponnt,'))/2'); 
        case 3
new_str = strcat('(i*exp(-i*',exponnt,') - i*exp(+i*',exponnt,'))/(exp(i*',exponnt,') + exp(-i*',exponnt,'))');           
    end
    end_pos = lp2;
    out = strcat(out(1:L(lp1)-1),new_str,out(end_pos+1:end));
    L(lp1+1:end) = L(lp1+1:end)+length(new_str)+L(lp1)-end_pos-1;
    %correct for the length of the new string
end
end