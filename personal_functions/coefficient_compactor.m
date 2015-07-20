function coeffs = coefficient_compactor(S,t)

%reduces an expression to just a sum of exponentials A_k exp(B_k * t)
%assuming it can be reduced into this way, i.e. no exp(cos(t)) shit, note
%the t is the variable to find, which must be character

kk = char(sym(S));
for cses = 1:4
    
    switch cses
        case 1
L = strfind(kk,'sin(');
        case 2
L = strfind(kk,'cos(');
        case 3
L = strfind(kk,'exp(');            
    end

for lp1 = 1:length(L)
    
    lft_cnt = 1; rgt_cnt=0;%number of brackets counted
    
    for lp2 = L(lp1)+4:length(kk)
        if strcmp(kk(lp2),'(')
            lft_cnt = lft_cnt+1;
        elseif strcmp(kk(lp2),')')
            rgt_cnt = rgt_cnt+1;
            if rgt_cnt == lft_cnt %end of expression
                break
            end
        end
    end
    exponnt = kk(L(lp1)+4:(lp2-1));
    %find all the t's
    k2 = strfind(exponnt,strcat('*',t,'*'));
    
    k2 = strfind(exponnt,strcat('*',t));
    
    k2 = strfind(exponnt,strcat(t,'*'));    
    
    switch cses
        case 1
 new_str = strcat('(exp(i*',exponnt,') - exp(-i*',exponnt,'))/(2i)'); 
        case 2
new_str = strcat('(exp(i*',exponnt,') + exp(-i*',exponnt,'))/2'); 
        case 3
new_str = strcat('(i*exp(-i*',exponnt,') - i*exp(+i*',exponnt,'))/(exp(i*',exponnt,') + exp(-i*',exponnt,'))');           
    end
    end_pos = lp2;
    kk = strcat(kk(1:L(lp1)-1),new_str,kk(end_pos+1:end));
    L(lp1+1:end) = L(lp1+1:end)+length(new_str)+L(lp1)-end_pos-1;
    %correct for the length of the new string
end
end