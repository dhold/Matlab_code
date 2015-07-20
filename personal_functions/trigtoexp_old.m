function out = trigtoexp(s)

k = char(sym(s));

L = strfind(k,'sin(');

for lp = 1:length(L)
    
    lft_cnt = 1; rgt_cnt=0;%number of brackets counted
    
    for lp2 = L(lp1+4):length(k)
        if strcmp(k(lp2),'(')
            lft_cnt = lft_cnt+1;
        elseif strcmp(k(lp2),')')
            rgt_cnt = rgt_cnt+1;
            if rgt_cnt == lft_cnt
                break
            end
        end
    end
    
    exponnt{lp}= k((lp1+4):(lp2-1));
    end_pos_sin(lp) = lp2;
end

L2 = strfind(k,'cos(');

for lp = 1:length(L2)
    
    lft_cnt = 1; rgt_cnt=0;%number of brackets counted
    
    for lp2 = L(lp1+4):length(k)
        if strcmp(k(lp2),'(')
            lft_cnt = lft_cnt+1;
        elseif strcmp(k(lp2),')')
            rgt_cnt = rgt_cnt+1;
            if rgt_cnt == lft_cnt
                break
            end
        end
    end
    
    exponnt2{lp}= k((lp1+4):(lp2-1));
    end_pos_cos(lp) = lp2;
end
    out = k;
for lp1 = 1:length(L)
    new_str = '(exp(i*exponnt) - exp(-i*exponnt))/(2i)'; 
    old_rng = L(lp1):end_pos_sin(lp1);
    out = [out(1:L(lp1)-1),new_str,end_pos_sin(lp1)+1:end];
    k(L:end_pos_sin(lp1)) = 