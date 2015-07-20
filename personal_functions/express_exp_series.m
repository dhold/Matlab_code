function [coeffmatfinal,coeffmat,Sout,Sdiff] = express_exp_series(S,t,remove_dupe)

%
%reduces an expression to just a sum of exponentials A_k exp(B_k * t)
%assuming it can be reduced into this way, exp(-t^2) etc, will still work
%but not things like exp(cos(t))
%the t is the variable to find, which must be sym or character
%coeffmat is a 2XN array of c_k and b_k in S = sum_{k=0}^N c_k exp(b_k t)
% INPUT ARGS
% S is the input symbolic function (although this should also work for
% character inputs)
% t is whatever symbolic variable you are doing the expansion in
% remove_dupe is a logical value which determines whether or not you remove
% duplicates and collate into a single term, true by default
% OUTPUT ARGS
% coeffmatsym is set of symbolic c_k (row 1) and b_k*t (row 2)
% coeffmat is similar but character and not trimmed
% Sout is the series expressed in the new way
% Sdiff is a check that the output series is equal to the input one
%
%Author David Holdaway

%tic
 S = sym(S); %just in case
 S = rewrite(S, 'exp');
 S = expand(S);
  
 if nargin == 2
      remove_dupe = true;
 end
 
kk = char(S);

if S == 0 
    coeffmatsym{2,1} = []; return
end

L = strfind(kk,'exp(');  %find strings with exponentials
    
if isempty(L)
    coeffmatsym{1,1} = S;
    coeffmatsym{2,1} = 0;
    return
end


% L2 = strfind(kk,' ');  
% 
% LL = strfind(kk,'('); 
% LR = strfind(kk,')');
% bpairs = zeros(2,length(LL)); %brackets and which they pair too
% lft_cnt = 0; rgt_cnt=0;
% 
% LLUP = []; %unpaired brackets
% for lp = 1:length(kk)
%     if lft_cnt<length(LL) && lp == LL(lft_cnt+1) %found left bracket
%         
%         lft_cnt = lft_cnt+1;
%         LLUP = [LLUP,lp];
% 
%     elseif lp == LR(rgt_cnt+1) %found right
%         
%             rgt_cnt = rgt_cnt + 1;     
%             bpairs(:,rgt_cnt) = [LR(rgt_cnt);LLUP(end)];
%             LLUP = LLUP(1:end-1);
%         
%     end
%     if rgt_cnt == length(LL)
%         break
%     end
% end

cnt = 0; cnt2 = 0;  %number in chain and chain number
coeffmat{2,length(L)} = []; %trim later


end_ps = []; newchain = true; coeff_end = 1;
for lp1 = 1:length(L)
    
    if newchain
       cnt2 = cnt2 + 1;     
    %first identify the prefactor, unfortunately matlab does pre AND post
    %multiplies based on what it thinks looks nice...
    %lft_cnt = 0; rgt_cnt=0;%number of brackets counted    
    for lp2 = L(lp1)-1:-1:1
        if strcmp(kk(lp2),' ')
               coeff_end = lp2+1;
               %as a rule, because of how this is expressed a space at the
               %start will always indicate the end of the prefactor
               sgn{cnt2} = kk(coeff_end-2); %plus or minus
              break
        elseif lp2 == 1
            coeff_end = 1;
            sgn{cnt2} = ('+');
        end
    end
    if L(lp1) == 1 %previous loop is missed out if first exp has no prefact
        coeff_end = 1;
        sgn{cnt2} = ('+');
    end
    
    
    coeffmat{1,cnt2} = kk(coeff_end:L(lp1)-2);
    if isempty( kk(1:L(1)-2))
        coeffmat{1,cnt2} = char(sym(1));
    end
    if cnt2> 1 %add post factors and add them to the prefactor 
        %also with the appropriate sgn
        tmp_str = kk(end_ps(cnt)+1:coeff_end-4); %all values from the end
        %bracket of last exp up to the space before the next sgn value
        if strcmp(sgn{cnt2-1},'+')
            %+sign is basically implicit
    coeffmat{1,cnt2-1} = strcat(coeffmat{1,cnt2-1},tmp_str );           
        elseif strcmp(sgn{cnt2-1},'-')
    coeffmat{1,cnt2-1} = strcat(sgn{cnt2-1},coeffmat{1,cnt2-1},tmp_str );
        else
            warning('I think it broke, this should be +/-')
            sgn{cnt2-1}
        end
    end
    
    end
    
    lft_cnt = 1; rgt_cnt=0;%number of brackets counted
    for lp3 = L(lp1)+4:length(kk)
        if strcmp(kk(lp3),'(')
            lft_cnt = lft_cnt+1;
        elseif strcmp(kk(lp3),')')
            rgt_cnt = rgt_cnt+1;
            if rgt_cnt == lft_cnt %end of expression
                cnt = cnt + 1;
                end_ps(cnt) = lp3;
                break
            end
        end
    end
    
    if newchain
    coeffmat{2,cnt2} =kk(L(lp1)+4:end_ps(cnt)-1);    
    else
    coeffmat{2,cnt2} = strcat(coeffmat{2,cnt2},'+',kk(L(lp1)+4:end_ps(cnt)-1));
    end
    
    if end_ps(cnt) < length(kk)-4
       %attempt to find if there are other exponentials multipled together
       %Matlab does this as exp(asdasd)*exp(ggfgf)*exp...
       %after the rewriting that has been done already

       if strcmp(kk(end_ps(cnt)+1:end_ps(cnt)+4),'*exp')            
          newchain = false; %continue as part of the chain
       else
           newchain = true;          
       end                                 
    end
end

%add post factor to the last value, 
     tmp_str = kk(end_ps(cnt)+1:length(kk)); %all values from the end
        %bracket of last exp up to the space before the next sgn value
        if strcmp(sgn{cnt2},'+')
            %+sign is basically implicit
    coeffmat{1,cnt2} = strcat(coeffmat{1,cnt2},tmp_str );           
        elseif strcmp(sgn{cnt2},'-')
    coeffmat{1,cnt2} = strcat(sgn{cnt2},coeffmat{1,cnt2},tmp_str );
        else
            warning('I think it broke, this should be +/-')
            sgn{cnt2}
        end
%toc
%tic
%coeffmat{:,:}
nullval = false(1,cnt2);
for lp = 1:cnt2
%covert strings to symbolic variables.. hopefully
%sym(coeffmat{2,cnt2})^2
%sym(coeffmat{1,cnt2})^2

    coeffmatsym{2,lp} = sym(coeffmat{2,lp});
    try
    coeffmatsym{1,lp} = sym(coeffmat{1,lp});
    catch ME
        ME
        coeffmat{1,lp}
        warning('unable to interpret string')
    end

    %subs(coeffmatsym{2,lp},t,0)
    %o1=taylor(coeffmatsym{2,lp},t,'order',1);
    %o2=taylor(coeffmatsym{2,lp},t,'order',2);
    %o3=taylor(coeffmatsym{2,lp},t,'order',3);
    o1=subs(coeffmatsym{2,lp},t,0);    %subtract zero order in t
    tp1 = o1==0;
    if ~tp1
       coeffmatsym{1,lp} = coeffmatsym{1,lp}*exp(o1);
       coeffmatsym{2,lp} = coeffmatsym{2,lp}-o1;
    else
       coeffmatsym{2,lp} = coeffmatsym{2,lp}; 
    end
    %if ~o3==o2
    %    warning('exponentials are not all linear in t')
    %    lp
    %end
    %finally collate coefficients of equal exponentials, this step is
    %rather slow regretably
    if remove_dupe
        if coeffmatsym{1,lp} == 0 %don't save with 0 coefficient anyway
            nullval(lp) = true;
        else
    for lp2 = 1:lp-1
       if coeffmatsym{2,lp}==coeffmatsym{2,lp2}
          %duplicate found, remove this one and collate coefficients
          
          coeffmatsym{1,lp2} = coeffmatsym{1,lp2} + coeffmatsym{1,lp};
           nullval(lp) = true;
           coeffmatsym{2,lp} = 0;
           coeffmatsym{1,lp} = 0;
           
       end
    end
        end    
    end
end

%toc
%tic
cntvar = 0; %trim repeated values in the final one
for k = 1:cnt2
    if ~nullval(k)
        cntvar = cntvar+1;
   coeffmatfinal{1,cntvar} = coeffmatsym{1,k};
   coeffmatfinal{2,cntvar} = coeffmatsym{2,k};
    end
end
if nargout >= 3 
Sout = sym(0);
for k = 1:cntvar 
    Sout = Sout + coeffmatfinal{1,k}*exp(coeffmatfinal{2,k});
end
%toc
%test this is equal to the input value
if nargout==4 
    %tic
    Sdiff = simplify(Sout-S,'IgnoreAnalyticConstraints',true,'steps',100);
    if Sdiff~=0
   warning('function output may not be equal to input') 
    end
    %toc
end
end

end