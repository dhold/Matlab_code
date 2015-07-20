function [R_red_op_full,ham_rs,R_red_sc] = gen_redfield_Liouville(R_red,sz1,sz2,rescale)
if nargin == 3
    rescale = false;
end
 %take out imaginary elements of R_red(a,b,a,b) which are just rescales of
 %the transition frequencies omega_{ab}, due to redundencies need only
 %rescale omega_{0 e} and omega_{0 f} etc etc as otheres can be expressed
 %in terms of these
 R_red_sc  = R_red; %scaled redfield
ham_rs = zeros(1,sz2);  %will rescale Hamiltonian

%clean up values from numerical errors
 for k = 1:length(R_red_sc)
    R_red_sc(k,k,k,k) =  real(R_red(k,k,k,k)); %remove eps level imag parts
 end
 % Now also remove all the eps level miss matches that mean R(a,b,a,b) ~=
 % conj(R(b,a,b,a)) by averaging
 for k = 1:length(R_red_sc)
     for j= (k+1):length(R_red_sc)

    R_red_sc(k,j,k,j) = (R_red(k,j,k,j) + conj(R_red(j,k,j,k)))/2; 
    R_red_sc(j,k,j,k) = conj(R_red_sc(k,j,k,j)); 
         
     end
 end 

    
 if rescale
 for a=2:N+1
        R_red_sc(1,a,1,a) = R_red_sc(1,a,1,a) - 1i*imag(R_red(1,a,1,a));    
        R_red_sc(a,1,a,1) = R_red_sc(a,1,a,1) - 1i*imag(R_red(a,1,a,1)); 
        
       for b=N+1+(1:(N-1)*N/2)
        R_red_sc(1,b,1,b) = R_red_sc(1,b,1,b) - 1i*imag(R_red(1,b,1,b));    
        R_red_sc(b,1,b,1) = R_red_sc(b,1,b,1) - 1i*imag(R_red(b,1,b,1));           
%          actually don't rescale these, just ground->whatever
%         R_red_sc(a,b,a,b) = R_red_sc(a,b,a,b) - 1i*imag(R_red(a,b,a,b));
%         R_red_sc(b,a,b,a) = R_red_sc(b,a,b,a) - 1i*imag(R_red(b,a,b,a));        
       end                    
 end

for b= 2:size(R_red,1)
            ham_rs = [ham_rs,ones(1,sz2)*imag(R_red(1,b,1,b))]; 
end
else
    ham_rs = zeros(1,sz1*sz2);
end 

R_red_op_full = zeros(sz1^2*sz2^2); 
%pad this out to the full thing
temp = zeros(sz1*sz2); 

cnt1vib = 0;  cnt2vib = 1; cnt1 = 1; cnt2 = 1;
%tic
for lp = 1:sz1^2*sz2^2
    
    %iterate along one vibrational lvl
    cnt1vib = cnt1vib+1;
    if cnt1vib > sz2 
        cnt1vib=1; cnt1 = cnt1+1;
        if cnt1 > sz1
            cnt1 = 1;  cnt2vib = cnt2vib+1;
            if cnt2vib > sz2
                cnt2vib = 1; cnt2 = cnt2+1;
            end
        end
    end
    temp = temp*0;
   for a = 1:sz1 
        for b = 1:sz1 
            temp((a-1)*sz2+cnt1vib,(b-1)*sz2+cnt2vib) = R_red_sc(cnt1,cnt2,a,b);
        end
   end

    R_red_op_full(lp,:) = reshape(temp,1,sz1^2*sz2^2 ).';

end
 R_red_op_full = sparse( R_red_op_full);