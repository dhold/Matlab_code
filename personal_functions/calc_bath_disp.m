function [a_n,phi_n] = calc_bath_disp(nn,A_0,N,rho_HEOM,scalefct,x)
%calculates the total bath displacement distribution
% N is the number of system of degrees of freedom
%if scalefct is given then this scaling will be applied to counter any
%previous scaling



n2 = sum(nn,2);
max_tier = max(n2);
tL = size(rho_HEOM,2);

a_n = zeros(N^2,tL,max_tier+1);

fct_tab = [1,cumprod(1:max_tier)];

for n=0:max_tier

    lg = n2 == n; n_sec = nn(lg,:);    
weight_fct = (-1)^n/A_0^(n/2) *sqrt(fct_tab(n+1))./ prod( fct_tab(n_sec+1),2);
          
    lg2 = repmat(lg,[N*N,1]); %rep to size for HEOM
    rho_n = rho_HEOM(lg2,:); cnt = 0;      
if exist('scalefct','var')
    scalefct_tmp = scalefct(lg);
    for lp = 1:sum(lg)
        rho_sec = rho_n(cnt+1:cnt+N^2,:) ;
        a_n(:,:,n+1) = a_n(:,:,n+1) + scalefct_tmp(lp)*weight_fct(lp)*rho_sec;
        cnt = cnt+N^2;
    end
       
else   
    
    for lp = 1:sum(lg)
        rho_sec = rho_n(cnt+1:cnt+N^2,:) ;
        a_n(:,:,n+1) = a_n(:,:,n+1) + weight_fct(lp)*rho_sec;
        cnt = cnt+N^2;
    end
    
end
end

%if x range is given then also generate Hermite polynomials
if exist('x','var')
    x2 = x/sqrt(2*A_0); x2 = reshape(x2,1,length(x2));
    phi_n = zeros(max_tier+1,length(x));
    phi_n(1,:) = exp(-x2.^2)/sqrt(2*pi);
    phi_n(2,:) = sqrt(2)*phi_n(1,:).*x2;
    for n = 1:max_tier-1
        phi_n(n+2,:) = sqrt(2)*phi_n(n+1,:).*x2/sqrt(n+1) -...
                        sqrt(n/(n+1))*phi_n(n,:);
    end
else
    phi_n = [];
end