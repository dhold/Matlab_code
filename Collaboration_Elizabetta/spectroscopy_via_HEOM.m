function [OD,CD,SSden,chk,sup_op] = spectroscopy_via_HEOM(Kap2 ,om_0,reorg_E,...
                    drude_gam,freq_range,MM,Temp_kel,dipole_vec,R_disp) ;
                
                   
%generic units for scaling to cm^-1
    light_speed = 299792458; %units m s^{-1}, NOT cm 
    length_unit = 100; % units m
    hbar = 1.05457173*10^(-34); %units Js
    boltz_const = 1.3806488 * 10^(-23); %units J K^{-1}
                         
if nargin <=1 %use standard set of arguments
    warning('no inputs, using default');
    %first position relates to ground state
reorg_E = [25,25,25,25,25,25,25,25];
drude_gam = [35,35,35,35,35,35,35,35]; MM = -20;
om_0 = 12590;  Temp_kel = 77;
freq_range = linspace(om_0-150+25,om_0+150+25,31); %in wavenumber
%only used if any reorg energies differ
dipole_vec = ...
[ -0.823072771604003,	-0.5246284653639453,-0.21752284012024087;...
-0.21108138821716999,	-0.9529594677893477,	-0.2175152875063172;...
0.5246284653639459,	-0.8230727716040029,	-0.21752284012024128;...
0.9529594677893479,	-0.211081388217168,	-0.21751528750631696;...
0.8230727716040026,	0.524628465363946,	-0.2175228401202412;...
0.21108138821716746,	0.9529594677893481,	-0.2175152875063173;...
-0.5246284653639459,	0.8230727716040029,	-0.21752284012024128;...
-0.9529594677893481,	0.21108138821716843,	-0.2175152875063174];
%magnitude of dipole_vec isn't important here, just constant of
%proportionality, as it is not used to calculate the coupling between sites
R_disp = [72.062	,-12.081,90.095; 72.913	,9.924	,90.122;...
19.538	,11.081	,90.095; 18.687,	-10.924	,90.122;...
57.881	,25.262	,90.095;  35.876	,26.113	,90.122;...
33.719,	-26.262	,90.095; 55.724	,-27.113,90.122];
    if nargin == 0
        Kap2 = 2; 
    end
end
 N = length(reorg_E);
%these operators relate to calculation of absorption and CD
ddotd = zeros(size(R_disp,1)); 
R_jl = zeros(size(R_disp,1),size(R_disp,1),3); dcrossd = R_jl;
for k = 1:size(R_disp,1)
    for j=1:size(R_disp,1)
        
        ddotd(k,j) = dot(dipole_vec(k,:) , dipole_vec(j,:));

        R_jl(k,j,:) = R_disp(k,:) - R_disp(j,:);
        dcrossd(k,j,:) = cross(dipole_vec(k,:) , dipole_vec(j,:));
    end
end
CDop = R_jl(:,:,1).*dcrossd(:,:,1)+R_jl(:,:,2).*dcrossd(:,:,2)...
            +R_jl(:,:,3).*dcrossd(:,:,3);
Temp = Temp_kel*boltz_const/ (2 * pi * hbar * light_speed * length_unit);   
%% construct Heirarchy for elements relating to the ground state and coherences only, 
%trivial here but have the general code
if Kap2>0  %else trivial

numwithn = zeros(1,Kap2+1); %size(cc2,2) is number of poles in J(omega)

for kk = 0:Kap2
numwithn(kk+1) = nchoosek(N-1+kk,kk);
%number of density operators with each order, equal to size of 
%irreducible nth power rep of sym group U(N*(Kappa+size(cc2,2))) 
end
    
nn = zeros(sum(numwithn),N); 
count1=0; tmp3 = 0;
for k =2:Kap2+1

        count2 = count1 + size(tmp3,1); 
        
        tmp3 = zeros(length(count1+1:count2)*size(nn,2),size(nn,2));
       for j = count1+1:count2
        tmp = repmat(nn(j,:),size(nn,2),1); %select element from tier below
        tmp2 = tmp+eye(size(nn,2)); %add one to every position
        tmp3(1+(j-(count1+1))*size(nn,2):(j-count1)*size(nn,2) ,:) = tmp2;       
       end
       %remove any elements of tmp3 which fail to satisfy the high
       %frequency cut off condition

       %now remove duplicate entries
       jj=1;
       while jj < size(tmp3,1)-1
           lg = any(tmp3(jj+1:end,:)  ~= repmat(tmp3(jj,:),size(tmp3,1)-jj,1),2);
           tmp3 = tmp3([true(jj,1);lg],:);
           jj=jj+1;
       end
       count1 = count2;
       nn(count1+1:count1 + size(tmp3,1),:) = tmp3;
end
else
   nn=0; numwithn=1; %trivial case
end

tmp = zeros(N*size(nn,1),1);
tmp(1:N) = (rand(N,1) + 1i*rand(N,1))*10^(-7);  %init cond
for lp2=1:N
    rho_eles{lp2} = tmp;
end
if Kap2>0
    const_fct = -1i*nn*(drude_gam.');
else
    const_fct = 0; 
end

%% Calculate super operator! Excluded freq dependent and dipole dep
   
M = MM*(diag(ones(N-1,1),1) + diag(ones(N-1,1),-1));
%periodic boundary conditions mean that first and last are also coupled
%not that M(1,:) relates to the ground state only
M(1,end) = MM;  M(end,1) = MM;

sup_op = kron(sparse(1:sum(numwithn),1:sum(numwithn),...
                ones(sum(numwithn),1)),sparse(M)); %no cc needed

%% Calculate coupling to different tiers

 tmp = cumsum(numwithn);
down_coup_mat = zeros(tmp(end)); cnt1 =0; rng = 1;
for k = 2:Kap2+1
    
    cnt2 = cnt1 + numwithn(k);
    down_coup_mat(cnt1+1:cnt2,rng) = k;
    rng = cnt1+1:cnt2; %previous range
    cnt1= cnt2;
end
rnj = 1:length(reorg_E);
lambda_mat = sparse(rnj, rnj ,reorg_E);
gamma_mat = sparse(rnj, rnj , -2*Temp + 1i*drude_gam);

up_coup_mat = down_coup_mat; up_coup_mat=up_coup_mat';
up_coup_mat(up_coup_mat~=0) = -1; 

sup_op = sup_op + kron(sparse(up_coup_mat),lambda_mat);
sup_op = sup_op + kron(sparse(down_coup_mat),gamma_mat);
clear down_coup_mat up_coup_mat

    temp = kron(const_fct,ones(N,1));
    %identical constant factor 
    sup_op = sup_op - sparse(1:length(temp),1:length(temp),temp);
    
%% Loop over specific diagonalisations


rnj = 1:sum(numwithn); padmat = sparse(rnj,rnj,ones(size(rnj)));
%if any(reorg_E ~= reorg_E(1))
    
    SSden = zeros(N,N,length(freq_range)); chk = zeros(N,length(freq_range));
    OD = zeros(size(freq_range)); CD = OD; 
    
for lp1 = 1:length(freq_range)
    
    omega = freq_range(lp1);
    freq_fact = (om_0 + reorg_E -omega); 
    rnj = 1:length(freq_fact); 
    sup_op_freq = kron(padmat,sparse(rnj,rnj,freq_fact));
    HEOM_sub(-1i*(sup_op+sup_op_freq),[]); %pass this to the function
    for lp2 = 1:N         
        opts = optimset('Diagnostics','off', 'Display','off','tolX',1e-9);
      [tmp1,tmp2,tmp3] = fsolve(@(x) HEOM_sub(x,1),rho_eles{lp2},opts); %find zero value
      tmp3
      %tmp3 is exit flag
      rho_eles{lp2} = tmp2; %use as guess for next frequency
      chk(lp2,lp1) = norm(tmp2); %should be near to zero, stationary state
%        if lp1>=2
%            rho_eles = SSden(:,lp2,lp1-1); %use as new initial condition
%        end
      %occurs when eigenvalue is zero as this state doesn't change in time
      SSden(:,lp2,lp1) = tmp1(1:N,1);              
      
      OD(lp1) = OD(lp1) + ddotd(lp2,:)*tmp1(1:N,1);  
      %just imag of this relates to absorb but I might as well save it all
      CD(lp1) = CD(lp1) + CDop(lp2,:)*tmp1(1:N,1); 
    end
end
OD = OD.*freq_range;  CD = CD.*freq_range; 
% else
%     SSden = zeros(N,N); chk = zeros(N);
%     OD = zeros(N,1); CD = OD; 
%     for lp2 = 1:N   
%         
%       tmp = sparse(N,N); tmp(lp2,1) = -1; %#ok<SPRIX>
%       selector = kron(padmat,tmp);
%         opts.v0 = rho_eles;
% 
%       [tmp1,tmp2] = eigs(-1i*(sup_op+selector),1,0,opts);
%       
%       chk(lp2) = tmp2(1,1); %should be near to zero, stationary state
%       %occurs when eigenvalue is zero as this state doesn't change in time
%       SSden(:,lp2) = tmp1(1:N,1);              
%       
%       OD(lp2) = ddotd(lp2,:)*tmp1(1:N,1);  
%       %just imag of this relates to absorb but I might as well save it all
%       CD(lp2) = CDop(lp2,:)*tmp1(1:N,1); 
%     end
%     
%     
% end

end
function out = HEOM_sub(x,jj)

    persistent bigmat
        if isempty(jj)
           bigmat = x; out =[]; return 
        end

        out = bigmat*x;
        out(jj) = out(jj)+1i; 
end
