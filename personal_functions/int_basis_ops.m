function [V_ge,V2_ge,V_ef,V2_ef,V_ge_ex,V2_ge_ex,V_ef_ex,V2_ef_ex,...
    V_ge_full,V2_ge_full,V_ef_full,V2_ef_full,V_ge_ex_full,V2_ge_ex_full,V_ef_ex_full,V2_ef_ex_full]=...
        int_basis_ops(output_type,fock_space_rep,mu,mdip,M_e,M_f,totvib,M_e_full,M_f_full);
    
%this function calculates the vector operator that projects the interaction
%operator components in the exciton basis and also padded out with
%vibrational levels, and in the exciton vibrational basis
%if output_type is given as 'amp_angle' this will convert the outputs to 
%spherical polars, which can help with the averaging
N = size(mu,1);

V_ge = mu; V2_ge = mdip;
V_ge = reshape(V_ge,N,1,3); V2_ge = reshape(V2_ge,N,1,3);
if strcmp(output_type,'amp_angle') %output as amplitude theta, phi
[az,elev,r] = cart2sph(V_ge(:,:,1),V_ge(:,:,2),V_ge(:,:,3));
V_ge = cat(3,r,az,elev); %convert to spherical polars
[az,elev,r] = cart2sph(V2_ge(:,:,1),V2_ge(:,:,2),V2_ge(:,:,3));
V2_ge = cat(3,r,az,elev);
end
% x = r .* cos(elevation (dim 3)) .* cos(azimuth (dim 2))
% y = r .* cos(elevation) .* sin(azimuth)
% z = r .* sin(elevation)

V_ef = zeros(N*(N-1)/2,N,3);  V2_ef = zeros(N*(N-1)/2,N,3); 

for lp =1:N
   
    V = zeros(N*(N-1)/2,N); 
    
    lg = fock_space_rep(N+2:end,lp)==1; %elements with excitation at lp
    lg2 = fock_space_rep(N+2:end,:); lg2(:,lp) = 0; %[~,lg2] = find(lg2(lg,:)); 
    V(lg,1:N) = lg2(lg,:);

    for j = 1:3
    V_ef(:,:,j)  = V_ef(:,:,j) + mu(lp,j)*V;
    V2_ef(:,:,j) = V2_ef(:,:,j) + mdip(lp,j)*V;
    end
end

if strcmp(output_type,'amp_angle') %output as amplitude theta, phi
[az,elev,r] = cart2sph(V_ef(:,:,1),V_ef(:,:,2),V_ef(:,:,3));
V_ef = cat(3,r,az,elev); %convert to spherical polars
[az,elev,r] = cart2sph(V2_ef(:,:,1),V2_ef(:,:,2),V2_ef(:,:,3));
V2_ef = cat(3,r,az,elev);
end


% now with vibrational levels added
if nargin >=7
    
    V_ge_full = zeros(N*totvib,totvib,3); V2_ge_full = V_ge_full; 
    V_ef_full = zeros(N*totvib*(N-1)/2,N*totvib,3);
    V2_ef_full = zeros(N*totvib*(N-1)/2,N*totvib,3);
    for j = 1:3
        V_ge_full(:,:,j) = kron(V_ge(:,:,j),eye(totvib));
        V2_ge_full(:,:,j) = kron(V2_ge(:,:,j),eye(totvib));
        V_ef_full(:,:,j) = kron(V_ef(:,:,j),eye(totvib));
        V2_ef_full(:,:,j) = kron(V2_ef(:,:,j),eye(totvib));    
    end
    
end

if strcmp(output_type,'amp_angle') %output as amplitude theta, phi
[az,elev,r] = cart2sph(V_ge_full(:,:,1),V_ge_full(:,:,2),V_ge_full(:,:,3));
V_ge_full = cat(3,r,az,elev); %convert to spherical polars
[az,elev,r] = cart2sph(V2_ge_full(:,:,1),V2_ge_full(:,:,2),V2_ge_full(:,:,3));
V2_ge_full = cat(3,r,az,elev);
[az,elev,r] = cart2sph(V_ef_full(:,:,1),V_ef_full(:,:,2),V_ef_full(:,:,3));
V_ef_full = cat(3,r,az,elev); %convert to spherical polars
[az,elev,r] = cart2sph(V2_ef_full(:,:,1),V2_ef_full(:,:,2),V2_ef_full(:,:,3));
V2_ef_full = cat(3,r,az,elev);
end

% in exciton basis
if ~isempty(M_e)
try %try using mtimesx for speed
    V_ge_ex = mtimesx(M_e,'C',V_ge); V2_ge_ex = mtimesx(M_e,'C',V2_ge);
    V_ef_ex = mtimesx(mtimesx(M_f,'C',V_ef),M_e); 
    V2_ef_ex = mtimesx(mtimesx(M_f,'C',V2_ef),M_e);
    if nargin >=6
    V_ge_ex_full = mtimesx(M_e_full,'C',V_ge_full); 
    V2_ge_ex_full = mtimesx(M_e_full,'C',V2_ge_full);
    V_ef_ex_full = mtimesx(mtimesx(M_f_full,'C',V_ef_full),M_e_full); 
    V2_ef_ex_full = mtimesx(mtimesx(M_f_full,'C',V2_ef_full),M_e_full);  
    end
catch ME %probably mtimesx isn't compiled or something
    ME
    V_ge_ex = zeros(size(V_ge));  V2_ge_ex = V_ge_ex;
    V_ef_ex =  V2_ge_ex ;  V2_ef_ex =  V2_ge_ex ; 
    for j = 1:3
    V_ge_ex(:,:,j)  = M_e'*V_ge(:,:,j);
    V2_ge_ex(:,:,j) = M_e'*V2_ge(:,:,j);
    V_ef_ex(:,:,j)  = M_f'*V_ef(:,:,j)*M_e;
    V2_ef_ex(:,:,j) = M_f'*V2_ef(:,:,j)*M_e;
    if  nargin >=7
    V_ge_ex_full(:,:,j)  = M_e_full'*V_ge_full(:,:,j) ; 
    V2_ge_ex_full(:,:,j)  = M_e_full'*V2_ge_full(:,:,j) ;
    V_ef_ex_full(:,:,j)  = M_f_full'*V_ef_full(:,:,j) *M_e_full; 
    V2_ef_ex_full(:,:,j)  = M_f_full'*V2_ef_full(:,:,j) *M_e_full;         
    end
    end
end

if strcmp(output_type,'amp_angle') %output as amplitude theta, phi
    
[az,elev,r] = cart2sph(V_ge_ex(:,:,1),V_ge_ex(:,:,2),V_ge_ex(:,:,3));
V_ge_ex = cat(3,r,az,elev); %convert to spherical polars
[az,elev,r] = cart2sph(V2_ge_ex(:,:,1),V2_ge_ex(:,:,2),V2_ge_ex(:,:,3));
V2_ge_ex = cat(3,r,az,elev);
[az,elev,r] = cart2sph(V_ef_ex(:,:,1),V_ef_ex(:,:,2),V_ef_ex(:,:,3));
V_ef_ex = cat(3,r,az,elev); %convert to spherical polars
[az,elev,r] = cart2sph(V2_ef_ex(:,:,1),V2_ef_ex(:,:,2),V2_ef_ex(:,:,3));
V2_ef_ex = cat(3,r,az,elev);      
    
[az,elev,r] = cart2sph(V_ge_ex_full(:,:,1),V_ge_ex_full(:,:,2),V_ge_ex_full(:,:,3));
V_ge_ex_full = cat(3,r,az,elev); %convert to spherical polars
[az,elev,r] = cart2sph(V2_ge_ex_full(:,:,1),V2_ge_full(:,:,2),V2_ge_ex_full(:,:,3));
V2_ge_ex_full = cat(3,r,az,elev);
[az,elev,r] = cart2sph(V_ef_ex_full(:,:,1),V_ef_ex_full(:,:,2),V_ef_ex_full(:,:,3));
V_ef_ex_full = cat(3,r,az,elev); %convert to spherical polars
[az,elev,r] = cart2sph(V2_ef_ex_full(:,:,1),V2_ef_ex_full(:,:,2),V2_ef_ex_full(:,:,3));
V2_ef_ex_full = cat(3,r,az,elev);
end

end

%can average this in the usual way by selecting the appropriate
%coefficients etc etc