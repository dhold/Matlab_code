function [dipole_av_diag,dipole_av_diag_f,mag_av_diag,mag_av_diag_f,...
             dipole_av_coh ,  mag_av_coh]= ori_av_fn_w_vib(N,Nv,N2, ...
        mu_ex,m_ex,mu_ex2,m_ex2,conj_lg,pol_set,mag_set,calc_OD,calc_coh)
%calculates dipole averages which will be needed for third order
%spectroscopy with additional magnetic and well as electric dipole
%interactions.

% N is the number of 1 exciton states
% N2 is number of 2 exciton states
% Nv is number of vibrational states

if ~exist('calc_OD','var') %calculate mixings of the form V_j1 rho V_j2
    calc_OD = false;
end
if ~exist('cal_coh','var') %include averages relevant for coherence contributions
    calc_coh = false;
end

mpol_comb{1} =  cond_conj([mag_set(1,:);pol_set(2:4,:)],conj_lg);
mpol_comb{2} =  cond_conj([pol_set(1,:);mag_set(2,:);pol_set(3:4,:)],conj_lg);
mpol_comb{3} =  cond_conj([pol_set(1:2,:);mag_set(3,:);pol_set(4,:)],conj_lg);
mpol_comb{4} =  cond_conj([pol_set(1:3,:);mag_set(4,:)],conj_lg);
pol_set = cond_conj(pol_set,conj_lg);

if ~calc_OD
    
    dipole_av_diag = zeros(Nv,N*Nv,N*Nv); dipole_av_diag_f =  zeros(Nv,N*Nv,N*Nv,N2*Nv);
    mag_av_diag = zeros(Nv,N*Nv,N*Nv); mag_av_diag_f =  zeros(Nv,N*Nv,N*Nv,N2*Nv);
     if calc_coh
         dipole_av_coh = zeros(Nv,N*Nv,N*(N-1)*Nv/2);
         mag_av_coh = zeros(Nv,N*Nv,N*(N-1)*Nv/2);
     else
         dipole_av_coh = [];
         mag_av_coh = [];         
     end



for g_lp = 1:Nv %first loop parameter goes over the possible ground state
% initial conditions, these are assumed to be vib eigenstates
for j = 1:N*Nv %excited state ex-vibrational states
    
    if calc_coh %calculate V_{g,j} V_{j,f} V_{f,j} V_{j,g'} shit       
        for f = 1:N*(N-1)*Nv/2
            
            dipole_av_coh(g_lp,j,f) = tensor_av(cond_conj([mu_ex(g_lp,j,:);mu_ex2(f,j,:)...
                                  ;mu_ex2(f,j,:);mu_ex(g_lp,j,:)],conj_lg), pol_set);
            mag_av_coh(g_lp,j,f) = tensor_av(cond_conj([m_ex(g_lp,j,:);mu_ex2(f,j,:)...
                                  ;mu_ex2(f,j,:);mu_ex(g_lp,j,:)],conj_lg),mpol_comb{1})+...
      tensor_av(cond_conj([mu_ex(g_lp,j,:);m_ex2(f,j,:);mu_ex2(f,j,:)...
                                  ;mu_ex(g_lp,j,:)],conj_lg),mpol_comb{2})+...
      tensor_av(cond_conj([mu_ex(g_lp,j,:);mu_ex2(f,j,:);m_ex2(f,j,:)...
                                  ;mu_ex(g_lp,j,:)],conj_lg), mpol_comb{3})+...
     tensor_av(cond_conj([mu_ex(g_lp,j,:);mu_ex2(f,j,:);mu_ex2(f,j,:)...
                                  ;m_ex(g_lp,j,:)],conj_lg),mpol_comb{4});       
        end
    end
    
    for j2 = 1:N*Nv
        
        mu_set = [mu_ex(g_lp,j,:);mu_ex(g_lp,j,:);mu_ex(g_lp,j2,:);mu_ex(g_lp,j2,:)];
        
dipole_av_diag(g_lp,j,j2) = tensor_av(cond_conj(mu_set,conj_lg),pol_set);

       m_set1 = [m_ex(g_lp,j,:);mu_ex(g_lp,j,:);mu_ex(g_lp,j2,:);mu_ex(g_lp,j2,:)]; 
       m_set2 = [mu_ex(g_lp,j,:);m_ex(g_lp,j,:);mu_ex(g_lp,j2,:);mu_ex(g_lp,j2,:)]; 
       m_set3 = [mu_ex(g_lp,j,:);mu_ex(g_lp,j,:);m_ex(g_lp,j2,:);mu_ex(g_lp,j2,:)]; 
       m_set4 = [mu_ex(g_lp,j,:);mu_ex(g_lp,j,:);mu_ex(g_lp,j2,:);m_ex(g_lp,j2,:)]; 
       
mag_av_diag(g_lp,j,j2) = tensor_av(cond_conj(m_set1, conj_lg),mpol_comb{1})+...
      tensor_av(cond_conj(m_set2,conj_lg),mpol_comb{2})+...
      tensor_av(cond_conj(m_set3,conj_lg),mpol_comb{3})+...
      tensor_av(cond_conj(m_set4,conj_lg),mpol_comb{4});
  
                    for f=1:N2*Nv
      mu_set = [mu_ex(g_lp,j,:);mu_ex(g_lp,j,:);mu_ex(f,j2,:);mu_ex(f,j2,:)];                  
                        
dipole_av_diag_f(g_lp,j,j2,f)  = tensor_av(cond_conj(mu_set,conj_lg),pol_set);

       m_set1 = [m_ex(g_lp,j,:);mu_ex(g_lp,j,:);mu_ex(f,j2,:);mu_ex(f,j2,:)]; 
       m_set2 = [mu_ex(g_lp,j,:);m_ex(g_lp,j,:);mu_ex(f,j2,:);mu_ex(f,j2,:)]; 
       m_set3 = [mu_ex(g_lp,j,:);mu_ex(g_lp,j,:);m_ex(f,j2,:);mu_ex(f,j2,:)]; 
       m_set4 = [mu_ex(g_lp,j,:);mu_ex(g_lp,j,:);mu_ex(f,j2,:);m_ex(f,j2,:)]; 

mag_av_diag_f(g_lp,j,j2,f) = tensor_av(cond_conj(m_set1, conj_lg),mpol_comb{1})+...
      tensor_av(cond_conj(m_set2,conj_lg),mpol_comb{2})+...
      tensor_av(cond_conj(m_set3,conj_lg),mpol_comb{3})+...
      tensor_av(cond_conj(m_set4,conj_lg),mpol_comb{4});                 
                    
                    end
    end
end
end
else
    warning('not yet finished')
   dipole_av_diag = zeros(Nv,N*Nv,N*Nv,N*Nv,N*Nv); 
   dipole_av_diag_f =  zeros(Nv,N*Nv,N*Nv,N*Nv,N*Nv,N*Nv*(N-1)/2);
mag_av_diag = zeros(Nv,N*Nv,N*Nv,N*Nv,N*Nv); 
mag_av_diag_f =  zeros(Nv,N*Nv,N*Nv,N*Nv,N*Nv,N*Nv*(N-1)/2);
    if calc_coh %calculate V_{g,j} V_{j,f} V_{f,j} V_{j,g} shit       

            dipole_av_coh = zeros(N*Nv,N*Nv*(N-1)/2,N*Nv*(N-1)/2,N*Nv);
            mag_av_coh = zeros(N*Nv,N*Nv*(N-1)/2,N*Nv*(N-1)/2,N*Nv);
  
    end

for j = 1:N
    for j2 = 1:N
        for j3=1:N
            for j4=1:N
dipole_av_diag(j,j2,j3,j4) = tensor_av(cond_conj([mu_ex(j,:);mu_ex(j2,:);mu_ex(j3,:);mu_ex(j4,:)],conj_lg),...
                        pol_set);
                    
mag_av_diag(j,j2,j3,j4) = tensor_av(cond_conj([m_ex(j,:);mu_ex(j2,:);mu_ex(j3,:);mu_ex(j4,:)],conj_lg),...
                        mpol_comb{1})+...
                    tensor_av(cond_conj([mu_ex(j,:);m_ex(j2,:);mu_ex(j3,:);mu_ex(j4,:)],conj_lg),...
                        mpol_comb{2})+...
                    tensor_av(cond_conj([mu_ex(j,:);mu_ex(j2,:);m_ex(j3,:);mu_ex(j4,:)],conj_lg),...
                        mpol_comb{3})+...
                    tensor_av(cond_conj([mu_ex(j,:);mu_ex(j2,:);mu_ex(j3,:);m_ex(j4,:)],conj_lg),...
                        mpol_comb{4});
                    
                    for f=1:N2
dipole_av_diag_f(j,j2,j3,j4,f)  = tensor_av(cond_conj([mu_ex(j,:);mu_ex(j2,:);...
                                 mu_ex2(j3,:,f);mu_ex2(j4,:,f)],conj_lg),pol_set);
                    
mag_av_diag_f(j,j2,j3,j4,f) = tensor_av(cond_conj([m_ex2(j,:);mu_ex(j2,:);mu_ex2(j3,:,f);mu_ex2(j4,:,f)],conj_lg),...
                        mpol_comb{1})+...
                    tensor_av(cond_conj([mu_ex(j,:);m_ex2(j2,:);mu_ex2(j3,:,f);mu_ex2(j4,:,f)],conj_lg),...
                        mpol_comb{2})+...
                    tensor_av(cond_conj([mu_ex(j,:);mu_ex2(j2,:);m_ex2(j3,:,f);mu_ex2(j4,:,f)],conj_lg),...
                        mpol_comb{3})+...
                    tensor_av(cond_conj([mu_ex(j,:);mu_ex2(j2,:);mu_ex2(j3,:,f);m_ex2(j4,:,f)],conj_lg),...
                        mpol_comb{4});                    
                    
                    end
            end
        end
    end
end    
end
end

function out=cond_conj(x_set,lg)

    out  = squeeze(x_set);
    out(lg,:) = conj(x_set(lg,:));

end