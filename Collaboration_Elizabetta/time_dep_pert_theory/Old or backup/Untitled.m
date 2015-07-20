%% Calculate first order polarization with HEOM

t_end_ps = 2;
t_end = convfact*t_end_ps; numpoints =1000;
t_lin = linspace(0,t_end,numpoints);

Kap1 = [];
freq_scale =1000;
H0scaled = H0 - eye(size(H0))*freq_scale;
%numerically easier if you take out a factor out of the
%Hamiltonian to reduce the oscillatory components, this is like adding a
%factor of exp(-1i*freq_scale*t) into it

%compute evolution of density matricies with just one element non zero,
%scaled to unity e.g. [0,0;1,0] to work stuff out later
clear evo_1 evo_2
for k = 1:N 
    
 tmp = zeros(N+1);   tmp(k+1,1) = 1; tmp(1,k+1) = 1; 
[out1,t_range1,out2,t_range2,ftharg] = reduced_HEOM(init_cond10{k},...
    t_end,H0scaled,QQ(2:end,:),{cc_com{2},cc_com{3}},...
    {cc_acom{2},cc_acom{3}},vv,Kap1,3,numpoints);
%reduced HEOM only have to propogate a system which is proportial to the 
%number of excited states times the size of the heirarchy where
%as normal would be N^2 times this
                            if isempty(Kap1)
                                Kap1=ftharg;
                            end
                            
out1 = out1(lg1,:); t_range1 = t_range1(lg1); 
out2 = out2(lg2,:); t_range2 = t_range2(lg2); 

evo_1{k} = interp1(t_range1,out1,t_lin,'pchip');                           
evo_2{k} = interp1(t_range2,out2,t_lin,'pchip'); 

end

FO_pol_HEOM = zeros(3,numpoints);

rho0 = [1,0,0;0,0,0;0,0,0];

syms mu_1 mu_2 mm_1 mm_2 f_11 f_21 f_12 f_22
syms g_11 g_21 g_12 g_22
test1 =mu_1*[0,0,0;f_11,0,0;f_21,0,0]+ mu_2*[0,0,0;f_12,0,0;f_22,0,0]+...
-mu_1*[0,g_11,g_21;0,0,0;0,0,0] -  mu_2*[0,g_21,g_22;0,0,0;0,0,0]  ;
test2 =mm_1*[0,0,0;f_11,0,0;f_21,0,0]+ mm_2*[0,0,0;f_12,0,0;f_22,0,0]+...
mm_1*[0,g_11,g_21;0,0,0;0,0,0] +  mm_2*[0,g_21,g_22;0,0,0;0,0,0]  ;
%signs on commutators cancels as mm is complex

mu_sym = [0,mu_1,mu_2;mu_1,0,0;mu_2,0,0];
mm_sym = [0,-mm_1,-mm_2;mm_1,0,0;mm_2,0,0];

test11 = trace(mu_sym*test1);
test11 = subs(test11,{f_11,f_12,f_21,f_22,g_11,g_12,g_21,g_22},...
    {evo_1{1}(:,1),evo_1{1}(:,2),evo_1{2}(:,1),evo_1{2}(:,2),...
    evo_2{1}(:,1),evo_2{1}(:,2),evo_2{2}(:,1),evo_2{2}(:,2)});
test22 = trace(mu_sym*test2);
test22 = subs(test22,{f_11,f_12,f_21,f_22,g_11,g_12,g_21,g_22},...
    {evo_1{1}(:,1),evo_1{1}(:,2),evo_1{2}(:,1),evo_1{2}(:,2),...
    evo_2{1}(:,1),evo_2{1}(:,2),evo_2{2}(:,1),evo_2{2}(:,2)});

for k = 1:3
   FO_pol_HEOM(3,:) = FO_pol_HEOM(3,:)+subs(test11,{mu_1,mu_2},{mu(1,k),mu(2,k)}).';
    FO_pol_HEOM(2,:) = FO_pol_HEOM(2,:)+subs(test22,{mu_1,mu_2,mm_1,mm_2},...
                        {mu(1,k),mu(2,k),mm_tot(1,k),mm_tot(2,k)}).';
end

temp = mu_sym*rho0 - rho0*mu_sym ; 
temp2 = mm_sym*rho0 - rho0*mm_sym ; 
for k =1:N
    tmp1{k} = temp(k+1,1)*evo_1{k}
end