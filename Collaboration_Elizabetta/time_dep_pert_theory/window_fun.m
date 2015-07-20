function  [tmp_gg_sav, tmp_ee_sav, tmp_ef_sav] = window_fun(...
          des_freq_r,H_exciton,R_red_sc,supop_ge ,V_ge,V_eg,tmpge,supop_ef,V_ef,V_fe,tmpef,sz2) 
%%  Calculate the window function
% Assume the pulse is short enough that no time evolution occurs during the
% interaction.
% Doesn't actually calculate the window function, just the components
% before the integral over the electric field in frequency space are
% performed

mid_freq_r = mean( des_freq_r );
om_sr = -1i*mid_freq_r ;
om_f = om_sr; %lack of a better idea for the scaling
N = length(V_ge);
if nargin < 11
sz2 = length(V_ge{1})/(1+N+N*(N-1)/2); 
%vibrational levels, if you change %the system 
    if abs(round(sz2)-sz2)<eps(length(V_ge{1}))
        sz2 = length(V_ge{1})/(1+N); %assume just single excited
    end
    if abs(round(sz2)-sz2)<eps(length(V_ge{1}))
        warning('size of vibrational manifold not clear from input')
    end

end


t_step = 2*pi/(des_freq_r(end)-des_freq_r(1)); %set by range
t_end = 2*pi/(des_freq_r(2)-des_freq_r(1)); %set by frequency spacing
t3_range = 0:t_step:t_end ;
fft_sc_e = repmat(exp(om_sr*t3_range),length(supop_ge),1);     %scale middle freq back to centre fft
fft_sc_f = repmat(exp(om_f*t3_range),length(supop_ef),1);     %scale middle freq back to centre fft

for e4 = 1:N
    
    Vge_L = reshape(V_ge{e4},numel(V_ge{e4}),1)'; %conjugated
    Vge_L = Vge_L(tmpge); %reduced to elements that it can be mixed to
    %prop in time as usual but with backwards acting super op

    om_shift = H_exciton(1+e4*sz2,1+e4*sz2)-H_exciton(1,1); 
    
    om_ge = -om_shift;  

    Gamma_e = real(R_red_sc(1,1+e4,1,1+e4)); %~coherence decay rate

    num_scaling_ge = exp((-1i*om_ge-Gamma_e)*t3_range); 

    
    %num_scaling = exp(-Gamma_s*t3_range); 
    %scale out oscillations and decay from electronic DOF
    supop_ge_scaled = supop_ge  + (+1i*(om_ge)+Gamma_e)*sparse(eye(length(supop_ge)));

    de_fn_bk_ge = @(t,v) mtimesx(v,'T',supop_ge_scaled).';       
    %ode45 only takes column so can't do left acting directly
    
    output_DE_fun(length(t3_range),Vge_L,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);
    ode45(de_fn_bk_ge ,[t3_range(1),t3_range(end)],Vge_L,options);
    
    [tmp1,tmp2]  =  output_DE_fun(length(t3_range),Vge_L,'get_data+clean');             
    tmp3 = (interp1(tmp1,tmp2,t3_range).') ;    
     tmp3 = tmp3  .*repmat(num_scaling_ge,size(tmp2,2),1);
     
    %initial window state, time dimension two
    %fourier transform to freq space, shift zero freq to centre, note this
    %is ofset by om_s due to scaling
    tmp3 = fftshift(fft(tmp3.*fft_sc_e,[],2),2)/length(t3_range);   
     tmp3 = reshape(tmp3,sz2,sz2*N,length(t3_range));

     for f = 1:N*(N-1)/2
     Vef_L = reshape(V_ef{e4,f},numel(V_ef{e4,f}),1)'; %conjugated
     Vef_L = Vef_L(tmpef); %reduced to elements that it can be mixed to        

        Gamma_f = real(R_red_sc(e4+1,N+1+f,e4+1,N+1+f))/2; %~coherence decay rate
     om_shift = H_exciton(1+(N+f)*sz2,1+(N+f)*sz2)-H_exciton(1+e4*sz2,1+e4*sz2); 
    om_ef = +om_shift;      
    
    num_scaling_ef = exp((-1i*om_ef-Gamma_f)*t3_range); 
 
     %I can scale these elements with a different om
        supop_ef_scaled = supop_ef  + (om_sr - 1i*om_ef+Gamma_f)*sparse(eye(length(supop_fe))); 
        
        de_fn_bk_ef =  @(t,v) mtimesx(v,'T',supop_ef_scaled).';          
        
            output_DE_fun(length(t3_range),Vef_L,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);
    ode45(de_fn_bk_ef ,[t3_range(1),t3_range(end)],Vef_L,options);
    
    [tmp1,tmp2]  =  output_DE_fun(length(t3_range),Vef_L,'get_data+clean');             
    tmp4 = (interp1(tmp1,tmp2,t3_range,'pchip').')  ;    
    tmp4 = tmp4.*repmat(num_scaling_ef,size(tmp2,2),1);
    
    tmp4 = fftshift(fft(tmp4.*fft_sc_f,[],2),2)/length(t3_range);
    tmp4 = reshape(tmp4,sz2*N,sz2*N*(N-1)/2,length(t3_range));    

    tmp5{f} = tmp4;
    
    end
    for e3 = 1:N  %apply the other operator, left acting commutator
        
        op_3 = V_eg{e3};    op_3 = op_3(sz2+1:sz2*(N+1),1:sz2); 
        op_3f = V_fe{e3};   op_3f = op_3f(sz2*(N+1)+1:end,sz2+1:sz2*(N+1));           
        
        tmp_gg_sav{e3,e4} = mtimesx(tmp3,op_3);  %need also
        tmp_gg_sav{e3,e4} = tmp_gg_sav{e3,e4} +...
            conj(permute(tmp_gg_sav{e3,e4},[2,1,3]));
        tmp_ee_sav{e3,e4} = mtimesx(op_3,tmp3);       
        tmp_ee_sav{e3,e4} = tmp_ee_sav{e3,e4} +...
            conj(permute(tmp_ee_sav{e3,e4},[2,1,3]));
        
        %excited states
        for f= 1:N*(N-1)/2
        tmp_ef_sav{e3,e4,f} = mtimesx(tmp5{f},op_3f);
        tmp_ef_sav{e3,e4,f} = tmp_ef_sav{e3,e4,f} +...
            conj(permute(tmp_ef_sav{e3,e4,f},[2,1,3]));     
        end
        
 
        
    end
    
end