function drho = HEOM_light_DE(tim,rho_vc) %Couples to same numwithn value

persistent tt viblvls QQ H1 H2 spmat decay_term U_tev coup_com_save coup_acom_save

        if iscell(rho_vc)
%note that symbolic functions must be in terms of symbolic parameter t
            tt = sym('t','real');

            H1 = rho_vc{1}{1};  %Time indep Hamiltonian, this will only 
            %matter if vibrations are included quantum mechanically
            %**write fn that doesn't have this for when I don't want to**
            H2 = rho_vc{1}{2}; %Time dependent Hamiltonian,                                    
            QQ = rho_vc{1}{3};  %params for truncation operator 
            U_tev = rho_vc{1}{4}; %time evolution operator
            viblvls= rho_vc{1}{5}; %vibrational lvls

            coup_com_save = rho_vc{2}{1}; %Heirarchy terms
            coup_acom_save = rho_vc{2}{2};
            spmat = rho_vc{2}{3}; %big sparse matrix
            decay_term = rho_vc{2}{4}; %does not need projection to interaction             
            
            drho =[];

            return           
        end
        %outtest = [t,max(max(abs(rho_vc)))]
        %construct propogation matricies at these times
 % tic

   H2sub = kron(double(subs(H2,tt,tim)),eye(viblvls));   %this is already in the interaction basis
 L2 = -1i*(kron(eye(length(H2sub)),H2sub)-kron(H2sub',eye(length(H2sub))));
 
tmp_var = double(subs(U_tev ,tt,tim)); %subs is slow, don't use often
tmp_var = kron(tmp_var,eye(viblvls)); %expand to full basis

 LL = -1i*(kron(eye(length(H1)),tmp_var'*H1*tmp_var)-...
        kron((tmp_var'*H1*tmp_var).',eye(length(H1))));

       total_prop_op  = kron(spmat,LL+L2) + decay_term;          
        trunc_op = zeros(size(LL)); %truncation operator
        eyeQ = sparse(eye(size(H2sub)));
    for j = 1 : length(coup_com_save) %now include upcoupling etc etc
        
        Qjsec = zeros(size(H2));  Qjsec(j,j) = 1;
        Qjsec = tmp_var'*kron(Qjsec,eye(viblvls))*tmp_var ; 
        tmp_var2 = -1i*sparse(kron(eyeQ,Qjsec) - kron(Qjsec.',eyeQ)); %commtator
        tmp_var3 =  sparse(kron(eyeQ,Qjsec) + kron(Qjsec.',eyeQ));    

        trunc_op = trunc_op + QQ(j,1).*tmp_var2*tmp_var2+ QQ(j,2).*tmp_var3*tmp_var2; 
        %check the complex coefficient of QQ(j,2) when I actually use this
    if ~isempty(coup_com_save{j})
    total_prop_op = total_prop_op + kron(sparse(coup_com_save{j}),tmp_var2);
    end
    if ~isempty(coup_acom_save{j})
    total_prop_op = total_prop_op + kron(sparse(coup_acom_save{j}),tmp_var3);
    end
    end
    
        total_prop_op = total_prop_op + kron(spmat,trunc_op); 
        %truncation correction
%toc
        %now actually... calculate 
    drho = total_prop_op*rho_vc;
end