      
tim = sym('t');
U_tev = expm(-1i*tim*(kron(eye(length(Htot)),Htot)-kron(Htot.',eye(length(Htot)))));
%rowpos,colpos,expression needed to construct sparse operator

bigK = length(coup_com_save{1});
    cnt=0;
    rowpos = zeros(bigK^2*N,1); colpos=rowpos; 
    valu = sym(colpos); %preallocate as symbolic, might help? 
    
    cnt2=0; %this count relates to the number of 
    L2 = -1i*(kron(eye(length(Htot)),V_int)-kron(V_int.',eye(length(Htot))));
    L2 = U_tev*L2*(U_tev'); %interaction picture
    
    rowpos2 = zeros(bigK^2*N,1); colpos2=rowpos2; 
    valu2 = sym(colpos); %preallocate as symbolic, might help?     
     
    
for j = 1:N    

    Qjsec = zeros(N);  Qjsec(j,j) = 1;
    Qjsec = kron(Qjsec,eye(viblvls)); eyeQ = eye(length(Qjsec)); 
    
    trunc_opj = QQ(j,1).*(kron(Qjsec.',eyeQ) + kron(eyeQ,Qjsec) - 2*kron(Qjsec.',Qjsec))...
                + QQ(j,2) .*(kron(Qjsec.',eyeQ) - kron(eyeQ,Qjsec));
    trunc_opj = U_tev*trunc_opj*(U_tev'); %interaction picture       
    
    Qj = -1i*(kron(eyeQ,Qjsec) - kron(Qjsec.',eyeQ)); %commtator
    Qj = U_tev*Qj*(U_tev'); %interaction picture
    %keep factor of -1i from -1i tilde(V)^X term
    Uj =  (kron(eyeQ,Qjsec) + kron(Qjsec.',eyeQ));     
    Uj = U_tev*Uj*(U_tev'); %note U_tev is symbolic
    

    for k = 1:bigK
        for kk = 1:bigK
            if coup_com_save{j}(k,kk)~=0  || coup_acom_save{j}(k,kk)~=0 
                %non zero element is coupling
                for j1 = 1:N^2
                    for j2 = 1:N^2
                        if Uj(j1,j2) ~=0 
                            %will give symbolic for non zero entries
                            %which is interpreted as "true" by if
                cnt = cnt+1;
              rowpos(cnt) = (kk-1)*N^2+j2;
              colpos(cnt) = (k-1)*N^2+j1;
              valu(cnt) = Qj(j1,j2)*coup_com_save{j}(k,kk);
              valu(cnt) = valu(cnt) + Uj(j1,j2)*coup_acom_save{j}(k,kk);
              
                        end
                    end
                end
            end
            if k ==kk %on diagonal 
                for j1 = 1:N^2
                    for j2 = 1:N^2
                        if trunc_opj(j1,j2) ~=0 

                cnt = cnt+1;
              rowpos(cnt) = (kk-1)*N^2+j2;
              colpos(cnt) = (k-1)*N^2+j1;
              valu(cnt) = trunc_opj(j1,j2);
              
                        end
                        if j==1 
                            if L2(j1,j2) ~=0 
                                cnt2 = cnt2 + 1;
              rowpos2(cnt2) = (kk-1)*N^2+j2;
              colpos2(cnt2) = (k-1)*N^2+j1;
              valu2(cnt2) = L2(j1,j2);
                            end
                        end    
                           
                    end
                end
            end
        end
    end
end
%trim only ones that were included
rowpos = rowpos(1:cnt); colpos = colpos(1:cnt); valu = valu(1:cnt);
rowpos2 = rowpos2(1:cnt2); colpos2 = colpos2(1:cnt2); valu2 = valu2(1:cnt2);
%% collate duplicates, there will be some after all
bignum = bigK*N^2+1;
lg = sum([rowpos,colpos*bignum],2);
[lg2,lg3] = sort(lg); 
lg4 = diff(lg3)==0; %identical position entries
lg5 = true(size(rowpos)); %entires that will be included in final
jlast = -1; 

for j = find(lg4)
%note lg2 maps j back to original positions 
    if j-jlast == 1 %part of same block
        jlast = j; 
        valu(lg2(jblk_start)) = valu(lg2(jblk_start)) + valu(lg2(j+1));
        lg5(lg2(j+1)) = false; %this element is not needed
    else %start of new block
        if jlast >0 %not for the very first
           lg5(lg2(jlast+1)) = false;
        end
        jblk_start = j; jlast = j;
        valu(lg2(j)) = valu(lg2(j)) + valu(lg2(j+1)); %collate with next value
    end
end

rowpos = rowpos(lg5); colpos = colpos(lg5); valu = valu(lg5);

% to get the super operator sub in for time value and do
% sup_op = sparse(rowpos,colpos,double(sub(valu,time,t)))
% for the interaction with light one must take  
% L_plus = sparse(rowpos2,colpos2,double(sub(exp(-1i*omega*time)*valu,time,t)))
% L_minus = sparse(rowpos2,colpos2,double(sub(exp(+1i*omega*time)*valu,time,t)))
% the reason for plus and minus seeming counter intuitive is these raise /
% lower the index related to the factor of exp(i*l*dot(k,r))
