%% all real, PP setup, x-pump
mu_sym = sym('mu',[4,3]); mu_sym = sym(mu_sym,'real');

%E_vec1 = [[1,0,0];[1,0,0];[1,1i,0]/sqrt(2);[1,-1i,0]/sqrt(2)];
E_vec1 = [[1,0,0];[1,0,0];[1,1i,0]/sqrt(2);[1,1i,0]/sqrt(2)];
%E_vec2 = [[1,0,0];[1,0,0];[1,-1i,0]/sqrt(2);[1,1i,0]/sqrt(2)];
E_vec2 = [[1,0,0];[1,0,0];[1,-1i,0]/sqrt(2);[1,-1i,0]/sqrt(2)];

tmp1 = tensor_av(mu_sym,E_vec1);
tmp2 = tensor_av(mu_sym,E_vec2);

tmp3 = (tmp1+tmp2)/2; tmp3x = simplify(tmp3);
tmp4 = (tmp1-tmp2); tmp4x = simplify(tmp4);

%% 5th order
R_sym = sym('R',[4,3]);  R_sym = sym(R_sym,'real');
ku_sym = sym('ku',[1,3]);  ku_sym = sym(ku_sym,'real');
%kr_sym = sym('kr',[1,3]);  kr_sym = sym(kr_sym,'real');
%k_sym = [ku_sym;-ku_sym;kr_sym;-kr_sym];
k_sym = [ku_sym;-ku_sym;[0,0,1];-[0,0,1]];

temp1 = sym(0); temp2 = sym(0);
for j = 1:4
temp1 =temp1 +  tensor_av([mu_sym;R_sym(j,:)],[E_vec1;1i*k_sym(j,:)]);
temp2 =temp2 +  tensor_av([mu_sym;R_sym(j,:)],[E_vec2;1i*k_sym(j,:)]);
end

temp3 = (temp1+temp2)/2; temp3x = simplify(temp3);
temp4 = (temp1-temp2); temp4x = simplify(temp4);


%% all real, PP setup w/ y comp
mu_sym = sym('mu',[4,3]); mu_sym = sym(mu_sym,'real');
E_sym = sym('E',[1,3]);  E_sym = sym(E_sym,'real');
E_vec1 = [E_sym;conj(E_sym);[1,1i,0]/sqrt(2);[1,-1i,0]/sqrt(2)];
E_vec2 = [E_sym;conj(E_sym);[1,-1i,0]/sqrt(2);[1,1i,0]/sqrt(2)];

tmp1 = tensor_av(mu_sym,E_vec1);
tmp2 = tensor_av(mu_sym,E_vec2);

tmp3 = (tmp1+tmp2)/2;
tmp4 = (tmp1-tmp2);

%% 5th order
R_sym = sym('R',[4,3]);  R_sym = sym(R_sym,'real');
ku_sym = sym('ku',[1,3]);  ku_sym = sym(ku_sym,'real');
kr_sym = sym('kr',[1,3]);  kr_sym = sym(kr_sym,'real');
k_sym = [ku_sym;-ku_sym;kr_sym;-kr_sym];

temp1 = sym(0); temp2 = sym(0);
for j = 1:4
temp1 =temp1 +  tensor_av([mu_sym;R_sym(j,:)],[E_vec1;1i*k_sym(j,:)]);
temp2 =temp2 +  tensor_av([mu_sym;R_sym(j,:)],[E_vec2;1i*k_sym(j,:)]);
end

temp3 = (temp1+temp2)/2; temp3 = simplify(temp3);
temp4 = (temp1-temp2); temp4 = simplify(temp4);

%%
pretty(tmp3x)
%%
pretty(temp4x)

%%
pretty(tmp4x)
%%
pretty(temp3x)

%% 
pretty(tmp3)
%% 
pretty(tmp4)
