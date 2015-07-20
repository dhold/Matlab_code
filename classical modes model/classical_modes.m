function [time_units,rho_vec,convfact,tu2,rv2,exvibmix] = classical_modes
%Hel= [17113.9, 319.374, 30.4857, 7.5811;319.374, 17033.3, 20.3238, -43.8736;...
% 30.4857, 20.3238, 15807.4, 3.3873;7.5811, -43.8736, 3.3873, 16371.9];
%om_0 = {1108,1108,1108,1108};
%om_0 is the values of brownian modes included
%lambda ={44.32,44.32,44.32,44.32}; %reorganisation energy
%gamma = {5.30884,5.30884,5.30884,5.30884}; %damping of modes
%in general each of these can be a vector of different lengths, even empty
%omegavib = {1108,1108,1108,1108};%[1108,938,1111,1450];

    Temp =300; %units Kelvin
    light_speed = 299792458; %units m s^{-1}, NOT cm 
    length_unit = 100; % units m, this is 1/[1(cm)^{-1}]
    hbar = 1.05457173*10^(-34); %units Js
    boltz_const = 1.3806488 * 10^(-23); 

 E1 = 1042; %Energy of initial state (compared to excitation in other chromopore)
 V = 92; %Coupling
 H = [0,V;V,E1];
 [eigvec,~] = eig(H);
 omegavib = [1111;1111]; %Frequency of vibrations / level spacing,
 coupg = [267.1;267.1]; %coupling term between electronic and vibrational DOF
 gamma = [5.30884;5.30884];
 tendps = 2; %end time in pico seconds
 convfact = 2 * pi * light_speed * length_unit * 10^(-12);
 tend = tendps * convfact;
 gamma_dru = {};%{100,100}; %drude decay constant
 lambda_dru = {};%{60,60}; 

xi = zeros(length(omegavib),1); %position of oscillators
xi = [0.1;0.2];
yi = 0*xi;  %velocity
%choose some thermal distribution for these.  I.e. random expectation
%values of centre of mass wavefunction and phase gradients with a mean of
%zero but a varience of what would be expected

 coupled_sys_solve({1:numel(H),numel(H)+1:length(xi)+numel(H),...
            numel(H)+length(xi)+1:numel(H)+length(xi)+length(yi)},...
            {coupg,H,gamma,omegavib });
 rho0 = eigvec'*[0,0;0,1]*eigvec;
 rho_init = [reshape(rho0,4,1);xi;yi];  
 [time_units,rho_vec] =ode45(@coupled_sys_solve,[0,tend],rho_init);
%  if nargout >3  %cba making this work
%       coupled_sys_solve_2({},{coupg,H,gamma,omegavib,xi,yi });
%     [tu2,rv2] =ode45(@coupled_sys_solve_2,[0,tend],reshape(rho0,4,1)); 
%     xx2=coupled_sys_solve_2([],{});
%  end

%now solve fully quantum
    if nargout >3  
viblvls =4; %seems a reasonable number

occmat = zeros(viblvls^2,2);  Emat = zeros(viblvls^2,1);
n1=0; n2=0;  cnt = 1;
while n2 <= viblvls-1
    
    occmat(cnt,:) = [n1,n2];
    Emat(cnt) = dot([n1;n2]+1/2,omegavib);
    if n1==viblvls-1
        n1=0; n2=n2+1;
    else
        n1=n1+1;
    end
    cnt=cnt+1;

end

Hvibtot = diag(Emat);
%calculate mixing for b_1 and b_2
exvibmix = {zeros(viblvls^2,viblvls^2),zeros(viblvls^2,viblvls^2)};

for j = 1:viblvls^2
    for jj = 1:2
    tmp = occmat(j,:); fct = sqrt(tmp(jj)); tmp(jj) = tmp(jj)-1;
    lg = all(abs(occmat-repmat(tmp,size(occmat,1),1)) < eps(viblvls),2);
    exvibmix{jj}(j,lg) = fct;
    end
end

Htot = kron(H,eye(length(Hvibtot)))+kron(eye(length(H)),Hvibtot);

Htot(1:length(Htot)/2,1:length(Htot)/2) = coupg(1)*exvibmix{1}'+...
    Htot(1:length(Htot)/2,1:length(Htot)/2)+ coupg(1)*exvibmix{1};
Htot(length(Htot)/2+1:end,length(Htot)/2+1:end) = coupg(2)*exvibmix{2}'+...
    Htot(length(Htot)/2+1:end,length(Htot)/2+1:end)+ coupg(2)*exvibmix{2};

LL = sparse(-1i*(kron(eye(length(Htot)),Htot)-kron((Htot ).',eye(length(Htot)))));
% size(LL)
       coupled_sys_solve_quantum({1},LL);
       
       tmp = 0*Hvibtot; tmp(1,1)=1;
       rho02 = kron(rho0,tmp);
       
       
     [tu2,rv2] =ode45(@coupled_sys_solve_quantum,[0,tend],reshape(rho02,numel(rho02),1)); 

    end
end
 
function  drho=coupled_sys_solve_quantum(t,rho)

persistent  L 
        %first time give t = {density matrix relating to the electronic
        %occupation, exp oscillator displacements, exp oscillator
        %velocities} pass rho as Louivillian
        
        if iscell(t)
           
             L = rho; drho = [];   return
            
        end

        drho = L*rho;
        
 
end

function  drho=coupled_sys_solve(t,rho)

persistent sections H cpg om gam diagpic 
        %first time give t = {density matrix relating to the electronic
        %occupation, exp oscillator displacements, exp oscillator
        %velocities} pass rho as Louivillian
        
        if iscell(t)
           
            cpg = rho{1};  H = rho{2}; gam = rho{3}; om = rho{4};
            sections = t;
            diagpic = eye(length(H)); 
            diagpic = sections{1}(logical(diagpic));
                        drho = [];   return
            
        end

        drho = rho*0; 
        tmp = diag(cpg.*rho(sections{2}));
        tmp2 = sparse(-1i*(kron(eye(length(H)),H+tmp)-kron((H+tmp).',eye(length(H)))));
        drho(sections{1}) = tmp2*rho(sections{1});
        drho(sections{2}) = rho(sections{3}); %gradient of x is exactly velocity
        drho(sections{3}) = -om.*(om.*rho(sections{2}) + 2*cpg.*rho(diagpic)) ...
                            -gam.*rho(sections{3});
        
 
end

function  drho=coupled_sys_solve_2(t,rho)
%this version assumes that the density matrix doesn't change significantly
%during compared to the oscillator postion and can be assumed constant
%during a timestep.  Might need to write my own ODE solver to get this to
%actually work
% here rho is only the electronic degrees of freedom
persistent H cpg om gam diagpic xi xx tt
        %first time give t = {density matrix relating to the electronic
        %occupation, exp oscillator displacements, exp oscillator
        %velocities} pass rho as Louivillian
        
        if iscell(t)
           
            cpg = rho{1};  H = rho{2}; gam = rho{3}; om = rho{4};
             xx = rho{5}; tt =rho{6}; %initial values
            xi = sqrt(gam.^2/4-om.^2); %will be complex for underdamped
            diagpic = reshape(logical(eye(length(H))),numel(H),1); 
                        drho = [];   return
           
        end
        if iscell(rho)
            %pass out xx and tt
            drho = {xx,tt}; return
        end
        tstep = t-tt(end);
        xstep = exp(-gam/2-xi).*expm1((gam/2+xi)*tstep).*rho(diagpic) + ...
        exp(-gam/2+xi).*expm1((gam/2-xi)*tstep).*rho(diagpic);
        
        tmp = diag(cpg.*rho(diagpic).*(xx(:,end)+xstep));
        xx(:,end+1) = xx(:,end)+xstep; 
        tt = [tt;t];
        tmp2 = sparse(-1i*(kron(eye(length(H)),H+tmp)-kron((H+tmp).',eye(length(H)))));


        
        drho  = tmp2*rho;

        
 
end
 