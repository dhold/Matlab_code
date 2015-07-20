function [tsteps2,denatt2,tsteps,denatt]=test_efield_int_solve
%test to solve optical equations for simple system
% system iteracts with a Gaussian pulse of freq omega
%Assume units cm^-1
 %width


sigma = 10;
rho0 = [1,0;0,0];
rho0 = reshape(rho0,numel(rho0),1);
%peak of pulse hits at t=0, go 6 sigma before and after
t0 = -6*sigma;
tend = 6*sigma;

tic
[tsteps2,denatt2]=ode45(@pulse_with_system_2,[t0,tend],rho0);
toc
if nargout >2
tic
[tsteps,denatt]=ode45(@pulse_with_system_3,[t0,tend],rho0);
toc
end
end
function drho = pulse_with_system(t,rho)

    sigma = 10;
    Emax = 5/sqrt(2*pi*sigma); %field at peak
    omega = 10000;
    mu = [0,1;1,0];
    H0 = [0,0;0,10000];
    Htot = H0 + mu*Emax*exp(-(t/2/sigma)^2)*cos(omega*t);
    %Htot = H0 + mu*(eye(size(mu))*Emax*exp(-(t/2/sigma)^2)*exp(-1i*omega*t)/2);

    L = -1i*(kron(eye(length(Htot)),Htot)-kron(Htot.',eye(length(Htot))));

    drho = L*rho;

end
function drho = pulse_with_system_2(t,rho)

persistent L0 L1 sigma omega t0 L0proj gamma

            if isempty(L0)
                t0 = t;
                sigma = 10;
                Emax = 0.7/sqrt(2*pi*sigma); %field at peak
                omega = 9999.9;
                H0 = [0,0;0,10000];
                L0 = -1i*(kron(eye(length(H0)),H0)-kron(H0.',eye(length(H0))));
                % U(t,t_0) is matrix exponential of this fn, in order to
                % make computing this easier diagonalise it

                mu = [0,1;1,0];
                L1 = -1i*Emax*(kron(eye(length(mu)),mu)-kron(mu.',eye(length(mu))));
                gam = 10^(-8);
                gamma = [0,0,0,gam;0,-gam/2,0,0;0,0,-gam/2,0;0,0,0,-gam];
            end


            %LI = expm(-L0*(t-t0)) *( L1*(exp(-(t/2/sigma)^2)*cos(omega*t)) )* expm(L0*(t-t0));
            Utt0 = expm(L0*(t-t0));
            LI = Utt0' *( L1*(exp(-(t/2/sigma)^2)*cos(omega*t)) )* Utt0;
            %LI =  L1*(exp(-(t/2/sigma)^2) )/2;
            drho = (LI+gamma*(t-t0))*rho ; %rho in interaction picture
            
end

function drho = pulse_with_system_3(t,rho)
%makes RWA
persistent L1 sigma om0 delta t0 gamma

            if isempty(L1)
                t0 = t;
                sigma = 10;
                Emax = 0.7/sqrt(2*pi*sigma); %field at peak
                om0 = 10000; %resonant frequency 
                delta = 0.1;  %detuning
               
                mu = [0,1;1,0];
                L1 = -1i*Emax*(kron(eye(length(mu)),mu)-kron(mu.',eye(length(mu))));
                gam = 10^(-8);
                gamma = [0,0,0,gam;0,-gam/2,0,0;0,0,-gam/2,0;0,0,0,-gam];
            end
          
            u0RW = (diag(exp([0,1i*delta*(t-t0),-1i*delta*(t-t0),0]))...
                    + diag(exp([0,-inf,-inf,0])))/2;
           % u0RW = (diag(exp([0,1i*delta*(t-t0),-1i*delta*(t-t0),0]))...
           %   + diag(exp([0,-1i*(2*om0-delta)*(t-t0),1i*(2*om0-delta)*(t-t0),0])))/2;
            LI = u0RW*L1*(exp(-(t/2/sigma)^2) )*u0RW';
            drho = (LI+gamma*(t-t0))*rho ; %rho in interaction picture
            
end