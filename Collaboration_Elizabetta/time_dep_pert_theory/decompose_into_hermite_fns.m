
%unfortunately GH quadrature cannot be used iteratively, increase
%quad_order until convergence is reached
  quad_order= 20;   tau = tau_u;
   
  [omm, weig] = GaussHermite(quad_order,eps(20));
  %[omm1,omm2] = meshgrid(omm,omm);

  weights = weig*weig'; %quad_order by quad_order 
  %tabulate hermite ppolynomials
  %nmax is max order of poly
  nmax = 20;
  
  hpol = zeros(length(omm),nmax+1);
  hpol(:,1) = omm.^(0)/sqrt(sqrt(pi));
  hpol(:,2) = -sqrt(2).*omm.*hpol(:,1);
  
  for k = 2:nmax
      
      hpol(:,k+1) = -sqrt(2/k).*omm.*hpol(:,k) + sqrt((k-1)/k).*hpol(:,k-1);
      
  end
  
  sym_to_fn('temp_fn.m',vpa(S3tilde2,16),[om,om1,om2])
  
  %make function which sums over all quadrature points
    omm = omm/tau; %will always be called like this from now on
  S_poly_1 = zeros(length(om_plot_rng),nmax,nmax); % w_1 + w_pu
  S_poly_2 = zeros(length(om_plot_rng),nmax,nmax); % w_1 - w_pu
  
  omu = om_u_range(1);
  
  for n1 = 1:nmax
     for n2 = 1:nmax
         prefct = (hpol(:,n1)*(hpol(:,n2)')).*weights;
  %find coefficients of polynomial by manipulating orthogonality
  temp1 = om_plot_rng*0; temp2 = temp1;
      for lp2 = 1:quad_order
          for lp3 = 1:quad_order
    
temp1 = temp1 + prefct(lp2,lp3)*temp_fn(om_plot_rng,omm(lp2)-omu,omm(lp3)+omu);
temp2 = temp2 + prefct(lp2,lp3)*temp_fn(om_plot_rng,omm(lp2)+omu,omm(lp3)-omu);              
                                  
          end
      end
S_poly_1(:,n1,n2)  = temp1;      S_poly_2(:,n1,n2)  = temp2;    
     end
  end