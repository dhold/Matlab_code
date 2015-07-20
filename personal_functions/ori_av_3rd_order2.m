function [av_set1,av_set2 ,av_set3] = ori_av_3rd_order2(mu,m,Q,k,lg,save_sep)

% System parameters:
% mu are the 4 exciton dipole moments, 4 by 3 by n
% m are the 4 exciton intrinsic magnetic moments, can also include mR in m 
% as 1i sum_{ma} psi_{xi,ma} mu_{ma} cross R_{m} as this works out the same
% Q is the rank 2 tensor component  for exciton the quadrupole moment and
% (also the position dependent terms mu_{ma}^(v1) r_m^(v2) if required)
% Note Q should be 3 by 3 by n by 4 array
% Beam parameters:
% k is the set of wavevectors (recommended to pass unit and deal with the
% amplitude later)
% k is the set of wavevectors of the three beams followed by the 
% polarisations of each beam, and the 4th Heterodyne beam
% lg is logical length 3 determining which configurations should be included
% in the averages (rephasing nonrephasing and coherence).

if length(lg) ~=3 
    error('lg should be vector of length 3 and logical')
end
if ~exist('save_sep','var')
   save_sep = false; %combine all contributions with different 
   %proportionality to k so this can be added in later
end

%permute to put final dimension first, this just works out easier for array
%operations but for other purposes the order is more sensible
v = permute(v,[3,2,1]); k = permute(k,[3,2,1]);
n_beam = size(v,1);  
%calculate lowest order terms

tmp = [-1,-1,-1];
prefct = -prod(tmp(lg),2); %prefactor determining whether abs or emission

M4 = [4,-1,-1;-1,4,-1;-1,-1,4]/30;
M5 = [3,-1,-1,1,1,0;-1,3,-1,-1,0,1;-1,-1,3,0,-1,-1;...
      1,-1,0,3,-1,1; 1,0,-1,-1,3,-1;0,1,-1,1,-1,3 ]/30;
  
%pass a cell mu = {mu_ge, me_ef} to get this to calculate 
      %all possible relevant dipole averages for each configuration

      %mu_{ge} mu_{e f}  mu_{f e}  ->  mu_{e g}, mu_{g e}
      mu_ge = mu{1};  me_ef = mu{2};
      if lg(3) %coherence pathways
          
              mu_set_coh = zeros(N^2*(N-1)/2,3,4);  
              for lp = 1:N
              tmp = repmat(mu_ge(lp,:),N*(N-1)/2,1);
              mu_set_coh(cnt+1:N*(N-1)/2,:,1) = tmp;
              mu_set_coh(cnt+1:N*(N-1)/2,:,4) = tmp;
              cnt = 0;

                  tmp = repmat(mu_ef,N,1);
                    mu_set_coh(cnt+(1:N*(N-1)/2),:,2) = squeeze(mu_ef(:,lp,:));
                    mu_set_coh(cnt+(1:N*(N-1)/2),:,3) = squeeze(mu_ef(:,lp,:));
                    cnt = cnt + N*(N-1)/2;
              end
          
      end
            %mu_{ge} mu_{e' g} -> mu_{e'' g} mu_{g e''' } ,mu_{e'' f} mu_{f e'''}
      if lg(2) || lg(1)
          mu_set = zeros(
          
          
      end
  

if istilde(1) %check if I actually need to output this
    av_set1 = [];
else
    %permute to put final dimension first
    mu = permute(mu,[3,2,1]);
    
        V4 = [(mu(:,:,4)*mu(:,:,3).' ) .* (mu(:,:,2)*mu(:,:,1).' ),...
              (mu(:,:,4)*mu(:,:,2).' ) .* (mu(:,:,3)*mu(:,:,1).' ),...
              (mu(:,:,4)*mu(:,:,1).' ) .* (mu(:,:,3)*mu(:,:,2).' )];  
          
        F4 = [(v(:,:,4)*v(:,:,3).' ) .* (v(:,:,2)*v(:,:,1).' ),...
              (v(:,:,4)*v(:,:,2).' ) .* (v(:,:,3)*v(:,:,1).' ),...
              (v(:,:,4)*v(:,:,1).' ) .* (v(:,:,3)*v(:,:,2).' )];

    av_set1 = prefct.*(F4.' * M4 * V4); 

end

if istilde(2) || isempty(m)
    av_set2 = []; %don't bother to calculate
else
    %permute to put final dimension first
    m = permute(m,[3,2,1]);
    vm = cross(v,k); %magnetic field direction
    
    %calculate the 4 different terms with the mag interaction at different
    %positions
        V44 = [(m(:,:,4)*mu(:,:,3).' ) .* (mu(:,:,2)*mu(:,:,1).' ),...
              (m(:,:,4)*mu(:,:,2).' ) .* (mu(:,:,3)*mu(:,:,1).' ),...
              (m(:,:,4)*mu(:,:,1).' ) .* (mu(:,:,3)*mu(:,:,2).' )];  
        V43 = [(mu(:,:,4)*m(:,:,3).' ) .* (mu(:,:,2)*mu(:,:,1).' ),...
              (mu(:,:,4)*mu(:,:,2).' ) .* (m(:,:,3)*mu(:,:,1).' ),...
              (mu(:,:,4)*mu(:,:,1).' ) .* (m(:,:,3)*mu(:,:,2).' )];  
        V42 = [(mu(:,:,4)*mu(:,:,3).' ) .* (m(:,:,2)*mu(:,:,1).' ),...
              (mu(:,:,4)*m(:,:,2).' ) .* (mu(:,:,3)*mu(:,:,1).' ),...
              (mu(:,:,4)*mu(:,:,1).' ) .* (mu(:,:,3)*m(:,:,2).' )];  
        V41 = [(mu(:,:,4)*mu(:,:,3).' ) .* (mu(:,:,2)*m(:,:,1).' ),...
              (mu(:,:,4)*mu(:,:,2).' ) .* (mu(:,:,3)*m(:,:,1).' ),...
              (mu(:,:,4)*m(:,:,1).' ) .* (mu(:,:,3)*mu(:,:,2).' )];   
               
          
        F44 = [(vm(:,:,4)*v(:,:,3).' ) .* (v(:,:,2)*v(:,:,1).' ),...
              (vm(:,:,4)*v(:,:,2).' ) .* (v(:,:,3)*v(:,:,1).' ),...
              (vm(:,:,4)*v(:,:,1).' ) .* (v(:,:,3)*v(:,:,2).' )];
        F43 = [(v(:,:,4)*vm(:,:,3).' ) .* (v(:,:,2)*v(:,:,1).' ),...
              (v(:,:,4)*v(:,:,2).' ) .* (vm(:,:,3)*v(:,:,1).' ),...
              (v(:,:,4)*v(:,:,1).' ) .* (vm(:,:,3)*v(:,:,2).' )];
        F42 = [(v(:,:,4)*v(:,:,3).' ) .* (vm(:,:,2)*v(:,:,1).' ),...
              (v(:,:,4)*vm(:,:,2).' ) .* (v(:,:,3)*v(:,:,1).' ),...
              (v(:,:,4)*v(:,:,1).' ) .* (v(:,:,3)*vm(:,:,2).' )];
        F41 = [(v(:,:,4)*v(:,:,3).' ) .* (v(:,:,2)*vm(:,:,1).' ),...
              (v(:,:,4)*v(:,:,2).' ) .* (v(:,:,3)*vm(:,:,1).' ),...
              (v(:,:,4)*vm(:,:,1).' ) .* (v(:,:,3)*v(:,:,2).' )];

          if save_sep %this allows us to reintroduce the proportionality to |k|
              %well calculating every frequency of interaction we can
              %introduce back this magnitude
         av_set2(:,:,1)= -1i*prefct.*F44.' * M4 * V44;
         av_set2(:,:,2)= 1i*prefct.* s(3,:).*(F43.' * M4 * V43);
         av_set2(:,:,3)= 1i*prefct.* s(2,:).*(F42.' * M4 * V42);
         av_set2(:,:,4)= 1i*prefct.* s(1,:).*(F41.' * M4 * V41);
          else
%otherwise just weight the mag dipole moment with the exciton res freq           
    av_set2= 1i*(-F44.' * M4 * V44 + s(3,:).*(F43.' * M4 * V43)+...
                 s(2,:).*(F42.' * M4 * V42) + s(1,:).*(F41.' * M4 * V41)) ; 
    av_set2= prefct.*av_set2;
          end
end


if istilde(3) || isempty(Q) %rank 5 terms, can include the oligomer 
    %structure terms in either and it should give the same answer.
    av_set3 = []; %don't bother to calculate
else    
    warning('this is not yet finished')
    av_set3 =zeros(size(Q,4),n_beam);
    %calculate the 4 different terms with the mag interaction at different
    %positions
    unit_vecs = eye(3);
    mu = permute(mu,[2,4,1,3]) ; %reshape mu to a more useful shape to match Q
    
    for lp = 1:3
    
         tmp = v; tmp(:,:,5) = repmat(unit_vecs(:,lp),[size(tmp,1),1]);
         F_tmp =F5(tmp);
    
        V54 = [(m(:,:,4)*mu(:,:,3).' ) .* (mu(:,:,2)*mu(:,:,1).' ),...
              (m(:,:,4)*mu(:,:,2).' ) .* (mu(:,:,3)*mu(:,:,1).' ),...
              (m(:,:,4)*mu(:,:,1).' ) .* (mu(:,:,3)*mu(:,:,2).' )];  
          
          
          
          
        V53 = [(mu(:,:,4)*m(:,:,3).' ) .* (mu(:,:,2)*mu(:,:,1).' ),...
              (mu(:,:,4)*mu(:,:,2).' ) .* (m(:,:,3)*mu(:,:,1).' ),...
              (mu(:,:,4)*mu(:,:,1).' ) .* (m(:,:,3)*mu(:,:,2).' )];  
        V52 = [(mu(:,:,4)*mu(:,:,3).' ) .* (m(:,:,2)*mu(:,:,1).' ),...
              (mu(:,:,4)*m(:,:,2).' ) .* (mu(:,:,3)*mu(:,:,1).' ),...
              (mu(:,:,4)*mu(:,:,1).' ) .* (mu(:,:,3)*m(:,:,2).' )];  
        V51 = [(mu(:,:,4)*mu(:,:,3).' ) .* (mu(:,:,2)*m(:,:,1).' ),...
              (mu(:,:,4)*mu(:,:,2).' ) .* (mu(:,:,3)*m(:,:,1).' ),...
              (mu(:,:,4)*m(:,:,1).' ) .* (mu(:,:,3)*mu(:,:,2).' )];   
          
        av_set3 =  av_set3 + F_tmp*M5*V54;      
    end     
end
end

function  out = F5(vv)  %calculates the 5th rank beam dep params

     out = [vv(:,:,4)*cross(vv(:,:,3),vv(:,:,2)).' .* (vv(:,:,1)*vv(:,:,5).'),...
      vv(:,:,4)*cross(vv(:,:,3),vv(:,:,1)).' .* (vv(:,:,2)*vv(:,:,5).'),...
      vv(:,:,4)*cross(vv(:,:,3),vv(:,:,5)).' .* (vv(:,:,2)*vv(:,:,1).'),...
      vv(:,:,4)*cross(vv(:,:,2),vv(:,:,1)).' .* (vv(:,:,3)*vv(:,:,5).'),...
      vv(:,:,4)*cross(vv(:,:,2),vv(:,:,1)).' .* (vv(:,:,3)*vv(:,:,1).'),...
      vv(:,:,4)*cross(vv(:,:,1),vv(:,:,1)).' .* (vv(:,:,3)*vv(:,:,2).')];
    out = out/30;
end

function  out = rank2dot(d,tt)  
    %same as d dot t1 cross t2 if tt^{a,b} = t1^{a} t2^{b}

     out = d(3,1,:).*(tt(1,2,:)-tt(2,1,:))-...
           d(2,1,:).*(tt(1,3,:)-tt(3,1,:)) +...
           d(1,1,:).*(tt(2,3,:)-tt(3,2,:));
end

function  [out1,out2,out3,out4] = V5(dd,t)  %calculates the 5th rank beam dep params

    out4 = [mtimesx(mtimesx(dd(:,1,:,1),'T',t(:,:,:,4)), cross(dd(:,1,:,3),dd(:,1,:,2))),...
      mtimesx(mtimesx(dd(:,1,:,2),'T',t(:,:,:,4)), cross(dd(:,1,:,3),dd(:,1,:,1))),...
      rank2dot(dd(:,1,:,3),t(:,:,:,4)) .* mtimesx(dd(:,1,:,1),'T',dd(:,1,:,2)) ,...
       mtimesx(mtimesx(dd(:,1,:,3),'T',t(:,:,:,4)), cross(dd(:,:,:,2),dd(:,1,:,1))),...
      rank2dot(dd(:,1,:,2),t(:,:,:,4)) .* mtimesx(dd(:,1,:,3),'T',dd(:,1,:,1)) ,...
      rank2dot(dd(:,1,:,1),t(:,:,:,4)) .* mtimesx(dd(:,1,:,3),'T',dd(:,1,:,2)) ];

    out3 = [mtimesx(mtimesx(dd(:,1,:,1),'T',t(:,:,:,3)), cross(dd(:,1,:,2),dd(:,1,:,4))),...
      mtimesx(mtimesx(dd(:,1,:,2),'T',t(:,:,:,3)), cross(dd(:,1,:,1),dd(:,1,:,4))),...
      -rank2dot(dd(:,1,:,4),t(:,:,:,3)) .* mtimesx(dd(:,1,:,2),'T',dd(:,1,:,1)) ,...
       (t(1,1,:,3)+t(2,2,:,3)+t(3,3,:,3)).*mtimesx(dd(:,1,:,4),'T',cross(dd(:,1,:,2),dd(:,1,:,1))),...
    
end