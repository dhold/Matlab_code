function [av_set1,av_set2 ,av_set3,prefct] = ori_av_3rd_order(mu,m,Q,k,v,s,save_sep)

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
% v is the set of polarisations of each beam
% s is length in dim2 which is plus or minus 1, with plus denoting a
% conjugate, if it has any length in dim1 this gets summed over

if size(s,2) ~=3 || any(abs(s)~=1)
    error('s matrix of length 3 in dim 2 of plus of minus values')
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

prefct = -prod(s,2); %prefactor determining whether abs or emission

M4 = [4,-1,-1;-1,4,-1;-1,-1,4]/30;
M5 = [3,-1,-1,1,1,0;-1,3,-1,-1,0,1;-1,-1,3,0,-1,-1;...
      1,-1,0,3,-1,1; 1,0,-1,-1,3,-1;0,1,-1,1,-1,3 ]/30;

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

   % in hind sight the +/- 1 factor is pretty trivial to add later
   %       av_set1 = zeros(n_beam,size(mu,1),length(prefct));
   %       for j = 1:length(prefct)
   % av_set1(:,:,j) = prefct(j).*(F4.' * M4 * V4); 
   %       end
    av_set1 = (F4.' * M4 * V4);

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

          for j = 1:length(prefct) %loop over the different configurations considered
          if save_sep %this allows us to reintroduce the proportionality to |k|
              %well calculating every frequency of interaction we can
              %introduce back this magnitude
         av_set2(:,:,j,1)= -1i*prefct(j).*F44.' * M4 * V44;
         av_set2(:,:,j,2)= 1i*prefct(j).* s(3,:).*(F43.' * M4 * V43);
         av_set2(:,:,j,3)= 1i*prefct(j).* s(2,:).*(F42.' * M4 * V42);
         av_set2(:,:,j,4)= 1i*prefct(j).* s(1,:).*(F41.' * M4 * V41);
          else
%otherwise just weight the mag dipole moment with the exciton res freq           
    av_set2(:,:,j)= 1i*(-F44.' * M4 * V44 + s(j,3).*(F43.' * M4 * V43)+...
                 s(j,2).*(F42.' * M4 * V42) + s(j,1).*(F41.' * M4 * V41)) ; 
    av_set2(:,:,j)= prefct(j).*av_set2(:,:,j);
          end
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

warning('unfinished function')

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