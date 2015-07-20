function F=ori_F_calc(v)
%this generates the elements for F^(4) and F^(5)

if size(v,1) ~= 4 &&  size(v,1) ~= 5
    error('input v must be 4 or 5 by 3 by N')
end

if size(v,1) == 4
F = zeros(3,size(v,3));
else
F = zeros(6,size(v,3));    
end
if strcmp(class(v),'sym')
    F=sym(F);
end
    
for lp =1:size(v,3)
    
    vv = v(:,:,lp);
    
if size(v,1) == 4
    
    F(:,lp) = [vv(4,:)*vv(3,:).' .* vv(2,:)*vv(1,:).',...
               vv(4,:)*vv(2,:).' .* vv(3,:)*vv(1,:).',...
               vv(4,:)*vv(1,:).' .* vv(3,:)*vv(2,:).'];
    
else

    F(:,lp) = [vv(4,:)*cross(vv(3,:),vv(2,:)).' .* (vv(1,:)*vv(5,:).'),...
               vv(4,:)*cross(vv(3,:),vv(1,:)).' .* (vv(2,:)*vv(5,:).'),...
               vv(4,:)*cross(vv(3,:),vv(5,:)).' .* (vv(2,:)*vv(1,:).'),...
               vv(4,:)*cross(vv(2,:),vv(1,:)).' .* (vv(3,:)*vv(5,:).'),...
               vv(4,:)*cross(vv(2,:),vv(5,:)).' .* (vv(3,:)*vv(1,:).'),...
               vv(4,:)*cross(vv(1,:),vv(5,:)).' .* (vv(3,:)*vv(2,:).')];
    
end

end