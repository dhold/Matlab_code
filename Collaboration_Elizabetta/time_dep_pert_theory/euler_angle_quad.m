%euler angle quadrature rule (by me should look this up really)
maxorder = 4;
intlist = intpartgen2(maxorder,3);
for k =1:maxorder
    temp = intlist{k};
   if size(temp,2)~=3
      temp = [temp,zeros(size(temp,1),3-size(temp,2))]; 
   end
   temp2 = zeros(0,3); %empty to start
    for j = 1:size(temp,1)
    temp2 = [temp2;perms(temp(j,:))];
    end
    permlist{k} = temp2;
end

phirng{1} = [0,pi];
thetarng{1} = pi/2;
xirng{1} = [0,pi];

eangles{1} = [0,pi/2,0;pi,pi/2,0;0,pi/2,pi;pi,pi/2,pi];
weights{1} = [1/4,1/4,1/4,1/4];

%interpolate by taking half way points

for k = 2:maxorder
    
    tmp = eangles{k-1};
    
	phirng{k} = 2*pi*((2^(-k)):2^(-k+1):(1-2^(-k))) ; %#ok<*SAGROW>
    thetarng{k} = pi*((2^(-k)):2^(-k+1):(1-2^(-k))) ;
    xirng{k} = 2*pi*((2^(-k+1)):2^(-k+2):(2-2^(-k+1)));
    sizeoftier = length(phirng{k})*length(thetarng{k-1})*length(xirng{k-1})
    
    for j=1:size(temp,1)
    eangles{k} = 
    end
end



