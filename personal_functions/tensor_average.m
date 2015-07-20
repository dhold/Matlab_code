function out = tensor_average(vctors,a,b)
% assumes first N-2 components are in the x direction
% a, b, can be 1,2 or 3 respectively, 1 is x, 2,y,3,z etc etc
%input vctors as a 3XN matrix
if nargin == 1
   a=1; b=1; %assume x x x x,..,x average 
end
if size(vctors,2)==1
    out=0; %always averages to zero
elseif  size(vctors,2)==2
    if a==b
    out = dot(vctors(:,1),vctors(:,2))/3;
    else
       out = 0; 
    end
elseif  size(vctors,2)==3
    if a==b || a==1 || b==1
        out=0;
    else
        if a==2 && b==3
            tmp = cross(vctors(:,2),vctors(:,3));
        elseif  b==2 && a==3
            tmp = cross(vctors(:,3),vctors(:,2));
        end
        out = dot(vctors(:,1),tmp)/6; 
    end
elseif  size(vctors,2)==4
    if a~=b
        out=0;
    elseif a==b && a==1
     out = (dot(vctors(:,1),vctors(:,2))*dot(vctors(:,3),vctors(:,4))...
           +dot(vctors(:,1),vctors(:,3))*dot(vctors(:,2),vctors(:,4))...
           +dot(vctors(:,1),vctors(:,4))*dot(vctors(:,2),vctors(:,3)))/15;     
    elseif a==b && a~=1  
     out = (4*dot(vctors(:,1),vctors(:,2))*dot(vctors(:,3),vctors(:,4))...
           -dot(vctors(:,1),vctors(:,3))*dot(vctors(:,2),vctors(:,4))...
           -dot(vctors(:,1),vctors(:,4))*dot(vctors(:,2),vctors(:,3)))/30;                    
    end
elseif  size(vctors,2)==5
    if a==2 && b==3
        tmp = cross(vctors(:,4),vctors(:,5));
    elseif b==2 && a==3
        tmp = cross(vctors(:,5),vctors(:,4));
    else
        tmp=[0;0;0];
    end
     out = (dot(vctors(:,1),vctors(:,2))*dot(vctors(:,3),tmp)...
           +dot(vctors(:,1),vctors(:,3))*dot(vctors(:,2),tmp)...
           +dot(vctors(:,1),tmp)*dot(vctors(:,2),vctors(:,3)))/15;       
else
    out ='I havent writen averages to this rank yet'
end
end               