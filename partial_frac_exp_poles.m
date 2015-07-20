N=60;
matr = diag(2*(1:N-1).*(2*(1:N-1)+1),1);
matr(N,:) = -2*N*(2*N+1);

%     matr2 = matr*0;
% for m = 1:N
%     for n=1:N
%         if m+1==n
%             matr2(m,n) = 2*m*(2*m+1);
%         elseif m==N
%             matr2(m,n) = -2*N*(2*N+1);
%         end
%     end
% end
% 
% all(all(matr2==matr))

ev = eig(matr);
xj = 2*[sqrt(ev);-sqrt(ev)];
%xj =2*sqrt(ev);
if 1==0
plot(real(xj),imag(xj),'Marker','+','LineStyle','none');
hold on
plot(0*xj,2*pi*([-(1:length(xj)/2),1:length(xj)/2]),'Marker','X','LineStyle','none');
% Create xlabel
xlabel('Real');

% Create ylabel
ylabel('Imag');
end

