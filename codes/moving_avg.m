% Feb. 1, 2018, Y.C.
% apply moving window to the data
% Input:
% X: 2D matrix with dimension of NX*NY
% N: length of the window
% type: 'gauss' or 'constant'
% dim: 1 (column) or 2 (row)
function [Y,w]=moving_avg(X,N,type,dim)
if nargin==3
    dim=1;
elseif nargin==2
    dim=1;
    type='constant';
end
if dim==2
   X=X'; 
end
switch type
    case 'gauss'
        w=gausswin(N,2.5);
    case 'constant'
        w=ones(1,N)/N;
end
if mod(N+1,2)~=0
    warning('N needs to be an odd number!');
    return;
end
[nx,ny]=size(X);
Y=zeros(nx,ny);
% padding the element to the edge
X=padarray(X,floor(N/2),'symmetric','both');
for i=1:ny
    % convolve with the gaussian window
    d=X(:,i);
    dnew=conv(d,w);
    d=dnew(N:end-N+1);
    Y(:,i)=d;
end
if dim==2
   Y=Y'; 
end
