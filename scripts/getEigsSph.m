function [eiglist,func] = getEigsSph(a,b,nmax,typebounds,dx,cutoff)
% get eigenvalues for cylindrically symmetric diffusion problem
% in annular region between radius a and b
% nmax is the number of eigenvalues to calculate
% typebounds = for inner, outer boundary, 0 for absorbing; 1 for reflecting
% dx is the step in initial x values for finding zeros
% default is 1

if (~exist('dx','var'))
    dx = 0.1/b;
end
if (~exist('cutoff','var'))
    cutoff = 2*pi*5/(b-a);
end

if (typebounds(1)==0 & typebounds(2)==1)
    % inner absorbing, outer reflecting
    error('not set up')
elseif(typebounds(1)==1 & typebounds(2)==0)
    % inner reflecting, outer absorbing
    func = @(x) a*x.*cos(x*(b-a)) + sin(x*(b-a));
elseif (typebounds(1)==0 & typebounds(2)==0)
    % both bounds absorbing
     error('not set up')      
elseif (typebounds(1)==1 & typebounds(2) == 1)
    % both bounds reflecting
    func = @(x) sin(x*(b-a)).* (x.^2 + 1/a/b)- cos(x*(b-a)).*x.*(1/a-1/b);
else
    error('invalid typebounds. Must be 0 for abs, 1 for ref')
end

eiglist = [];
ct = 0;

x0 = [dx/2,3*dx/2];

while length(eiglist)<nmax
    if (func(x0(1))*func(x0(2)) > 0)        
        % same sign, just use this starting point
        val =  fzero(func,mean(x0));
    else % different sign, search in this interval
        val = fzero(func, x0);
    end
    [x0,val];
    if (isempty(eiglist))
        ct=ct+1;
        eiglist(ct) = val;
        if (val>cutoff)
            x0 = val+pi/(b-a);
            x0 = [x0-dx,x0+dx];
        end
    elseif (val > eiglist(end)+1e-10/b)
        ct =ct+1;
        eiglist(ct) = val;
        if (val>cutoff)
            x0 = val+pi/(b-a);
            x0 = [x0-dx,x0+dx]; % try an interval
        end
    else
        x0 = x0 + dx;
    end
   
end


if (typebounds(1)==1 & typebounds(2)==1)
    % include 0 eigenvalue
    eiglist = [0 eiglist(eiglist>1e-7)];
else
    % if there are absorbing bounds, get rid of 0 eigenvalue
    eiglist(eiglist<1e-10) = [];
end


%plot(xlist,flist,eiglist,zeros(size(eiglist)),'*')
end