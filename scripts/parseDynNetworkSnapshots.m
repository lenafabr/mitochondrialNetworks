function [networks,allinfo] = parseDynNetworkSnapshots(filename,nsnapshots)
% parse a file containing dynamic network snapshots into an array of
% network objects

%%
data = dlmread(filename);%,'',[0 0 64000 499]);

if(exist("nsnapshots","var"))
    networks = NetworkObj.empty(nsnapshots,1,0);
    allinfo = cell(nsnapshots,1);
else
    nsnapshots = inf;
end

nc = 0; % current network
startline = 1; % initial line for current network

% keep going until read whole file or reached highest desired snapshot
while (startline< size(data,1) && nc < nsnapshots) 
    nc = nc+1;
    networks(nc,1,1) = [NetworkObj()];
    
    infoline = data(startline,:);
    dim = infoline(1);
    nnode = infoline(2);
    nedge = infoline(3);
    nnv = infoline(4); % number extra node values (beyond pos)
    nev = infoline(5); % number extra edge values (beyond n1,n2)
    ninfo = infoline(6); % number extra info values
    info = infoline(7:6+ninfo);
   
    % check for integer values
    if (mod(nnv,1)~=0 || mod(nev,1)~=0 || mod(dim,1)~=0)
        error('something is wrong with infoline')
    end

    allinfo{nc} = info;
    
    NT = networks(nc);
    NT.dim = dim;
    NT.nnode = nnode;
    NT.nedge = nedge;
    lc = startline;
    % node coordinates
    NT.nodepos = data(lc+1:lc+dim,1:nnode)';
    
    % extra node values;
    lc = lc+dim;
    NT.nodevals = data(lc+1:lc+nnv,1:nnode)';
    
    % edge connectivity
    lc = lc+nnv;
    NT.edgenodes = data(lc+1:lc+2,1:nedge)';
    
    lc = lc+2;
    % edge lengths
    if (nev>=1) % assume first value is edge length
        NT.edgelens = data(lc+1,1:nedge)';
    end
    % other edge values
    NT.edgevals = data(lc+2:lc+nev,1:nedge)';
    
    % setup network
    %NT.setupNetwork()
    
    startline = lc+nev+1;    
end