function mitodiff = get_mean_cluster_diffusivity(networks,lagtime,simstepspersnap,delt)
    % networks: a 1D array of NetworkObj
    % lagtime: the time interval over which to measure instantaneous MSD
    % simstepspersnap: how frequently are network snapshots saved, in units
    % of simulation steps
    % delt: how long is one simulation step, in simulation time units
    nedges = networks(end).nedge;
    nsnaps = numel(networks);
    ndt = round(lagtime/(simstepspersnap*delt));
    epos = zeros(nsnaps,nedges,3);
    esizes = zeros(nsnaps,nedges);
        
    for sc = 1:nsnaps
        NT = networks(sc);
        NT.setupNetwork;
    
        edgelist = squeeze(NT.edgeedges(:,2,:));
        A = zeros(NT.nedge,NT.nedge);
        for ec = 1:NT.nedge
            neighbs = edgelist(ec,:);
            A(ec,edgelist(ec,neighbs > 0)) = 1;
        end
        G = graph(A);
        edgeclusters = conncomp(G);
        sizes = histcounts(edgeclusters, 1:max(edgeclusters)+1);
        edgeclustsizes = sizes(edgeclusters);
    
        epos_curr = 0.5 * (NT.nodepos(NT.edgenodes(:,1),:) + NT.nodepos(NT.edgenodes(:,2),:));
        epos(sc,:,:) = epos_curr;
        esizes(sc,:) = edgeclustsizes;
    end
    
    delta = epos(1+ndt:end,:,:) - epos(1:end-ndt,:,:);
    edgeD = sum(delta.*delta,3)/(6*ndt*simstepspersnap*delt);
    sizes_save = esizes(1:end-ndt,:);
    
    % calculate from edgeD and sizes
    [DofN, Nvals, Ncounts]=groupsummary(edgeD(:),sizes_save(:),"mean");
    mitodiff = sum(Ncounts./Nvals.*DofN)/sum(Ncounts./Nvals);
end