function [meancs, largestcs] = get_mean_cluster_size_single(NT,use_edge_lengths)
    % return in terms of edge lengths or edges, depending on
    % use_edge_lengths
    if(numel(NT.edgeedges) == 0)
        NT.setupNetwork;
    end
    if ~use_edge_lengths
        edgelens = ones(NT.nedge,1);
    else
        edgelens = NT.edgelens;
    end
    edgelist = squeeze(NT.edgeedges(:,2,:));
    A = zeros(NT.nedge,NT.nedge);
    for ec = 1:NT.nedge
        neighbs = edgelist(ec,:);
        A(ec,edgelist(ec,neighbs > 0)) = 1;
    end
    NTgraph = graph(A);
    bins = conncomp(NTgraph);
    allsz = zeros(max(bins),1);
    num = 0;
    for bc = 1:max(bins)
        sz = sum(edgelens(bins==bc));
        allsz(bc) = sz;
        num = num+sz^2;
    end
    meancs = num/sum(edgelens);
    largestcs = max(allsz);
end