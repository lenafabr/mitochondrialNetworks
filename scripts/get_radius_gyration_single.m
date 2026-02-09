function [rg, x1, x2] = get_radius_gyration_single(NT,PERCLUSTER,rs,mitolen,use_edge_lengths)
% rg = radius of gyration per cluster (return in terms of microns)
% x1 = typical number of tips per cluster
% x2 = typical number of deg-2 nodes per cluster

% rs = steric radius. Enter zero if using real rather than simulation networks
    nodepos = NT.nodepos;
    if use_edge_lengths
        edgelens = NT.edgelens;
    else
        edgelens = mitolen*ones(NT.nedge,1);
    end
    if rs > 0
        tipmask = 1:NT.nnode;
        for nc = tipmask(NT.degrees==1)
            n2p = nodepos(NT.nodenodes(nc,1),:);
            ncp = nodepos(nc,:);
            vec = (ncp-n2p)/sqrt(sum((ncp-n2p).^2));
            nodepos(nc,:) = ncp + rs*vec;
        end
    end

    if ~PERCLUSTER
        com = [0 0 0];
        comsqr = 0;
        edgelocs = zeros(NT.nedge,3);
        tot = sum(edgelens);
        for ec = 1:NT.nedge
            loc = 0.5*(nodepos(NT.edgenodes(ec,1),:) + nodepos(NT.edgenodes(ec,2),:));
            com = com + loc*edgelens(ec)/tot;
            comsqr = comsqr + sum(loc.^2)*edgelens(ec)/tot;
            edgelocs(ec,:) = loc;
        end
        currrgsqr = comsqr+sum(com.^2);
        for ec = 1:NT.nedge
            currrgsqr = currrgsqr - 2*sum(edgelocs(ec,:).*com)*edgelens(ec)/tot;
        end
        rg = sqrt(currrgsqr);
        degs = NT.degrees;
        x1 = sum(degs==1);
        nxg2 = sum(degs(degs>2));
        x2 = sum(edgelens)/mitolen - 0.5*(x1+nxg2);
    else
        if(numel(NT.edgeedges) == 0)
            NT.setupNetwork;
        end
        
        bins = conncomp(NT.makeGraph);
        if ~use_edge_lengths
            allweights = mitolen * NT.degrees/2;
        else
            allweights = zeros(NT.nnode,1);
            for nc = 1:NT.nnode
                deg = NT.degrees(nc);
                allweights(nc) = sum(edgelens(NT.nodeedges(nc,1:deg)))/2;
            end
        end
        weightedavg = 0;
        weightedx1 = 0;
        weightedx2 = 0;
        for bc = 1:max(bins)
            weights = allweights(bins==bc);
            pos = nodepos(bins==bc,:);
            sz = sum(weights(:,1),'all');
            com = sum(pos.*weights,1)/sz;
            comsqr = sum(pos.^2.*weights,'all')/sz;
            currrgsqr = comsqr+sum(com.^2);
            for nc = 1:sum(bins==bc)
                currrgsqr = currrgsqr - 2*sum(pos(nc,:).*weights(nc,:).*com)/sz;
            end
            weightedavg = weightedavg + sqrt(currrgsqr)*sz;
            degs = NT.degrees(bins==bc);
            nx1 = sum(degs==1);
            if use_edge_lengths
                nxg2 = sum(degs(degs>2));
                nx2 = sum(edgelens)/mitolen - 0.5*(nx1+nxg2);
            else
                nx2 = sum(degs==2);
            end
            weightedx1 = weightedx1 + nx1*sz;
            weightedx2 = weightedx2 + nx2*sz;
        end
        tot = sum(allweights);
        rg = weightedavg/tot;
        x1 = weightedx1/tot;
        x2 = weightedx2/tot;
    end
end