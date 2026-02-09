function [df, gdavg, gdmax] = get_fractal_dimension_graph_distance(NT,mitolen,PLOT)
    if(~exist("PLOT","var"))
        PLOT = false;
    end
    % find the fractal dimension of a network
    % also get the mean graph distance (gdavg)
    % also get the max graph distance (gdmax)
    nedges = NT.nedge;
    
    if(isscalar(mitolen))
        edgelens = mitolen.*ones(nedges,1);
    else
        edgelens = mitolen;
        mitolen = mean(edgelens);
    end
    G = graph(NT.edgenodes(:,1),NT.edgenodes(:,2),edgelens);
    
    graphDistances = distances(G);
    graphDistances = triu(graphDistances);
    graphDistances = graphDistances(graphDistances>0 & ~isinf(graphDistances));

    % use the mean distance between edges to set the cutoff distance
    mask = ~isinf(graphDistances);
    gdavg = mean(graphDistances(mask));
    gdmax = max(graphDistances(mask));
    degs = degree(G);
    seglenavg = 2*sum(G.Edges.Weight) / sum(degs(degs~=2));

    if(gdavg <= seglenavg)
        df = 1;
        return
    end
    
    rvals = log10(mitolen/2:mitolen:max(graphDistances(~isinf(graphDistances))));    
    LofR = cumsum(histcounts(graphDistances(~isinf(graphDistances)), [-inf; 10.^rvals(:); inf]));
    LofR = LofR(1:end-1);   % final result: counts(i) = # distances <= radii(i)

    idx1 = find(10.^rvals>seglenavg,1);
    if(isempty(idx1))
        df=1;
        return
    end
    idx2 = find(10.^rvals>gdavg,1);
    if(isempty(idx2))
        idx2 = length(rvals);
    end
    if(idx1 >= idx2)
        df = 1;
        return
    end
    f2 = fit(rvals(idx1:idx2)',log10(LofR(idx1:idx2))','a*x+b','start',[1.5 0]);
    df = max([f2.a 1]);

    if(PLOT)
        figure
        scatter(10.^rvals',LofR','filled','HandleVisibility','off');
        hold on
        plot(10.^[rvals(idx1) rvals(idx2)], (10.^f2.b)*(10.^([rvals(idx1) rvals(idx2)])).^f2.a,'LineWidth',2,'DisplayName','slope = $'+string(f2.a)+'$');
        
        set(gca,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',20)
        xlabel("$r_g$ (graph distance)")
        ylabel("$n$ (number of nodes)")
        xlim(10.^[min(rvals) max(rvals)])
        ylim([min(LofR) max(LofR)])
    
        set(gca,'XScale','log')
        set(gca,'YScale','log')

        leg = legend('Location','northwest');
        leg.Interpreter = "latex";
        leg.FontSize = 15;

        xticks([1 10])
        yticks([1 10 100])
        
        set(gcf,'color','w');
        grid on
    end
end