% this is a matlab file to be used in analysis of example transport
% simulations. The output will be a plot similar to one of the curves from
% figure 4a in our manuscript.

% NOTE: you must also go to github and add:
% github.com/lenafabr/networktools/NetworkObj.m
% to your path for this script to work correctly

%% set the relevant parameters
signature = "constconc"; % the name of your param files/out files
mitochondrialNetworksPath = '../'; % replace this with the path to mitochondrialNetworks on your computer
simloc = mitochondrialNetworksPath+"param_files/transport_example/";
snaploc = simloc+signature+"_";
ffloc = simloc+signature+"_";
imax = 1; 
jmax = 6; % simulation index - corresponds to decay rate for this example
kmax = 1;

savesnapstart = 2000000;
bdsteps = [18810000 4580000 3160000 3020000 3000000 3000000];
simstepspersnap = [16810 2580 1160 1020 1000 1000];
nsnaps = 1001;
beginavgsnap = round((nsnaps-1)/2);
delt=1e-4;
nspecies = 9;
nedges = 250;
mitolen = 0.5;

ku1 = 300;
ku2 = 100;
kfiss = 1;
partdiff = 4800;
kdecay = logspace(-2.5,2.5,jmax);

rs = 0.15;
rc = 0.2;

%% measure structural parameters and material content
clustersize = zeros(kmax*imax,1); % mean cluster size
networkdim = zeros(kmax*imax,1); % network fractal dimension
graphdiststats = zeros(kmax*imax,2); % graph distance statistics
clustertipstats = zeros(kmax*imax,2); % number of tips, deg.2 nodes per cluster
edgeconnections = zeros(kmax*imax,1); % mean number of neighbors for each unit
radg = zeros(kmax*imax,1); % cluster radius of gyration
meanConcRaw = zeros(kmax*imax*jmax,nspecies,nsnaps); % mean material per unit

parfor idx = 1:imax*kmax*jmax
    [k, i, j] = ind2sub([kmax imax jmax],idx);
    disp("running idx = "+string(idx)+" of "+string(kmax*imax))
    networks = parseDynNetworkSnapshots(snaploc+string(k-1)+"_"+string(i-1)+"_"+string(j-1)+".snap.out");
    networks = networks(2:end); %ignore the first snapshot
    edgevals = {networks.edgevals};
    edgevals = cell2mat(reshape(edgevals,1,1,[]));

    meanConcRaw(idx,:,:) = squeeze(mean(edgevals,1));

    weight = 1/(nsnaps-beginavgsnap+1);
    cs=0; econn=0; rgyr=0; x1x2=[0 0]; df=0; gdmax=0; gdavg=0; cvhr=0;
    for sc = beginavgsnap:nsnaps
        NT = networks(sc);
        [df1, gda1, gdm1] = get_fractal_dimension_graph_distance(NT,mitolen);
        df = df + df1;
        gdmax = gdmax + gdm1;
        gdavg = gdavg + gda1;
        cs = cs + get_mean_cluster_size_single(NT,false);
        econn = econn + 2*(mean(NT.degrees(NT.degrees>1))-1);
        [rg, x1, x2] = get_radius_gyration_single(NT,true,rs,mitolen,false);
        rgyr = rgyr+rg;
        x1x2 = x1x2 + [x1 x2];
    end
    clustersize(idx) = cs*weight;
    edgeconnections(idx) = econn*weight;
    radg(idx) = rgyr*weight;
    clustertipstats(idx,:) = x1x2*weight;
    networkdim(idx) = df*weight;
    graphdiststats(idx,:) = [gdavg gdmax]*weight;
end

meanclustersize = mean(reshape(clustersize, [kmax imax jmax]),[1 3]);
meandimension = mean(reshape(networkdim, [kmax imax jmax]),[1 3]);
meanedgeconn = mean(reshape(edgeconnections,[kmax imax jmax]),[1 3]);
meangraphdist = reshape(mean(reshape(graphdiststats, [kmax imax jmax 2]),[1 3]), [imax 2]);
meanradg = mean(reshape(radg, [kmax imax jmax]),[1 3]);
meanclustercomp = reshape(mean(reshape(clustertipstats, [kmax imax jmax 2]),[1 3]), [imax 2]);

% mean material in the network: average over trials and species (species are just replicates)
meanConc = reshape(mean(reshape(meanConcRaw, [kmax imax jmax nspecies nsnaps]),[1 4]), [imax jmax nsnaps]);
% find steady state network material
meanConcss = mean(meanConc(:,:,beginavgsnap:end),3);

%% check that total material reaches a steady-state in each simulation
figure
hold on
for j = 1:jmax
    data = meanConc(1,j,:)*nedges - 1;
    timeseries = (1:nsnaps)*delt*simstepspersnap(j);
    plot(timeseries/max(timeseries),data(:),'LineWidth',2)
end
set(gca,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',20)
xlabel("$t$ (normalized simulation time)")
ylabel("S (material in the network)")
set(gcf,'color','w')
%set(gca,'xscale','log')
set(gca,'yscale','log')


%% measure mitochondrial mobility
lagtime = max(simstepspersnap*delt);
mitodiff = zeros(imax*jmax*kmax,1);

parfor idx = 1:imax*jmax*kmax
    [k,i,j] = ind2sub([kmax imax jmax], idx);
    disp("running idx = "+string(idx)+" of "+string(kmax*imax*jmax))
    ndt = round(lagtime/(simstepspersnap(j)*delt));
    epos = zeros(nsnaps-beginavgsnap,nedges,3);
    esizes = zeros(nsnaps-beginavgsnap,nedges);
    
    snapfile = snaploc+string(k-1)+"_"+string(i-1)+"_"+string(j-1)+".snap.out";
    networks = parseDynNetworkSnapshots(snapfile);
    
    for sc = beginavgsnap:nsnaps
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
        epos(sc-beginavgsnap+1,:,:) = epos_curr;
        esizes(sc-beginavgsnap+1,:) = edgeclustsizes;
    end
    
    delta = epos(1+ndt:end,:,:) - epos(1:end-ndt,:,:);
    edgeD = sum(delta.*delta,3)/(6*ndt*simstepspersnap(j)*delt);
    sizes_save = esizes(1:end-ndt,:);
    
    % calculate from edgeD and sizes
    [DofN, Nvals, Ncounts]=groupsummary(edgeD(:),sizes_save(:),"mean");
    mitodiff(idx) = sum(Ncounts./Nvals.*DofN)/sum(Ncounts./Nvals)
end

meanmitodiff = mean(reshape(mitodiff, [kmax imax jmax]), [1 3]);

%% calculate the static limit contribution to the analytic model
kdecaymath = logspace(-4,4,101);
lambda = sqrt(partdiff./kdecaymath);
totstuffssstatic = zeros(imax,length(kdecaymath));

parfor i = 1:imax
    df = meandimension(i);
    mgd = meangraphdist(i,1);
    mcl = meanclustersize(i);
    mec = meanedgeconn(i);
    beta = 2^(df+1)*gamma(df+1)^3/((df+1)*gamma(df/2+1/2)^2*gamma(2*df+1));
    R = (2*df+1)/(2*df*beta)*mgd;
    R1 = 1/(df+1)*R;
    R2 = (2*R^df-R1^df)^(1/df);
    nu = 1-df/2;
    a = R/mcl^(1/df);
    x1 = a*sqrt(df)./lambda;
    x2_1 = R1*sqrt(df)./lambda;
    x2_2 = R2*sqrt(df)./lambda;
    BbyA_1 = besseli(nu-1,x2_1)./besselk(nu-1,x2_1);    
    sdf_1 = df./x1.*(-besseli(nu-1,x1)+besselk(nu-1,x1).*BbyA_1)./(besseli(nu,x1)+besselk(nu,x1).*BbyA_1);
    sdf_1(isinf(BbyA_1)) = df./x1(isinf(BbyA_1)).*besselk(nu-1,x1(isinf(BbyA_1)))./besselk(nu,x1(isinf(BbyA_1)));
    BbyA_2 = besseli(nu-1,x2_2)./besselk(nu-1,x2_2);    
    sdf_2 = df./x1.*(-besseli(nu-1,x1)+besselk(nu-1,x1).*BbyA_2)./(besseli(nu,x1)+besselk(nu,x1).*BbyA_2);
    sdf_2(isinf(BbyA_2)) = df./x1(isinf(BbyA_2)).*besselk(nu-1,x1(isinf(BbyA_2)))./besselk(nu,x1(isinf(BbyA_2)));

    sdf = (sdf_1+sdf_2)/2;
    totstuffssstatic(i,:) = sdf;
end


%% calculate the social network contribution to the analytic model
R = 5;
kT = 1; mu = 1;
B1 = 2; alpha1 = B1/mitolen/kT; 
B2 = 6; alpha2 = B2/mitolen/kT;

Nfree = nedges - meanclustersize;
vcnodes = meanclustercomp(:,1)'*(0.5*4*pi/3*(((2*rc)^2 - (2*rs)^2)^(3/2) + (2*rc)^3 - (2*rs)^3)) + ...
    meanclustercomp(:,2)'*(4*pi/3*((2*rc)^2 - (2*rs)^2)^(3/2));
a = 2*meanradg;
b = (a.^3 + 3*vcnodes/(4*pi)).^(1/3);

b(b > R) = R;
a(a > R-2*(rc-rs)) = R-2*(rc-rs);
vc = 4/3*pi*(b.^3-a.^3);

b = b'*ones(size(kdecaymath));
a = a'*ones(size(kdecaymath));
vc = vc'*ones(size(kdecaymath));

Dmito = 2*meanmitodiff; % relative diffusivity of mito pairs
lambda = sqrt(Dmito'*(1./kdecaymath));

expfact = exp(2*(R-b)./lambda);
z1 = -(b-lambda).*(R+lambda);
z2 = (b+lambda).*(R-lambda).*expfact;
top = z1 + z2;
z3 = (R+lambda);
z4 = (R-lambda).*expfact;
zmat = 4*pi*b.*lambda.*top./(z3 + z4);
zmat(isinf(expfact)) = 4*pi*b(isinf(expfact)).*lambda(isinf(expfact)).*(lambda(isinf(expfact))+b(isinf(expfact)));

frac = sum(meanclustercomp,2)';
FUDGE2 = (meanclustercomp(:,1).^2)'./frac;
FUDGE3 = 2*(meanclustercomp(:,1).*meanclustercomp(:,2))'./frac;

ku = FUDGE2.*predict_pesc_deg2_public3(ku1,B1,alpha1,mitolen,rs,rc,kT,mu)/...
    (0.5*4*pi/3*(((2*rc)^2 - (2*rs)^2)^(3/2) + (2*rc)^3 - (2*rs)^3)) + ...
    FUDGE3.*predict_pesc_deg3_public2(ku2,B1,B2,alpha2,mitolen,rs,rc,kT,mu)/...
    (4*pi/3*((2*rc)^2 - (2*rs)^2)^(3/2));

% solve for static soln at kd = kf
skf = zeros(1,imax);
for i = 1:imax
    [~,kfidx] = min(abs(log10(kdecaymath)-log10(kfiss)));
    skf(i) = totstuffssstatic(i,kfidx);
end

kumat = ku'*ones(1,length(kdecaymath));
kdmat = ones(imax,1)*kdecaymath;
skfmat = ((skf'+1)./(meanclustersize)')*ones(1,length(kdecaymath));
srcconcmat = (totstuffssstatic+1)./((meanclustersize)'*ones(1,length(kdecaymath)));

kubykd = kumat./kdmat.*skfmat;
hc = (kubykd.*srcconcmat)./(1+kubykd+zmat./vc);
mathmeanConcss = hc.*(zmat + vc)./(4*pi/3*(R^3-a.^3));
totstuffssfrag = totstuffssstatic + mathmeanConcss.*(Nfree'*ones(size(kdecaymath)));

%% Finally, plot the result
figure
hold on
color = [0.206085 0.717553 0.473947 1];

data = (meanConcss(1,:)*nedges-1);
plot(1./kdecay,data,'-','color',color,'Linewidth',2)
plot(1./kdecaymath,totstuffssfrag,'--','color',color,'Linewidth',2);

set(gca,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',20)
xlabel("$\tau_d$ (decay time)")
ylabel("S (steady-state content)")
set(gca,'xscale','log')
set(gca,'yscale','log')

legend("Simulation","Analytic Model","Location","southeast","Interpreter","latex")

xticks(logspace(-4,4,3))
xlim([1e-4 1e4])
ylim([0.05 nedges])
set(gcf,'color','w');
