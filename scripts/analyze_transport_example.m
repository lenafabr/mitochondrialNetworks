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
jmax = 6; % how many decay rates were run?
kmax = 1;

% simulation parameters
savesnapstart = 2000000;
bdsteps = [18810000 4580000 3160000 3020000 3000000 3000000];
simstepspersnap = [16810 2580 1160 1020 1000 1000];
nsnaps = 1001;
beginavgsnap = round((nsnaps-1)/2);
delt=1e-4;
% unit size
rs = 0.15;
rc = 0.2;
mitolen = 0.5;
% simulation dynamic parameters
ku1 = 300;
ku2 = 100;
kfiss = 1;
partdiff = 4800;
kdecay = logspace(-2.5,2.5,jmax);
% more simulation parameters
nspecies = 9;
nedges = 250;
R = 5;
kT = 1; mu = 1;
B1 = 2; alpha1 = B1/mitolen/kT; 
B2 = 6; alpha2 = B2/mitolen/kT;

% for measuring mito diffusivities
lagtime = max(simstepspersnap*delt); 



%% measure structural and dynamic parameters and material content
clustersize = zeros(kmax*imax,1); % mean cluster size
networkdim = zeros(kmax*imax,1); % network fractal dimension
graphdiststats = zeros(kmax*imax,2); % graph distance statistics
clustertipstats = zeros(kmax*imax,2); % number of tips, deg.2 nodes per cluster
edgeconnections = zeros(kmax*imax,1); % mean number of neighbors for each unit
radg = zeros(kmax*imax,1); % cluster radius of gyration
meanConcRaw = zeros(kmax*imax*jmax,nspecies,nsnaps); % mean material per unit
mitodiff = zeros(imax*jmax*kmax,1); % mean mitochondrial diffusivity by cluster

parfor idx = 1:imax*kmax*jmax
    [k, i, j] = ind2sub([kmax imax jmax],idx);
    disp("running idx = "+string(idx)+" of "+string(kmax*imax*jmax))
    networks = parseDynNetworkSnapshots(snaploc+string(k-1)+"_"+string(i-1)+"_"+string(j-1)+".snap.out");
    networks = networks(2:end); %ignore the first snapshot

    % get the average material per unit
    edgevals = {networks.edgevals};
    edgevals = cell2mat(reshape(edgevals,1,1,[]));
    meanConcRaw(idx,:,:) = squeeze(mean(edgevals,1));

    % get the average cluster diffusivity
    mitodiff(idx) = get_mean_cluster_diffusivity(networks(beginavgsnap:end),lagtime,simstepspersnap(j),delt);

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
meanmitodiff = mean(reshape(mitodiff, [kmax imax jmax]), [1 3]);

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
set(gca,'yscale','log')


%% calculate the analytic solution based on the structural and dynamic parameters
kdecaymath = logspace(-4,4,101);
lambda = sqrt(partdiff./kdecaymath);
totstuffssstatic = zeros(imax,length(kdecaymath));

% static network contribution
parfor i = 1:imax
    df = meandimension(i);
    mgd = meangraphdist(i,1);
    mcl = meanclustersize(i);
    totstuffssstatic(i,:) = get_static_limit_material(lambda,df,mgd,mcl);
end

% social network contribution
totstuffssfrag = totstuffssstatic + get_social_limit_material(kdecaymath,totstuffssstatic, ...
    meanclustersize,meanclustercomp,meanradg,nedges,meanmitodiff,kfiss,R,ku1,ku2,mitolen, ...
    rc,rs,kT,mu,B1,B2,alpha1,alpha2);


%% Finally, plot the result
figure
hold on
color = [0.206085 0.717553 0.473947 1];

data = (meanConcss(1,:)*nedges-1);
plot(1./kdecay,data,'-','color',color,'Linewidth',2)
plot(1./kdecaymath,totstuffssfrag(1,:),'--','color',color,'Linewidth',2);

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
