% This is a matlab file to be used in analysis of 10 example simulations (example_0_0_j)
% this will show you how to visualize the network structures over time and
% create a plot similar to that shown in Fig. 3B of our manuscript

% NOTE: you must also go to github and add:
% github.com/lenafabr/networktools/NetworkObj.m
% to your path for this script to work correctly

%% file location parameters
mitochondrialNetworksPath = '../'; % replace this with the path to mitochondrialNetworks on your computer
simloc = mitochondrialNetworksPath+"param_files/";
snaploc = simloc+"example_";
ffloc = simloc+"example_";

%% read in a particular simulation file to visualize the dynamic network evolution
% simulation index (0-9)
j = 0;
networks = parseDynNetworkSnapshots(snaploc+"0_0_"+string(j)+".snap.out");

%% play the network evolution as a movie in a figure
% cell size 
hw = 5;
% number of snapshots to skip between viewings
viewinterval = 100;

figure
for sc = 1:viewinterval:length(networks)
    pause(1)
    hold off
    NT = networks(sc); 
    NT.plotNetwork();
    xlim([-hw hw]);
    ylim([-hw hw]);
    zlim([-hw hw]);
    view(45,30);
    xlabel("x")
    ylabel("y")
    zlabel("z")
    title("snapshot "+string(sc))
    grid on
    hold on
    drawnow;
end

%% simulation parameters to make C2 plot
nedges = 250;
mitolen = 0.5;
R = 5;
rs = 0.1;
rc = 0.15;
ku1 = 1500;
ku2 = logspace(1,5,10);
kfiss = 1;
bendmod2 = 6;
bendmod1 = 2;
alpha2 = 12;
alpha1 = 4;
kt = 1;
mu = 1;
%length of the simulation
nsteps = 200e4;
stepsize = 1e-4;
%how often are snapshot files output by the simulation
snapinterval = 1000;
nsnaps = round(nsteps/snapinterval)+1;
%what fraction of the simulation to omit when calculating averages. Use
%this to make sure you have reached steady state
frac = 0.7;

start = round(frac*nsteps/snapinterval)+1;
finish = nsnaps;

%% setup output data matrices
imax = length(ku1);
jmax = length(ku2);

% these will store the mean number of degrees of each type and number of
% fusions of each type for each simulation
mean_ndeg = zeros(imax,jmax,3);
mean_kfiss_deg3 = zeros(imax,jmax);
mean_kfuse_deg3 = zeros(imax,jmax);
mean_kfiss_deg2 = zeros(imax,jmax);
mean_kfuse_deg2 = zeros(imax,jmax);

%% run through the set of simulations, count fusions, fissions, and nodes
% if you do not have the parallel processing package for MATLAB, you can
% use a regular for loop instead of parfor
parfor j = 1:jmax
    for i = 1:imax
        plotstr = 'i='+string(i-1)+'/'+string(imax-1)+', j='+string(j-1)+'/'+string(jmax-1);
        disp('Running '+plotstr)
        % load the networks for a particular simulation
        networks = parseDynNetworkSnapshots(snaploc+"0_"+string(i-1)+"_"+string(j-1)+".snap.out");
        
        % find the average number of nodes of each degree over the course
        % of the simulation
        degrees = [0 0 0];
        for sc = start:finish
            NT = networks(sc);
            NT.setupNetwork();
            degrees = degrees + [sum(NT.degrees==1) sum(NT.degrees==2) sum(NT.degrees==3)]/(finish-start);
        end
        mean_ndeg(i,j,:) = degrees;

        % load the recorded fusion and fission events and count the frequency of each 
        eventsfile = ffloc+"0_"+string(i-1)+"_"+string(j-1)+".ffevents.out";
        events = readmatrix(eventsfile,'FileType','text');
        if(size(events,2) < 4)
            disp("no events found")
            mean_kfiss_deg3(i,j) = 0; mean_kfuse_deg3(i,j) = 0;
            mean_kfiss_deg2(i,j) = 0; mean_kfuse_deg2(i,j) = 0;
        else
            first = size(events,1);
            last = first;
            while(events(first,4) > frac*nsteps && first > 1)
                first = first-1;
            end
            while(events(last,4) > nsteps && last > 1)
                last = last-1;
            end
            first = first+1;

            fuse3idxs = events(:,1) == -3;
            fiss3idxs = events(:,1) == 3;
            fuse2idxs = events(:,1) == -2;
            fiss2idxs = events(:,1) == 2;
        
            tottime = nsteps*(1-frac)*stepsize;
            
            mean_kfiss_deg3(i,j) = sum(fiss3idxs)/tottime; 
            mean_kfuse_deg3(i,j) = sum(fuse3idxs)/tottime;
            mean_kfiss_deg2(i,j) = sum(fiss2idxs)/tottime;
            mean_kfuse_deg2(i,j) = sum(fuse2idxs)/tottime;
        end
    end
end

%% calculate C2 from the simulations using both the fusion rate and the node counts
c2fromnodes = 1.5*mean_ndeg(:,:,3)./(mean_ndeg(:,:,1).*mean_ndeg(:,:,2));

a2 = mean_kfuse_deg3./(mean_ndeg(:,:,1).*mean_ndeg(:,:,2));
b2 = mean_kfiss_deg3./mean_ndeg(:,:,3);
c2fromrates = 1.5*a2./b2;

%% estimate C2 using the analytic model
c2theory = zeros(imax,jmax);
psim2 = zeros(imax,jmax);
peq2 = zeros(imax,jmax);

Vtot = 4/3*pi*(R-mitolen)^3 +17/12*pi*mitolen^3 - 4*mitolen^2*pi*R + 3*mitolen*pi*R^2;

for i = 1:imax
    disp("i="+string(i))
    [c2theory(i,:),psim2(i,:),peq2(i,:)] = predict_pesc_deg3_public2(ku2,bendmod1,bendmod2,alpha2,mitolen,rs,rc,kt,mu);
end

%% make the c2 vs ku2 plot
i = 1;
figure
color = [0.123463,0.581687,0.547445];
plot(ku2,c2fromnodes(i,:)*nedges,'-o','Color',color,'MarkerFaceColor',color,'LineWidth',2)
hold on
plot(ku2,c2theory(i,:)*nedges/Vtot,'--k','LineWidth',3)
plot(ku2,c2theory(i,:).*peq2(i,:)./psim2(i,:)*nedges/Vtot,':k','LineWidth',3)

set(gca,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',20)
xlabel("$k_{u2}/k_f$")
set(gca,'xscale','log')
set(gca,'yscale','log')
ylabel("$\hat{\rho}C_2$")
set(gcf,'color','w');
grid on

xticks([1e1 1e3 1e5])
xlim([min(ku2(:)) max(ku2(:))])
ylim([1e-2 1e4])
