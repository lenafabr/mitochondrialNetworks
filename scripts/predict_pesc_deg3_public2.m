function [c2meanfield,pescavg,pescavgeq,vnear,ctheta1list, weight2full,weight2,Efuse,ctheta1map,drhophi] = predict_pesc_deg3_public(kulist,B1,B2,alpha2,L0,rs,rcontact,kT,mu)
% This function computes the degree-3 association constant C2 from the
% mean-field theory
% as well as the correction factors pesc and pesceq
% Inputs:
% kulist = list of fusion rate prefactors (k_{u,2})
% B1, B2 = degree 2 and degree 3 bending moduli
% alpha2 = angle sensitivity for degree 3 fusion. Set to B2/(L0*kT) by
% default.
% note: kT and mu are assumed to be equal to 1
% L0 = ground-state segment length
% D = relative diffusivity of two nodes
% tau = timescale for mechanical relaxation after fission
% rs = steric radius
% rcontact = contact radius
% kT = effective thermal energy; should be 1 by default
% mu = friction coefficient for an individual bead. Should be 1 by default.
% -----
% Outputs:
% ------
% c2meanfield = mean fields predicted deg3 association constant, for all
% fusion rates in kulist
% pescavg = escape probability for newly fissed nodes
% pescavgeq = escape probability for an equilibrium system (nodes fiss to the
% same configuration they fused at)

%% set relative diffusivity and timescale for reorientation for newly separated segments
tau_mech = mu*L0^3/B2;
tau_diff = mu*L0^2/(4*kT);
if tau_mech<tau_diff
    tau = tau_diff;
    D = 1.61*kT/mu;
else
    tau = tau_mech;
    D = 1.52*kT/mu;
end
%% assumed starting distance between newly fissed nodes, *beyond* the 2r_s
% steric contact.
% Used for estimating sterically excluded zone
%del = (rcontact-rs);

% use the average separation distance
a = 2*rs; b = 2*rcontact;
del = (b^4-a^4)/(b^3 - a^3)*3/4 - 2*rs;
%% get energies and weights for different starting conditions
% integrate over cos(theta1) (between connected segments)
% cos(theta2), phi for orientation of the free segment 

L1 = L0-rs; % length of segment to inset node

nt2 = 50; np = 50; nt1=50;
rhovals = linspace(-1,1,nt2);

point1s = [0 0 -L1];
point1e = [0 0 0];

drho = rhovals(2)-rhovals(1);

%
ctheta1list = linspace(-1,1,nt1);
cc=0; kvals = []; weight = []; cc2=0;
 Evals = [];

% estimated range of phi angles to avoid steric overlap
angmax = asin(2*rs/(2*rs+del));
phivals = linspace(angmax,2*pi-angmax,np);
dphi = phivals(2)-phivals(1);

drhophi = drho*dphi;

clear Evals weighteq weight weightA weightj
for t1c = 1:nt1 % integral over cos(theta1)
    ctheta1 = ctheta1list(t1c);
    
    point2s = [0 0 0];
    point2e = L1*[0 sqrt(1-ctheta1^2) ctheta1];
    
    for t2c = 1:nt2 % integrate  over cos(theta2)
        ct = rhovals(t2c);
        st = sqrt(1-ct.^2);                
        
        for pc = 1:np % integrate over phi
            cc = cc+1;   

            phi = phivals(pc);
            
            point3s = [-(del+2*rs) 0 0];
            point3e = (L1-rs)*[st*cos(phi),st*sin(phi),ct] + point3s;
            tip = point3s- rs*[st*cos(phi),st*sin(phi),ct];
            outer = rs * tip/sqrt(sum(tip.^2));
            rj = 0.5*(tip+outer);
            %rj_old = point3s/2; % position of putative junction
            % calculate energy at the putative 3-way junction
            % Efuse is the energy component that contributes to the fusion
            % rate
            % Ebend is the mechanical bending energy
            [Efusetmp,Ebend] = junc3energy_r(rj,point1s,point2e,point3e);
            Efuse(cc) = Efusetmp*alpha2; Ebend = Ebend*B2/L0;
                           
            weightFuse(cc) = exp(-Efuse(cc));
            weight3(cc) = exp(-Ebend/kT);
            
            % weight based on angle of connected segments only
            weight2(cc) = exp(B1/L0*ctheta1/kT);    
            ctheta1map(cc) = ctheta1;
        end
    end
    
    % weight for connected junction, regardless of the other angles
    % includes orientations that are sterically excluded as well
    weight2full(t1c) = exp(B1/L0*ctheta1/kT);       
end
% total number of configurations sampled in the integrals
ncc = cc;


% equilibrium fusion energy at long times
Elong = sum(Efuse.*weight2)/sum(weight2);

%% get eigenvalues with absorbing outer boundary
a = 2*rs; b = 2*rcontact;
E0list = linspace(0,3*max(alpha2,B2/L0)/L0,100);

typebounds = [1,0];
nmax = 300;

[eiglist,func] = getEigsSph(a,b,nmax,typebounds);
eiglist = eiglist';

% normalization factor
Nvals = 2*(eiglist.^2*a^2 + 1)./((b-a)*(eiglist.^2*a^2+1) + a);

% starting at a=2*rs
ba = eiglist*(b-a);
coeff = Nvals.*sin(ba)/a;    
Rfactint = b./eiglist.*coeff;

%% Get escape probability (post-fission and equilibrated), as a function of starting energy and ku0

clear pesc pesceq
for kuc = 1:length(kulist)
    ku0 = kulist(kuc);
    
    tlist = logspace(log10(min(1e-4,1/ku0*1e-3)),log10(min(1/ku0*1e2,0.5)),1e3);
    dt = diff(tlist);

    
    for ec = 1:length(E0list)
        E0 = E0list(ec);   

        %rxn rate integral over time
        kintegral = tau*ku0*exp(-Elong)*(-expint(E0-Elong) + expint((E0-Elong)*exp(-tlist/tau)));
      
        %%
        Dbt = D*tlist'*eiglist'.^2;   
        % timing factor t by b
        timefact = exp(-Dbt - kintegral');
        Grint = timefact*Rfactint;     
        
        % energy over time
        Etime = Elong+(E0-Elong)*exp(-tlist/tau);        
        kvals = ku0*exp(-Etime);
        integ = kvals.*Grint';
        integavg = (integ(2:end)+integ(1:end-1))/2;   
        pesc(ec,kuc) = 1-sum(integavg.*dt);
                
        % starting at uniform distribution        
        Rfactinteq = 3*b^2/(b^3-a^3)* Nvals./eiglist.^2;
        Grinteq = timefact*Rfactinteq;         
        integ = kvals.*Grinteq';
        integavg = (integ(2:end)+integ(1:end-1))/2;   
        pesceq(ec,kuc) = 1-sum(integavg.*dt);
    end        
end

%% Average escape probabilities over starting configurations

pescavg = zeros(1,length(kulist));
pescavgeq = pescavg;
kuavgeq = pescavg;


for kuc = 1:length(kulist)
    ku0=kulist(kuc);
    
    %overall fusion rate, averaged over orientations
    % weighted by the mechanical energy of the remaining fused segments    
    %kuavgeq(kuc) =sum(ku0*exp(-Efuse).*weight2)/sum(weight2);
    kuavgeq(kuc) =sum(ku0*exp(-Efuse).*weight2)/sum(weight2full);
    
    % escape probability after fission
    % weight the starting energy by mechanical Boltzman factor
    pescvals = interp1(E0list,pesc(:,kuc),Efuse); 
    pescavg(kuc) = sum(pescvals.*weight3)/sum(weight3);

    % escape probability from average radius at fusion
    % weight the starting energy by the fusion probability factor
    pesceqvals = interp1(E0list,pesceq(:,kuc),Efuse);              
    pescavgeq(kuc) = sum(pesceqvals.*weightFuse.*weight2)/sum(weightFuse.*weight2);    
    %pescavgeq(kuc) = sum(pesceqvals.*weightFuse)/sum(weightFuse);    

end


% volume within contact radius
vnear = 4*pi/3*((2*rcontact)^2 - (2*rs)^2)^(3/2);
%c2meanfield = vnear*kuavgeq;
c2meanfield = vnear*kuavgeq*drho*dphi/(4*pi);

end