function [c1meanfield,pescavg,pescavgeq] = predict_pesc_deg2_public3(kulist,B1,alpha1,L0,rs,rcontact,kT,mu)
% This function computes the degree-2 association constant C1 from the
% mean-field theory
% as well as the correction factors pesc and pesceq
% Inputs:
% kulist = list of fusion rate prefactors (k_{u,1})
% B1 = degree 2 bending modulus
% alpha1 = angle sensitivity for degree 2 fusion. Set to B1/(L0*kT) by
% default
% L0 = ground-state edge length. Should be 0.5um by default
% rs = steric radius
% rcontact = contact radius
% kT = effective thermal energy; should be 1 by default
% mu = friction coefficient for an individual bead. Should be 1 by default.
% -----
% Outputs:
% ------
% c1meanfield = mean fields predicted deg2 association constant, for all
% fusion rates in kulist
% pescavg = escape probability for newly fissed nodes
% pescavgeq = escape probability for an equilibrium system (nodes fiss to the
% same configuration they fused at)

%% set relative diffusivity and timescale for reorientation for newly separated segments
D = 5/3*kT/mu;
tau = mu/kT*L0^2/8;

%% estimate fraction of contact volume
% volume within contact radius
a = 2*rs; b = 2*rcontact;
vnear = 0.5*4*pi/3*((b^2 - a^2)^(3/2) + b^3 - a^3);

%% compute weights for all possible orientations

Ebend = [];

% rho1 = cos(angle) from first segment to inter-node distance
% rho2 = cos(angle), phi for orientation of second segment
Lfrag = L0-2*rs;
nrho1 = 50; nrho2 = 400; nphi = 50;
rho1min = -sqrt(1-(rs/rcontact)^2);
rho1list = linspace(rho1min,1,nrho1);
point1s = [0 0 -Lfrag];
point1e = [0 0 0];
v1 = (point1e-point1s); v1 = v1/norm(v1);

philist = linspace(0,2*pi,nphi+1);
philist = philist(1:end-1); % so as not to double-count

% limit how bent the junction can be to possibly fuse based on alpha1
rho2list = linspace(max(1-40/alpha1,-1),1,nrho2);

ct = 0;
for cc1 = 1:nrho1
    rho1 = rho1list(cc1);
    % this does not account for sterically forbidden orientations
    if(rho1 > 0)
        del = 3/4*(b^4-a^4)/(b^3-a^3) - 2*rs;
    else
        rmin = a/sqrt(1-rho1.^2);
        del = 3/4*(b^4-rmin^4)/(b^3-rmin^3) - 2*rs;
    end

    for cc2 = 1:nrho2
        rho2 = rho2list(cc2);
        st2 = sqrt(1-rho2.^2);

        point2s = (2*rs+del)*[sqrt(1-rho1^2) 0 rho1];
        for pc = 1:nphi
            phi = philist(pc);
            point2e = point2s + Lfrag*[cos(phi)*st2 sin(phi)*st2 rho2];

            % junction point
            v2 = point2e - point2s;
            v2 = v2/norm(v2);
            ptj = (point1e + rs*v1 + point2s-rs*v2)/2;
            % cosine angle at the junction
            d1 = ptj-point1s; d2 = point2e-ptj;
            cang = (d1)*(d2)'/(norm(d1)*norm(d2));

            ct = ct+1;
            % mechanical bending energy
            Ebend(ct) = B1/L0*(1-cang);
            % weighting factor for fusion
            Efuse(ct) = alpha1*(1-cang);
        end
    end
end

da = (rho1list(2)-rho1list(1))*(rho2list(2)-rho2list(1))*(philist(2)-philist(1));
% weighting factor for fused junction orientation
weight1 = exp(-Ebend/kT);
% weighting factor for which orientation fusion occurs in
weightFuse = exp(-Efuse);

% average energy factor for fusion at long times, when orientations are
% uniform
Elong = alpha1;

%% get eigenvalues with absorbing outer boundary
typebounds = [1,0];
nmax = 100;
E0list = linspace(0,10*B1/L0,100);

[eiglist,func] = getEigsSph(a,b,nmax,typebounds);
eiglist = eiglist';

% normalization factor
Nvals = 2*(eiglist.^2*a^2 + 1)./((b-a)*(eiglist.^2*a^2+1) + a);

% starting at 2*rs
rstart = a;
ba = eiglist*(b-rstart);
coeff = Nvals.*sin(ba)/rstart;    
Rfactint = b./eiglist.*coeff;

%% Get escape probability (post-fission and equilibrated), as a function of starting energy and ku0


clear pesc pesceq
for kuc = 1:length(kulist)
    ku0 = kulist(kuc);
    
    tlist = logspace(log10(min(1e-4,1/ku0*1e-3)),log10(min(1/ku0*1e2,0.5)),1e3);
    dt = diff(tlist);  
    
    for ec = 1:length(E0list)
        % initial energy
        E0 = E0list(ec);        
        
        % get escape probability for starting at a, varying energy from E0
        % to Elong
        % fusion rates integrated over time
        kintegral = tau*ku0*exp(-Elong)*(-expint(E0-Elong) + expint((E0-Elong)*exp(-tlist/tau)));        
        %
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
                

        %get escape probability for starting at uniform distrib, energy E0            
        Rfactinteq = 3*b^2/(b^3-a^3)* Nvals./eiglist.^2;
        Grinteq = timefact*Rfactinteq;         
        integ = kvals.*Grinteq';
        integavg = (integ(2:end)+integ(1:end-1))/2;
        pesceq(ec,kuc) = 1-sum(integavg.*dt);
    end   
end

%% Average escape probabilities over starting configurations

Ebend(Ebend<0) = 0;
Efuse(Efuse<0) = 0;


pescavg = zeros(1,length(kulist));
pescavgeq = pescavg;
kuavgeq = pescavg;

%dca=calphalist(2)-calphalist(1);
for kuc = 1:length(kulist)
    ku0=kulist(kuc);
    
    % escape probability after fusion, average weighted by bending
    % Boltzmann factor
    pescvals = interp1(E0list,pesc(:,kuc),Ebend);          
    pescavg(kuc) = sum(pescvals.*weight1)/sum(weight1);

    % escape probability from the point where fusion first occurred
    % average weighted by fusion rate factor
    pesceqvals = interp1(E0list,pesceq(:,kuc),Efuse);    
    pescavgeq(kuc) = sum(pesceqvals.*weightFuse)/sum(weightFuse);  

    % average fusion rate, over all orientations
    kuavgeq(kuc) =sum(ku0*exp(-Efuse))*da/(4*pi*(1-rho1min));
end

c1meanfield = vnear*kuavgeq;
%c1meanfield = vnear*kulist*0.5*(1-exp(-2*alpha1))/alpha1;

end