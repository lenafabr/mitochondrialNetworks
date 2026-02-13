function totstuffsocial = get_social_limit_material(kdecaymath,totstuffssstatic,meanclustersize,meanclustercomp,meanradg,nedges,meanmitodiff,kfiss,R,ku1,ku2,mitolen,rc,rs,kT,mu,B1,B2,alpha1,alpha2)
    % kdecaymath [1 x N]: array of decay rates for which to calculate the steady-state material content
    % totstuffssstatic [imax x N]: matrix containing results from the
    % static model for each network type and decay rate
    % meanclustersize [imax x 1]: the average cluster size (number of
    % units) for each network type
    % meanclustercomp [imax x 2]: matrix containing the average number of tips and degree-2 nodes
    % per cluster for each network type
    % meanradg [imax x 1]: the average radius of gyration per cluster for
    % each network type
    % nedges: the total number of units in the network
    % meanmitodiff [imax x 1]: the mean cluster diffusivity for each
    % network type
    % kfiss: the fission rate per node
    % R: the size of the simulation cell boundary
    % ku1: the local tip-tip fusion rate for nearby nodes
    % ku2: the local tip-side fusion rate for nearby nodes ...
    
    nnettypes = size(meanclustercomp,1);

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
    tiptipfactor = (meanclustercomp(:,1).^2)'./frac;
    tipsidefactor = 2*(meanclustercomp(:,1).*meanclustercomp(:,2))'./frac;
    
    ku = tiptipfactor.*predict_pesc_deg2_public3(ku1,B1,alpha1,mitolen,rs,rc,kT,mu)/...
        (0.5*4*pi/3*(((2*rc)^2 - (2*rs)^2)^(3/2) + (2*rc)^3 - (2*rs)^3)) + ...
        tipsidefactor.*predict_pesc_deg3_public2(ku2,B1,B2,alpha2,mitolen,rs,rc,kT,mu)/...
        (4*pi/3*((2*rc)^2 - (2*rs)^2)^(3/2));
    
    % solve for static soln at kd = kf
    [~,kfidx] = min(abs(log10(kdecaymath)-log10(kfiss)));
    skf = totstuffssstatic(:,kfidx);
    
    kumat = ku'*ones(1,length(kdecaymath));
    kdmat = ones(nnettypes,1)*kdecaymath;
    skfmat = ((skf+1)./(meanclustersize)')*ones(1,length(kdecaymath));
    srcconcmat = (totstuffssstatic+1)./((meanclustersize)'*ones(1,length(kdecaymath)));
    
    kubykd = kumat./kdmat.*skfmat;
    hc = (kubykd.*srcconcmat)./(1+kubykd+zmat./vc);
    mathmeanConcss = hc.*(zmat + vc)./(4*pi/3*(R^3-a.^3));
    Nfree = nedges - meanclustersize;
    totstuffsocial = mathmeanConcss.*(Nfree'*ones(size(kdecaymath)));
end