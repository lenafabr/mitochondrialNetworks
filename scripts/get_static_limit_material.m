function [totstuffssstatic, totstuffss1d] = get_static_limit_material(lambda,df,mgd,mcl,mec,mitolen)
    % lambda: array of diffusive lengthscales for which to calculate the steady-state
    % material content
    % df: the intrinsic (graph distance) fractal dimension of the network
    % mgd: the mean graph distance between nodes in the network (in um)
    % mcl: the mean cluster size (number of units) of the network
    % mec: the average number of units connected to each unit. Needed for
    % 1d spoke calculation
    % mitolen: the length of a single unit in um
    if(~exist("mec","var") || ~exist("mitolen","var"))
        mec = inf;
        mitolen = inf;
    end

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

    % get the d-dimensional continuum result
    sdf = (sdf_1+sdf_2)/2;
    totstuffssstatic = sdf;

    % optionally, calculate the 1d result
    totstuffss1d = mec*lambda/mitolen.*tanh(mcl*mitolen/mec./lambda);
end