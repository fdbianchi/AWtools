function [Taw,Kaw] = computeCtrl(G,K,H,sysInfo)

% COMPUTECTRL computes the controller and Aw compensator given the plant
% G, the controller K and the gain H

% fbianchi - 2021-07-27

% problem dimensions and options
ns = sysInfo.ns;
ny = sysInfo.ny;
nu = sysInfo.nu;
nv = sysInfo.nv;
nk = sysInfo.nk;
nr = sysInfo.nr;
mth = sysInfo.method;
F1 = sysInfo.F1;
F2 = sysInfo.F2;

if strcmp(mth(1),'f')
    % case full order

    % plant data:
    [A,B,C,D] = ssdata(G);
    
    % anti-windup compensator
    Aaw  = zeros(ns,ns,nv);
    Baw  = B;
    Cawu = H;
    Dawu = zeros(nu,nu,nv);
    Cawy = zeros(ny,ns,nv);
    Dawy = D;
    for ii = 1:nv
        Aaw(:,:,ii)  = A(:,:,ii) + B(:,:,ii)*H(:,:,ii);
        Cawy(:,:,ii) = C(:,:,ii) + D(:,:,ii)*H(:,:,ii);
    end
    Taw = ss(Aaw,Baw,[Cawu;Cawy],[Dawu;Dawy]);
    
    % controller w/ anti-windup
    % first, check if Kaw would be polytopic
    nD = 0;
    if (nv > 1)
        nD = norm(D(:,:,1));
        for ii = 2:nv
            nD = max(nD,norm(D(:,:,ii)));
        end
        if (nD > 0)
            warning('Cannot compute Kaw for an LPV plant with nonzero D')
        end
    end
    if (nD == 0)
        [Ac,Bc,Cc,Dc] = ssdata(K);
        Bcr = Bc(:,1:nr,:);
        Bcy = Bc(:,nr+1:end,:);
        Dcr = Dc(:,1:nr,:);
        Dcy = Dc(:,nr+1:end,:);
        
        Au = zeros(ns+nk,ns+nk,nv);
        Bu = zeros(ns+nk,nr+ny+nu,nv);
        Cu = zeros(nu,ns+nk,nv);
        Du = zeros(nu,nr+ny+nu,nv);
        if (nr == 0)
            Cawy = -Cawy;
            Dawy = -Dawy;
        end
        for ii = 1:nv
            Au(:,:,ii) = [Ac(:,:,ii),  Bcy(:,:,ii)*Cawy(:,:,ii);
                zeros(ns,nk) Aaw(:,:,ii)];
            Bu(:,:,ii) = [Bcr(:,:,ii), Bcy(:,:,ii), Bcy(:,:,ii)*Dawy(:,:,ii);
                zeros(ns,nr), zeros(ns,ny), Baw(:,:,ii)];
            Cu(:,:,ii) = [Cc(:,:,ii), (Dcy(:,:,ii)*Cawy(:,:,ii)-Cawu(:,:,ii))];
            Du(:,:,ii) = [Dcr(:,:,ii), Dcy(:,:,ii), (Dcy(:,:,ii)*Dawy(:,:,ii)-Dawu(:,:,ii))];
        end
        Kaw = ss(Au,Bu,Cu,Du);
    else
        Kaw = [];
    end
    
elseif strcmp(mth(1),'p')
    % case low order

    % anti-windup compensator
    [A1,B1,C1,D1] = ssdata(F1); n1 = order(F1);
    [A2,B2,C2,D2] = ssdata(F2); n2 = order(F2);
    Aaw  = blkdiag(A1,A2);
    Cawu = [C1 zeros(nu,n2)];
    Cawy = [zeros(ny,n1) C2];
    for ii = 1:nv
        Hu = H(1:nu,:,ii);
        Hy = H(nu+1:end,:,ii);
        Baw(:,:,ii)  = [B1*Hu; B2*Hy];
        Dawu(:,:,ii) = D1*Hu;
        Dawy(:,:,ii) = D2*Hy;
    end
    Taw = ss(Aaw,Baw,[Cawu;Cawy],[Dawu;Dawy]);
    na = n1 + n2;
    
    % controller w/ anti-windup
    [Ac,Bc,Cc,Dc]  = ssdata(K);
    Bcr = Bc(:,1:nr,:);
    Bcy = Bc(:,nr+1:end,:);
    Dcr = Dc(:,1:nr,:);
    Dcy = Dc(:,nr+1:end,:);
    
    Au = zeros(na+nk,na+nk,nv);
    Bu = zeros(na+nk,nr+ny+nu,nv);
    Cu = zeros(nu,na+nk,nv);
    Du = zeros(nu,nr+ny+nu,nv);
    if (nr == 0)
        Cawy = -Cawy;
        Dawy = -Dawy;
    end
    for ii = 1:nv
        Au(:,:,ii) = [Ac(:,:,ii),  Bcy(:,:,ii)*Cawy; 
                      zeros(na,nk) Aaw];
        Bu(:,:,ii) = [Bcr(:,:,ii), Bcy(:,:,ii), Bcy(:,:,ii)*Dawy(:,:,ii);
                      zeros(na,nr), zeros(na,ny), Baw(:,:,ii)];
        Cu(:,:,ii) = [Cc(:,:,ii), (Dcy(:,:,ii)*Cawy-Cawu)];
        Du(:,:,ii) = [Dcr(:,:,ii), Dcy(:,:,ii), (Dcy(:,:,ii)*Dawy(:,:,ii)-Dawu(:,:,ii))];
    end
    Kaw = ss(Au,Bu,Cu,Du);
    

elseif strcmp(mth(1),'s')
    % case static
    
    % anti-windup compensator
    Dawu = H(1:nu,:,:);
    Dawy = H(nu+1:end,:,:);
    Taw = ss(H);
    
    % controller w/ anti-windup
    [Ac,Bc,Cc,Dc]  = ssdata(K);
    Bcr = Bc(:,1:nr,:);
    Bcy = Bc(:,nr+1:end,:);
    Dcr = Dc(:,1:nr,:);
    Dcy = Dc(:,nr+1:end,:);
    
    Bu = zeros(nk,nr+ny+nu,nv);
    Du = zeros(nu,nr+ny+nu,nv);
    if (nr == 0)
        Dawy = -Dawy;
    end
    for ii = 1:nv
        Bu(:,:,ii) = [Bcr(:,:,ii), Bcy(:,:,ii), Bcy(:,:,ii)*Dawy(:,:,ii)];
        Du(:,:,ii) = [Dcr(:,:,ii), Dcy(:,:,ii), (Dcy(:,:,ii)*Dawy(:,:,ii)-Dawu(:,:,ii))];
    end
    Kaw = ss(Ac,Bu,Cc,Du);
    
else
    error('COMPUTECTRL:inputError','Invalid AW method')

end

