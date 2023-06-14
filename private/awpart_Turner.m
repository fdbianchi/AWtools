function [H,g] = awpart_Turner(G,K,sysInfo,opts)

% AWPART_TURNER computes the partial order compensator based on sector 
% boundedness as proposed in: 
%
% M. C. Turner and I. Postlethwaite, “A new perspective on static and
%   low order anti-windup synthesis,” Int. J. Control, vol. 77, no. 1,
%   pp. 27–44, Jan. 2004.

% fbianchi - 2021-07-27

% problem dimensions
ns = sysInfo.ns;
ny = sysInfo.ny;
nu = sysInfo.nu;
nv = sysInfo.nv;
nk = sysInfo.nk;
nr = sysInfo.nr;
F1 = sysInfo.F1;
F2 = sysInfo.F2;

% solver options
etol = opts.etol;
solver_opts = opts.solver_opts;
minU = opts.MinBndU;
wY = opts.weigthY;
wU = opts.weigthU;

% plant data:
[A,B,C,D] = ssdata(G);
% controller data:
[Ac,Bc,Cc,Dc] = ssdata(K);
% filter data:
n1 = order(F1);
[A1,B1,C1,D1] = ssdata(F1);
n2 = order(F2);
[A2,B2,C2,D2] = ssdata(F2);

% checking B, C, D or Dcy
dB = 0; dC = 0;   dD = 0;  
nD = 0; nDcy = 0;
flagLdp = 0;
if (nv > 1)
    B0 = B(:,:,1);
    C0 = C(:,:,1);
    D0 = D(:,:,1);
    for ii = 1:nv
        dB   = max(dB,norm(B(:,:,ii) - B0));
        dC   = max(dC,norm(C(:,:,ii) - C0));
        dD   = max(dD,norm(D(:,:,ii) - D0));
        nD   = max(nD,norm(D(:,:,ii)));
        nDcy = max(nDcy,norm(Dc(:,(1:ny)+nr,ii)));
    end
    if (nD > 0) && (nDcy > 0)
        error('In case of LPV systems, D or Dcy must be zero')
    end
    if (dB > 0) && (dC > 0) && (dD > 0)
        error('In case of LPV systems, B, C and D must be parameter independent')
    end
    % check if V can be parameter dependant
    flagLdp = (nD == 0);
    D = D0;
end

% optimization variables
Y = sdpvar(ns+nk+n1+n2);
if flagLdp
    Lu = sdpvar(nu,nu,nv,'full');
    auxLy = sdpvar(nu,nu,'full');
    Ly = repmat(auxLy,1,1,nv);
    L  = [Lu;Ly];
else
    auxL = sdpvar(nu+ny,nu,'full');
    L = repmat(auxL,1,1,nv);
end
Ud = sdpvar(nu,1);
U = diag(Ud);
g = sdpvar(1);

lmis = [];
for ii=1:nv
    
    if (nr == 0)
        Bcy = -Bc(:,:,ii);
        Dcy = -Dc(:,:,ii);
    else
        Bcy = Bc(:,nr+1:end,ii);
        Dcy = Dc(:,nr+1:end,ii);
    end
    
    % augmented plant
    J   = eye(ny)/(eye(ny) - D*Dcy);
    Jt  = eye(nu)/(eye(nu) - Dcy*D);
    Ab  = [A(:,:,ii) + B(:,:,ii)*Jt*Dcy*C(:,:,ii),...
           B(:,:,ii)*Jt*Cc(:,:,ii),...
           B(:,:,ii)*Jt*C1, -B(:,:,ii)*Jt*Dcy*C2;
           Bcy*J*C(:,:,ii),...
           Ac(:,:,ii) + Bcy*J*D*Cc(:,:,ii),...
           Bcy*J*D*C1, -Bcy*J*C2;
           zeros(n1,ns+nk), A1, zeros(n1,n2);
           zeros(n2,ns+nk+n1), A2];
    Bo  = [B(:,:,ii)*Jt;        Bcy*J*D; zeros(n1+n2,nu)];
    Bb  = [B(:,:,ii)*Jt*D1,    -B(:,:,ii)*Jt*Dcy*D2; 
           Bcy*J*D*D1, -Bcy*J*D2;
           B1, zeros(n1,ny);
           zeros(n2,nu), B2];
    C1b = [Jt*Dcy*C(:,:,ii), Jt*Cc(:,:,ii), Jt*C1, -Jt*Dcy*C2];
    D01 = Jt*Dcy*D;
    D1b = [Jt*D1, -Jt*Dcy*D2];
    C2b = [J*C(:,:,ii), J*D*Cc(:,:,ii),...
           J*D*C1,     -J*D*Dcy*C2];
    D02 = J*D;
    D2b = [J*D*D1 -J*D*Dcy*D2];
    
    % constraints
    const  = blkvar;
    const(1,1) = (Ab*Y) + (Ab*Y)';
    const(1,2) = Bo*U + Bb*L(:,:,ii) - Y*C1b';
    aux = U + D01*U + D1b*L(:,:,ii);
    const(2,2) = -(aux + aux');
    const(1,3) = zeros(ns+nk+n1+n2,nu);
    const(2,3) = eye(nu);
    const(3,3) = -g*eye(nu);
    const(1,4) = Y*C2b';
    const(2,4) = (D02*U + D2b*L(:,:,ii))';
    const(3,4) = zeros(nu,ny);
    const(4,4) = -g*eye(ny);
    const = sdpvar(const);
    
    lmis = [lmis, const <= -etol*eye(size(const))];
end
lmis = [lmis, Y >= etol*eye(size(Y))];
lmis = [lmis, Ud >= minU];

% objective function
obj = g + wY*trace(Y) + wU*trace(U);

% solver
diagnostics = optimize(lmis,obj,solver_opts);
if diagnostics.problem == 1
    % Infeasible
    H = [];
    g = inf;
    fprintf('\n')
    fprintf('The design is infeasible with the settings:\n')
    fprintf('\t MinBndU = %d\n',opts.MinBndU)
    fprintf('\t weigthY = %d\n',opts.weigthY)
    fprintf('\t weigthU = %d\n',opts.weigthU)
    fprintf('Try to use the opts argument to relax these constraints\n')
else
    % Feasible
    U = value(U); 
    su = svd(U);
    fprintf('min Sigma(U) = %d, max Sigma(U) = %d\n',min(su),max(su))
    L = value(L);
    H = zeros(nu+ny,nu,nv);
    for ii = 1:nv
        H(:,:,ii) = L(:,:,ii)/U;
    end
    g = value(g);
end
