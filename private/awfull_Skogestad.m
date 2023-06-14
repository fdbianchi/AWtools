function [H,g] = awfull_Skogestad(G,sysInfo,opts)

% AWFULL_SKOGESTAD computes the full order compensator based on coprime
% factors and sector boundeness as proposed in
%
% S. Skogestad and I. Postlethwaite, Multivariable feedback control:
%   analysis and design. Chichester, UK: John Wiley and Sons Ltd, 2005.
%   Multivariable Feedback (page 489)

% fbianchi - 2021-07-27

% problem dimensions
ns = sysInfo.ns;
ny = sysInfo.ny;
nu = sysInfo.nu;
nv = sysInfo.nv;

% solver options
etol = opts.etol;
solver_opts = opts.solver_opts;
wY = opts.weigthY;
wU = opts.weigthU;

% plant data:
[A,B,C,D] = ssdata(G);

% optimization variables
Y = sdpvar(ns);
U = diag(sdpvar(nu,1));
g = sdpvar(1);

% checking if H can be parameter dependent
flagVdp = 1;
if (nv > 1)
    B0 = B(:,:,1);
    D0 = D(:,:,1);
    for ii = 1:nv
        if (norm(B(:,:,ii)-B0) > 1e-6) || (norm(D(:,:,ii)-D0) > 1e-6)
            Warning('B or D is parameter dependent, then Taw will be constant')
            flagVdp = 0;
        end
    end
end
if flagVdp
    L = sdpvar(nu,ns,nv,'full');
else
    Laux = sdpvar(nu,ns,'full');
    L = repmat(Laux,1,1,nv);
end

lmis = [];
for ii = 1:nv
    
    % constraints on:
    %   u_check -> [ud,yd]
    const = blkvar;
    AclY = A(:,:,ii)*Y + B(:,:,ii)*L(:,:,ii);
    const(1,1) = (AclY + AclY');
    const(1,2) = B(:,:,ii)*U - L(:,:,ii)';
    const(2,2) = -2*U;
    const(1,3) = zeros(ns,nu);
    const(2,3) = eye(nu);
    const(3,3) = -g*eye(nu);
    const(1,4) = (C(:,:,ii)*Y + D(:,:,ii)*L(:,:,ii))';
    const(2,4) = (D(:,:,ii)*U)';
    const(3,4) = zeros(nu,ny);
    const(4,4) = -g*eye(ny);
    const = sdpvar(const);
    lmis = [lmis, const <= -etol*eye(size(const))];
    
    % pole location
    if (opts.MinDamping > 0)
        beta = atan(sqrt(1-opts.MinDamping^2)/opts.MinDamping);
        Mpoles = [sin(beta) -cos(beta); cos(beta) sin(beta)];        
        constPoles = kron(Mpoles,AclY) + kron(Mpoles',AclY');
        lmis = [lmis, constPoles <= -etol*eye(size(constPoles))];
    end
    if ~isinf(opts.MaxFreq)
        Lpoles = -opts.MaxFreq*eye(2);
        Mpoles = [0 0;1 0];
        constPoles = kron(Lpoles,Y) + kron(Mpoles,AclY) + kron(Mpoles',AclY');
        lmis = [lmis, constPoles <= -etol*eye(size(constPoles))];
    end
    if (opts.MinDecay > 0)
        Lpoles = 2*opts.MinDecay;
        constPoles = Lpoles*Y + AclY + AclY';
        lmis = [lmis, constPoles <= -etol*eye(size(constPoles))];
    end
    
end
lmis = [lmis, Y >= etol*eye(ns)];

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
    fprintf('\t MinDecay   = %6.2f\n',opts.MinDecay)
    fprintf('\t MinDamping = %6.2f\n',opts.MinDamping)
    fprintf('\t MaxFreq    = %d\n',opts.MaxFreq)
    fprintf('\t weigthY    = %d\n',opts.weigthY)
    fprintf('\t weigthU    = %d\n',opts.weigthU)
    fprintf('Try to use the OPTS argument to relax these constraints\n')
else
    % Feasible
    L = value(L);
    Y = value(Y);
    H = zeros(nu,ns,nv);
    for ii = 1:nv
        H(:,:,ii) = L(:,:,ii)/Y;
    end
    g = value(g);
end

