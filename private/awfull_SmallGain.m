function [H,g] = awfull_SmallGain(G,sysInfo,opts)

% AWFULL_SMALLGAIN computes the full order compensator based on coprime
% factors and Small Gain Theorem

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
wG = opts.weigthGamma;

% optimization variables
Y = sdpvar(ns);
if (wG == 1)
    gu = 1;    
else
    gu = sdpvar(1);
end
g  = sdpvar(1);

% checking if H can be parameter dependent
flagVdp = 1;
if (nv > 1)
    B0 = G.B(:,:,1);
    D0 = G.D(:,:,1);
    for ii = 1:nv
        Bi = G.B(:,:,ii);
        Di = G.D(:,:,ii);
        if (norm(Bi-B0) > 1e-6) || (norm(Di-D0) > 1e-6)
            warning('B or D is parameter dependent, then Taw will be constant')
            flagVdp = 0;
        end
    end
end
if flagVdp
    V = sdpvar(nu,ns,nv,'full');
else
    Vaux = sdpvar(nu,ns,'full');
    V = repmat(Vaux,1,1,nv);
end

% plant data:
[A,B,C,D] = ssdata(G);

lmis = [];
for ii = 1:nv
    
    % constraints on:
    %   u_check -> ud
    AclY = A(:,:,ii)*Y + B(:,:,ii)*V(:,:,ii);
    const  = blkvar;
    const(1,1) = (AclY+AclY')';
    const(1,2) = B(:,:,ii);
    const(2,2) = -eye(nu);
    const(1,3) = V(:,:,ii)';
    const(2,3) = zeros(nu);
    const(3,3) = -gu*eye(nu);
    const = sdpvar(const);
    lmis = [lmis, const <= -etol*eye(size(const))];

    %   u_check -> yd
    const  = blkvar;
    const(1,1) = (AclY+AclY')';
    const(1,2) = B(:,:,ii);
    const(2,2) = -eye(nu);
    const(1,3) = (C(:,:,ii)*Y+D(:,:,ii)*V(:,:,ii))';
    const(2,3) = D(:,:,ii)';
    const(3,3) = -g*eye(ny);
    const = sdpvar(const);
    lmis = [lmis, const <= -etol*eye(size(const))];

    %   pole location
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
if (wG == 1)
    obj = g + wY*trace(Y);
else
    obj = gu + wG*g + wY*trace(Y);
end    

% solver
diagnostics = optimize(lmis,obj,solver_opts);
if diagnostics.problem == 1
    % Infeasible
    H = [];
    g = inf;
    fprintf('\n')
    if (wG == 1)
        fprintf('The design is infeasible, ')
        fprintf('try to set opts.weigthGammause < 1\n')
        fprintf('this will substitute the hard constraint ||M-I|| < 1 for ')
        fprintf('||M-I|| < gamma_u\n\n')
    else
        fprintf('The design is infeasible with the settings:\n')
        fprintf('\t MinDecay    = %6.2f\n',opts.MinDamping)
        fprintf('\t MinDamping  = %6.2f\n',opts.MinDamping)
        fprintf('\t MaxFreq     = %d\n',opts.MaxFreq)
        fprintf('\t weigthY     = %d\n',opts.weigthY)
        fprintf('\t weigthGamma = %d\n',opts.weigthGamma)
        fprintf('Try to use the opts argument to relax these constraints\n')
    end
else
    % Feasible
    Y = value(Y);
    V = value(V);
    H = zeros(nu,ns,nv);
    for ii = 1:nv
        H(:,:,ii) = V(:,:,ii)/Y;
        
        % to check
        Mpoles = ss(A(:,:,ii)+B(:,:,ii)*H(:,:,ii),B(:,:,ii),H(:,:,ii),eye(nu));
        nM = norm(Mpoles-eye(nu),inf);
        if (nM >= 1.001)
            fprintf('||M-1||=%6.3f > 1 \n',nM);
        end
    end
    g = sqrt(value(g));
end

