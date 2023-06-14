function [Kaw,Taw,g] = awsyn(G,K,strMethod,nr,opts,F1,F2)

% Synthesis of anti-windup compensators based on coprime factors
%
% The synthesis relies on the works:
%   [1] M. C. Turner and I. Postlethwaite, “A new perspective on static and
%   low order anti-windup synthesis,” Int. J. Control, vol. 77, no. 1,
%   pp. 27–44, Jan. 2004.
%   [2] S. Skogestad and I. Postlethwaite, Multivariable feedback control:
%   analysis and design. Chichester, UK: John Wiley and Sons Ltd, 2005.
%   Multivariable Feedback (page 489)
%
% Use:
%   [Kaw,Taw,g] = AWSYN(G,K,[method],[nr],[options],[F1],[F2])
%
% Outputs:
% - Kaw: controller + AW compensator
%            u = Kaw*[r; y; u-sat(u)]
% - Taw: anti-windup compensator
%            [ud: yd] = Taw*(u-sat(u))
% - g: is the norm of the operator (u-sat(u)) -> yd
%
% Inputs:
% - G: plant to be controlled (lti, lmitool objects, lpvtools or psys)
% - K: stabilizing controller (lti, lmitool objects, lpvtools or psys)
% - method:
%    'fsec': full order compensator according to Skogestad et al. [2] using
%            sector boundedness. The AW compensator Taw has the same order 
%            than the plant G 
%    'fsgt': full order compensator based on small gain theorem. The
%            compensator Taw has the same order than the plant G
%    'ssec': static compensator according to Tunner et al. [1] using sector
%            boundedness.
%    'psec': low-order compensator according to Tunner et al. [1] using 
%            sector boundedness + filters F1 and F2. The compesator has the
%            sum of the orders of F1 and F2.
% - nr: number of reference inputs, by default nr = 0, e.i., the controller
%            input is e = r - y.
% - options: this is a struct with the settings for the Taw design and
%            depends on the particular method:
%
%                               default     fsec    fsgt    ssec    psec
%            opts.MinDecay        (1)         x       x
%            opts.MinDamping      (1)         x       x
%            opts.MaxFreq         (1)         x       x
%            opts.weigthY        1e-8         x       x       x       x
%            opts.weigthU        1e-1         x               x       x  
%            opts.MinBndU        1                            x       x 
%            opts.weigthGamma    0.01                 x
%   
%            (1) The first three parameters allow imposing the location of 
%            the poles of Taw, and the default values are determinated 
%            according to the poles of the plant.
% - F1, F2: are transfer functions as described in [1]
%
% See the documentation in the folder DOC for more details.

% fbianchi - 19/03/2015
% fbianchi - 2021-07-22 - rev

% =========================================================================
% default input arguments:
if (nargin < 2)
    error('AWSYN:InputError','Not Enough Input Arguments')
elseif (nargin < 3)
    strMethod = 'fsec'; 
    nr = 0;
    opts = struct();
    F1 = 0;
    F2 = 0;    
elseif (nargin < 4)
    nr = 0;
    opts = struct();
    F1 = 0;
    F2 = 0;    
elseif (nargin < 5)
    opts = struct();
    F1 = 0;
    F2 = 0;    
elseif (nargin < 6)
    F1 = 1;
    F2 = 1;    
elseif (nargin < 7)
    F2 = 1;    
end
if isempty(nr)
    nr = 0;
end

% =========================================================================
% Plant data
[G,sysInfo]  = standardizeSys(G);
[K,sysInfoK] = standardizeSys(K);
sysInfo.nk = sysInfoK.ns;
if (sysInfo.nv ~= sysInfoK.nv)
    error('AWSYN:InputError',...
        'Plant and controller must have the same number of elements');
end
if isnumeric(nr) && isscalar(nr) && (nr >= 0)
    sysInfo.nr = nr;
else 
    error('AWSYN:InputError',...
        'nr must be a positive numeric scalar');
end    
if ischar(strMethod)
    sysInfo.method = strMethod;
else
    error('AWSYN:InputError',...
        'Method must be a string with the desired design method');
end
if isnumeric(F1)
    if isinf(F1(end))
        sysInfo.F1 = mat2lti(F1);
    else
        sysInfo.F1 = ss(F1);
    end
elseif ismember(class(F1),{'ss','tf','zpk'})
    sysInfo.F1 = F1;
else
    error('AWSYN:InputError',...
        'F1 must be a valid LTI model');
end
if isnumeric(F2)
    if isinf(F2(end))
        sysInfo.F2 = mat2lti(F2);
    else
        sysInfo.F2 = ss(F2);
    end
elseif ismember(class(F2),{'ss','tf','zpk'})
    sysInfo.F2 = F2;
else
    error('AWSYN:InputError',...
        'F2 must be a valid LTI model');
end

% options
opts.etol = 1e-8;
opts.solver_opts = sdpsettings('solver','sedumi','verbose',1);
% pole location options
wnMx = 0; epsMn = inf;
for ii = 1:sysInfoK.nv
    [wn,z] = damp(K(:,:,ii));
    wnMx  = max(wnMx,max(wn));
    epsMn = min(wnMx,min(z));
end
if ~isfield(opts,'MinDecay')
    opts.MinDecay = 0;
end
if ~isfield(opts,'MinDamping')
    opts.MinDamping = epsMn;
end
if ~isfield(opts,'MaxFreq')
    opts.MaxFreq = wnMx;
end
% weigth for Lyapunov function 
if ~isfield(opts,'weigthY')
    opts.weigthY = 1e-8;
end
% weigth for diag matrux U (method:fsec/psec)
if ~isfield(opts,'weigthU')
    opts.weigthU = 1e-1;
end
% lower bound for U in static AW
if ~isfield(opts,'MinBndU')
    opts.MinBndU = 1;
end
% weigth on the u_chk -> yd
if ~isfield(opts,'weigthGamma')
    opts.weigthGamma = max(0,min(1,0.01));
end
 
% =========================================================================
% cases:
switch strMethod
    case 'fsgt'
        % full order based on small gain theorem
        [H,g] = awfull_SmallGain(G,sysInfo,opts);
        
    case 'fsec'
        % full order according to Skogestad et al.
        [H,g] = awfull_Skogestad(G,sysInfo,opts);

%   ToDo: method ssgt needs more evaluation            
%     case 'ssgt'
%         % static based on small gain theorem
%         H = awstat_SmallGain(G,K,sysInfo,opts);
        
    case 'ssec'
        % static according to Tunner et al.
        [H,g] = awstat_Turner(G,K,sysInfo,opts);
    
    case 'psec'
        % static according to Tunner et al.
        [H,g] = awpart_Turner(G,K,sysInfo,opts);

    otherwise
        error('AWSYN:inputError','Invalid design method')
end

% returns the compensator and controller + AW
if isempty(H)
    Taw = []; Kaw = [];
else
    [Taw,Kaw] = computeCtrl(G,K,H,sysInfo);
    Taw = unstandardizeSys(Taw,sysInfo);
    Kaw = unstandardizeSys(Kaw,sysInfo);
end

