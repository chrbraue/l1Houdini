function [x, t, signAXMinusB, ip, jp, numelJP, ipMinusID, jdMinusJP, eHat] = primalActiveSet(x, signAXMinusB, ip, jp, numelJP, aTY, ipMinusID, jdMinusJP, dHat, delta, A, aT, b, Delta, m, n, supTol, conTol, opt)
%PRIMALACTIVESET - Active set method for primal subproblem of l1-HOUDINI
%
% Syntax:  [x, t, signAXMinusB, ip, jp, numelJP, ipMinusID, jdMinusJP, eHat] = primalActiveSet(x, signAXMinusB, ip, jp, numelJP, aTY, ipMinusID, jdMinusJP, dHat, delta, A, aT, b, Delta, m, n, supTol, conTol, opt)
%
% Inputs:
%    x            - last primal solution
%    signAXMinusB - sign(A * x - b) on ip / A * x - b on complement of ip
%    ip           - last primal active set
%    jp           - last primal support
%    numelJP      - size of last primal support
%    aTY          - A' * y for last dual certificate y
%    ipMinusID    - ip \ id
%    jdMinusJP    - jd \ jp
%    dHat         - initial direction
%    delta        - last homotopy parameter
%    A            - mxn matrix (from original problem)
%    aT           - transpose of A
%    b            - m vector (from original problem)
%    Delta        - target homotopy parameter (from original problem)
%    m            - row dimension of A
%    n            - column dimension of A
%    supTol       - numerical tolerance for support estimation
%    conTol       - numerical tolerance for active set estimation
%    opt          - linprog options
%
% Outputs:
%    x            - new primal solution
%    t            - decrease of homotopy parameter
%    signAXMinusB - sign(A * x - b) on ip / A * x - b on complement of ip
%    ip           - new primal active set
%    jp           - new primal support
%    numelJP      - size of new primal support
%    ipMinusID    - ip \ id
%    jdMinusJP    - jd \ jp
%    eHat         - final Lagrange multipliers
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: l1HOUDINI,  DUALACTIVESET

% Authors: Christoph Brauer, Dirk A. Lorenz (TU Braunschweig)
% and Andreas M. Tillmann (TU Darmstadt)
% contact email address: ch.brauer@tu-braunschweig.de
% June 2016; Last revision: 14-February-2018

% initialize decrease of homotopy parameter, ...
t = 0;
% ... Lagrange multipliers, ...
eHat = zeros(m, 1);
mu = zeros(m, 1);
nu = zeros(n, 1);
% ... direction, ...
d = zeros(n, 1);
% ... A * d, ...
aD = zeros(m, 1);
% ... auxiliary variables, ...
actSetUpperCand = false(m, 1);
actSetLowerCand = false(m, 1);
actSetUpperAlpha = inf(m, 1);
actSetLowerAlpha = inf(m, 1);
actSetUpperAdd = false(m, 1);
actSetLowerAdd = false(m, 1);
actSetRem = false(m, 1);
suppCand = false(n, 1);
suppAlpha = inf(n, 1);
suppAdd = false(n, 1);
suppRem = false(n, 1);
% ... step size, ...
alpha = 1;
% ... and counter
l = 0;

% adapt support based on initial direction (cf. equation (35))
suppAdd(jdMinusJP) = abs(dHat(jdMinusJP)) > 0;
numelSuppAdd = sum(suppAdd);
jp(suppAdd) = true;
numelJP = numelJP + numelSuppAdd;
jdMinusJP(suppAdd) = false;

% adapt active set based on initial direction (cf. equation (35))
% note that the primal support has been adapted above so that we can use
% dHat(jp) instead of dHat(jd) here
actSetRem(ipMinusID) = abs(A(ipMinusID, jp) * dHat(jp) - signAXMinusB(ipMinusID, 1)) > 0;
numelActSetRem = sum(actSetRem);
ip(actSetRem) = false;
ipMinusID(actSetRem) = false;

% update signAXMinusB (on the complement of the active set
% we do not need the sign but A * x - b itself)
signAXMinusB(actSetRem) = A(actSetRem, jp) * x(jp) - b(actSetRem, 1);

% iterate
while true
    
    % determine d
    if l == 0 % (cf. subsection 3.3; note that the sign of dHat is flipped)     
        d(jp) = -dHat(jp);
    else % (cf. equation (16))
        d(jp) = linsolve(-A(ip, jp), signAXMinusB(ip));
        if contains(lastwarn, 'singular')
            % troubleshoot in case matrix is close to singular and linsolve
            % fails (a solution may exist, though)
            try
                d(jp) = linprog(zeros(sum(jp), 1), [], [], -A(ip, jp), signAXMinusB(ip), [], [], [], opt);
            catch
                [d(jp), flag] = pcg(aT(jp, ip) * A(ip, jp), -aT(jp, ip) * signAXMinusB(ip), 1e-4, 100); %#ok<*ASGLU>
            end
            lastwarn('');
        end
    end
    
    % calculate error
    conErr = norm(A(ip, jp) * d(jp) + signAXMinusB(ip), inf);
    
    % update aD only on IP^c (cf. equation (18))
    aD(~ip) = A(~ip, jp) * d(jp);

    % determine alpha (cf. equations (17)-(19))
    alphaOld = alpha;
    actSetUpperCand(~ip) = aD(~ip) > -1 + conTol;
    actSetUpperAlpha(actSetUpperCand) = (delta - t - signAXMinusB(actSetUpperCand)) ./ (aD(actSetUpperCand) + 1);
    actSetLowerCand(~ip) = aD(~ip) < 1 - conTol;
    actSetLowerAlpha(actSetLowerCand) = (delta - t + signAXMinusB(actSetLowerCand)) ./ (-aD(actSetLowerCand) + 1);
    suppCand(jp) = d(jp) .* aTY(jp) > supTol;
    suppAlpha(suppCand) = -x(suppCand) ./ d(suppCand);
    alpha = min([actSetUpperAlpha; actSetLowerAlpha; suppAlpha; delta - t - Delta]);
    
    % check whether the step alpha * d is feasible
    if conErr <= 1e-4 && alpha * conErr <= conTol && alphaOld >= 0

        if alpha > 0
            
            % update iterates
            x(jp) = x(jp) + alpha * d(jp);
            t = t + alpha;

            % update signAXMinusB (on the complement of the active set
            % we do not need the sign but A * x - b itself)
            signAXMinusB(~ip) = A(~ip, jp) * x(jp) - b(~ip, 1);

            % update active set (cf. equation (21))
            actSetUpperAdd(actSetUpperCand) = abs(delta - t - signAXMinusB(actSetUpperCand)) <= conTol;
            actSetUpperAdd(actSetRem) = abs(delta - t - signAXMinusB(actSetRem)) <= conTol;
            % (cf. section A.7 and the case alpha > 0)            
            actSetLowerAdd(actSetLowerCand) = abs(-delta + t - signAXMinusB(actSetLowerCand)) <= conTol;
            actSetLowerAdd(actSetRem) = abs(-delta + t - signAXMinusB(actSetRem)) <= conTol;
            % (cf. section A.7 and the case alpha > 0)            
            ip(actSetUpperAdd) = true;
            ip(actSetLowerAdd) = true;
            ipMinusID(actSetUpperAdd) = true;
            ipMinusID(actSetLowerAdd) = true;
            
            % update support
            suppRem(suppCand) = abs(x(suppCand)) <= supTol;
            suppRem(suppAdd) = abs(x(suppAdd)) <= supTol;
            % (cf. section A.7 and the case alpha > 0)            
            jp(suppRem) = false;
            numelJP = numelJP - sum(suppRem);
            jdMinusJP(suppRem) = true;

            % update signAXMinusB (on the active set we need the sign)
            signAXMinusB(actSetUpperAdd) = sign(signAXMinusB(actSetUpperAdd));
            signAXMinusB(actSetLowerAdd) = sign(signAXMinusB(actSetLowerAdd));
            
            % check optimality
            if delta - t - Delta <= 1e-6
                return;
            end

            % reset auxiliary variables
            actSetRem(actSetRem) = false;
            numelActSetRem = 0;
            suppAdd(suppAdd) = false;
            numelSuppAdd = 0;            
        else
            % (cf. section A.7 and the case alpha = 0)

            % update active set
            actSetUpperAdd(actSetRem) = actSetUpperAlpha(actSetRem) <= 0;
            actSetLowerAdd(actSetRem) = actSetLowerAlpha(actSetRem) <= 0;
            ip(actSetUpperAdd) = true;
            ip(actSetLowerAdd) = true;
            ipMinusID(actSetUpperAdd) = true;
            ipMinusID(actSetLowerAdd) = true;
            actSetRem(actSetUpperAdd) = false;
            actSetRem(actSetLowerAdd) = false;
            numelActSetRem = numelActSetRem - sum(actSetUpperAdd) - sum(actSetLowerAdd);
            
            % update support
            suppRem(suppAdd) = suppAlpha(suppAdd) <= 0;
            jp(suppRem) = false;
            numelJP = numelJP - sum(suppRem);
            jdMinusJP(suppRem) = true;
            suppAdd(suppRem) = false;
            numelSuppAdd = numelSuppAdd - sum(suppRem);            
        end
        
        % reset auxiliary variables
        actSetUpperAdd(actSetUpperAdd) = false;
        actSetLowerAdd(actSetLowerAdd) = false;
        suppRem(suppRem) = false;    
    else
        
        % determine index of smallest Lagrange multiplier
        eHat(ip) = linsolve([aT(jp, ip); signAXMinusB(ip)'], [zeros(numelJP, 1); 1]);
        % (cf. equation (22))
        if contains(lastwarn, 'singular')
            % troubleshoot in case matrix is close to singular and linsolve
            % fails (a solution may exist, though)
            try
                eHat(ip) = linprog(zeros(sum(ip), 1), [], [], [aT(jp, ip); signAXMinusB(ip)'], [zeros(numelJP, 1); 1], [], [], [], opt);
            catch
                [eHat(ip), flag] = pcg([A(ip, jp), signAXMinusB(ip)] * [aT(jp, ip); signAXMinusB(ip)'], [A(ip, jp), signAXMinusB(ip)] * [zeros(numelJP, 1); 1], 1e-4, 100);
            end
            lastwarn('');
        end
        if any(jdMinusJP) % (cf. equation (24))
           nu(jdMinusJP) = -aTY(jdMinusJP) .* (aT(jdMinusJP, ip) * eHat(ip));
        end
        if any(ipMinusID) % (cf. equation (23))
           mu(ipMinusID) = signAXMinusB(ipMinusID) .* eHat(ipMinusID);
        end
        [~, i] = min(mu);
        [~, j] = min(nu);
        
        % check optimality / update active set or support
        if min(mu(i), nu(j)) > -supTol
            
            % update variables and return
            ip(actSetRem) = true;
            ipMinusID(actSetRem) = true;
            jp(suppAdd) = false;
            numelJP = numelJP - sum(suppAdd);
            jdMinusJP(suppAdd) = true;
            signAXMinusB(actSetRem) = sign(signAXMinusB(actSetRem));
            return;            
            
        elseif mu(i) < nu(j)
            
            % update active set
            ip(i) = false;
            ipMinusID(i) = false;
            actSetRem(i) = true;
            numelActSetRem = numelActSetRem + 1;
            
            % update signAXMinusB (on the complement of the active set
            % we do not need the sign but A * x - b itself)
            signAXMinusB(i) = signAXMinusB(i) * (delta - t);
            
            % reset i-th multipliers
            eHat(i) = 0;
            mu(i) = 0;
        else   
            
            % update support
            jp(j) = true;
            numelJP = numelJP + 1;
            jdMinusJP(j) = false;
            suppAdd(j) = true;
            numelSuppAdd = numelSuppAdd + 1;            
            
            % reset j-th multiplier
            nu(j) = 0;
        end
    end
    
    % reset auxiliary variables
    actSetUpperAlpha(actSetUpperCand) = inf;
    actSetLowerAlpha(actSetLowerCand) = inf;
    actSetUpperCand(actSetUpperCand) = false;
    actSetLowerCand(actSetLowerCand) = false;
    suppAlpha(suppCand) = inf;
    suppCand(suppCand) = false;        
    
    % increase counter   
    l = l + 1;
end