function [y, aTY, id, jd, numelJD, ipMinusID, jdMinusJP, dHat] = dualActiveSet(y, aTY, id, jd, numelJD, signAXMinusB, ipMinusID, jdMinusJP, eHat, A, aT, m, n, supTol, conTol, opt)
%DUALACTIVESET - Active set method for dual subproblem of l1-HOUDINI
%
% Syntax:  [y, aTY, id, jd, numelJD, ipMinusID, jdMinusJP, dHat] = dualActiveSet(y, aTY, id, jd, numelJD, signAXMinusB, ipMinusID, jdMinusJP, eHat, A, aT, m, n, supTol, conTol, opt)
%
% Inputs:
%    y            - last dual certificate
%    aTY          - A' * y
%    id           - last dual support
%    jd           - last dual active set
%    numelJD      - size of last dual active set
%    signAXMinusB - sign(A * x - b) for last primal iterate x
%    ipMinusID    - ip \ id
%    jdMinusJP    - jd \ jp
%    eHat         - initial direction
%    A            - mxn matrix (from original problem)
%    aT           - transpose of A
%    m            - row dimension of A
%    n            - column dimension of A
%    supTol       - numerical tolerance for support estimation
%    conTol       - numerical tolerance for active set estimation
%    opt          - linprog options
%
% Outputs:
%    y            - new dual certificate
%    aTY          - A' * y
%    id           - new dual support
%    jd           - new dual active set
%    numelJD      - size of new dual active set
%    ipMinusID    - ip \ id
%    jdMinusJP    - jdMinusJP
%    dHat         - final Lagrange multipliers
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: l1HOUDINI,  PRIMALACTIVESET

% Authors: Christoph Brauer, Dirk A. Lorenz (TU Braunschweig)
% and Andreas M. Tillmann (TU Darmstadt)
% contact email address: ch.brauer@tu-braunschweig.de
% June 2016; Last revision: 14-February-2018

% initialize Lagrange multipliers, ...
dHat = zeros(m, 1);
mu = zeros(n, 1);
nu = zeros(m, 1);
% ... direction, ...
e = zeros(m, 1);
% ... A' * e, ...
aTE = zeros(n, 1);
% ... auxiliary variables, ...
actSetCand = false(n, 1);
actSetAlpha = inf(n, 1);
actSetAdd = false(n, 1);
actSetRem = false(n, 1);
suppCand = false(m, 1);
suppAlpha = inf(m, 1);
suppAdd = false(m, 1);
suppRem = false(m, 1);
% ... step size, ...
alpha = 1;
% ... and counter
l = 0;

% adapt support based on initial direction (cf. equation (34))
suppAdd(ipMinusID) = abs(eHat(ipMinusID)) > 0;
numelSuppAdd = sum(suppAdd);
id(suppAdd) = true;
ipMinusID(suppAdd) = false;

% adapt active set based on initial direction (cf. equation (34))
% note that the dual support has been adapted above so that we can use
% eHat(id) instead of eHat(ip) here
actSetRem(jdMinusJP) = abs(aT(jdMinusJP, id) * eHat(id)) > 0;
numelActSetRem = sum(actSetRem);
jd(actSetRem) = false;
numelJD = numelJD - numelActSetRem;
jdMinusJP(actSetRem) = false;

% iterate
while true    
    
    % determine e
    if l == 0 % (cf. subsection 3.3)
        e(id) = eHat(id);   
    else % (cf. equation (26))
        e(id) = linsolve([aT(jd, id); signAXMinusB(id)'], [zeros(numelJD, 1); 1]);
        if contains(lastwarn, 'singular')
            % troubleshoot in case matrix is close to singular and linsolve
            % fails (a solution may exist, though)
            try
                e(id) = linprog(zeros(sum(id), 1), [], [], [aT(jd, id); signAXMinusB(id)'], [zeros(numelJD, 1); 1], [], [], [], opt);                
            catch
                [e(id), flag] = pcg([A(id, jd), signAXMinusB(id)] * [aT(jd, id); signAXMinusB(id)'], [A(id, jd), signAXMinusB(id)] * [zeros(numelJD, 1); 1], 1e-4, 100); %#ok<*ASGLU>
            end
            lastwarn('');
        end
    end
    
    % calculate error
    conErr = norm(aT(jd, id) * e(id), inf);
    objErr = abs(signAXMinusB(id)' * e(id) - 1);
    
    % update aTE only on JD^c (cf. equation (28))
    aTE(~jd) = aT(~jd, id) * e(id);

    % determine alpha (cf. equations (27)-(29))
    alphaOld = alpha;
    actSetCand(~jd) = abs(aTE(~jd)) > conTol;
    actSetAlpha(actSetCand) = (1 - sign(aTE(actSetCand)) .* aTY(actSetCand)) ./ abs(aTE(actSetCand));
    suppCand(id) = e(id) .* signAXMinusB(id) < -supTol;
    suppAlpha(suppCand) = -y(suppCand) ./ e(suppCand);
    alpha = min([actSetAlpha; suppAlpha]);
    
    % check whether the step alpha * e is feasible
    if max(conErr, objErr) <= 1e-6 && alpha * conErr <= conTol && alphaOld >= 0
        
        if alpha > 0

            % update iterate
            y(id) = y(id) + alpha * e(id);
            
            % update aTY only on JD^c (cf. equation (30))
            aTY(~jd) = aT(~jd, id) * y(id);
            
            % update active set (cf. equation (30))
            actSetAdd(actSetCand) = abs(abs(aTY(actSetCand)) - 1) <= conTol;
            actSetAdd(actSetRem) = abs(abs(aTY(actSetRem)) - 1) <= conTol;
            % (cf. section A.7 and the case alpha > 0)
            jd(actSetAdd) = true;
            numelJD = numelJD + sum(actSetAdd);
            jdMinusJP(actSetAdd) = true;
            
            % update support (cf. equation (30))
            suppRem(suppCand) = abs(y(suppCand)) <= supTol;
            suppRem(suppAdd) = abs(y(suppAdd)) <= supTol;
            % (cf. section A.7 and the case alpha > 0)
            id(suppRem) = false;
            ipMinusID(suppRem) = true;

            % update aTY (the following is already approximately true by
            % the definition of actSetAdd above; however, we assign the
            % exact sign vector to aTY for reasons of numerical stability)
            aTY(actSetAdd) = sign(aTY(actSetAdd));
            
            % reset auxiliary variables
            actSetRem(actSetRem) = false;
            numelActSetRem = 0;
            suppAdd(suppAdd) = false;
            numelSuppAdd = 0;            
        else
            % (cf. section A.7 and the case alpha = 0)
            
            % update active set
            actSetAdd(actSetRem) = actSetAlpha(actSetRem) <= 0;
            jd(actSetAdd) = true;
            numelJD = numelJD + sum(actSetAdd);
            jdMinusJP(actSetAdd) = true;
            actSetRem(actSetAdd) = false;
            numelActSetRem = numelActSetRem - sum(actSetAdd);
            
            % update support
            suppRem(suppAdd) = suppAlpha(suppAdd) <= 0;            
            id(suppRem) = false;
            ipMinusID(suppRem) = true;
            suppAdd(suppRem) = false;
            numelSuppAdd = numelSuppAdd - sum(suppRem);            
        end        
        
        % reset auxiliary variables
        actSetAdd(actSetAdd) = false;
        suppRem(suppRem) = false;   
    else
        
        % determine index of smallest Lagrange multiplier
        dHat(jd) = linsolve(A(id, jd), signAXMinusB(id));
        % (cf. equation (31) with flipped sign)
        if contains(lastwarn, 'singular')
            % troubleshoot in case matrix is close to singular and linsolve
            % fails (a solution may exist, though)
            if l > 0
                try
                    dHat(jd) = linprog(zeros(sum(jd), 1), [], [], A(id, jd), signAXMinusB(id), [], [], [], opt);
                catch
                    [dHat(jd), flag] = pcg(aT(jd, id) * A(id, jd), aT(jd, id) * signAXMinusB(id), 1e-4, 100);
                end
            end
            lastwarn('');
        end
        if any(ipMinusID) % (cf. equation (33) with flipped sign)
            nu(ipMinusID) = signAXMinusB(ipMinusID) .* (A(ipMinusID, jd) * dHat(jd)) - 1;
        end
        if any(jdMinusJP) % (cf. equation (32) with flipped sign)
            mu(jdMinusJP) = aTY(jdMinusJP) .* dHat(jdMinusJP);
        end
        [~, j] = min(mu);
        [~, i] = min(nu);
        
        % check optimality / update active set or support
        if min(mu(j), nu(i)) >= - supTol
            
            % update variables and return
            jd(actSetRem) = true;
            numelJD = numelJD + sum(actSetRem);
            jdMinusJP(actSetRem) = true;
            id(suppAdd) = false;
            ipMinusID(suppAdd) = true;
            return;
            
        elseif mu(j) < nu(i)
            
            % update active set
            jd(j) = false;
            numelJD = numelJD - 1;
            jdMinusJP(j) = false;
            actSetRem(j) = true;
            numelActSetRem = numelActSetRem + 1;
            
            % reset j-th multipliers
            dHat(j) = 0;
            mu(j) = 0;
        else
            
            % update support
            id(i) = true;
            ipMinusID(i) = false;
            suppAdd(i) = true;
            numelSuppAdd = numelSuppAdd + 1;
            
            % reset i-th multiplier
            nu(i) = 0;
        end
    end
        
    % reset auxiliary variables
    actSetAlpha(actSetCand) = inf;
    actSetCand(actSetCand) = false;
    suppAlpha(suppCand) = inf;
    suppCand(suppCand) = false;
        
    % increase counter
    l = l + 1;
end