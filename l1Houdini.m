function [x, y, k] = l1Houdini(A, b, Delta)
%L1HOUDINI - l1-minimization under linf-constraints.
%l1Houdini(A, b, Delta) solves the problem
%min ||x||_1 s.t. ||Ax - b||_inf <= Delta
%using a homotopy method.
%
% Syntax:  [x, y, k] = l1Houdini(A, b, Delta)
%
% Input:
%    A - matrix of size mxn
%    b - vector of size m
%    Delta - positive scalar
%    
% Output:
%    x - primal solution
%    y - dual certificate
%    k - number of performed iterations
%
% Example: 
%    A = randn(10, 15);
%    b = randn(10, 1);
%    Delta = 0.1;
%    [x, y, k] = l1Houdini(A, b, Delta);
%
% Other m-files required: primalActiveSet.m,  dualActiveSet.m
% Subfunctions: none
% MAT-files required: none
%
% See also: PRIMALACTIVESET,  DUALACTIVESET

% Authors: Christoph Brauer, Dirk A. Lorenz (TU Braunschweig)
% and Andreas M. Tillmann (TU Darmstadt)
% contact email address: ch.brauer@tu-braunschweig.de
% November 2015; Last revision: 14-February-2017

% check whether x = 0 is optimal
if norm(b, inf) <= Delta
    x = zeros(size(A, 1), 1);
    y = zeros(size(A, 2), 1);
    k = 0;
    return;
end

% initialize numerical tolerances, ...
supTol = 1e-9;
conTol = 1e-7;
optTol = 1e-6;
% (if the algorithm gets stuck at some point, it may help to change these
% tolerances)
% ... output frequency, ...
dispFreq = 100;
% ... problem size, ...
[m, n] = size(A);
% ... homotopy parameter, ...
delta = norm(b, inf);
% ... primal iterate, support and active set, ...
x = zeros(n, 1);
ip = abs(abs(b) - delta) <= conTol;
jp = false(n, 1);
numelJP = 0;
% ... sign(A * x - b) (on the complement of the active set it is actually
% not the sign that we need), ...
signAXMinusB = zeros(m, 1);
signAXMinusB(ip) = -sign(b(ip));
signAXMinusB(~ip) = -b(~ip);
% ... dual iterate, direction, support and active set, ...
y = zeros(m, 1);
eHat = zeros(m, 1);
eHat(ip) = linsolve(-sign(b(ip))', 1);
id = false(m, 1);
jd = false(n, 1);
numelJD = 0;
% ... A' * y, ...
aTY = zeros(n, 1);
% ... complements, ...
ipMinusID = ip & ~id;
jdMinusJP = jd & ~jp;
% ... transposed matrix, ...
aT = A';
% ... options for linprog (which we use as an lp solver in some cases), ...
opt = optimoptions('linprog', 'Display', 'off');
% ... and counter.
k = 1;

% main loop
while true

    % dual update
    [y, aTY, id, jd, numelJD, ipMinusID, jdMinusJP, dHat] = dualActiveSet(y, aTY, id, jd, numelJD, signAXMinusB, ipMinusID, jdMinusJP, eHat, A, aT, m, n, supTol, conTol, opt);

    % primal update
    [x, t, signAXMinusB, ip, jp, numelJP, ipMinusID, jdMinusJP, eHat] = primalActiveSet(x, signAXMinusB, ip, jp, numelJP, aTY, ipMinusID, jdMinusJP, dHat, delta, A, aT, b, Delta, m, n, supTol, conTol, opt);
    delta = delta - t;

    % optimality check (original problem)
    if abs(delta - Delta) < optTol || t == 0
        if dispFreq < inf
            fprintf('%d\t%d\t%d\t%d\t%d\t%d\n', k, delta, sum(jp), sum(ip), sum(id), sum(jd));
            fprintf('solved\n');
        end
        return;
    end
    % (if t is zero, then the original problem is not feasible and the
    % function returns the solutions associated to the smallest delta such
    % that the problem is feasible)
    
    % display iteration, homotopy parameter and support size
    if mod(k, dispFreq) == 0
        fprintf('%d\t%d\t%d\t%d\t%d\t%d\n', k, delta, sum(jp), sum(ip), sum(id), sum(jd));
    end
    
    % increase iteration counter
    k = k + 1;
    
end