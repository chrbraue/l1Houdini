# l1-HOUDINI

l1-HOUDINI is a homotopy algorithm for l1-minimization with maximum norm constraints.

The call [x, y, k] = l1Houdini(A, b, Delta) computes a primal solution x and a dual solution y for the problem min ||x||_1  s.t. ||Ax - b||_inf <= Delta. k is the number of iterations.