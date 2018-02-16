# l1-HOUDINI

l1-HOUDINI is a homotopy algorithm for l1-minimization with maximum norm constraints.

[x, y] = l1Houdini(A, b, Delta) computes an optimal pair (x, y) for the problem

min ||x||_1  s.t. ||Ax - b||_inf <= Delta.

[pathX, pathDelta, pathY] = l1HoudiniPath(A, b, Delta) returns the entire primal and dual solution paths as well as the corresponding parameter path.