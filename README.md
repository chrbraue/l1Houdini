# l1-HOUDINI

l1-HOUDINI is a homotopy algorithm for l1-minimization with maximum norm constraints.

![equation](http://www.sciweavers.org/tex2img.php?eq=1%2Bsin%28mc%5E2%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=)

[x, y] = l1Houdini(A, b, Delta) computes an optimal pair (x, y) for the problem

min ||x||_1  s.t. ||Ax - b||_inf <= Delta.

[pathX, pathDelta, pathY] = l1HoudiniPath(A, b, Delta) returns the entire primal and dual solution paths as well as the corresponding parameter path.