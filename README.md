![equation](http://www.sciweavers.org/tex2img.php?eq=%5Cell_1&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)-HOUDINI is a homotopy algorithm sparse recovery with maximum norm constraints.

**[x, y] = l1Houdini(A, b, Delta)** computes an optimal pair (x, y) for the problem

![equation](http://www.sciweavers.org/tex2img.php?eq=%5Cmin_x%20%5C%20%5CVert%20x%5CVert_1%20%5Cquad%20%5Cmathrm%7Bs.t.%7D%20%5C%20%5CVert%20Ax%20-%20b%5CVert_%7B%5Cinfty%7D%20%5Cleq%20%5Cdelta&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

**[pathX, pathDelta, pathY] = l1HoudiniPath(A, b, Delta)** returns the entire primal and dual solution paths as well as the corresponding parameter path.