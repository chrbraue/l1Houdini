l1-HOUDINI is a homotopy algorithm sparse recovery with maximum norm constraints.


**[x, y] = l1Houdini(A, b, Delta)** computes an optimal pair (x, y) for the problem


![p delta](https://redaktionssystem.tu-braunschweig.de/default/NPS/preview/b_38/Medien-DB/iaa/brauer/p_delta.jpg)


**[pathX, pathDelta, pathY] = l1HoudiniPath(A, b, Delta)** returns the entire primal and dual solution paths as well as the corresponding parameter path.