![l_1](https://github.com/chrbraue/l1Houdini/blob/master/images/ell_1.jpg)-HOUDINI is a homotopy algorithm for sparse recovery with maximum norm constraints.



```matlab
[x, y] = l1Houdini(A, b, Delta)
```
computes an optimal pair (x, y) for the following problem:<br><br>


> ![p_delta](https://github.com/chrbraue/l1Houdini/blob/master/images/p_delta.jpg)


```matlab
[pathX, pathDelta, pathY] = l1HoudiniPath(A, b, Delta)
```
returns the entire primal and dual solution paths as well as the corresponding parameter path.