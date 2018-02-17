```matlab
[x, y] = l1Houdini(A, b, Delta)
```
computes an optimal pair (x, y) for the problem<br><br>


> ![p_delta](https://github.com/chrbraue/l1Houdini/blob/master/images/p_delta.jpg).


```matlab
[pathX, pathDelta, pathY] = l1HoudiniPath(A, b, Delta)
```
returns the entire primal and dual solution paths as well as the corresponding parameter path.