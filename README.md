- to compute a primal-dual optimal pair (x, y) for the problem ![p_delta](https://github.com/chrbraue/l1Houdini/blob/master/images/p_delta.jpg) use
```matlab
[x, y] = l1Houdini(A, b, Delta)
```


-
```matlab
[pathX, pathDelta, pathY] = l1HoudiniPath(A, b, Delta)
```
returns the entire primal and dual solution paths as well as the corresponding parameter path