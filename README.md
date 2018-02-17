![ell_1_big](https://github.com/chrbraue/l1Houdini/blob/master/aux/ell_1_big.jpg)-HOUDINI
>![p_delta](https://github.com/chrbraue/l1Houdini/blob/master/aux/p_delta.jpg):

- compute a primal-dual optimal pair (x, y):
```matlab
[x, y] = l1Houdini(A, b, Delta)
```

- compute the entire primal and dual solution paths as well as the corresponding parameter path:
```matlab
[pathX, pathDelta, pathY] = l1HoudiniPath(A, b, Delta)
```