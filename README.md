## ![ell_1_big](https://github.com/chrbraue/l1Houdini/blob/master/aux/ell_1_big.jpg)-HOUDINI

![ell_1](https://github.com/chrbraue/l1Houdini/blob/master/aux/ell_1_big.jpg)-HOUDINI is a primal-dual homotopy algorithm for the problem
>![p_delta](https://github.com/chrbraue/l1Houdini/blob/master/aux/p_delta.jpg)

- compute a single primal-dual optimal pair:
```matlab
[x, y] = l1Houdini(A, b, delta)
```

- compute the entire primal and dual solution paths together with the corresponding parameters:
```matlab
[pathX, pathDelta, pathY] = l1HoudiniPath(A, b, Delta)
```