# ![ell_1_big](https://github.com/chrbraue/l1Houdini/blob/master/aux/ell_1_big.jpg)-HOUDINI
>![p_delta](https://github.com/chrbraue/l1Houdini/blob/master/aux/p_delta.jpg)

- compute a single primal-dual optimal pair (x, y) associated with ![delta_geq_0](https://github.com/chrbraue/l1Houdini/blob/master/aux/delta_geq_0.jpg):
```matlab
[x, y] = l1Houdini(A, b, delta)
```

- compute the entire primal and dual solution paths together with the corresponding parameters:
```matlab
[pathX, pathDelta, pathY] = l1HoudiniPath(A, b, Delta)
```