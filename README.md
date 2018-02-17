## ![ell_1_big](https://github.com/chrbraue/l1Houdini/blob/master/aux/ell_1_big.jpg)-HOUDINI

![ell_1](https://github.com/chrbraue/l1Houdini/blob/master/aux/ell_1.jpg)-HOUDINI is a primal-dual homotopy algorithm for the problem
>![p_delta](https://github.com/chrbraue/l1Houdini/blob/master/aux/p_delta.jpg)

- Compute a single primal-dual optimal pair:
```matlab
[x, y] = l1Houdini(A, b, delta)
```

- Compute the entire primal and dual solution paths together with the corresponding parameters:
```matlab
[pathX, pathDelta, pathY] = l1HoudiniPath(A, b, Delta)
```

This repository includes a MATLAB implementation of the algorithm and some illustrating examples. A detailed derivation and discussion can be found in

>Brauer, Christoph, Lorenz, Dirk A. and Tillmann, Andreas M.  A primal-dual homotopy algorithm for ![ell_1](https://github.com/chrbraue/l1Houdini/blob/master/aux/ell_1.jpg)-minimization with ![ell_1](https://github.com/chrbraue/l1Houdini/blob/master/aux/ell_infty.jpg)-constraints.  *Computational Optimization and Applications*, February 2018.