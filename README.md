![ell_1_big](https://github.com/chrbraue/l1Houdini/blob/master/aux/ell_1_big.jpg)-HOUDINI

![ell_1](https://github.com/chrbraue/l1Houdini/blob/master/aux/ell_1.jpg)-HOUDINI is a primal-dual homotopy algorithm for the problem

![p_delta](https://github.com/chrbraue/l1Houdini/blob/master/aux/p_delta.jpg)

This repository includes a MATLAB implementation of the algorithm and some illustrating examples.

- Compute a single primal-dual optimal pair:
```matlab
[x, y] = l1Houdini(A, b, delta)
```

- Compute the entire primal and dual solution paths together with the corresponding parameters:
```matlab
[pathX, pathDelta, pathY] = l1HoudiniPath(A, b, 0)
```

Further reading:

- Brauer, Christoph, Lorenz, Dirk A. and Tillmann, Andreas M. <html>&nbsp;</html> [A primal-dual homotopy algorithm for ![ell_1](https://github.com/chrbraue/l1Houdini/blob/master/aux/ell_1.jpg)-minimization with ![ell_1](https://github.com/chrbraue/l1Houdini/blob/master/aux/ell_infty.jpg)-constraints.](http://em.rdcu.be/wf/click?upn=KP7O1RED-2BlD0F9LDqGVeSIkdqD3za-2Fu8mgBnf6P3LnA-3D_Rzh4RH5OrDsKsdHGLAwaROgUZ0o-2Bcob5mkrZAHb-2B7Y2F3kWNskI0wbS0BGpHCU7T8B2D0ndSPWhLuTsfxAB9fV8sCLhX34ViYh9Vs562uRbtkGAQeLpaa8wJrRav1os7UnQJ-2FNhC3wT7VWZ73HquXQu6wg4hzdB43w2QvOiqN2yR4G6aXwuLoG3WjXucNXl0hqTSAxztZSRhJE0x6bpn33p5Qx6jEMxsAnq8ddDKNNZhNpi8NYmBKxanqgdGrDyRHelF1LPEbErmKorhJobfXPXg5Y7814t0jxywtyp26RY-3D) *Computational Optimization and Applications*, February 2018.

<img src="https://github.com/chrbraue/l1Houdini/blob/master/aux/houdini.jpg" width="200">