### Exact Evaluation of Catmull-Clark Subdivision limit surface ###

This repo is a simple implementation of the algorithm implemented by
Jos Stam in his paper:


**Exact Evaluation of Catmull-Clark Subdivision Surfaces at Arbitrary Parameter Values**,
Jos Stam, SIGGRAPH 98 Conference Proceedings, Annual Conference
Series, July 1998

[http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/sig98.pdf](http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/sig98.pdf)

#### Usage ####
This is implemented as a C# Class Library project, which when compiled
produces a DLL. The evaluation of the limit surface patches can be
done by calling the method ```EigenEvalStructure.EvaluateSurface```.

```
Parameters:
* double[] controlPointCoords: 
Array of control points has the length 2*N+8 where N is the valence of
this patch, sorted in the order described in the Stam's paper. The
coordinates of the control points should be packed into the flat array
in the form - [x1, x2, x3, x4, ... y1, y2, y3, y4, ... z1, z2, z3, z4 ...]
* double[][] uvParams:
This should be a nested array, representing the u,v coordinates of the
points to be evalutated. The inner arrays are expected to be of length 2.
* EvaluationType evalType:
This enumeration tells the function what property to evaluate at the
given parameters. See the documentation of this enumeration for more details.
* double uvTolerance:
This will be used as the minimum value of u and v since the evaluation
at an extraordinary points is not possible (see the paper). If not
supplied, then the default value is used defined as a const member
variable (1e-12).

Return Value:
If the evalType is set to evaluate the position (or tangent or
normal), then a flat array of coordinates of the evaluated points (or
normal vectors or tangent vectors) is returned in the form - [x1, y1,
z1, x2, y2, z2 ...].
```

#### Acknowledgements ####
All thanks to Jos Stam, obviously for the original paper, but also for
helping me with some details of this implementation.
