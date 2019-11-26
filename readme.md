This is the `sagemath` code used for computations when preparing
[arXiv:bla](https://arxiv.org/abs/bla).

Here is a sample session:

```
sage: %attach chromatic_symmetric_homology.py
sage: g = graphs.CompleteGraph(5)
sage: CC = ChromaticSymmetricChainComplex(g, verbose=True, name = "Complete graph on 5 vertices")
2019-11-26 11:11:35.151003: Complete graph on 5 vertices: Setting up computations
sage: CC.homology(j=0)
2019-11-26 11:11:44.415276: Complete graph on 5 vertices: Computing differential (0, 0)
2019-11-26 11:11:44.534713: Complete graph on 5 vertices: Computing differential (1, 0)
2019-11-26 11:11:44.916352: Complete graph on 5 vertices: Computing differential (2, 0)
2019-11-26 11:11:46.391448: Complete graph on 5 vertices: Computing differential (3, 0)
2019-11-26 11:11:51.476128: Complete graph on 5 vertices: Computing differential (4, 0)
2019-11-26 11:12:04.606883: Complete graph on 5 vertices: Reducing differentials using AMT
2019-11-26 11:12:05.885761: Complete graph on 5 vertices: Computing homology
2019-11-26 11:12:05.933610: Complete graph on 5 vertices: Done
{(0, 0): Z, (1, 0): Z^24 x C2^5, (2, 0): Z^91, (3, 0): Z^24}
```
