from datetime import datetime 

class ChromaticSymmetricChainComplex(SageObject):

    def __init__(self, G, verbose=False, name=None):
        r"""
        Sets up a class infrastructure to compute the Chromatic Symmetric Homology of ``G``.

        INPUT:

        - ``G`` -- a graph
        - ``verbose`` bool (default ``False``); whether to print out status messages while performing computations
        - ``name`` string (default ``G._repr_()``): the label to use when being verbose
        """
        self._G = G.relabel(range(1,G.num_verts()+1), inplace=False)
        self._verbose = verbose
        self._name = name if name else self._G._repr_()
        if verbose:
            print("%s: %s: Setting up computations"%(datetime.now(),self._name))
        # we use the broken circuit model to cut down the size of the chain complex we need to consider
        self._broken_circuits = Matroid(self._G).broken_circuits()
        self._S = SymmetricGroup(self._G.num_verts())
        self._V = GroupAlgebra(self._S)
        self._B = self._V.basis()

    @cached_method
    def _edges_module(self, edges):
        r"""
        Returns the module associated to a given subset of edges
        
        INPUT:

        - ``edges`` frozenset: the set of edges
        """
        # if the set of edges contains a broken circuit return the empty module
        if any( bc.issubset(edges) for bc in self._broken_circuits):
            return self._V.submodule(self._V(0))
        transpositions = list(map(self._S, edges))
        H = self._S.subgroup(transpositions)
        w = sum(self._B[h] for h in H)
        W = self._V.submodule([b*w for b in self._B])
        return W
    
    def _degree_1_map(self, source, target):
        r"""
        Return the map in between the module associated to ``source`` and the one associated to ``target``.

        INPUT:

        - ``source`` a set of k edges
        - ``target`` a set of k-1 edges

        WARNING: 

        No check is performed on the cardinalities of source and target.
        """
        source = frozenset(source)
        target = frozenset(target)
        U = self._edges_module(source)
        u = U.dimension()
        W = self._edges_module(target)
        w = W.dimension()
        extras =  [e for e in source if e not in target]
        e = extras.pop()
        if not extras:
            phi = W.retract.pre_compose(U.lift)
            M = matrix(u, w, map(lambda x: phi(x).to_vector(), U.basis()), sparse=True).transpose()
            sign = (-1)**sum(1 for f in source if f <= e)
            return sign*M
        return matrix(w, u, sparse=True)

    @cached_method
    def differential(self, i, j=0):
        r"""
        Return the ``i``-th differential in q-degree ``j``.

        INPUT:

        - ``i`` a non-negative integer
        - ``j`` a non-negative integer

        WARNING:

        The algorithm is only implemented in q-degree 0
        """
        if self._verbose:
            print("%s: %s: Computing differential (%d, %d)"%(datetime.now(),self._name,i,j))
        if j != 0:
            raise NotImplementedError("At the moment we can only compute in degree 0.")
        if i == 0:
            return matrix(0, self._V.dimension(), sparse=True)
        sources = Subsets(self._G.edges(labels=False), i)
        targets = Subsets(self._G.edges(labels=False), i-1)
        blocks = [ [ self._degree_1_map(s,t) for s in sources] for t in targets]
        return block_matrix(blocks,sparse=True)

    def homology(self, i=None, j=None, base_ring=ZZ, algorithm='auto'):
        r"""
        Compute the homology group(s)

        INPUT:

        - ``i`` a non-negative integer or ``None`` (default). 
            If an integer is provided compute only the homology group in homological degree ``i``. 
            Otherwise compute the homology groups in all degrees.

        - ``j`` a non-negative integer or ``None`` (default).
            If an integer is provided compute only the homology in q-degree ``j``.
            Otherwise compute the homology in all q-degrees.

            WARNING: currently only j=0 is implemented.

        - ``base_ring`` a ring (default ``ZZ``). The ring over which the computation will happen.
            Due to several limitations sprinkled here and there this can only be ``ZZ``, ``QQ``, or ``Zmod(p)`` for ``p`` prime.

        - ``algorithm`` string (default 'auto'). Which algorithm to use to compute homology (e.g. 'chomp')
        """
        if j !=0:
            raise NotImplementedError("Only implemented for q-degree 0")
        diff = { k: self.differential(k, 0) for k in ([i, i+1] if i is not None else range(self._G.num_verts())) }
        if self._verbose:
            print("%s: %s: Reducing differentials using AMT"%(datetime.now(),self._name))
        diff = amt_reduction(diff)
        if self._verbose:
            print("%s: %s: Computing homology"%(datetime.now(),self._name))
        CC = ChainComplex(diff, degree=-1, base_ring=base_ring)
        if i is not None:
            homology = {i: CC.homology(i, algorithm=algorithm)}
        else:
            homology = CC.homology(algorithm=algorithm) 
        if self._verbose:
            print("%s: %s: Done"%(datetime.now(),self._name))
        return { (k,0): homology[k] for k in homology }

def amt_reduction(diff):
    r"""
    Apply Algebraic Morse Theory to reduce the differentials in a chain complex.

    INPUT:

    - ``diff`` a dictionary ``{i:D}`` where ``D`` is the matrix of the ``i``-th differential 

    ALGORITHM:

    This is a rough implementation of the algorithm in arXiv:1903.00783 without the reordering steps
    """
    
    diff = {k:copy(D) for (k,D) in diff.items()}
   
    def matching():
        return { k:[ (u,v) for u in range(D.nrows()) for v in D.nonzero_positions_in_row(u)[:1] if D[u,v].is_unit() and all(w <= u for w in D.nonzero_positions_in_column(v)) ] for (k,D) in diff.items() }

    M = matching()
    
    while sum(M.values(),[]):
        Ip = { k:[v for _,v in N] for (k,N) in M.items() } 
        Im = { k-1:[u for u,_ in N] for (k,N) in M.items() }

        for k in diff:
            for v in Im.get(k,[]):
                diff[k][:,v] = 0
            for u in Ip.get(k-1,[]):
                diff[k][u,:] = 0

            for up in range(diff[k].nrows()):
                if up not in Ip.get(k-1,[])+Im.get(k-1,[]):
                    v = next( iter( filter( lambda x: x in Ip.get(k,[]), diff[k].row(up).nonzero_positions() )), None)
                    while v != None:
                        u = next( iter( filter(lambda x: x[1]==v, M[k])))[0]
                        diff[k][up] -= diff[k][u]*diff[k][up,v]/diff[k][u,v]
                        v = next( iter( filter( lambda x: x in Ip.get(k,[]), diff[k].row(up).nonzero_positions() )), None)
            diff[k] = diff[k].delete_rows(Im.get(k-1,[])+Ip.get(k-1,[])).delete_columns(Ip.get(k,[])+Im.get(k,[]))

        M = matching()

    return diff
