"""A module for analyzing plumbed manifolds.

"""
from copy import copy, deepcopy
from itertools import product, groupby
import math
import sys

import igraph as ig


class Plumbing:
    """A class for analyzing plumbed manifolds.

    Parameters
    ----------
    vertices_dict : dict
        A dictionary of the form {a:b} where a is the index of a vertex of the
        plumbing and b is its corresponding weight.
    edges : array_like
        A list of the form [(a,b)] where (a,b) represents an edge between the
        vertices of indicies a and b.

    Attributes
    ----------
    vertex_count
    edge_count
    weight_vector
    degree_vector
    intersection_form
    is_tree
    definiteness_type
    bad_vertices
    artin_fcycle
    is_weakly_elliptic
    is_rational
    is_almost_rational
    homology


    Methods
    -------


    Example
    -------
    >>> P = Plumbing({0:-3,1:-10}, [(0,1)])

    Here P is the plumbing consisting of 2 vertices joined by a single
    edge. One vertex has weight -3 and the other has weight -10.

    """


    def __init__(self, vertices_dict, edges):
        try:
            self._vertices_dict = vertices_dict
            self._edges = edges

            self._vertex_count = len(vertices_dict)
            self._edge_count = len(edges)

            self._graph = ig.Graph()
            self._graph.add_vertices(self._vertex_count)

            if self._edge_count != 0:
                self._graph.add_edges(edges)

            self._graph.vs["v_weights"] = list(vertices_dict.values())
            vertex_labels = [("v" + str(i), list(vertices_dict.values())[i])
                             for i in range(0, self._vertex_count)]
            self._graph.vs["label"] = vertex_labels

            self._weight_vector = Matrix(list(vertices_dict.values())).T
            self._degree_vector = Matrix(self._graph.degree()).T

            self._lay = self._graph.layout_reingold_tilford()
            self._plot = ig.plot(self._graph, layout=self._lay,
                                 vertex_size=int(10),
                                 vertex_label_dist=int(2),
                                 vertex_label_angle=int(0),margin=int(90),
                                 vertex_color="black",
                                 vertex_label_color="black")

            self._intersection_form = None
            self._intersection_smith_form = None
            self._is_tree = None
            self._definiteness_type = None
            self._bad_vertices = None
            self._artin_fcycle = None
            self._is_weakly_elliptic = None
            self._is_rational = None
            self._is_almost_rational = None
            self._homology = None

        except:
            print("Error: Plumbing entered incorrectly. Please check"
                  " formatting.")

    @property
    def vertex_count(self):
        """int: The number of vertices in the plumbing."""
        return self._vertex_count

    @property
    def edge_count(self):
        """int: The number of edges in the plumbing"""
        return self._edge_count

    @property
    def weight_vector(self):
        """Matrix: An sx1 matrix of the form [[m_1],...,[m_s]] where m_i is the
        weight of vertex i."""
        return self._weight_vector

    @property
    def degree_vector(self):
        """Matrix: An sx1 matrix of the form [[d_1],...,[d_s]] where d_i is the
        degree (or valence) of vertex i."""
        return self._degree_vector

    @property
    def intersection_form(self):
        """Matrix: A matrix representing the intersection form of the
        plumbing."""
        if self._intersection_form is None:
            intersection_form = Matrix(list(self._graph.get_adjacency()))
            for i in range(0, self._vertex_count):
                intersection_form[i, i] = self._weight_vector[i,0]
            self._intersection_form = intersection_form
        return self._intersection_form

    @property
    def intersection_smith_form(self):
        if self._intersection_smith_form is None:
            self._intersection_smith_form = self.intersection_form.smith_form()
        return self._intersection_smith_form

    @property
    def is_tree(self):
        """bool: True if the plumbing diagram is a finite tree, False
        otherwise."""
        if self._is_tree is None:
            self._is_tree = self._graph.is_connected()\
                            and (self._edge_count == self._vertex_count - 1)
        return self._is_tree

    @property
    def definiteness_type(self):
        """str: The definiteness type of the intersection form of the plumbing.

        Warnings
        --------
        Since the eigenvalues are computed numerically, they may contain small
        error terms. Therefore, to check the sign of an eigenvalue, we have
        chosen a small error threshold (1e-8). This potentially could lead to
        incorrect answers in some edge cases when the true eigenvalues are very
        close to zero, but non-zero.

        """
        if self._definiteness_type is None:
            eigenvalues = self.intersection_form.eigenvalues()
            if all(i < -1e-8 for i in eigenvalues):
                self._definiteness_type = "negative definite"
            elif all(i > 1e-8 for i in eigenvalues):
                self._definiteness_type = "positive definite"
            elif all(i == 0 for i in eigenvalues):
                self._definiteness_type = "zero matrix"
            elif all(i <= 1e-8 for i in eigenvalues):
                self._definiteness_type = "negative semi-definite"
            elif all(i >= -1e-8 for i in eigenvalues):
                self._definiteness_type = "positive semi-definite"
            else:
                return "positive and negative eigenvalues"
        return self._definiteness_type

    @property
    def bad_vertices(self):
        """tuple: A tuple of the form (bv, bv_count) where bv is a string
        formatted for LaTeX listing the bad vertices, and bv_count is the
        number of bad vertices.
        """
        if self._bad_vertices is None:
            bv_count = 0
            bv = ''
            test = False
            for i in range(0, self._vertex_count):
                if test and self._weight_vector[i,0] > -self._degree_vector[i,0]:
                    bv = bv + ", v_{"+str(i)+"}"
                    bv_count += 1
                else:
                    if self._weight_vector[i,0] > -self._degree_vector[i,0]:
                        bv = ": v_{"+str(i)+"}"
                        bv_count += 1
                        test = True
            if bv_count == 0:
                bv = '0 bad vertices.'
            elif bv_count == 1:
                bv = "1 bad vertex" + bv + "."
            else:
                bv = str(bv_count) + " bad vertices" + bv + "."
            self._bad_vertices = bv, bv_count
        return self._bad_vertices

    @property
    def artin_fcycle(self):
        """tuple: A tuple of the form (x, comp_seq) where x is the Artin
        fundamental cycle of the plumbing and comp_seq is the associated
        computation sequence used to compute x. The Artin fundamental cycle is
        used to determine the rationality of the plumbing graph.
        """
        if self._artin_fcycle is None and self.definiteness_type ==\
                "negative definite" and self.is_tree:
            x = [0] * self.vertex_count
            x[0] = 1
            x = Matrix(x)
            z = x*self.intersection_form
            comp_seq = [deepcopy(x)]
            while any(i > 0 for i in z.row(0)):
                j = 0
                while z[0,j] <= 0:
                    j = j + 1
                x[0,j] = x[0,j] + 1
                comp_seq.append(deepcopy(x))
                z = x * self.intersection_form
            self._artin_fcycle = x, comp_seq
        return self._artin_fcycle

    @property
    def is_weakly_elliptic(self):
        """bool: True if the plumbing is weakly elliptic, False otherwise."""
        if self._is_weakly_elliptic is None:
            if self.is_tree and self.definiteness_type == "negative definite":
                K = [-i-2 for i in self._graph.vs["v_weights"]]
                K = Matrix(K)
                m = -(K * self.artin_fcycle[0].T
                      + self.artin_fcycle[0]
                      * self.intersection_form
                      * self.artin_fcycle[0].T)[0,0] / 2
                if m == 0:
                    self._is_weakly_elliptic = True
                    self._is_rational = False
                else:
                    self._is_weakly_elliptic = False
            else:
                self._is_weakly_elliptic = False
        return self._is_weakly_elliptic

    @property
    def is_rational(self):
        """bool: True if the plumbing is rational, False otherwise."""
        if self._is_rational is None:
            if self.is_tree and self.definiteness_type == "negative definite":
                K = [-i-2 for i in self._graph.vs["v_weights"]]
                K = Matrix(K)
                m = -(K * self.artin_fcycle[0].T
                      + self.artin_fcycle[0]
                      * self.intersection_form
                      * self.artin_fcycle[0].T)[0,0] / 2
                if m == 1:
                    self._is_rational = True
                    self._is_weakly_elliptic = False
                else:
                    self._is_rational = False
            else:
                self._is_rational = False
        return self._is_rational

    @property
    def is_almost_rational(self):
        """bool/str: True if the plumbing is almost rational. Possibly
        inconclusive depending on the test threshold. This is the amount by
        which framings are decreased to test for rationality.
        """
        if self._is_almost_rational is None:
            if self.is_tree and self.definiteness_type == "negative definite":
                if self.bad_vertices[1] < 2:
                    self._is_almost_rational = True
                elif self._is_rational or self._is_weakly_elliptic:
                    self._is_almost_rational = True
                else:
                    very_bad_vert_count = 0
                    for i in range(0, self._vertex_count):
                        if -self._weight_vector[i,0] \
                                <= self._degree_vector[i,0] - 2:
                            very_bad_vert_count = very_bad_vert_count + 1
                    if very_bad_vert_count > 1:
                        self._is_almost_rational = False
                    else:
                        self._is_almost_rational = "inconclusive, using test" \
                                                   " threshold of 20"
                        counter = 1
                        while counter <= 20:
                            for i in range(0, self._vertex_count):
                                v = deepcopy(self._vertices_dict)
                                v[i] = v[i] - counter
                                plumb = Plumbing(v, self._edges)
                                K = [-j-2 for j in v.values()]
                                K = Matrix(K)
                                m = -(K * plumb.artin_fcycle[0].T
                                      + plumb.artin_fcycle[0]
                                      * self.intersection_form
                                      * plumb.artin_fcycle[0].T)[0,0] / 2
                                if m == 1:
                                    self._is_almost_rational = True
                                    break
                            counter = counter + 1
        return self._is_almost_rational

    @property
    def homology(self):
        """tuple: A tuple of the form (homology_group, homology_generators,
        rank) where homology_group is the first homology of the plumbed
        3-manifold, homology generators are the corresponding generators of
        homology_group, and rank is the Z-rank of the homology.
        """
        if self._homology is None:
            smith = self.intersection_smith_form
            D = smith[0]
            U = smith[1]
            U_inv = U.inverse()
            s = self.vertex_count
            rank = D.diagonal().count(0)
            num_of_pivots = s - rank
            invariant_factors = [D[i,i] for i in range(0, num_of_pivots)]
            p = invariant_factors.count(1)
            invariant_factors = invariant_factors[p:]
            finite_ord_coker_gens = [U_inv[:, i] for i in range(p,
                                                                num_of_pivots)]
            infinite_ord_coker_gens = [U_inv[:, i] for i in range(num_of_pivots,
                                                                  s)]
            homology_generators = []
            if rank == 0:
                if len(invariant_factors) == 0:
                    homology_group = "0"
                    homology_generators.append("N/A")
                else:
                    homology_group = "Z_{" + str(invariant_factors[0]) + "}"
                    homology_generators.append(finite_ord_coker_gens[0])
                    for i in range(1, len(invariant_factors)):
                        homology_group = homology_group + " + Z_{"\
                                         + str(invariant_factors[i]) + "}"
                        homology_generators.append(finite_ord_coker_gens[i])
            elif rank == 1:
                homology_group = "Z"
                homology_generators.append(infinite_ord_coker_gens[0])
                for i in range(0, len(invariant_factors)):
                    homology_group = homology_group + " + Z_{"\
                                     + str(invariant_factors[i]) + "}"
                    homology_generators.append(finite_ord_coker_gens[i])
            else:
                homology_group = "Z^{" + str(rank) + "}"
                for i in range(0, rank):
                    homology_generators = homology_generators.append(infinite_ord_coker_gens[i])
                for i in range(0, len(invariant_factors)):
                    homology_group = homology_group + " + Z_{"\
                                     + str(invariant_factors[i]) + "}"
                    homology_generators.append(finite_ord_coker_gens[i])

            self._homology = homology_group, homology_generators, rank

        return self._homology

    def display(self):
        """Displays a plot of the plumbing diagram."""
        self._plot.show()

    def is_in_integer_image(self, K):
        """Given a vector K, check if it is in the integer image of the
        intersection form.

        Parameters
        ----------
        K: list
            A list of integers of length = self.vertex_count.

        Returns
        -------
        bool
            True if K is in the integer image of the intersection form, False
            otherwise.

        """
        K = Matrix(K).T
        if self.definiteness_type == "negative definite":
            H = self.intersection_form.inverse()*K
            for x in H.column(0):
                if float(x) % 1 != 0:
                    return False
            return True
        else:
            smith = self.intersection_smith_form
            D = smith[0]
            U = smith[1]
            num_of_pivots = self.vertex_count - D.diagonal().count(0)
            J = U * K
            for i in range(0, num_of_pivots):
                if float(J[i, 0]) % float(D[i, i]) != 0:
                    return False
            for i in range(num_of_pivots, self.vertex_count):
                if float(J[i, 0]) != 0:
                    return False
            return True

    def equiv_spinc_reps(self, K1, K2):
        """Given two characteristic vectors, check if they represent the same
        spinc structure

        Parameters
        ----------
        K1: list
            A list of integers [k_1, ..., k_s] where s = self.vertex_count.

        K2: list
            A list of integers [k_1, ..., k_s] where s = self.vertex_count.

        Returns
        -------
        bool
            True if K1 and K2 represent the same spinc structure on the
            plumbed 3-manifold, False otherwise.

        """

        try:
            K1 = Matrix(K1).T
            K2 = Matrix(K2).T
            for i in range(0, self.vertex_count):
                if (float(K1[i, 0])-float(self.weight_vector[i, 0])) % 2 != 0:
                    raise Exception
                if (float(K2[i, 0])-float(self.weight_vector[i, 0])) % 2 != 0:
                    raise Exception
            K = (1/2)*(K1-K2)
            K = K.column(0)
            return self.is_in_integer_image(K)
        except:
            print("Error: one or more of the inputs are not a characteristic "
                  "vector.")

    def char_vector_properties(self, K):
        """Given a characteristic vector K, compute some basic properties.

        Parameters
        ----------
        K: list
            A list of integers [k_1, ..., k_s] where s = self.vertex_count.

        Returns
        -------
        tuple
            (a,b,c) where: a is a string which says if the associated spinc
            structure on the plumbed 3-manifold is torsion or non-torsion, b is
            the order of the 1st Chern class of the associated spinc structure
            on the plumbed 3-manifold, c is the square of the 1st Chern class of
            the associated spinc structure on the plumbed 4-manifold (in
            other words, c = K^2).

        """
        try:
            K = Matrix(K).T
            for i in range(0, self.vertex_count):
                if (float(K[i, 0])-float(self.weight_vector[i, 0])) % 2 != 0:
                    raise Exception

            if self.definiteness_type == "negative definite":
                H = self.intersection_form.inverse()*K
                denominators_of_H_entries = [x.denominator() for x in
                                             H.column(0)]
                order_of_chern_class = abs(lcm(denominators_of_H_entries))
                square = (K.T * H)[0, 0]
                return "Torsion", order_of_chern_class, square
            else:
                smith = self.intersection_smith_form
                D = smith[0]
                U = smith[1]
                V = smith[2]
                num_of_pivots = self.vertex_count - D.diagonal().count(0)
                J = U * K

                for i in range(num_of_pivots, self.vertex_count):
                    if J[i, 0] != 0:
                        return "Non-Torsion", "N/A", "N/A"

                H = self.vertex_count*[0]

                for i in range(0, num_of_pivots):
                    H[i] = J[i, 0]/D[i, i]

                denoms_of_non_zero_H_entries = [H[i].denominator() for i in
                                                range(0, num_of_pivots)]

                order_of_chern_class = abs(lcm(denoms_of_non_zero_H_entries))
                H = V * Matrix(H).T
                square = (K.T * H)[0,0]
                return "Torsion", order_of_chern_class, square
        except:
            print("Error: input is not a characteristic vector.")

    def chi(self, K, x):
        """
        Given a vector K and a lattice point x (represented as a vector),
        compute chi_K(x) = -1/2(K(x) + (x,x)).

        Parameters
        ----------
        K: list
            A list of integers [k_1, ..., k_s] where s = self.vertex_count.

        x: list
            A list of integers [x_1, ..., x_s] where s = self.vertex_count.

        Returns
        -------
        sage constant
            The value of chi_K(x)

        """

        K = Matrix(K)
        x = Matrix(x).T
        return -(1/2)*(K * x + x.T * self.intersection_form * x)[0,0]

    def chi_min(self, K):
        """
        Given a vector K, computes the minimum of the function chi_K on
        Euclidean space and computes the vector which achieves this minimum.
        Note this vector, in general, need not be integral.

        Parameters
        ----------
        K: list
            A list of integers [k_1, ..., k_s] where s = self.vertex_count.

        Returns
        -------
        tuple
            (a,b) where: a is the minimum value of Chi_K over R^s and b is a
            list representing the unique vector which achieves this minimum.

        """

        chi_min = self.char_vector_properties(K)[2]/8
        K = Matrix(K).T
        chi_min_vector = -(1/2) * self.intersection_form.inverse() * K
        return chi_min, list(chi_min_vector.column(0))

    def F(self, K, x):
        """
        Given a vector K and a lattice element x, computes F(x).

        Parameters
        ----------
        K: list
            A list of integers [k_1, ..., k_s] where s = self.vertex_count.
        x: list
            A list of integers [x_1, ..., x_s] where s = self.vertex_count.

        Returns
        -------
        int
            The value F(x).

        """
        K = Matrix(K).T
        x = Matrix(x).T
        y = 2*self.intersection_form*x + K - self.weight_vector\
                                           - self.degree_vector
        F = 1
        for i in range(0, self.vertex_count):
            if self.degree_vector[i,0] == 0:
                if y[i, 0] == 0:
                    F = -2*F
                elif y[i, 0] != 2 and y[i, 0]!= -2:
                    F = 0
                    return F
            elif self.degree_vector[i,0] == 1:
                if y[i, 0] == 1:
                    F = -F
                elif y[i, 0] != -1:
                    F = 0
                    return F
            elif self.degree_vector[i, 0] == 2:
                if y[i, 0] != 0:
                    F = 0
                    return F
            else:
                if abs(y[i, 0]) >= self.degree_vector[i, 0]-2:
                    F = F*(1/2)*sign(y[i,0])^(self.degree_vector[i, 0])
                    F = F*binomial((self.degree_vector[i, 0]
                                    + abs(y[i, 0]))/2-2,
                                    self.degree_vector[i, 0] -3)
                else:
                    F = 0
                    return F
        return F

    def chi_local_min_bounds(self, K):
        """
        Given a vector K, computes two lists [-chi_K(-e_1), ..., -chi(-e_s)]
        and [chi_K(e_1), ..., chi_K(e_s)] where e_i = (0, ..., 0, 1, 0, ..., 0)
        is the ith standard basis vector.

        Parameters
        ----------
        K: list
            A list of integers [k_1, ..., k_s] where s = self.vertex_count.

        Returns
        -------
        tuple
            (a,b) where: a = [-chi_K(-e_1), ..., -chi(-e_s)] and
            b = [chi_K(e_1), ..., chi_K(e_s)]

        """
        I = Matrix.identity(self.vertex_count)
        negative_I = -I
        positive_basis = [I.row(i) for i in range(0, self.vertex_count)]
        negative_basis = [negative_I.row(i) for i in
                          range(0, self.vertex_count)]
        chi_upper = [self.chi(K, x) for x in positive_basis]
        chi_lower = [-self.chi(K, x) for x in negative_basis]
        return chi_lower, chi_upper

    def chi_local_min_set(self, K):
        """
        Given a vector K, computes the set of lattice points at
        which chi_K achieves a local min, when restricted to the lattice. In
        other words, it computes the lattice points x such that
        chi_K(x) <= chi_K(x +/- e_i) for all i where
        e_i = (0, ..., 0, 1, 0, ..., 0)  is the ith standard basis vector. Note
        chi_K(x +/- e_i) = chi_K(x)+ chi_K(+/- e_i) -/+ (x, e_i). Hence, x is
        in the min set iff -chi_K(-e_i) <= (x, e_i) <= chi_K(e_i) for all i.
        This explains the reason for the helper function chi_local_min_bounds.

        Parameters
        ----------
        K: list
            A list of integers [k_1, ..., k_s] where s = self.vertex_count.
        Returns
        -------
        lists
            Each element of the output list is a tuple (a, b, c) where a is an
            element of the local min set, b is chi_K(a), and c = F(a).
        """
        if self.definiteness_type == "negative definite":
            bounds = self.chi_local_min_bounds(K)
            M_inv = self.intersection_form.inverse()
            iterator = [range(bounds[0][i], bounds[1][i]+1) for i in
                        range(0, self.vertex_count)]
            iterator = product(*iterator)
            lms = []
            for x in iterator:
                y = M_inv*Matrix(x).T
                if y in MatrixSpace(ZZ, self.vertex_count, 1):
                    u = tuple(y.column(0))
                    lms.append((u,self.chi(K, u), self.F(K, u)))
            lms.sort(key = lambda x:x[1])
            return lms
        else:
            print("not yet implemented for non-negative definite plumbings")


    def chi_sublevels(self, K, n):
        """
        Given a characteristic vector K and a positive integer n, computes the
        lattice elements in each of the first n non-empty sublevel sets of
        chi_K. Also, computes the value of chi_K and F on each lattice point in
        each sublevel set.

        Parameters
        ----------
        K: list
            A list of integers [k_1, ..., k_s] where s = self.vertex_count.

        n: int
            A positive integer.

        Returns
        -------
        list
            A list of the form [S_1, ..., S_n] where S_i is the ith non-empty
            sublevel set. Each S_i is a set whose elements are tuples of the
            form (a, b, c) where a is a lattice point in S_i, b is the value
            of chi_K on that element, and c is the value of F on that element.

        """
        if self.definiteness_type == "negative definite":
            lms = self.chi_local_min_set(K)

            groups = groupby(lms, operator.itemgetter(1))
            lms_partition = [tuple(group) for key, group in groups]

            min_level = lms_partition[0][0][1]
            sublevels = [set(lms_partition[0])]

            for i in range(1, n):
                sublevel_height = i + min_level
                sublevel_temp1 = copy(sublevels[-1])
                sublevel_temp2 = copy(sublevels[-1])
                for x in sublevel_temp1:
                    for j in range(0, self.vertex_count):
                        y = list(x[0])
                        z = list(x[0])
                        y[j] = y[j]-1
                        z[j] = z[j]+1
                        if self.chi(K, y) == sublevel_height:
                            sublevel_temp2.add((tuple(y), sublevel_height,
                                                self.F(K, tuple(y))))
                        if self.chi(K, z) == sublevel_height:
                            sublevel_temp2.add((tuple(z), sublevel_height,
                                                self.F(K, tuple(z))))
                for u in lms_partition:
                    if u[0][1] == sublevel_height:
                        sublevel_temp2 = sublevel_temp2.union(set(u))
                        break
                sublevels.append(sublevel_temp2)

            return sublevels

        else:
            print("Not yet implemented for non-negative definite plumbings.")

    def chi_sublevels_with_graphs(self, K, n):
        if self.definiteness_type == "negative definite":
            sublevels = self.chi_sublevels(K, n)
            graphs = []
            basis_vectors = []

            big_graph = ig.Graph()
            big_sublevel_list = list(sublevels[-1])
            big_sublevel_list.sort()
            vertices = [list(w[0]) for w in big_sublevel_list]
            v_labels = []
            for x in big_sublevel_list:
                if x[2] == 0:
                    v_labels.append(None)
                else:
                    v_labels.append(x[2])
            num_of_vertices = len(vertices)
            big_graph.add_vertices(int(num_of_vertices))
            big_graph.vs["label"] = v_labels

            edges = []
            for i in range(1, num_of_vertices):
                for j in range(0, self.vertex_count):
                    x = copy(vertices[i])
                    x[j] = x[j] - 1
                    y = copy(vertices[i])
                    y[j] = y[j] + 1
                    if x in vertices[:i]:
                        edges.append((int(vertices[:i].index(x)), int(i)))
                    if y in vertices[:i]:
                        edges.append((int(vertices[:i].index(y)), int(i)))

            big_graph.add_edges(edges)

            graphs = []
            for sl in sublevels[:-1]:
                v_list = [big_sublevel_list.index(x) for x in sl]
                graphs.append(big_graph.subgraph(v_list))

            graphs.append(big_graph)

            plots = [ig.plot(g, bbox = (2000, 2000), vertex_size = int(5),
                     vertex_color = "black", vertex_label_size = 20,
                     vertex_label_color="red", vertex_label_dist = int(2))
                     for g in graphs]

            return sublevels, graphs, plots

        else:
            print("Not yet implemented for non-negative definite plumbings.")

    def zhat(self, a, n):
        a = Matrix(a).T
        M_inv = self.intersection_form.inverse()
        min_level = ((a.T*M_inv*a)[0,0])/4
        bounding_box = []
        for i in range(0, self.vertex_count):
            if min_level % 1 == 0:
                x = 2*math.sqrt(-(self.weight_vector[i, 0])*(n-1))
            else:
                x = 2*math.sqrt(-(self.weight_vector[i, 0])*n)
            bounding_box.append([-x, x])

        F_supp = []
        for i in range(0, self.vertex_count):
            if self.degree_vector[i, 0] == 0:
                if bounding_box[i][0]<= -2 and bounding_box[i][1]>= 2:
                    F_supp.append([-2, 0, 2])
                elif bounding_box[i][0]<= -2 and bounding_box[i][1]< 2:
                    F_supp.append([-2, 0])
                elif bounding_box[i][0]> -2 and bounding_box[i][1]>= 2:
                    F_supp.append([0, 2])
                else:
                    F_supp.append([0])
            elif self.degree_vector[i, 0] == 1:
                if bounding_box[i][0]<= -1 and bounding_box[i][1]>= 1:
                    F_supp.append([-1, 1])
                elif bounding_box[i][0]<= -1 and bounding_box[i][1]<1:
                    F_supp.append([-1])
                elif bounding_box[i][0]> -1 and bounding_box[i][1]>=1:
                    F_supp.append([1])
                else:
                    return 0
            elif self.degree_vector[i, 0] == 2:
                F_supp.append([0])
            else:
                r = self.degree_vector[i, 0]-2
                values = []
                if bounding_box[i][0] <=-r:
                    values.append(-r)
                    for j in range(1, floor((-r-bounding_box[i][0])/2)+1):
                        values.append(-r - 2*j)
                if bounding_box[i][1] >=r:
                    values.append(r)
                    for j in range(1, floor((bounding_box[i][1]-r)/2)+1):
                        values.append(r + 2*j)
                if len(values)==0:
                    return 0
                F_supp.append(copy(values))
        iterator = product(*F_supp)

        def F_hat(y):
            F = 1
            for i in range(0, self.vertex_count):
                if self.degree_vector[i, 0] == 0:
                    if y[i, 0] == 0:
                        F = -2*F
                    elif y[i, 0] != 2 and y[i, 0]!= -2:
                        F = 0
                        return F
                elif self.degree_vector[i, 0] == 1:
                    if y[i, 0] == 1:
                        F = -F
                    elif y[i, 0] != -1:
                        F = 0
                        return F
                elif self.degree_vector[i, 0] == 2:
                    if y[i, 0] != 0:
                        F = 0
                        return F
                else:
                    if abs(y[i, 0]) >= self.degree_vector[i, 0]-2:
                        F = F*(1/2)*sign(y[i,0])^(self.degree_vector[i, 0])
                        F = F*binomial((self.degree_vector[i, 0]
                                        + abs(y[i, 0]))/2-2,
                                        self.degree_vector[i, 0] -3)
                    else:
                        F = 0
                        return F
            return F

        exponents = [ceil(min_level) + i for i in range(0, n)]
        coefficients = n*[0]

        for y in iterator:
            y = Matrix(y).T
            c = -((y.T*M_inv*y)[0,0])/4
            if frac(min_level) == 0:
                if c <= n-1:
                    x = (1/2)*M_inv*(y-a)
                    if x in MatrixSpace(ZZ, self.vertex_count, 1):
                        ind = c
                        coefficients[ind] = coefficients[ind] + F_hat(y)
            else:
                if c <= n:
                    x = (1/2)*M_inv*(y-a)
                    if x in MatrixSpace(ZZ, self.vertex_count, 1):
                        ind = floor(c)
                        coefficients[ind] = coefficients[ind] + F_hat(y)


        u = Matrix(self.vertex_count*[1])
        normalization_term = -min_level -(3*self.vertex_count +
                                         (u*self.weight_vector)[0,0])/4
        zhat_list = [(coefficients[i], exponents[i] + normalization_term)
                     for i in range(0, n)]

        R.<q> = PuiseuxSeriesRing(QQ)
        truncated_zhat = 0
        for x in zhat_list:
            truncated_zhat = truncated_zhat + x[0]*q^(x[1])
        return truncated_zhat
