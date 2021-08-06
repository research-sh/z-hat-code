"""A sage package for analyzing plumbed manifolds.

"""

from copy import copy, deepcopy
from itertools import product, groupby

from sage.graphs.graph_plot import GraphPlot

import math
import sys

class CustomVertex(tuple):
    """A class to allow for non-unique vertex labels for a sage Graph() object.

    A vertex will be specified by a tuple v whose last entry is it's label. The
    subtuple v[:-1] must uniquely specify the vertex.
    """
    def __init__(self, v):
        self.vertex = v

    def __str__(self):
        return str(self.vertex[-1])

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
    intersection_smith_form
    is_intersection_form_non_singular
    is_tree
    definiteness_type
    bad_vertices
    artin_fcycle
    is_weakly_elliptic
    is_rational
    is_almost_rational
    homology

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

            self._sage_formatted_graph_data = {}
            for i in range(0, self._vertex_count):
                adjacent_verts = []
                for x in edges:
                    if x[0] == i:
                        adjacent_verts.append('$v_{' + str(x[1]) + '}\hspace{2} '
                                              + str(vertices_dict[x[1]])
                                              + '$')
                    elif x[1] == i:
                        adjacent_verts.append('$v_{' + str(x[0]) + '}\hspace{2} '
                                              + str(vertices_dict[x[0]])
                                              + '$')
                self._sage_formatted_graph_data['$v_{' + str(i) + '}\hspace{2} '
                                                + str(vertices_dict[i])
                                                + '$']\
                                                = copy(adjacent_verts)

            self._graph = Graph(self._sage_formatted_graph_data)

            self._plot_options = options = {'vertex_color': 'black',
                                            'vertex_size': 20,
                                            'layout': 'tree'}

            self._graph_plot = GraphPlot(self._graph, self._plot_options)

            self._weight_vector = Matrix(list(vertices_dict.values())).T
            self._degree_vector = Matrix(self._graph.degree()).T

            self._intersection_form = None
            self._intersection_smith_form = None
            self._is_intersection_form_non_singular = None
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
                  " input.")

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
            intersection_form = self._graph.adjacency_matrix()
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
    def is_intersection_form_non_singular(self):
        if self._is_intersection_form_non_singular is None:
            d = self.intersection_form.det()
            if d == 0:
                self._is_intersection_form_non_singular = False
            else:
                self._is_intersection_form_non_singular = True
        return self._is_intersection_form_non_singular

    @property
    def is_tree(self):
        """bool: True if the plumbing diagram is a finite tree, False
        otherwise."""
        if self._is_tree is None:
            self._is_tree = self._graph.is_tree()
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
        else:
            self._artin_fcycle = "Not applicable; plumbing is not a negative\
                                  definite tree."
        return self._artin_fcycle

    @property
    def is_weakly_elliptic(self):
        """bool: True if the plumbing is weakly elliptic, False or N/A
        otherwise.
        """
        if self._is_weakly_elliptic is None:
            if self.is_tree and self.definiteness_type == "negative definite":
                k = [-i-2 for i in self._graph.vs["v_weights"]]
                k = Matrix(k)
                m = -(k * self.artin_fcycle[0].T
                      + self.artin_fcycle[0]
                      * self.intersection_form
                      * self.artin_fcycle[0].T)[0,0] / 2
                if m == 0:
                    self._is_weakly_elliptic = True
                    self._is_rational = False
                else:
                    self._is_weakly_elliptic = False
            else:
                self._is_weakly_elliptic = "Not applicable; plumbing is not a\
                                            negative definite tree."
        return self._is_weakly_elliptic

    @property
    def is_rational(self):
        """bool: True if the plumbing is rational, False or N/A otherwise."""
        if self._is_rational is None:
            if self.is_tree and self.definiteness_type == "negative definite":
                k = [-i-2 for i in self._graph.vs["v_weights"]]
                k = Matrix(k)
                m = -(k * self.artin_fcycle[0].T
                      + self.artin_fcycle[0]
                      * self.intersection_form
                      * self.artin_fcycle[0].T)[0,0] / 2
                if m == 1:
                    self._is_rational = True
                    self._is_weakly_elliptic = False
                else:
                    self._is_rational = False
            else:
                self._is_rational = "Not applicable; plumbing is not a negative\
                                     definite tree."
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
                                k = [-j-2 for j in v.values()]
                                k = Matrix(k)
                                m = -(k * plumb.artin_fcycle[0].T
                                      + plumb.artin_fcycle[0]
                                      * self.intersection_form
                                      * plumb.artin_fcycle[0].T)[0,0] / 2
                                if m == 1:
                                    self._is_almost_rational = True
                                    break
                            counter = counter + 1
            else:
                self._is_almost_rational = "Not applicable; plumbing is not a\
                                            negative definite tree."
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
        self._graph_plot.show()

    def is_in_integer_image(self, k):
        """Given a vector k, check if it is in the integer image of the
        intersection form.

        Parameters
        ----------
        k: list
            A list of integers of length = self.vertex_count.

        Returns
        -------
        bool
            True if k is in the integer image of the intersection form, False
            otherwise.

        """
        k = Matrix(k).T
        if self.is_intersection_form_non_singular:
            h = self.intersection_form.inverse()*k
            for x in h.column(0):
                if float(x) % 1 != 0:
                    return False
            return True
        else:
            smith = self.intersection_smith_form
            D = smith[0]
            U = smith[1]
            num_of_pivots = self.vertex_count - D.diagonal().count(0)
            j = U * k
            for i in range(0, num_of_pivots):
                if float(j[i, 0]) % float(D[i, i]) != 0:
                    return False
            for i in range(num_of_pivots, self.vertex_count):
                if float(j[i, 0]) != 0:
                    return False
            return True

    def equiv_spinc_reps(self, k1, k2):
        """Given two characteristic vectors, check if they represent the same
        spinc structure

        Parameters
        ----------
        k1: list
            A list of integers [x_1, ..., x_s] where s = self.vertex_count.

        k2: list
            A list of integers [y_1, ..., y_s] where s = self.vertex_count.

        Returns
        -------
        bool
            True if k1 and k2 represent the same spinc structure on the
            plumbed 3-manifold, False otherwise.

        """

        try:
            k1 = Matrix(k1).T
            k2 = Matrix(k2).T
            for i in range(0, self.vertex_count):
                if (float(k1[i, 0])-float(self.weight_vector[i, 0])) % 2 != 0:
                    raise Exception
                if (float(k2[i, 0])-float(self.weight_vector[i, 0])) % 2 != 0:
                    raise Exception
            k = (1/2)*(k1-k2)
            k = k.column(0)
            return self.is_in_integer_image(k)
        except:
            print("Error: one or more of the inputs are not a characteristic "
                  "vector.")

    def char_vector_properties(self, k):
        """Given a characteristic vector k, compute some basic properties.

        Parameters
        ----------
        k: list
            A list of integers [x_1, ..., x_s] where s = self.vertex_count.

        Returns
        -------
        tuple
            (a,b,c) where: a is a string which says if the associated spinc
            structure on the plumbed 3-manifold is torsion or non-torsion, b is
            the order of the 1st Chern class of the associated spinc structure
            on the plumbed 3-manifold, c is the square of the 1st Chern class of
            the associated spinc structure on the plumbed 4-manifold (in
            other words, c = k^2).

        """
        try:
            k = Matrix(k).T
            for i in range(0, self.vertex_count):
                if (float(k[i, 0])-float(self.weight_vector[i, 0])) % 2 != 0:
                    raise Exception

            if self.is_intersection_form_non_singular:
                h = self.intersection_form.inverse()*k
                denominators_of_h_entries = [x.denominator() for x in
                                             h.column(0)]
                order_of_chern_class = abs(lcm(denominators_of_h_entries))
                square = (k.T * h)[0, 0]
                return "Torsion", order_of_chern_class, square
            else:
                smith = self.intersection_smith_form
                D = smith[0]
                U = smith[1]
                V = smith[2]
                num_of_pivots = self.vertex_count - D.diagonal().count(0)
                j = U * k

                for i in range(num_of_pivots, self.vertex_count):
                    if j[i, 0] != 0:
                        return "Non-Torsion", "N/A", "N/A"

                h = self.vertex_count*[0]

                for i in range(0, num_of_pivots):
                    h[i] = j[i, 0]/D[i, i]

                denoms_of_non_zero_h_entries = [h[i].denominator() for i in
                                                range(0, num_of_pivots)]

                order_of_chern_class = abs(lcm(denoms_of_non_zero_h_entries))
                h = V * Matrix(h).T
                square = (k.T * h)[0,0]
                return "Torsion", order_of_chern_class, square
        except:
            print("Error: input is not a characteristic vector.")

    def chi(self, k, x):
        """
        Given a vector k and a lattice point x (represented as a vector),
        compute chi_k(x) = -1/2(k(x) + (x,x)).

        Parameters
        ----------
        k: list
            A list of integers [a_1, ..., a_s] where s = self.vertex_count.

        x: list
            A list of integers [x_1, ..., x_s] where s = self.vertex_count.

        Returns
        -------
        sage constant
            The value of chi_k(x)

        """

        k = Matrix(k)
        x = Matrix(x).T
        return -(1/2)*(k * x + x.T * self.intersection_form * x)[0,0]

    def chi_min(self, k):
        """
        Given a vector k, computes the minimum of the function chi_k on
        Euclidean space and computes the vector which achieves this minimum.
        Note this vector, in general, need not be integral.

        Parameters
        ----------
        k: list
            A list of integers [a_1, ..., a_s] where s = self.vertex_count.

        Returns
        -------
        tuple
            (a,b) where: a is the minimum value of chi_k over R^s and b is a
            list representing the unique vector which achieves this minimum.

        """
        if self.definiteness_type == "negative definite":
            chi_min = self.char_vector_properties(k)[2]/8
            k = Matrix(k).T
            chi_min_vector = -(1/2) * self.intersection_form.inverse() * k
            return chi_min, list(chi_min_vector.column(0))
        else:
            return "Only implemented for negative definite plumbings."

    def F(self, k, x):
        """
        Given a vector k and a lattice element x, computes F(x).

        Parameters
        ----------
        k: list
            A list of integers [a_1, ..., a_s] where s = self.vertex_count.
        x: list
            A list of integers [x_1, ..., x_s] where s = self.vertex_count.

        Returns
        -------
        int
            The value F(x).

        """
        k = Matrix(k).T
        x = Matrix(x).T
        y = 2*self.intersection_form*x + k - self.weight_vector\
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

    def chi_local_min_bounds(self, k):
        """
        Given a vector k, computes two lists [-chi_k(-e_1), ..., -chi_k(-e_s)]
        and [chi_k(e_1), ..., chi_k(e_s)] where e_i = (0, ..., 0, 1, 0, ..., 0)
        is the ith standard basis vector. For the purpose of this function, see
        the function chi_local_min_set.

        Parameters
        ----------
        k: list
            A list of integers [x_1, ..., x_s] where s = self.vertex_count.

        Returns
        -------
        tuple
            (a,b) where: a = [-chi_k(-e_1), ..., -chi_k(-e_s)] and
            b = [chi_k(e_1), ..., chi_k(e_s)]

        """
        I = Matrix.identity(self.vertex_count)
        negative_I = -I
        positive_basis = [I.row(i) for i in range(0, self.vertex_count)]
        negative_basis = [negative_I.row(i) for i in
                          range(0, self.vertex_count)]
        chi_upper = [self.chi(k, x) for x in positive_basis]
        chi_lower = [-self.chi(k, x) for x in negative_basis]
        return chi_lower, chi_upper

    def chi_local_min_set(self, k):
        """
        Given a vector k, computes the set of lattice points at
        which chi_k achieves a local min, when restricted to the lattice. In
        other words, it computes the lattice points x such that
        chi_k(x) <= chi_k(x +/- e_i) for all i where
        e_i = (0, ..., 0, 1, 0, ..., 0)  is the ith standard basis vector. Note
        chi_k(x +/- e_i) = chi_k(x)+ chi_k(+/- e_i) -/+ (x, e_i). Hence, x is
        in the min set iff -chi_k(-e_i) <= (x, e_i) <= chi_k(e_i) for all i.
        This explains the reason for the helper function chi_local_min_bounds.

        Parameters
        ----------
        k: list
            A list of integers [x_1, ..., x_s] where s = self.vertex_count.
        Returns
        -------
        lists
            Each element of the output list is a tuple (a, b, c, d) where a is an
            element of the local min set, b is chi_k(a), c = F(a),
            d = a dot (weight_vector + degree_vector).
        """
        if self.definiteness_type == "negative definite":
            bounds = self.chi_local_min_bounds(k)
            M_inv = self.intersection_form.inverse()
            iterator = [range(bounds[0][i], bounds[1][i]+1) for i in
                        range(0, self.vertex_count)]
            iterator = product(*iterator)
            lms = []
            for x in iterator:
                y = M_inv*Matrix(x).T
                if y in MatrixSpace(ZZ, self.vertex_count, 1):
                    u = tuple(y.column(0))
                    pairing = (Matrix(u)*(self.weight_vector
                               + self.degree_vector))[0,0]
                    lms.append((u,self.chi(k, u), self.F(k, u), pairing))
            lms.sort(key = lambda x:x[1])
            return lms
        else:
            print("Only implemented for negative definite plumbings")


    def chi_sublevels(self, k, n):
        """
        Given a characteristic vector k and a positive integer n, computes the
        lattice points in each of the first n non-empty sublevel sets of
        chi_k. Also, computes chi_k(x), F(x), and,
        x dot (weight_vector + degree_vector) associated to each lattice
        point x in each sublevel set.

        Parameters
        ----------
        k: list
            A list of integers [x_1, ..., x_s] where s = self.vertex_count.

        n: int
            A positive integer.

        Returns
        -------
        list
            A list of the form [S_1, ..., S_n] where S_i is the ith non-empty
            sublevel set. Each S_i is a set whose elements are tuples of the
            form (a, b, c, d) where a is a lattice point in S_i, b = chi_k(a),
            c = F(a), and d = a dot (weight_vector + degree_vector).

        """
        if self.definiteness_type == "negative definite":
            lms = self.chi_local_min_set(k)

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
                        if self.chi(k, y) == sublevel_height:
                            pairing = (Matrix(y)*(self.weight_vector
                                       + self.degree_vector))[0,0]
                            sublevel_temp2.add((tuple(y), sublevel_height,
                                                self.F(k, tuple(y)), pairing))
                        if self.chi(k, z) == sublevel_height:
                            pairing = (Matrix(z)*(self.weight_vector
                                       + self.degree_vector))[0,0]
                            sublevel_temp2.add((tuple(z), sublevel_height,
                                                self.F(k, tuple(z)), pairing))
                for u in lms_partition:
                    if u[0][1] == sublevel_height:
                        sublevel_temp2 = sublevel_temp2.union(set(u))
                        break
                sublevels.append(sublevel_temp2)

            return sublevels

        else:
            print("Only implemented for negative definite plumbings.")

    def bigraded_root(self, k, n):
        """
        Given a characteristic vector k and a positive integer n, computes
        the first n levels of the bigraded root.
        """
        if self.definiteness_type == "negative definite":
            sublevels = self.chi_sublevels(k, n)

            k_sqaured = self.char_vector_properties(k)[2]

            for element in sublevels[0]:
                break
            min_chi_level = element[1]
            d_inv = -2*(min_chi_level) + self.vertex_count/4\
                    + k_sqaured/4
            normalization_term = -(k_sqaured-3*self.vertex_count
                                   + 2*sum(self.weight_vector)[0])/4 + sum(k)/2\
                                   - sum(self.degree_vector)[0]/4

            top_sublevel = list(sublevels[-1])
            top_sublevel.sort()
            vertices = [list(w[0]) for w in top_sublevel]
            num_of_vertices = len(vertices)

            top_sublevel_graph = Graph(num_of_vertices)

            edges = []
            for i in range(1, num_of_vertices):
                for j in range(0, self.vertex_count):
                    x = copy(vertices[i])
                    x[j] = x[j] - 1
                    y = copy(vertices[j])
                    y[j] = y[j] + 1
                    if x in vertices[:i]:
                        edges.append((vertices[:i].index(x), i))
                    if y in vertices[:i]:
                        edges.append((vertices[:i].index(y), i))

            top_sublevel_graph.add_edges(edges)

            sublevel_graphs = []
            for sl in sublevels[:-1]:
                v_list = [top_sublevel.index(x) for x in sl]
                sublevel_graphs.append(top_sublevel_graph.subgraph(v_list))

            sublevel_graphs.append(top_sublevel_graph)

            connected_components = []
            for graph in sublevel_graphs:
                connected_components.append(graph.connected_components())

            bgr_vertices = []
            bgr_vertex_two_variable_weights = []
            index = 0
            h_index = 0
            v_index = 0
            h_index_dictionary = {}
            Q.<q> = PuiseuxSeriesRing(QQ)
            T.<t> = LaurentPolynomialRing(QQ)
            for x in connected_components:
                for y in x:
                    two_vpoly = 0
                    for z in y:
                        w = top_sublevel[z]
                        two_vpoly = two_vpoly\
                                    + (w[2])*q^(2*(w[1]) + w[3])*t^(w[3])
                    bgr_vertex_two_variable_weights.append(two_vpoly)
                    bgr_vertices.append((h_index, v_index, y[0], '$\hspace{4}'
                                         + str(latex(d_inv + 2*v_index))
                                         + ': P_{'
                                         + str(index) + '}$'))

                    h_index = h_index + 1
                    index = index + 1
                h_index_dictionary[v_index] = h_index
                h_index = 0
                v_index = v_index + 1

            bgr = Graph()
            bgr.add_vertices(bgr_vertices)

            bgr_edges = []
            for x in bgr_vertices:
                if x[1]< v_index-1:
                    for i in range(0, h_index_dictionary[x[1]+1]):
                        if x[2] in connected_components[x[1]+1][i]:
                            ind = sum([h_index_dictionary[i] for
                                          i in range(0,x[1]+1)]) + i
                            bgr_edges.append((x, bgr_vertices[ind]))
                            break
            bgr.add_edges(bgr_edges)

            bgr = Graph([(CustomVertex(a), CustomVertex(b)) for (a, b) in
                         bgr.edges(labels=False)])

            options = {'layout': 'forest',
                       'forest_roots': bgr_vertices[-h_index_dictionary[v_index-1]:],
                       'vertex_color': 'black', 'vertex_size': 20}

            bgr_plot = GraphPlot(bgr, options)

            for i in range(0, index):
                bgr_vertex_two_variable_weights[i] = q^(normalization_term)*bgr_vertex_two_variable_weights[i]

            bgr_plot.show()

            for i in range(0, index):
                print('P_{' + str(i) + '} = '
                               + str(bgr_vertex_two_variable_weights[i]))
                print('\n')

            return bgr_plot, bgr_vertex_two_variable_weights

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