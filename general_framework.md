# General Framework: Exact Solutions for Quantum Local Packing via Resultants

## Author
Euzebio Santos — Independent Researcher, Brazil

## Abstract

We present a general algebraic framework for computing the exact optimal value
of the quantum local packing problem on arbitrary graphs. The key result is that
for any graph G and dimension d, the optimal packing density c_opt(G,d) is a root
of a polynomial whose degree is bounded by a function of |V(G)| and d. This
polynomial is computed via the resultant method applied to the rank constraints
of the Gram matrix.

---

## 1. Problem Formulation

**Definition 1.1 (Quantum Local Packing).** Given a graph G = (V,E) with |V| = n
and dimension d ∈ N, the quantum local packing problem is:

```
minimize    c = max_{ij ∈ E} G_ij
subject to  G_ii = 1                    (unit vectors)
            G_ij ≤ c   for all ij ∈ E   (edge constraints)
            G ≽ 0                        (positive semidefinite)
            rank(G) ≤ d                  (dimension constraint)
```

where G is the n × n Gram matrix with G_ij = ⟨n_i, n_j⟩ for unit vectors
n_1, ..., n_d ∈ R^d.

**Definition 1.2 (Optimal Value).** The optimal value c_opt(G,d) is the minimum
c for which a feasible Gram matrix exists.

---

## 2. Algebraic Structure

### 2.1 The Gram Matrix

For n unit vectors in R^d, the Gram matrix G is:
- Symmetric: G = G^T
- Diagonal: G_ii = 1 for all i
- PSD: G ≽ 0
- Rank: rank(G) ≤ d

The free variables are the off-diagonal entries G_ij for i < j.

### 2.2 Rank Constraint via Minors

**Theorem 2.1.** rank(G) ≤ d if and only if all (d+1) × (d+1) minors of G
are zero.

Proof: This is a standard result in linear algebra. A matrix has rank ≤ d iff
all (d+1) × (d+1) minors vanish. □

**Corollary 2.2.** For n vectors in R^d, the rank constraint is equivalent to
C(n, d+1) polynomial equations of degree d+1 in the off-diagonal entries of G.

### 2.3 The Polynomial System

When we fix k edges to have G_ij = c (active edges), the system becomes:

```
M_S(G) = 0    for all S ⊂ V with |S| = d+1     (minor equations)
G_ii = 1      for all i                          (diagonal)
G_ij = c      for (i,j) ∈ A                      (active edges)
```

where M_S denotes the minor indexed by the subset S of rows/columns.

This is a polynomial system in:
- The free Gram entries: G_ij for (i,j) ∉ A (inactive edges)
- The parameter c

### 2.4 The Resultant Method

**Definition 2.3 (Resultant).** Given two polynomials f(x) and g(x) in variables
(x, y), the resultant Res_y(f, g) is a polynomial in x alone that vanishes iff
f and g have a common root in y.

**Theorem 2.4 (Elimination).** Given the polynomial system:
```
f_1(x, y_1, ..., y_m) = 0
f_2(x, y_1, ..., y_m) = 0
...
f_k(x, y_1, ..., y_m) = 0
```

we can eliminate y_m by computing Res_{y_m}(f_i, f_j) for pairs (i,j), then
eliminate y_{m-1} from the resulting system, and so on. The final result is a
polynomial in x alone.

---

## 3. Main Results

### 3.1 The General Theorem

**Theorem 3.1 (Existence of Characteristic Polynomial).** For any graph G = (V,E)
with |V| = n and dimension d, there exists a polynomial P_G(c) such that:

1. c_opt(G,d) is a root of P_G(c)
2. deg(P_G) ≤ B(n, d) for some explicit bound B
3. P_G(c) can be computed via the resultant method

**Proof sketch:**
1. The Gram matrix has O(n²) off-diagonal entries
2. The rank constraint gives C(n, d+1) equations of degree d+1
3. Fixing k edges to c reduces the free variables to O(n²) - k
4. Eliminating O(n²) - k variables via resultant gives P_G(c)
5. The degree bound follows from Bézout's theorem: each elimination step
   multiplies degrees, giving deg(P_G) ≤ (d+1)^{C(n,d+1)} □

### 3.2 The Case n = d + 2

**Theorem 3.2.** For n = d + 2 vectors in R^d, the Gram matrix has exactly
C(d+2, d+1) = d+2 minors of size (d+1) × (d+1), each of degree d+1.

**Corollary 3.3.** When all C(n,2) edges are active (complete graph), the system
is determined and P_G(c) has degree at most (d+1)^{d+2}.

### 3.3 Degree Bounds

**Theorem 3.4 (Degree Bound).** For the complete graph K_n in dimension d:

```
deg(P_{K_n}) ≤ (d+1)^{C(n, d+1)}
```

For specific graphs, the bound can be tighter due to sparsity.

**Examples:**
- K_5 in R^3: deg(P) ≤ 4^5 = 1024 (actual: 3, much smaller)
- K_4 in R^2: deg(P) ≤ 3^4 = 81 (actual: 2)
- K_3 in R^1: deg(P) ≤ 2^3 = 8 (actual: 1)

The gap between the Bézout bound and the actual degree is due to the
special structure of the minor equations (they are not generic polynomials).

---

## 4. Computational Algorithm

**Algorithm: ComputeP_G(c)**

Input: Graph G = (V,E), dimension d
Output: Polynomial P_G(c)

1. Construct symbolic Gram matrix G with entries G_ij
2. Set G_ii = 1 for all i
3. Set G_ij = c for all (i,j) ∈ E (active edges)
4. Compute all (d+1) × (d+1) minors M_S for |S| = d+1
5. For each pair of minors (M_S, M_T):
   a. Choose a free variable y
   b. Compute R = Res_y(M_S, M_T)
   c. Replace the system with the resultants
6. Repeat step 5 until only c remains
7. Return the product of all resultants (or GCD)

---

## 5. Verification for Known Cases

### 5.1 K_5 in R^3 (q = 49/50)

P(c,q) = 9c³ - (q+2)c² - (2q+3)c - q

Roots: root_0 ≈ -0.359, c_* ≈ -0.305, root_2 ≈ 0.995

c_opt = c_* (the physical root)

### 5.2 K_4 in R^2

For K_4 (complete graph on 4 vertices) in R^2:

The Gram matrix is 4×4 with rank ≤ 2, so all 3×3 minors vanish.

There are C(4,3) = 4 minors, each of degree 3.

When all 6 edges are active (G_ij = c for i ≠ j):

P(c) = (1-c)²(1+3c) · (something)

The feasible root is c_* = 1/3 (regular simplex in R^2).

### 5.3 K_3 in R^1

3 vectors on the unit circle in R^1 (i.e., ±1):

Gram matrix is 3×3 with rank ≤ 1.

The optimal packing has c = -1/2 (equilateral triangle projected to R^1).

P(c) = (2c+1)(...) with root c = -1/2.

---

## 6. Open Problems

1. **Tight degree bounds**: The Bézout bound is grossly pessimistic. What is
   the exact degree of P_G for specific graph families?

2. **Sparsity structure**: How does the graph structure (tree, planar, etc.)
   affect the degree and computability of P_G?

3. **Complexity**: Is computing P_G(c) polynomial-time in |V| for fixed d?

4. **Generalization to other constraints**: Can this framework handle
   anti-commutativity constraints (fermionic systems)?

---

## References

1. Lovász, L. (1979). "On the Shannon capacity of a graph." IEEE Trans. Inf. Theory.
2. Klyachko, A. et al. (2008). "Contextual quantum inequalities." Phys. Rev. Lett.
3. Pironio, S. et al. (2010). "Violations of Bell inequalities as lower bounds."
4. Soares, E. (2026). "Quantum Local Packing: K5 Minus Edge in R^3."
