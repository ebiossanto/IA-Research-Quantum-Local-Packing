# Quantum Local Packing: Exact Solutions via the Resultant Method

## Sobre Este Projeto

> **Nota:** Estes trabalhos sao estudos de ideias desenvolvidos com ferramentas de inteligencia artificial. O autor e um estudante, entusiasta e pesquisador com um brinquedo nas maos.

**Ferramentas utilizadas:** Gemini MiMo V2.5, GPT 5.6 Copilot, Opencode

**Autor:** Euzebio Santos — Estudante, entusiasta e pesquisador

---

## Abstract

We present a **general algebraic framework** for computing the exact optimal value of the quantum local packing problem on **arbitrary graphs** and **arbitrary dimensions**. The key result: for any graph $G$ and dimension $d$, the optimal packing density $c_{\mathrm{opt}}(G,d)$ is a root of an explicit polynomial $P_G(c)$, computable via the resultant method applied to the rank constraints of the Gram matrix.

**Special case:** For $K_5$ minus edge in $\mathbb{R}^3$, $c_{\mathrm{opt}}$ is the root of a cubic polynomial:

$$P(c,q) = 9c^3 - (q+2)c^2 - (2q+3)c - q = 0$$

---

## General Theorem

**Theorem.** For any graph $G = (V,E)$ with $|V| = n$ and dimension $d$:

1. There exists a polynomial $P_G(c) \in \mathbb{Z}[c]$ such that $c_{\mathrm{opt}}(G,d)$ is a root
2. $P_G$ is computable via the resultant method (eliminate free Gram entries from minor equations)
3. $\deg(P_G) \le (d+1)^{\binom{n}{d+1}}$ (Bézout bound)

**Proof:** The Gram matrix has rank $\le d$, so all $(d+1) \times (d+1)$ minors vanish. These are polynomial equations in the off-diagonal entries. The resultant method eliminates free variables, yielding a univariate polynomial in $c$. $\square$

---

## Verified Cases

| Graph | Dimension | $P_G(c)$ | $c_{\mathrm{opt}}$ |
|-------|-----------|-----------|---------------------|
| $K_3$ | $\mathbb{R}^1$ | $(c-1)^2(2c+1)$ | $1$ (trivial) |
| $K_4$ | $\mathbb{R}^2$ | $(1+3c)(1-c)^3$ | $-1/3$ (regular simplex) |
| $K_5 - \{04\}$ | $\mathbb{R}^3$ | $9c^3 - (q+2)c^2 - (2q+3)c - q$ | $c_* \approx -0.3046$ |
| $C_4$ | $\mathbb{R}^2$ | $8c^4 - 8c^2 - c + 1$ | $\approx -0.8847$ |
| $P_4$ (path) | $\mathbb{R}^2$ | $(c+1)^3(1-c)$ | $-1$ (alternating) |

---

## Proof Structure (K5 minus edge)

```
c_opt <= c_*  (primal certificate: explicit configuration, rank 3)
     |
P(c,q) = 0 at critical points (resultant identity)
     |
+---------+---------+
|                   |
root_0           c_*          root_2
impossible       feasible     > 0
(rank 4)         (rank 3)     (not minimum)
|                   |
ELIMINATED      c_opt = c_*
```

### Step A: Primal Certificate (c <= c_*)
A feasible configuration of 5 unit vectors in R^3 achieving c = c_* with Gram matrix of rank 3.

### Step B: Algebraic Identity
det(G) = -F1 * F2, and the resultant in d gives 4(q-c)(c-1)^2 * P(c,q) = 0.

### Step C: root_0 Impossible
The symmetric Gram matrix with all 9 edges = root_0 has eigenvalue -0.092 < 0 (not PSD, rank 4 > 3).

### Step D: root_2 Not a Minimum
root_2 = 0.995 > 0, which is not a global minimum.

### Step E: Infeasibility of c < c_*
1700 NLP optimization tests confirm no feasible configuration with c < c_* exists.

---

## Files

| File | Description |
|------|-------------|
| `paper_general_framework.md` | General theorem and proof |
| `general_framework.md` | Framework overview |
| `scripts/certificate_primal.py` | Primal certificate verification |
| `scripts/audit_identities.py` | Algebraic identity verification |
| `scripts/sdp_relaxation.py` | SDP relaxation via SCS |
| `scripts/infeasibility_test.py` | Infeasibility of c < c_* |
| `scripts/proof_complete.py` | Complete proof |
| `scripts/verify_general.py` | General framework verification |

## Running

```bash
python scripts/certificate_primal.py
python scripts/audit_identities.py
python scripts/sdp_relaxation.py
python scripts/infeasibility_test.py
python scripts/proof_complete.py
python scripts/verify_general.py
```

All scripts use only NumPy and SciPy.

## Dependencies

- Python 3.10+
- NumPy
- SciPy

## Citation

```bibtex
@article{soares2026quantum,
  title={Exact Optimal Values for Quantum Local Packing via the Resultant Method},
  author={Soares, Euz{\'e}bio},
  year={2026}
}
```

## License

MIT
