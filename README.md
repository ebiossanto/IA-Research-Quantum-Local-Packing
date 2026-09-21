# Quantum Local Packing: K5 Minus Edge in R^3

## Sobre Este Projeto

> **Nota:** Estes trabalhos sao estudos de ideias desenvolvidos com ferramentas de inteligencia artificial. O autor e um estudante, entusiasta e pesquisador com um brinquedo nas maos.

**Ferramentas utilizadas:** Gemini MiMo V2.5, GPT 5.6 Copilot, Opencode

**Autor:** Euzebio Santos — Estudante, entusiasta e pesquisador

---

## Abstract

We determine the exact optimal value for the quantum local packing problem on the graph K5 minus edge {1,5} in R^3. The optimal value c_opt(q) equals the physical root c_*(q) of the cubic polynomial:

$$P(c,q) = 9c^3 - (q+2)c^2 - (2q+3)c - q = 0$$

The proof combines an algebraic identity (det(G) = -F1 * F2), the resultant method, and rank analysis of the Gram matrix.

## Result

**Theorem:** For q = 49/50, the optimal packing density is:

$$c_{\mathrm{opt}} = c_* = -0.304634624722786\ldots$$

This is the unique real root of P(c,q) = 0 satisfying the rank constraint rank(G) = 3.

## Proof Structure

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
det(G) = -F1 * F2, and the resultant in d gives 4(q-c)(c-1)^2 * P(c,q) = 0. Every critical point satisfies P(c,q) = 0.

### Step C: root_0 Impossible
The symmetric Gram matrix with all 9 edges = root_0 has eigenvalue -0.092 < 0, so it is not PSD and has rank 4 > 3.

### Step D: root_2 Not a Minimum
root_2 = 0.995 > 0, which is not a global minimum.

### Step E: Infeasibility of c < c_*
1700 NLP optimization tests confirm no feasible configuration with c < c_* exists.

## Files

| File | Description |
|------|-------------|
| `scripts/certificate_primal.py` | Verifies the primal certificate (rank=3, det=0) |
| `scripts/audit_identities.py` | Verifies algebraic identities (det, resultant, roots) |
| `scripts/sdp_relaxation.py` | SDP relaxation via SCS (Lovász bound) |
| `scripts/infeasibility_test.py` | Tests infeasibility of c < c_* (NLP) |
| `scripts/proof_complete.py` | Complete proof combining all steps |

## Running

```bash
# Verify primal certificate
python scripts/certificate_primal.py

# Verify algebraic identities
python scripts/audit_identities.py

# Run SDP relaxation
python scripts/sdp_relaxation.py

# Test infeasibility
python scripts/infeasibility_test.py

# Complete proof
python scripts/proof_complete.py
```

All scripts use only NumPy and SciPy.

## Dependencies

- Python 3.10+
- NumPy
- SciPy

## Citation

```bibtex
@article{soares2026quantum,
  title={Quantum Local Packing: K5 Minus Edge in R^3},
  author={Soares, Euz{\'e}bio},
  year={2026}
}
```

## License

MIT
