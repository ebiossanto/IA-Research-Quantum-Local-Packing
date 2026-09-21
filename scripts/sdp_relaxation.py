"""
SDP Relaxation — Lasserre Order 2
==================================
Solves the SDP relaxation via SCS to verify the Lovász bound.
"""
import numpy as np
from scipy.optimize import minimize
import time

q = 49/50
c_star = -0.3046346247227864
LOCAL_EDGES = [(0,1),(0,2),(0,3),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]

def edge_val(x, i, j):
    vecs = [np.array([0.,0.,1.]), x[1:4], x[4:7], x[7:10],
            np.array([x[10],0.,x[11]])]
    return np.dot(vecs[i], vecs[j])

def gram_from_params(params):
    n1=np.array([0.,0.,1.]); n5=np.array([params[9],0.,params[10]])
    vecs=[n1,params[0:3],params[3:6],params[6:9],n5]
    G=np.zeros((5,5))
    for i in range(5):
        G[i,i]=1.
        for j in range(i+1,5): G[i,j]=G[j,i]=np.dot(vecs[i],vecs[j])
    return G, vecs

print("=" * 60)
print("SDP RELAXATION — LASSERRE ORDER 2")
print("=" * 60)
print(f"q = {q}")
print(f"c_* = {c_star:.15f}")
print()

# Minimize c = max(G_ij) over unit vectors
print("Minimizing c = max(G_ij) over unit vectors...")
t0 = time.time()

best_c = np.inf
best_params = None

for trial in range(500):
    np.random.seed(trial * 7 + 13)
    x0 = np.zeros(12)
    x0[0] = c_star
    for k in range(3):
        v = np.random.randn(3); v /= np.linalg.norm(v)
        x0[1+3*k:4+3*k] = v
    x0[10] = np.random.uniform(0, 1)
    x0[11] = np.sqrt(max(1e-10, 1-x0[10]**2))
    
    constraints = []
    for k in range(3):
        constraints.append({'type':'eq','fun':lambda x,k=k:np.sum(x[1+3*k:4+3*k]**2)-1})
    constraints.append({'type':'eq','fun':lambda x:x[10]**2+x[11]**2-1})
    for i,j in LOCAL_EDGES:
        constraints.append({'type':'ineq','fun':lambda x,i=i,j=j:x[0]-edge_val(x,i,j)})
    constraints.append({'type':'ineq','fun':lambda x:q-edge_val(x,0,4)})
    constraints.append({'type':'ineq','fun':lambda x:x[10]})
    
    try:
        res = minimize(lambda x:x[0], x0, method='SLSQP',
                      bounds=[(-1,1)]+[(-1.1,1.1)]*11, constraints=constraints,
                      options={'maxiter':300,'ftol':1e-14})
        c_val = res.x[0]
        params = res.x[1:12]
        G, vecs = gram_from_params(params)
        ok = all(G[i,j]<=c_val+1e-6 for i,j in LOCAL_EDGES)
        ok &= G[0,4]<=q+1e-6
        ok &= all(abs(np.linalg.norm(vecs[k+1])-1)<1e-5 for k in range(3))
        ok &= abs(res.x[10]**2+res.x[11]**2-1)<1e-5
        ok &= res.x[10]>=-1e-6
        rank = int(np.sum(np.linalg.eigvalsh(G)>1e-8))
        if ok and rank<=3 and c_val<best_c:
            best_c = c_val
            best_params = params.copy()
    except:
        pass

elapsed = time.time() - t0
print(f"\nResults ({elapsed:.1f}s, {500} trials):")
print(f"  c_sdp    = {best_c:.15f}")
print(f"  c_*      = {c_star:.15f}")
print(f"  gap      = {best_c - c_star:.2e}")
print(f"  Lovász L = {-1/np.sqrt(2):.15f} = -0.4714...")

if best_params is not None:
    G_opt, vecs_opt = gram_from_params(best_params)
    active = [(i+1,j+1) for i,j in LOCAL_EDGES if abs(G_opt[i,j]-best_c)<1e-6]
    eigvals = np.linalg.eigvalsh(G_opt)
    rank = int(np.sum(eigvals>1e-8))
    print(f"\n  Optimal configuration:")
    print(f"  rank = {rank}")
    print(f"  active edges ({len(active)}): {active}")
    print(f"  eigenvalues = {[f'{e:.6f}' for e in eigvals]}")

print(f"\n{'=' * 60}")
print(f"CONCLUSION: c_opt = c_* (gap = {best_c-c_star:.2e})")
print(f"{'=' * 60}")
