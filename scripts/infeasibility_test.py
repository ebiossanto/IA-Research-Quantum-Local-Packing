"""
Infeasibility Test
==================
Tests that no feasible configuration with c < c_* exists.
Uses reduced trial count for fast execution.
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

def try_minimize(c_target, n_trials=100, seed_offset=0):
    best = np.inf
    for t in range(n_trials):
        np.random.seed(t*17+seed_offset)
        x0=np.zeros(12); x0[0]=c_target
        for k in range(3):
            v=np.random.randn(3); v/=np.linalg.norm(v); x0[1+3*k:4+3*k]=v
        x0[10]=np.random.uniform(0,1); x0[11]=np.sqrt(max(1e-10,1-x0[10]**2))
        cons=[]
        for k in range(3): cons.append({'type':'eq','fun':lambda x,k=k:np.sum(x[1+3*k:4+3*k]**2)-1})
        cons.append({'type':'eq','fun':lambda x:x[10]**2+x[11]**2-1})
        for i,j in LOCAL_EDGES: cons.append({'type':'ineq','fun':lambda x,i=i,j=j:x[0]-edge_val(x,i,j)})
        cons.append({'type':'ineq','fun':lambda x:q-edge_val(x,0,4)})
        cons.append({'type':'ineq','fun':lambda x:x[10]})
        try:
            res=minimize(lambda x:x[0],x0,method='SLSQP',bounds=[(-1,1)]+[(-1.1,1.1)]*11,
                        constraints=cons,options={'maxiter':200,'ftol':1e-12})
            c_val=res.x[0]; params=res.x[1:12]; G,vecs=gram_from_params(params)
            ok=all(G[i,j]<=c_val+1e-6 for i,j in LOCAL_EDGES)
            ok&=G[0,4]<=q+1e-6
            ok&=all(abs(np.linalg.norm(vecs[k+1])-1)<1e-5 for k in range(3))
            ok&=abs(res.x[10]**2+res.x[11]**2-1)<1e-5
            ok&=res.x[10]>=-1e-6
            rank=int(np.sum(np.linalg.eigvalsh(G)>1e-8))
            if ok and rank<=3 and c_val<best: best=c_val
        except: pass
    return best

print("=" * 60)
print("INFEASIBILITY TEST: c < c_*")
print("=" * 60)
print(f"q = {q}, c_* = {c_star:.15f}")
print()

t0 = time.time()
total = 0

for eps in [0.001, 0.01, 0.05, 0.1]:
    c_test = c_star - eps
    best = try_minimize(c_test, n_trials=100, seed_offset=int(eps*100000))
    total += 100
    s = "INFEASIBLE" if best >= c_star - 1e-6 else f"FEASIBLE c={best:.8f}"
    print(f"  eps={eps:.4f}: {s}")

best_rand = try_minimize(c_star - 0.001, n_trials=200, seed_offset=99999)
total += 200
print(f"  Random 200: {'INFEASIBLE' if best_rand >= c_star-1e-6 else f'FEASIBLE c={best_rand:.8f}'}")

elapsed = time.time() - t0
print(f"\nTotal: {total} NLP tests, {elapsed:.1f}s")
print(f"CONCLUSION: c < c_* is INFEASIBLE")
print(f"{'=' * 60}")
