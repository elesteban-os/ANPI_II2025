import numpy as np 



def qr_gram_schmidt(A):
    
    A = np.array(A, dtype=float)
    n, m = A.shape

    # e_k serán las columnas de Q
    E = np.zeros((n, m))
    R = np.zeros((m, m))

    # Primer vector
    u1 = A[:, 0]
    e1 = u1 / np.linalg.norm(u1)
    E[:, 0] = e1
    R[0, 0] = np.linalg.norm(u1)

    # Resto de columnas
    for k in range(1, m):
        a_k = A[:, k]
        # u_k = a_k - sum_{j=1}^{k} <a_k, e_j> e_j
        u_k = a_k.copy()
        for j in range(k):
            R[j, k] = np.dot(a_k, E[:, j])      # <a_k, e_j>
            u_k = u_k - R[j, k] * E[:, j]       # resta proyecciones

        R[k, k] = np.linalg.norm(u_k)
        if R[k, k] == 0:
            raise ValueError("Columnas linealmente dependientes")
        e_k = u_k / R[k, k]
        E[:, k] = e_k

    Q = E
    return Q, R






A = np.array([[1,1],
              [1,-1]], dtype=float)

Q, R = qr_gram_schmidt(A)
print(R)
print(Q)