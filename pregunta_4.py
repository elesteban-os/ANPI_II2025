import pandas as pd
from sympy import symbols
import time

# Importar los métodos desde la carpeta methods4
from methods4.biseccion import metodo_biseccion as bs
from methods4.newton_rhapson import new_rhapson_solve as nr
from methods4.steffensen import metodo_steffensen as stf
from methods4.secante import secante as sc
from methods4.falsa_posicion import metodo_falsa_posicion as fp
from methods4.NHO import hno

# Definición de variables
f='x*exp(-x)-5-(cos(x))/(x)'
x = symbols('x')
tol = 1e-8
max_iter = 10000

# Método de la bisección
start_bs = time.time()
result_bs = bs(f, -0.3, -0.1, tol=tol, max_iter=max_iter)
end_bs = time.time()
time_bs = end_bs - start_bs

# Método Newton-Raphson
start_nr = time.time()
result_nr = nr(f, -0.1, tol=tol, max_iter=max_iter)
end_nr = time.time()
time_nr = end_nr - start_nr

# Método de Steffensen
start_stf = time.time()
result_stf = stf(f, -0.1, tol=tol, max_iter=max_iter)
end_stf = time.time()
time_stf = end_stf - start_stf

# Método Secante
start_sc = time.time()
result_sc = sc(f, -0.1, -0.3, tol=tol, max_iter=max_iter)
end_sc = time.time()
time_sc = end_sc - start_sc

# Método de la Falsa Posición
start_fp = time.time()
result_fp = fp(f, -0.3, -0.1, tol=tol, max_iter=max_iter)
end_fp = time.time()
time_fp = end_fp - start_fp

# Método NHO
start_hno = time.time()
result_hno = hno(f, -0.1,tol=tol, max_iter=max_iter)
end_hno = time.time()
time_hno = end_hno - start_hno

# Representacion tabulada de los resultados
data = {
    "Método": ["Bisección", "Newton-Raphson", "Steffensen", "Secante", "Falsa Posición", "HNO"],
    "Resultado Aproximado (xk)": [result_bs[0], result_nr[0], result_stf[0], result_sc[0], result_fp[0], result_hno[0]],
    "Error (|f(xk)|)": [result_bs[1], result_nr[1], result_stf[1], result_sc[1], result_fp[1], result_hno[1]],
    "Iteraciones (k)": [result_bs[2], result_nr[2], result_stf[2], result_sc[2], result_fp[2], result_hno[2]],
    "Tiempo de Ejecución (s)": [time_bs, time_nr, time_stf, time_sc, time_fp, time_hno]
}

df = pd.DataFrame(data)
print(df)
