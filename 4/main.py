import pandas as pd
from newton_rhapson import new_rhapson_solve as nr
from secante import secante as sc
from NHO import hno
from sympy import symbols
import time

f='x*exp(-x)-5-(cos(x))/(x)'
x = symbols('x')
tol = 1e-8
max_iter = 100

start_nr = time.time()
result_nr = nr(f, -0.1, tol=tol, max_iter=max_iter)
end_nr = time.time()
time_nr = end_nr - start_nr

start_sc = time.time()
result_sc = sc(f, -0.1, -0.3, tol=tol, max_iter=max_iter)
end_sc = time.time()
time_sc = end_sc - start_sc

start_hno = time.time()
result_hno = hno(f, -0.1,tol=tol, max_iter=max_iter)
end_hno = time.time()
time_hno = end_hno - start_hno


data = {
    "Método": ["Newton-Raphson", "Secante", "HNO"],
    "Resultado Aproximado": [result_nr[0], result_sc[0], result_hno[0]],
    "Error": [result_nr[1], result_sc[1], result_hno[1]],
    "Iteraciones": [result_nr[2], result_sc[2], result_hno[2]],
    "Tiempo de Ejecución (s)": [time_nr, time_sc, time_hno]
}

df = pd.DataFrame(data)
print(df)
#df.to_csv("resultados_metodos_numericos.csv", index=False)