import numpy as np
import sympy as sp
from math import ceil

def cota_simpson_puntos(f, a, b, tol, samples=10001):
    """
    Devuelve (n, alpha_max) donde n es el número mínimo de subintervalos (par)
    para que la cota de error de la regla compuesta de Simpson sea < tol,
    y alpha_max = max_{x in [a,b]} |f^(4)(x)|.
    f: expresión de sympy o string que puede parsearse con sympy.sympify
    a, b: intervalos (float)
    tol: tolerancia deseada (float)
    samples: número de puntos para estimar el máximo de |f^(4)(x)|
    """
    
    if isinstance(f, str):
        x = sp.symbols('x')
        f_expr = sp.sympify(f)
    elif isinstance(f, sp.Expr):
        f_expr = f
        syms = list(f_expr.free_symbols)
        if len(syms) == 0:
            x = sp.symbols('x')
        else:
            x = syms[0]
    else:
        raise ValueError("f debe ser una expresión sympy o string representando la función en x")

    # cuarta derivada
    f4 = sp.diff(f_expr, x, 4)
    f4_abs = sp.Abs(f4)

    # lambdify para evaluar numéricamente
    f4_num = sp.lambdify(x, f4_abs, 'numpy')

    xs = np.linspace(a, b, samples)
    vals = f4_num(xs)
    vals = np.array(vals, dtype=float)
    # manejar posibles NaN/inf
    vals = np.nan_to_num(vals, nan=0.0, posinf=np.max(vals[np.isfinite(vals)]) if np.any(np.isfinite(vals)) else 0.0)

    alpha_max = float(np.max(vals))

    # fórmula: ((b-a)^5 * alpha_max) / (2880 * n^4) < tol  => n > [((b-a)^5 * alpha_max)/(2880*tol)]^(1/4)
    numer = (b - a) ** 5 * alpha_max
    if tol <= 0:
        raise ValueError("tol debe ser > 0")
    n_req = (numer / (2880.0 * tol)) ** 0.25
    n = int(ceil(n_req))
    # Simpson compuesto requiere n par
    if n % 2 == 1:
        n += 1
    return n, alpha_max

if __name__ == "__main__":
    # Ejemplo del enunciado: f(x) = exp(x)*(26 - 10*x + x**2) en [5, 5.5], tol = 1e-8
    x = sp.symbols('x')
    f = sp.exp(x) * (26 - 10*x + x**2)
    a, b, tol = 5.0, 5.5, 1e-8
    n, alpha = cota_simpson_puntos(f, a, b, tol)
    print("n requerido:", n)
    print("alpha_max estimado:", alpha)