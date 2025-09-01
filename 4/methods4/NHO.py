from sympy import symbols, diff, lambdify, sympify
import time
f='x*exp(-x)-5-(cos(x))/(x)'
x = symbols('x')




def hno(f, x0, tol=1e-8, max_iter=10000):
    """ 
    Método de la Hibridación de Newton y la Secante (HNO) para encontrar raíces de una función.

    Parameters:
    f : str
        La función para la cual se busca la raíz, expresada como una cadena.
    x0 : float
        El primer valor inicial para el método.
    tol : float, optional
        La tolerancia para la convergencia (default es 1e-8).
    max_iter : int, optional
        El número máximo de iteraciones permitidas (default es 10000).
    Returns:
    root : float""" 

    f_sym = sympify(f, evaluate=False)
    f_num = lambdify(x, f_sym, 'numpy')
    
   
    f_d = diff(f_sym, x)
    f_d_num = lambdify(x, f_d, 'numpy')

    xn = x0
    for i in range(max_iter):
        z_n = xn - f_num(xn) / f_d_num(xn)
        H_n = H(xn, z_n, f_d_num)
        xn1 = z_n - H_n * f_num(z_n) / f_d_num(z_n)
        error = abs(f_num(xn1))
        if error < tol:
            return xn1, error, i
        xn = xn1

def H(x,z,f_d_num): 
    return (f_d_num(x)-f_d_num(z))/(3*f_d_num(z)-f_d_num(x))


