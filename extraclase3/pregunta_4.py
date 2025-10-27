import numpy as np
import matplotlib.pyplot as plt

def edo2(p, q, r, h, a, b, y0, yn):
    """
    Resuelve la ecuación diferencial y''(x) + p(x)y'(x) + q(x)y(x) = r(x)
    usando diferencias finitas.
    
    Parámetros:
    -----------
    p, q, r : función
        Funciones de una variable real que definen la EDO
    h : float
        Tamaño de paso
    a, b : float
        Extremos del intervalo [a,b]
    y0, yn : float
        Valores de frontera y(a) = y0, y(b) = yn
    
    Retorna:
    --------
    x : np.array
        Vector de puntos x donde se aproxima la solución
    y : np.array
        Vector de valores aproximados de y en los puntos x
    """
    # Generar puntos x
    n = int((b - a) / h) + 1  # n puntos incluyendo extremos
    x = np.linspace(a, b, n)
    
    
    aj = np.zeros(n-2)  # diagonal principal
    bj = np.zeros(n-2)  # diagonal superior
    cj = np.zeros(n-2)  # diagonal inferior
    d = np.zeros(n-2)   # Vector
    
    
    for i in range(n-2):
        xi = x[i+1]  # x_i+1 porque i va de 1 a n-2 en la formulación matemática
        
        # Coeficientes según la discretización de y''(x) + p(x)y'(x) + q(x)y(x) = r(x)
        aj[i] = -2 + h*h*q(xi)        
        bj[i] = 1 + (h/2)*p(xi)       
        cj[i] = 1 - (h/2)*p(xi)       
        d[i] = h*h*r(xi)              
    
    # Modificar primera y última ecuación para incluir condiciones de frontera
    d[0] -= (1 - (h/2)*p(x[1])) * y0      
    d[-1] -= (1 + (h/2)*p(x[-2])) * yn   

   
    y_aux = thomas(aj[1:-1], bj[1:-1], cj[1:-1], d)
    #y = np.zeros(len(x))
    y = []
    y.append(y0)
    y.extend(y_aux)
    y.append(yn)
    
    
    # Necesitamos implementar/usar aquí el método de Thomas
    # para resolver el sistema tridiagonal
    
    return x, y  # Reemplazar None por la solución y


def thomas(a, b, c, d):

    n = len(d)
    p = np.zeros(n)
    q = np.zeros(n)
    x = np.zeros(n)
    
    p[0] = b[0] / a[0]
    q[0] = d[0] / a[0]
    
    for i in range(0, n-2):
        denom = a[i] - c[i] * p[i-1]
        p[i] = b[i] / denom if i < n-1 else 0
        q[i] = (d[i] - c[i] * q[i-1]) / denom
    
    #Sustitución hacia atrás
    x[n-1] = q[n-1]
    for i in range(n-2, -1, -1):
        x[i] = q[i] - p[i] * x[i+1]
        
    return x

# Ejemplo del problema dado
def p_ejemplo(x):
    return -1/x

def q_ejemplo(x):
    return 1/(4*x*x) - 1

def r_ejemplo(x):
    return 0

def y_exacta(x):
    return (np.sin(6-x)/(np.sin(5)*np.sqrt(x)))


h_valores = [1, 0.5, 0.2, 0.1, 0.01]
x_exacta = np.linspace(1, 6, 1000)
y_exacta_vals = y_exacta(x_exacta)

plt.figure(figsize=(10, 6))
plt.plot(x_exacta, y_exacta_vals, 'k-', label='Solución exacta', linewidth=2)

for h in h_valores:
    x, y = edo2(p_ejemplo, q_ejemplo, r_ejemplo, h, 1, 6, 1, 0)
    if y is not None:  # Graficar todas las funciones
        plt.plot(x, y, 'o-', label=f'h = {h}', markersize=4)

plt.title('Solución EDO y comparación con diferentes tamaños de paso')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.legend()
plt.show()
