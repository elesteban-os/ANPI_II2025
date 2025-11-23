import numpy as np
import sympy as sp
from sympy import symbols, expand, simplify
import matplotlib.pyplot as plt
import math


def newton_simbolico_expandido(x, y, n):
    """
    Calcula el polinomio interpolador usando diferencias divididas de Newton.
    Retorna la expresión simbólica del polinomio.
    
    Parámetros:
    -----------
    x : list o np.ndarray
        Vector de puntos x (n+1 puntos)
    y : list o np.ndarray
        Vector de valores y correspondientes (n+1 valores)
    n : int
        Grado del polinomio (número de puntos - 1)
    
    Retorna:
    --------
    tabla : np.ndarray
        Matriz de diferencias divididas
    polinomio_simbolico : sympy.Expr
        Expresión simbólica del polinomio interpolador
    polinomio_numerico : Callable
        Función que evalúa el polinomio en cualquier punto
    """
    
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)
    
    # Paso 1: Calcular diferencias divididas simbólicamente
    # Declarar tabla de (n+1) x (n+1)
    tabla = np.zeros((n + 1, n + 1))
    
    # Inicializar primera columna con valores y
    for i in range(n + 1):
        tabla[i, 0] = y[i]
    
    # Calcular diferencias divididas recursivamente
    for j in range(1, n + 1):
        for i in range(n - j + 1):
            tabla[i, j] = (tabla[i + 1, j - 1] - tabla[i, j - 1]) / (x[i + j] - x[i])
    
    # Paso 2: Construir polinomio simbólicamente
    X = symbols('X')
    polinomio_simbolico = tabla[0, 0]
    
    # Primer término: f[x0]
    termino_acumulativo = 1
    
    # Términos restantes
    for k in range(1, n + 1):
        termino_acumulativo = termino_acumulativo * (X - x[k - 1])
        polinomio_simbolico = polinomio_simbolico + tabla[0, k] * termino_acumulativo
    
    # Expandir el polinomio
    polinomio_expandido = expand(polinomio_simbolico)
    
    # Función numérica para evaluación rápida
    def evaluar_polinomio(x_eval):
        """Evalúa el polinomio en el punto x_eval"""
        resultado = tabla[0, 0]
        termino = 1
        
        for k in range(1, n + 1):
            termino = termino * (x_eval - x[k - 1])
            resultado = resultado + tabla[0, k] * termino
        
        return resultado
    
    return tabla, polinomio_expandido, evaluar_polinomio


def mostrar_tabla_diferencias(tabla, x, y):
    """
    Muestra la tabla de diferencias divididas de forma legible.
    """
    n = len(x) - 1
    print("\n" + "=" * 100)
    print("TABLA DE DIFERENCIAS DIVIDIDAS")
    print("=" * 100)
    print(f"\n{'i':>3} {'x[i]':>12} {'y[i]':>15}", end="")
    
    for j in range(1, n + 1):
        print(f" {'D^' + str(j):>18}", end="")
    print("\n" + "-" * 100)
    
    for i in range(n + 1):
        print(f"{i:3d} {x[i]:12.6f} {tabla[i, 0]:15.8f}", end="")
        for j in range(1, min(i + 2, n + 1)):
            print(f" {tabla[i, j]:18.8f}", end="")
        print()
    
    print("=" * 100 + "\n")


def main():
    """
    Ejemplo principal con los puntos dados: (1,2), (2,5), (3,10), (4,17), (5,26)
    """
    
    # Definir puntos de interpolación
    #puntos = [(1, 2), (2, 5), (3, 10), (4, 17), (5, 26)]
    # f(x) = log_x(arcsin(x)) = ln(arcsin(x)) / ln(x)
    f = lambda t: np.log(np.arcsin(t)) / np.log(t)
    puntos = [(0.1, f(0.1)), (0.2, f(0.2)), (0.3, f(0.3)), (0.4, f(0.4)), (0.5, f(0.5)), (0.6, f(0.6)), (0.7, f(0.7)), (0.8, f(0.8))]
    print(f"Puntos de interpolación: {puntos}\n")
    
    
    
    x = np.array([p[0] for p in puntos], dtype=float)
    y = np.array([p[1] for p in puntos], dtype=float)
    
    n = len(x) - 1
    
    print("\n" + "=" * 100)
    print("INTERPOLACIÓN POR DIFERENCIAS DIVIDIDAS DE NEWTON")
    print("=" * 100)
    print(f"\nPuntos de interpolación (n = {n}):")
    print("puntos ← [(0.1, f(0.1)), (0.2, f(0.2)), (0.3, f(0.3)), (0.4, f(0.4)), (0.5, f(0.5)), (0.6, f(0.6)), (0.7, f(0.7)), (0.8, f(0.8))]")
    print(f"\nx = {x}")
    print(f"y = {y}\n")
    
    # Calcular interpolación
    tabla, polinomio_simbolico, polinomio_numerico = newton_simbolico_expandido(x, y, n)
    
    # Mostrar tabla de diferencias
    mostrar_tabla_diferencias(tabla, x, y)
    
    # Mostrar coeficientes del polinomio
    print("Coeficientes en forma de Newton:")
    print("P(X) = a₀ + a₁·(X-x₀) + a₂·(X-x₀)(X-x₁) + ... + aₙ·(X-x₀)...(X-xₙ₋₁)")
    print("\nDonde:")
    for k in range(n + 1):
        if k == 0:
            print(f"a₀ = {tabla[0, k]:.8f}")
        else:
            termino_desc = "·".join([f"(X-{x[j]:.0f})" for j in range(k)])
            print(f"a₍{k}₎ = {tabla[0, k]:.8f}  ← {termino_desc}")
    
    # Mostrar polinomio expandido
    print("\n" + "=" * 100)
    print("POLINOMIO INTERPOLADOR EXPANDIDO")
    print("=" * 100)
    print(f"\nP(X) = {polinomio_simbolico}\n")
    
    
    # Verificar en los puntos de interpolación
    print("=" * 100)
    print("VERIFICACIÓN EN PUNTOS DE INTERPOLACIÓN")
    print("=" * 100)
    print(f"\n{'x':>12} {'y real':>15} {'P(x)':>15} {'f(x)':>15} {'Error':>15}")
    print("-" * 75)
    
    for i, x_val in enumerate(x):
        p_val = polinomio_numerico(x_val)
        expected = f(x_val)
        error = abs(y[i] - p_val)
        print(f"{x_val:12.4f} {y[i]:15.8f} {p_val:15.8f} {expected:15.8f} {error:15.2e}")
    
    # Evaluar en puntos intermedios
    print("\n" + "=" * 100)
    print("EVALUACIÓN EN PUNTOS INTERMEDIOS (INTERVALO [0.1,0.8])")
    print("=" * 100)
    print(f"\n{'x':>12} {'P(x)':>15} {'f(x)':>15} {'Diferencia':>15}")
    print("-" * 60)
    
    x_eval = np.array([0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75])
    for x_val in x_eval:
        p_val = polinomio_numerico(x_val)
        expected = f(x_val)
        diff = abs(p_val - expected)
        print(f"{x_val:12.4f} {p_val:15.8f} {expected:15.8f} {diff:15.2e}")
    
    # Graficar
    print("\n" + "=" * 100)
    print("GENERANDO GRÁFICA...")
    print("=" * 100 + "\n")
    
    # Graficar en el intervalo [0.1, 0.8]
    x_plot = np.linspace(0.1, 0.8, 300)
    y_plot = np.array([polinomio_numerico(xi) for xi in x_plot])
    y_real = f(x_plot)
    
    plt.figure(figsize=(12, 7))
    
    # Graficar función real f(x) = log_x(arcsin(x))
    plt.plot(x_plot, y_real, 'b-', linewidth=2, label=r'$f(x)=\log_{x}(\arcsin(x))$ (Real)')
    
    # Graficar polinomio interpolador
    plt.plot(x_plot, y_plot, 'r--', linewidth=2, label='$P(x)$ (Interpolador Newton)')
    
    # Puntos de interpolación
    plt.plot(x, y, 'go', markersize=10, label='Puntos de interpolación', zorder=5)
    
    # Configurar gráfica
    plt.xlabel('x', fontsize=12)
    plt.ylabel('y', fontsize=12)
    plt.title('Interpolación por Diferencias Divididas de Newton', fontsize=14, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.axhline(y=0, color='k', linewidth=0.5)
    plt.axvline(x=0, color='k', linewidth=0.5)
    
    # Anotaciones de los puntos
    for xi, yi in zip(x, y):
        plt.annotate(f'({xi:.2f}, {yi:.4f})', 
                    xy=(xi, yi), 
                    xytext=(5, 5), 
                    textcoords='offset points',
                    fontsize=9,
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))

    plt.show()


if __name__ == "__main__":
    main()
