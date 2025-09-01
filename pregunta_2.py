import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

#---Definir variable simbolica---
x = sp.symbols('x')

#---Definir funcion como texto---
f_text = '(x^3-3*x^2+3*x-1)/(x^2-2*x)'
#Convertir texto a simbolico
f_sym = sp.sympify(f_text)
#Convertir simbolico a numerico
f_num = sp.lambdify(x, f_sym, 'numpy')

#---Dominio de la funcion---
num, den = sp.fraction(f_sym)
restriction = sp.Eq(den, 0)
excluded_values = sp.solve(restriction, x)
Real = sp.S.Reals
print("Dominio:",str(Real) + " - " + str(sp.FiniteSet(*excluded_values)) )

#---Interseccion entre eje x y y---

try: # Intersección con eje Y
    y_intersection = f_sym(0)
    print("Intersección con eje Y:", y_intersection)
except:
    print("No interseca con el eje Y (función no definida en x=0)")

try: # Intersección con eje X
    x_intersection = sp.solve(f_sym, x)
    print("Interseccion eje x:", x_intersection)
except:
    print("No interseca con el eje X (función no definida en x=0)")


#---Asintotas--- 
##Verticales 
asymptotes_verticals = sp.solve(sp.denom(f_sym), x)
print("Asintotas verticales:", asymptotes_verticals)
##Horizontales
asymptotes_horizontal = sp.limit(f_sym, x, sp.oo)
print("Asintotas horizontales:", asymptotes_horizontal)
##Oblicuas
asymptotes_oblique = sp.limit(f_sym - x, x, sp.oo)
print("Asintotas oblicuas:", asymptotes_oblique)

#---Derivadas---
f_d = sp.diff(f_sym,x)
print("Derivada de la función:", f_d)

f_d_II = sp.diff(f_d,x)
print("Segunda derivada de la función:", f_d_II)

#---Graficar---
f_d_num = sp.lambdify(x, f_d, 'numpy')
f_d_II_num = sp.lambdify(x, f_d_II, 'numpy')

Index = 1000
x_left = np.linspace(-10, -0.01, Index) #Sección izquierda
x_mid= np.linspace(0.01, 1.99, Index) #Sección media
x_right = np.linspace(2.01, 10, Index) #Sección derecha

x_vals_left = np.unique(x_left) #Sección izquierda
y_vals_left = f_num(x_vals_left)

x_vals_mid = np.unique(x_mid) #Sección media
y_vals_mid = f_num(x_vals_mid)

x_vals_right = np.unique(x_right) #Sección derecha
y_vals_right = f_num(x_vals_right)

fig, axs = plt.subplots(3, 1, figsize=(24, 12))

#plt.figure(1, figsize=(12, 6))
axs[0].plot(x_vals_left, y_vals_left, color = 'green', label='f(x) = (x^3 - 3*x^2 + 3*x - 1)/(x^2 - 2*x)')
axs[0].plot(x_vals_mid, y_vals_mid,  color = 'green')
axs[0].plot(x_vals_right, y_vals_right,  color = 'green')
axs[0].axhline(0, color='black', linewidth=0.5)
axs[0].axvline(0, color='black', linewidth=0.5)
axs[0].grid()
axs[0].legend()


x_f = np.linspace(-10, 10, 1000)

regiones = [
    np.linspace(-10, -0.1, 500),    # Región izquierda
    np.linspace(0.1, 1.9, 500),     # Región central  
    np.linspace(2.1, 10, 500)       # Región derecha
]


for region in regiones:
    y_f_d = f_d_num(region)
    axs[1].plot(region, y_f_d, color='blue', label='f\'(x)' if region is regiones[0] else "")
axs[1].axhline(0, color='black', linewidth=0.5)
axs[1].axvline(0, color='black', linewidth=0.5)
axs[1].grid()
axs[1].legend()

for region in regiones:
    y_f_d_II = f_d_II_num(region)
    axs[2].plot(region, y_f_d_II, color='red', label='f \'\'(x)' if region is regiones[0] else "")
axs[2].axhline(0, color='black', linewidth=0.5)
axs[2].axvline(0, color='black', linewidth=0.5)
axs[2].grid()
axs[2].legend()

plt.tight_layout()
plt.show()


#---Intervalos de crecimiento y decrecimiento de funcion---
# Puntos críticos (donde f'(x) = 0)
critical_points = sp.solve(sp.Eq(f_d, 0), x)

# Puntos de discontinuidad de f'(x)
num_f, den_f = sp.fraction(f_d)
discontinuities_f = sp.solve(sp.Eq(den_f, 0), x)

# Todos los puntos importantes (ordenados y únicos)
all_points_f = sorted(set(critical_points + discontinuities_f))

# Construir intervalos CORRECTAMENTE
intervals_f = []
n = len(all_points_f)

# Intervalo desde -∞ hasta el primer punto
intervals_f.append(sp.Interval(-sp.oo, all_points_f[0]))

# Intervalos entre puntos
for i in range(n - 1):
    intervals_f.append(sp.Interval(all_points_f[i], all_points_f[i + 1]))

# Intervalo desde el último punto hasta ∞
intervals_f.append(sp.Interval(all_points_f[-1], sp.oo))

print("\nIntervalos a analizar:", intervals_f)

# Analizar signo en cada intervalo
increasing = []
decreasing = []

for interval in intervals_f:
    # Elegir punto de prueba adecuado
    if interval.start == -sp.oo:
        test_point = interval.end - 1  # Punto antes del final
    elif interval.end == sp.oo:
        test_point = interval.start + 1  # Punto después del inicio
    else:
        test_point = (interval.start + interval.end) / 2  # Punto medio
    
    # Evaluar el signo
    try:
        sign_value = sp.sign(f_d.subs(x, test_point))
        
        if sign_value > 0:
            increasing.append(interval)
            print(f"En {interval}: f'(x) > 0 → CRECIENTE")
        elif sign_value < 0:
            decreasing.append(interval)
            print(f"En {interval}: f'(x) < 0 → DECRECIENTE")
        else:
            print(f"En {interval}: f'(x) = 0")
            
    except:
        print(f"En {interval}: No se pudo evaluar")


print("Intervalos de crecimiento:", increasing)
print("Intervalos de decrecimiento:", decreasing)

#---Intervalos cuando la funcion es concava hacia arriba o hacia abajo---

inflection_points = [sol for sol in sp.solve(sp.Eq(f_d_II, 0), x) if sol.is_real]
print("Puntos de inflexión REALES:", inflection_points)

num_dd, den_dd = sp.fraction(f_d_II)
discontinuities_dd = [sol for sol in sp.solve(sp.Eq(den_dd, 0), x) if sol.is_real]
print("Discontinuidades REALES:", discontinuities_dd)

#Todos los puntos importantes (solo reales)
all_points_dd = []
if inflection_points:
    all_points_dd.extend(inflection_points)
if discontinuities_dd:
    all_points_dd.extend(discontinuities_dd)

# Eliminar duplicados y ordenar (ahora sí funcionará porque todos son reales)
all_points_dd = sorted(list(set(all_points_dd)))
print("Todos los puntos REALES:", all_points_dd)

intervals_dd = []

if len(all_points_dd) == 0:
    intervals_dd.append(sp.Interval(-sp.oo, sp.oo))
    print("Solo un intervalo: (-∞, ∞)")
else:
    intervals_dd.append(sp.Interval(-sp.oo, all_points_dd[0]))
    for i in range(len(all_points_dd) - 1):
        intervals_dd.append(sp.Interval(all_points_dd[i], all_points_dd[i + 1]))
    intervals_dd.append(sp.Interval(all_points_dd[-1], sp.oo))

print("Intervalos a analizar:", intervals_dd)

concave_up = []
concave_down = []

for interval in intervals_dd:
    if interval.start == -sp.oo and interval.end == sp.oo:

        test_point = 0
    elif interval.start == -sp.oo:
        test_point = interval.end - 1
    elif interval.end == sp.oo:
        test_point = interval.start + 1
    else:
        test_point = (interval.start + interval.end) / 2
    
    try:
        sign_val = sp.sign(f_d_II.subs(x, test_point))
        if sign_val > 0:
            concave_up.append(interval)
            print(f"En {interval}: f''(x) > 0 → Cóncava arriba (∪)")
        elif sign_val < 0:
            concave_down.append(interval)
            print(f"En {interval}: f''(x) < 0 → Cóncava abajo (∩)")
    except:
        print(f"En {interval}: Error al evaluar")

print("Cóncava hacia arriba (∪):", concave_up)
print("Cóncava hacia abajo (∩):", concave_down)
print("Puntos de inflexión:", inflection_points)
