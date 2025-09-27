
from math import copysign, sqrt, sin


def muller(f, x0, x1, x2, tol, iterMax):

    # Valores de a, b y c
    den = lambda x0, x1, x2: ((x0 - x1) * (x0 - x2) * (x1 - x2))
    a = lambda x0, x1, x2, deno: ((x1 - x2) * (f(x0) - f(x2)) - (x0 - x2) * (f(x1) - f(x2))) / deno
    b = lambda x0, x1, x2, deno: ((x0 - x2)**2 * (f(x1) - f(x2)) - (x1 - x2)**2 * (f(x0) - f(x2))) / deno
    c = lambda x2: f(x2)

    # Raiz para obtener x3
    r = lambda a, b, c, x2: x2 - (2 * c) / (b + copysign(1, b) * sqrt(b ** 2 - 4 * a * c))

    # Valores iniciales
    x0_ = x0
    x1_ = x1
    x2_ = x2
    x3_ = 0         # Aqui es donde se guarda la raiz
    a_, b_, c_ = 0, 0, 0
    iter = 0
    den_ = 0
    err = 0

    # Iterar
    for iter in range(1, iterMax + 1):
        den_ = den(x0_, x1_, x2_)
        a_ = a(x0_, x1_, x2_, den_)
        b_ = b(x0_, x1_, x2_, den_)
        c_ = c(x2_)
        x3_ = r(a_, b_, c_, x2_)

        # Error
        err = abs(f(x3_))
        if(err < tol):
            break     
        
        # Mover xs
        x0_ = x1_
        x1_ = x2_
        x2_ = x3_

    print("xk =", x3_)
    print("|f(xk)| =", err)
    print("Iteraciones:", iter)

# Ejecutar
if __name__ == "__main__":
    func = lambda x: (1 + x) * sin(x) - 1
    x0 = 2.5
    x1 = 2.75
    x2 = 3
    muller(func, x0, x1, x2, 10 ** (-10), 1000)



        

