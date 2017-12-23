from sympy import symbols, expand, diff, lambdify
from sympy.plotting import plot

def full_lagrange_interpolate(X, Y):
    # lagrange dimension decided by length of X and Y
    n = len(X)
    x = symbols('x')
    Ls = []
    # find Fj for each j in n
    for j in range(n):
        X_sub = X[:j]+X[j+1:]
        F = 1
        for xr in X_sub:
            F *= x-xr
        F_f = lambdify(x, F)
        L = expand(F/F_f(X[j]))
        Ls.append(L)
    y = 0
    for j, L in enumerate(Ls):
        y += Y[j]*L
    return expand(y)


def cubic_hermite(X,Y):
    n = len(X)
    x = symbols('x')
    Us = []
    Vs = []
    for j in range(n):
        X_sub = X[:j] + X[j + 1:]
        F = 1
        for xr in X_sub:
            F *= x - xr
        F_f = lambdify(x, F)
        L = F / F_f(X[j])
        L_d = lambdify(x, diff(L))
        U = (1 - 2 * L_d(X[j]) * (x - X[j]))*(L**2)
        V = (x - X[j])*(L**2)
        Us.append(U)
        Vs.append(V)
    y = 0
    # initial b[j] = y'(j)
    b = [(Y[j+1]-Y[j])/(X[j+1]-X[j]) for j in range(n-1)]
    b.append(Y[-1]/X[-1])
    #print(b)
    for j in range(n):
        y += Y[j]*Us[j] + b[j]*Vs[j]
    return expand(y)


B = [0.0, 0.2,  0.4,  0.6,  0.8,   1.0,   1.1,   1.2,   1.3,   1.4,    1.5,    1.6,    1.7,    1.8,     1.9]
H = [0.0, 14.7, 36.5, 71.7, 121.4, 197.4, 256.2, 348.7, 540.6, 1062.8, 2318.0, 4781.9, 8687.4, 13924.3, 22650.2]

x = symbols('x')
lag = full_lagrange_interpolate(B[:6],H[:6])
p1 = plot(lag,(x,0,2.5),xlabel='B',ylabel='H')
print('Full Lagrange Interpolation 1: ',lag, '\n')

lag2 = full_lagrange_interpolate(B[:1]+B[8:10]+B[12:],H[:1]+H[8:10]+H[12:])
p2 = plot(lag2, (x,0,2),xlabel='B',ylabel='H')
print('Full Lagrange Interpolation 2: ',lag2, '\n')

cub = cubic_hermite(B[:1]+B[8:10]+B[12:],H[:1]+H[8:10]+H[12:])
p2 = plot(cub, (x,0,100),xlabel='B',ylabel='H')
print('Cubic Hermite Interpolation: ',cub, '\n')
