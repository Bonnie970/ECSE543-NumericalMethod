from sympy import symbols, expand, diff, lambdify, Piecewise

x = symbols('x')
B = [0.0, 0.2,  0.4,  0.6,  0.8,   1.0,   1.1,   1.2,   1.3,   1.4,    1.5,    1.6,    1.7,    1.8,     1.9]
H = [0.0, 14.7, 36.5, 71.7, 121.4, 197.4, 256.2, 348.7, 540.6, 1062.8, 2318.0, 4781.9, 8687.4, 13924.3, 22650.2]


def find_piecewise(X,Y):
    f = 0
    for j in range(len(X)-1):
        x1 = X[j]
        x2 = X[j+1]
        y1 = Y[j]
        y2 = Y[j+1]
        k = (y2-y1)/(x2-x1)
        b = y2 - k*x2
        if j==0:
            g = Piecewise((0, x >= x2), (k * x + b, True))
        elif j==len(X)-2:
            g = Piecewise((0, x < x1), (k * x + b, True))
        else:
            g = Piecewise((0, x<x1),(0, x>=x2),(k*x+b, True))
        #print x1,x2,y1,y2
        #print g.subs(x,x2)
        f += g
    #f = lambdify(x, f)
    return f

def run_newton_raphson(B,H,err):
    #find piecewise curve of B and H
    h = find_piecewise(B, H)
    print('Piecewise linear function: ', h)
    #initial flux and iterration counter
    flux = 0
    i = 0
    #initial reduction function F and its derivative
    F = 39.788735e6 * x + 0.3 * h.subs(x, x / 1e-4) - 8000
    dF = diff(F)
    while True:
        f = F.subs(x,flux)
        df = dF.subs(x,flux)
        #print(i, flux, f, df)
        print('Interation: ',i, ' Flux: ',flux, ' f: ',f, ' df: ',df)
        if (abs(f/df)<err):
            break
        else:
            i += 1
            flux += -f/df


def run_substitution(B,H,err):
    #find piecewise curve of B and H
    h = find_piecewise(B, H)
    print('Piecewise linear function: ', h)
    #initial flux and iterration counter
    flux = 0
    i = 0
    #initial reduction function F and its derivative
    F = 39.788735e6 * x + 0.3 * h.subs(x, x / 1e-4) - 8000
    F /= 1e8
    while True:
        f = F.subs(x,flux)
        #print(i, flux, f, df)
        print('Interation: ',i, ' Flux: ',flux, ' f: ',f)
        if (abs(f)<err):
            break
        else:
            i += 1
            flux += -f

err = 1e-6

run_newton_raphson(B,H,err)

run_substitution(B,H,err)