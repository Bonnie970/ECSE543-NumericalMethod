from sympy import symbols, expand, diff, lambdify
from sympy.functions import exp

v1,v2 = symbols('v1 v2')


def inverse2b2(a,b,c,d):
    det = a*d-b*c
    #check if possible
    if det==0:
        print('Matrix not invertible.')
        return [None,None,None,None]
    else:
        return [d/det,-1*b/det,-1*c/det, a/det]


def run_newton_raphson(err):
    #initial node voltages and iterration counter
    v = [0,0]
    i = 0
    #initial reduction function F and its derivative
    f1 = 44.06e-5 - 0.002 * v1 - 6e-7 * exp(40 * (v1 - v2))
    f2 = 44.12e-5 - 0.002 * v1 - 12e-7 * exp(40 * v2)
    df11 = diff(f1,v1)
    df12 = diff(f1,v2)
    df21 = diff(f2,v1)
    df22 = diff(f2,v2)
    while True:
        f1t = f1.subs([(v1,v[0]),(v2,v[1])])
        f2t = f2.subs([(v1, v[0]), (v2, v[1])])
        df11t = df11.subs([(v1, v[0]), (v2, v[1])])
        df12t = df12.subs([(v1, v[0]), (v2, v[1])])
        df21t = df21.subs([(v1, v[0]), (v2, v[1])])
        df22t = df22.subs([(v1, v[0]), (v2, v[1])])
        inv = inverse2b2(df11t, df12t, df21t, df22t)
        r = [inv[0]*f1t+inv[1]*f2t, inv[2]*f1t+inv[3]*f2t]
        #print('Interation: ',i, ' v1: ',v[0], ' v2: ',v[1], ' f\'*f: ',r)
        #print(i, ',', v[0], ',', v[1], ',', r[0], ',', r[1])
        r2 = [abs(x) for x in r]
        print(max(r2))
        if (max(r2)<err):
            break
        else:
            i += 1
            v[0] += -r[0]
            v[1] += -r[1]

run_newton_raphson(1e-5)