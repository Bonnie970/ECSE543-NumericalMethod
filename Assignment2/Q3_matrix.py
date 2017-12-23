import time
from math import sqrt
from choleski import transpose, check_pd, m_product, m_substract, m_sum,cholesky

n = 6
m = 6
h = 0.02
voltage = 15
ground = 0

#generate node matrix
def gen_nodes():
    nodes = []
    for i in range(1,n+1):
        start = (i-1)*m+1
        nodes.append(list(range(start,start+m)))
    return nodes


'''        node14
             |
   node7---node8---node9
             |
           node2
'''
interior_uk = [8,9,10,11,14,15,16,17,20,21,22,23,26,27]
neuman_uk = [12,18,24,32,33]
uk = interior_uk+neuman_uk
def gen_Ax_b():
    nodes = gen_nodes()
    A1 = []
    A2 = []
    for j in range(n):
        for i in range(m):
            if nodes[j][i] not in uk:
                continue
            else: 
                row = [0]*34
                if nodes[j][i] in interior_uk:
                    row[nodes[j][i]-1] = -4
                    row[nodes[j][i-1]-1] = 1
                    row[nodes[j][i+1]-1] = 1
                    row[nodes[j-1][i]-1] = 1
                    row[nodes[j+1][i]-1] = 1
                elif nodes[j][i] in neuman_uk:
                    row[nodes[j][i]-1] = -4
                    if (i+1)<m:
                        row[nodes[j][i+1]-1] = 1
                        row[nodes[j-1][i]-1] = 2
                        row[nodes[j][i-1]-1] = 1
                    if (j+1)<n:
                        row[nodes[j+1][i]-1] = 1
                        row[nodes[j][i-1]-1] = 2
                        row[nodes[j-1][i]-1] = 1
                #move matrix columns so that free nodes and fixed nodes separated
                col1 = [row[i] for i in range(len(row)) if (i+1) in uk]
                col2 = [row[i] for i in range(len(row)) if (i+1) not in uk]       
                A1.append(col1)
                A2.append(col2)
    print('A_free = ')
    for r in A1:
        print(r)
    print('A_fixed = ')
    for r in A2:
        print(r)
    x = [[0] for i in range(19)]
    print('v_free = {}\n'.format([i[0] for i in x]))
    f = [0,0,0,0,0,0,0,0,0,0,15,15,15,0,15]
    print('v_fixed = {}\n'.format(f))
    f = [[-1*x] for x in f]
    b = m_product(A2, f)
    print('A = ')
    for r in A1:
        print(r)
    print('b = {}\n'.format([i[0] for i in b]))
    return A1,x,b

def twonorm(x):
    if isinstance(x[0],list):
        x = [i[0] for i in x]
    result = 0
    for item in x:
        result += item**2
    return sqrt(result)

def infinitynorm(x):
    if isinstance(x[0],list):
        x = [i[0] for i in x]
    result = max([abs(i) for i in x])
    return result
    
def conjugate_gradient(A, b, x):
    #run conjugate
    start_time = time.time()
    r = m_substract(b,m_product(A, x))
    p = r
    print('iteration,infinity norm,2-norm')
    for k in range(len(x)):
        print('%d,%f,%f'%(k+1, round(max(r)[0],6),round(twonorm(r),6)))
        Ap = m_product(A,p)
        alpha = m_product(transpose(p),r)[0][0]/m_product(transpose(p),Ap)[0][0]
        x = [[X[0]+alpha*P[0]] for X,P in zip(x,p)]
        r = m_substract(b, m_product(A,x))
        Ar = m_product(A,r)
        beta = -1*m_product(transpose(p),Ar)[0][0]/m_product(transpose(p),Ap)[0][0]
        p = [[R[0]+beta*P[0]] for R,P in zip(r,p)]
    solving_time = time.time()-start_time
    print('Conjugate gradient solving time: {}s'.format(solving_time))
    return x
    

def pd_convert(A, b):
    if not check_pd(A):
        A_t = transpose(A)
        A = m_product(A_t, A)
        b = m_product(A_t, b)
        check_pd(A)
    return A,b

#run cheloski
A,x,b = gen_Ax_b()
A,b = pd_convert(A,b)
print('-------Chelisky--------')
# ???????? A2 = copy.copy(A)
x1 = cholesky(A, b, hb=None)
for i,v in enumerate(x1):
    print('%d : %f'%(i,v))

    
#run conjugate gradient
A,x,b = gen_Ax_b()
A,b = pd_convert(A,b)
print('A = ')
for r in A:
    print(r)
print('b = {}'.format([i[0] for i in b]))
print('-------Conjugate Gradient--------')
x2 = conjugate_gradient(A, b, x)
for i,v in enumerate(x2):
    print('%d : %f'%(i,v[0]))

