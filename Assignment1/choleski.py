"""
This script solves Ax=b matrix problem using Cholesky Decomposition, where A is represented by L*transpose(L)
Author: Guanqing Hu, Oct 3rd, 2017
"""
import math
import numpy as np
import time

# check S.P.D
def check_symmetric(A):
    n = len(A)
    for i in range(n):
        for j in range(i + 1, n):
            if A[i][j] != A[j][i]:
                exit('ERROR: Input matrix must be Symmetric!')


# positive definite if determinant of A is positive
# determinant equals to the product of all the eigenvalues of A
def check_pd(A):
    det = np.linalg.det(A)
    print det
    if det <= 0:
        exit('ERROR: Input matrix must be Positive Definite!')

def read_matrix():
    # read matrix, last line is b
    with open('matrix.txt', 'r') as f:
        A = []
        for line in f:
            if line.strip() == '':
                continue
            A.append([float(x) for x in line.split(',')])
            # check if data complete
            # todo
        b = A[-1]
        A = A[:-1]
        print 'Matrix being solved: A = ', A, '\nb = ', b
        return A,b


def cholesky(A, b, hb=None):
    # Cholesky implementation
    n = len(A)
    print 'Half band = ', hb
    start_time = time.time()
    for j in range(n - 1):
        if A[j][j] <= 0:
            exit('ERROR: Input matrix must be Positive Definite!')
        A[j][j] = math.sqrt(A[j][j])
        b[j] = b[j] / A[j][j]

        i_range = range(j + 1, n)
        if (hb!=None):
            if ((j+hb+1)<n):
                i_range = range(j+1, (j+hb+1))

        for i in i_range: #range(j+1, n):
            #print 'updating element A{}{}'.format(i,j)
            A[i][j] = A[i][j] / A[j][j]
            b[i] = b[i] - A[i][j] * b[j]
            for k in range(j + 1, i + 1):
                A[i][k] = A[i][k] - A[i][j] * A[k][j]
    # Back Substitution: T(L)x = y
    x = [0.0 for i in range(n)]
    for i in range(n - 1, -1, -1):
        x[i] = round((b[i] - sum([A[j][i] * x[j] for j in range(i + 1, n)])) / A[i][i], 2)
    solving_time = time.time()-start_time
    print 'Cholesky solving time: {}s'.format(solving_time)
    print 'Cholesky solving matrix dimension: ', len(A)

    return x

def solve_matrix():
    A, b = read_matrix()
    check_symmetric(A)
    #check_pd(A)
    results = cholesky(A, b, None)
    print 'Solution found: ', results

def read_network():
    with open('networkBranch.txt', 'r') as f:
        A = []
        for line in f:
            if line.strip() == '':
                continue
            A.append([float(x) for x in line.split(',')])
        J, R, E = A[-3], A[-2], A[-1]
        A = A[:-3]
        return A, J, R, E


def m_product(A, B):
    # return A*B
    # check dimension
    if len(A[0]) != len(B):
        exit('ERROR: Matrix production failed due to wrong dimensions!')
    # create correct dimension output matrix
    #print 'product', A, B
    C = [[0 for col in B[0]] for row in A]
    for i in range(len(C)):
        for j in range(len(C[0])):
            C[i][j] = sum([A[i][k] * B[k][j] for k in range(len(B))])
    return C


def transpose(A):
    A2 = [[0 for i in A] for j in A[0]]
    for i in range(len(A)):
        for j in range(len(A[0])):
            A2[j][i] = A[i][j]
    return A2

def m_substract(A, B):
    # check dimension
    if (len(A)!=len(B))|(len(A[0])!=len(B[0])):
        exit('ERROR: Matrix dimension error in substraction!')
    C = []
    for i in range(len(A)):
        C.append([A[i][k]-B[i][k] for k in range(len(A[0]))])
    return C


# this method generates network for a regular N by 2N finite-difference mesh and replace each horizontal and vertical
# line by a 1 kW resistor
def generate_network(N, R, sJ, sR, sE):
    # calculate number of branches
    nbranch = (N+1)*2*N+N*(2*N+1)
    J = [0]*nbranch
    R = [R]*nbranch
    E = [0]*nbranch
    # add teasting branch: voltage source branch
    J.append(sJ)
    R.append(sR)
    E.append(sE)
    # calculate nodes
    nnode = (N+1)*(2*N+1)
    A = [[0]*(nbranch+1) for i in range(nnode)]
    #print len(A), len(A[0])
    for i in range(1, nnode+1):
        level = i/(2*N+1)
        offset = i%(2*N+1)
        if offset==0:
            level -= 1
            offset = 2*N+1
        branch_per_level = (2*N+(2*N+1))
        # calculate surrounding branch indices
        right = branch_per_level * level + offset
        left = right - 1
        top = 2*N*level + (2*N+1)*(level-1)+ offset
        bottom = top + branch_per_level
        #print i, 'r',right, 'l', left, 't',top, 'b', bottom

        i -= 1
        # node on top border
        if (level==0):
            # top left
            if (offset==1):
                A[i][right - 1], A[i][bottom - 1] = -1, 1
            # top right
            elif (offset==2*N+1):
                A[i][left - 1], A[i][bottom - 1] = 1, 1
            else:
                A[i][right - 1], A[i][left - 1], A[i][bottom - 1] = -1, 1, 1
        # node on bottom border
        elif (level == N):
            # bottom left
            if (offset==1):
                A[i][right - 1], A[i][top - 1] = -1, -1
            # bottom right
            elif (offset==2*N+1):
                A[i][top - 1], A[i][left - 1] = -1, 1
            else:
                A[i][right - 1], A[i][top - 1], A[i][left - 1] = -1, -1, 1
        # node on left border
        elif (offset == 1):
            A[i][right - 1], A[i][top - 1], A[i][bottom - 1] = -1, -1, 1
        # node on right border
        elif (offset == 2*N+1):
            A[i][top - 1], A[i][left - 1], A[i][bottom - 1] = -1, 1, 1
        # nodes in middle
        else:
            A[i][right - 1], A[i][top - 1], A[i][left - 1], A[i][bottom - 1] = -1, -1, 1, 1

        # setup the voltage source
        # right top corner
        if (level==0)&(offset==2*N+1):
            A[i][-1] = -1
        elif (level==N)&(offset==1):
            A[i][-1] = 1
    # set up ground node by removing that node from A
    # ground left-bottom corner node
    del A[-(2*N+1)]

    return A, J, R, E


def solve_network(A, J, R, E, hb=None):
    # (A * y * T(A)) * v = A * (J - y * E)

    y = [[0 for r in R] for r in R]
    for i in range(len(y)):
        y[i][i] = 1.0/R[i]

    A2 = m_product(m_product(A, y), transpose(A))
    check_symmetric(A2)
    #check_pd(A2)

    '''print '### start of matrix'
    for i in A2:
        print i
    '''
    b = m_product(A, m_substract([[j] for j in J], m_product(y, [[e] for e in E])))
    # example of b: [[0.5], [0.0]], flat it
    b = [i[0] for i in b]
    #print 'b', b, 'A2', A2
    v = cholesky(A2, b, hb)
    return v

def fit_curve(x,y):
    print 'Fitting curve ... '
    print 'x = ',x
    print 'y = ',y
    x = np.array(x)
    y = np.array(y)
    k, b = np.polyfit(np.log(x), y, 1)
    print 'Curve found: y = {}*ln(x) + {}'.format(k,b)
    r = []
    for n in x:
        r.append(k * math.log(n, math.e) + b)
    error = max([(a0-a)/a for (a0,a) in zip(r, y)])
    print 'Error of curve = ', error


if __name__ == '__main__':
    ### q1
    solve_matrix()

    ### q1
    A, J, R, E = read_network()
    print 'Reading network ...\nA = {}\nJ = {}\nR = {}\nE = {}'.format(A,J,R,E)
    v = solve_network(A, J, R, E)
    print 'Voltage found: ',v

    ### q2
    # J, R ,E of the source branch
    sJ = 0
    sR = 1000.0
    sE = 100.0
    for N in range(2,2):
        print '\nN = {}'.format(N)
        A, J, R, E = generate_network(N, 1000, sJ, sR, sE) #N, R, sJ, sR, sE
        # half band, not sure ...
        #hb = 2*N+2
        v = solve_network(A, J, R, E)
        # get the voltage at top right and bottom left
        v1 = v[2*N]
        v2 = 0#v[-1-2*N]
        # R_total = sR+R_mesh
        I  = ((v2+sE)-v1)/sR
        R_mesh = (v1-v2)/I
        #print 'Voltage of {} nodes found. '.format(len(v))#, v
        print 'R_eq = ', R_mesh

    # fit curve
    N = [2, 3, 4, 5, 6, 7, 8, 9, 10]
    R = [2057.17, 2497.73, 2828.48, 3089.98, 3308.49, 3494.38, 3659.83, 3803.07, 3933.4]
    fit_curve(N, R)

