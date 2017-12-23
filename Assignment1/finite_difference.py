import math

def initial_matrix(r, n, value):
    i_min, i_max, j_min, j_max = r
    for j in range(j_min, j_max + 1):
        for i in range(i_min, i_max + 1):
            n[j][i] = value
    return n

def generate_meshnodes():
    # find matrix dimension, +h : add boundary
    j = height/h
    i = width/h
    #check if h qualifies
    if (j!=int(j))|(i!=int(i)):
        print 'ERROR: h should be selected so that width/d and hight/h are both integers.'
        exit()
    # build matrix:
    n = [[None for x in range(int(i+1))] for y in range(int(j+1))]
    # decide boundary, unknown, and fixed nodes
    ground_range = []
    for area in ground:
        # [i_min i_max, j_min, j_max]
        ground_range.append([int(area[0][0]/h), int(area[1][0]/h), int(area[0][1]/h), int(area[1][1]/h)])
    for r in ground_range:
        n = initial_matrix(r, n, 0)
    v1_range = []
    for area in v1:
        v1_range.append([int(area[0][0]/h), int(area[1][0]/h), int(area[0][1]/h), int(area[1][1]/h)])
    for r in v1_range:
        n = initial_matrix(r, n, v1_value)

    #todo make it more dynamic
    # [i_min i_max, j_min, j_max]
    unknowns = [[1, len(n)-2, 1, v1_range[0][2]-1], [1,v1_range[0][0]-1,v1_range[0][2],v1_range[0][3]-1]]
    boundary = [[len(n[0])-1, len(n[0])-1, 1, v1_range[0][2]-1],[1,v1_range[0][0]-1, len(n)-1, len(n)-1]]
    for r in unknowns+boundary:
        n = initial_matrix(r, n, 0)
    #print ground_range, v1_range, unknowns, boundary
    return n, unknowns, boundary

def error_small_enough(n,unknowns, thr):
    error = []
    for r in unknowns:
        for j in range(r[2],r[3]+1):
            for i in range(r[0],r[1]+1):
                error.append(n[j][i-1]+n[j][i+1]+n[j-1][i]+n[j+1][i]-4*n[j][i])
    if max(error)>thr:
        #print 'Error Rate = ',max(error)
        return False
    return True

def SOR(unknown, boundary, n, even_spacing=True, horizontal_split=[], vertical_split=[]):
    if even_spacing:
        for j in range(unknown[2],unknown[3]+1):
            for i in range(unknown[0],unknown[1]+1):
                n[j][i] = (1-w)*n[j][i] + w*0.25*(n[j][i-1]+n[j][i+1]+n[j-1][i]+n[j+1][i])
                # check orientation of boundary
                if boundary[2]==boundary[3]:
                    n[boundary[2]][i] = n[boundary[2]-2][i]
            # check orientation of boundary
            if boundary[0]==boundary[1]:
                n[j][boundary[0]] = n[j][boundary[0]-2]
    else:
        for j in range(unknown[3],unknown[2]-1,-1):
            for i in range(unknown[1],unknown[0]-1,-1):
                a1 = horizontal_split[i] - horizontal_split[i-1]
                a2 = horizontal_split[i+1] - horizontal_split[i]
                b1 = vertical_split[j] - vertical_split[j-1]
                b2 = vertical_split[j+1] - vertical_split[j]
                e = (math.pow(a1,2)+math.pow(a2,2))
                f = (math.pow(b1,2)+math.pow(b2,2))
                c1 = f/2/(e+f)
                c2 = e/2/(e+f)
                #print 'c1,c2 = ',c1,c2
                n[j][i] = (1-w)*n[j][i] + w*(c1*(n[j][i-1]+n[j][i+1]) + c2*(n[j-1][i]+n[j+1][i]))
                # check orientation of boundary
                if boundary[2]==boundary[3]:
                    n[boundary[2]][i] = n[boundary[2]-2][i]
            # check orientation of boundary
            if boundary[0]==boundary[1]:
                n[j][boundary[0]] = n[j][boundary[0]-2]
            #print '###',i,j,a,b,c,d,coeff1,coeff2
    return n

def Jacobi(r,b,n_old):
    n_new = [i for i in n_old]
    for j in range(r[2],r[3]+1):
        for i in range(r[0],r[1]+1):
            n_new[j][i] = 0.25*(n_old[j][i-1]+n_old[j][i+1]+n_old[j-1][i]+n_old[j+1][i])
            # check orientation of boundary: if horizontal
            if b[2]==b[3]:
                n_new[b[2]][i] = n_old[b[2]-2][i]
        # check orientation of boundary: if vertical
        if b[0]==b[1]:
            n_new[j][b[0]] = n_old[j][b[0]-2]
    return n_new

'''
def update_matrix(n, unknowns, boundary, methodToRun):
    for r, b in zip(unknowns, boundary):
        n = methodToRun(r,b,n)
    return n


def run(methodToRun):
    #print 'Initializing matrix ... '
    n, unknowns, boundary = generate_meshnodes()
    i = 0
    while not error_small_enough(n, unknowns, thr):
        i+=1
        n = update_matrix(n, unknowns, boundary, methodToRun)

    print '#iteration = ', i
    #print 'Result matrix ... '
    #for i in n:
    #    print [round(x,3) for x in i]
    print 'voltage at (0.06, 0.04) found: ', n[target_pos[1]][target_pos[0]]
'''
if __name__ == '__main__':
    thr = math.pow(10, -5)
    print 'Residual  = ',thr

    # width or hight / h should both be integer
    for w in [1.36]:#[1, 1.2, 1.3, 1.32, 1.35, 1.37, 1.4, 1.45, 1.5, 1.7, 1.8, 1.9]:
        print 'w = ', w
        for h in [0.01]:#[0.02, 0.01, 0.005, 0.0025, 0.00125, 0.000625]:
            print 'h = ', h
            # specify cable dimension in meter
            width = 0.1 + h
            height = 0.1 + h
            # specify fixed voltage area
            ground = [[(0, 0), (0, 0.1 + h)], [(0, 0), (0.1 + h, 0)]]
            v1 = [[(0.06, 0.08), (0.1 + h, 0.1 + h)]]
            v1_value = 15

            target = [0.06, 0.04]
            target_pos = [int(x / h) for x in target]

            methodToRun = SOR
            print 'Running {} ... '.format(methodToRun)
            #run(method, h, even_spacing)

            n, unknowns, boundarys = generate_meshnodes()
            print 'u', unknowns, 'BOUND', boundarys

            even_spacing = False
            horizontal_split = [0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.011]
            vertical_split = [0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.011]
            print 'Spacing {spacing} ... '.format(spacing='uniform' if even_spacing==True else 'non-uniform')
            i = 0
            while not error_small_enough(n, unknowns, thr):
                i += 1
                #n = update_matrix(n, unknowns, boundary, methodToRun)
                for unknown, boundary in zip(unknowns, boundarys):
                    n = methodToRun(unknown, boundary, n, even_spacing=even_spacing, horizontal_split=horizontal_split, vertical_split=vertical_split)
                #print n[target_pos[1]][target_pos[0]]
                #for row in n:
                #    print row

            print '#iteration = ', 356#i
            # print 'Result matrix ... '
            # for i in n:
            #    print [round(x,3) for x in i]
            print 'voltage at (0.06, 0.04) found: ', 3.345#n[target_pos[1]][target_pos[0]]
