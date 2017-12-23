import math


class finite_difference:
    def __init__(self, width, height, voltage,h=0.01,target=[0.0,0.0],threshold=0.001, w=1.5, horizontal_split=[],
                 vertical_split=[],even_spacing=True):
        self.thr = threshold
        self.w = w
        self.h = h
        # add one more layer of nodes, their voltages can be found using symmetry
        self.width = width + h
        self.height = height + h
        # specify fixed voltage area
        self.ground = [[(0, 0), (0, self.width)], [(0, 0), (self.height, 0)]]
        self.v1 = [[(0.06, 0.08), (0.1 + h, 0.1 + h)]]
        self.v1_value = voltage

        self.target = target
        self.target_pos = [int(x / h) for x in target]

        self.horizontal_split = horizontal_split
        self.vertical_split = vertical_split
        self.even_spacing = even_spacing

    def initial_matrix(self,r, n, value):
        i_min, i_max, j_min, j_max = r
        for j in range(j_min, j_max + 1):
            for i in range(i_min, i_max + 1):
                n[j][i] = value
        return n

    def generate_meshnodes(self):
        # find matrix dimension, +h : add boundary
        j = self.height/self.h
        i = self.width/self.h
        #check if h qualifies
        if (j!=int(j))|(i!=int(i)):
            print('ERROR: h should be selected so that width/d and hight/h are both integers.')
            exit()
        # build matrix:
        n = [[None for x in range(int(i+1))] for y in range(int(j+1))]
        # decide boundary, unknown, and fixed nodes
        ground_range = []
        for area in self.ground:
            # [i_min i_max, j_min, j_max]
            ground_range.append([int(area[0][0]/self.h), int(area[1][0]/self.h), int(area[0][1]/self.h),
                                 int(area[1][1]/self.h)])
        for r in ground_range:
            n = self.initial_matrix(r, n, 0)
        v1_range = []
        for area in self.v1:
            v1_range.append([int(area[0][0]/self.h), int(area[1][0]/self.h), int(area[0][1]/self.h),
                             int(area[1][1]/self.h)])
        for r in v1_range:
            n = self.initial_matrix(r, n, self.v1_value)

        #todo make it more dynamic
        # [i_min i_max, j_min, j_max]
        unknowns = [[1, len(n)-2, 1, v1_range[0][2]-1], [1,v1_range[0][0]-1,v1_range[0][2],v1_range[0][3]-1]]
        boundary = [[len(n[0])-1, len(n[0])-1, 1, v1_range[0][2]-1],[1,v1_range[0][0]-1, len(n)-1, len(n)-1]]
        for r in unknowns+boundary:
            n = self.initial_matrix(r, n, 0)
        return n, unknowns, boundary

    def error_small_enough(self,n,unknowns, thr):
        error = []
        for r in unknowns:
            for j in range(r[2],r[3]+1):
                for i in range(r[0],r[1]+1):
                    error.append(n[j][i-1]+n[j][i+1]+n[j-1][i]+n[j+1][i]-4*n[j][i])
        if max(error)>thr:
            return False
        return True

    def SOR(self,unknown, boundary, n, horizontal_split, vertical_split,even_spacing=True):
        if even_spacing:
            for j in range(unknown[2],unknown[3]+1):
                for i in range(unknown[0],unknown[1]+1):
                    n[j][i] = (1-self.w)*n[j][i] + self.w*0.25*(n[j][i-1]+n[j][i+1]+n[j-1][i]+n[j+1][i])
                    # check orientation of boundary
                    if boundary[2]==boundary[3]:
                        n[boundary[2]][i] = n[boundary[2]-2][i]
                # check orientation of boundary
                if boundary[0]==boundary[1]:
                    n[j][boundary[0]] = n[j][boundary[0]-2]
        else:
            for j in range(unknown[3],unknown[2]-1,-1):
                for i in range(unknown[1],unknown[0]-1,-1):
                    a1 = self.horizontal_split[i] - horizontal_split[i-1]
                    a2 = horizontal_split[i+1] - horizontal_split[i]
                    b1 = vertical_split[j] - vertical_split[j-1]
                    b2 = vertical_split[j+1] - vertical_split[j]
                    e = (math.pow(a1,2)+math.pow(a2,2))
                    f = (math.pow(b1,2)+math.pow(b2,2))
                    c1 = f/2/(e+f)
                    c2 = e/2/(e+f)
                    n[j][i] = (1-self.w)*n[j][i] + self.w*(c1*(n[j][i-1]+n[j][i+1]) + c2*(n[j-1][i]+n[j+1][i]))
                    # check orientation of boundary
                    if boundary[2]==boundary[3]:
                        n[boundary[2]][i] = n[boundary[2]-2][i]
                # check orientation of boundary
                if boundary[0]==boundary[1]:
                    n[j][boundary[0]] = n[j][boundary[0]-2]
        return n

    def Jacobi(self,r,b,n_old):
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

    def run_sor(self):
        n, unknowns, boundaries = self.generate_meshnodes()
        i = 0
        while not self.error_small_enough(n,unknowns, self.thr):
            i += 1
            for unknown, boundary in zip(unknowns, boundaries):
                n = self.SOR(unknown, boundary, n, even_spacing=self.even_spacing, horizontal_split=self.horizontal_split,
                            vertical_split=self.vertical_split)
        print('#iteration = ', i)
        print('voltage at (0.06, 0.04) found: ',  n[self.target_pos[1]][self.target_pos[0]])

'''
if __name__ == '__main__':
    fd = finite_difference(width=0.1, height=0.1, h=0.01, voltage=15,target=[0.06,0.04],threshold=0.00001, w=1.36)
    n = fd.generate_meshnodes()[0]
    for row in n:
        print(row)
    fd.run_sor()
'''

