n = 6
m = 6
h = 0.02
voltage = 15
ground = 0
f = open('file.dat','w')

#generate node matrix
nodes = []
for i in range(1,n+1):
    start = (i-1)*m+1
    nodes.append(list(range(start,start+m)))

#part 1
for j in range(n):
    for i in range(m):
        x = i*h
        y = j*h
        if nodes[j][i] not in [35,36]:
            f.write('%d %.3f %.3f\n'%(nodes[j][i], x, y))

#part 2
'''node1---node2
     |       |
   node3---node4
'''
f.write('\n')
for j in range(n-1):
    for i in range(m-1):
        two_element_mesh = [nodes[j][i],nodes[j][i+1],nodes[j+1][i],nodes[j+1][i+1]]
        if (35 in two_element_mesh) or (36 in two_element_mesh):
            continue
        f.write('%d %d %d %.3f\n'%(two_element_mesh[0],two_element_mesh[2],two_element_mesh[3],0))
        f.write('%d %d %d %.3f\n' % (two_element_mesh[0], two_element_mesh[3], two_element_mesh[1], 0))

#part 3
f.write('\n')
for node in nodes[0]:
    f.write('%d %.3f\n'%(node, ground))
for row in nodes:
    f.write('%d %.3f\n' % (row[0], ground))
for node in [28,29,30,34]:
    f.write('%d %.3f\n' % (node, voltage))

f.close()


