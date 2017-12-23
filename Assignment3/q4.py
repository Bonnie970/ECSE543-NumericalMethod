import math
    
def integral(func1,a, b, even=True, width_factor=10):
    if even == True:
        w = float(b - a) / float(width_factor)
        widths = [w]*int(width_factor)
    else:
        relativeWidths = [width_factor**i for i in range (0,10)]
        a,b = float(a),float(b)
        scale = (b - a) / sum(relativeWidths)
        widths = [width * scale for width in relativeWidths]
    summation = 0
    for w in widths:
        lowLim = a
        a += w
        highLim = a  
        #try:
        #    height = (func1(lowLim) + func1(highLim)) / 2  
        #except:
        height = func1((lowLim + highLim) / 2)
        summation += (highLim - lowLim) * height   
    return summation
    
ground_truth = math.sin(1) + math.sin(0)
print ("Integral of cos(x): ", ground_truth)
for i in range (1,21): 
    result =  integral (math.cos,0,1,True,i)
    #print ("i = ", i,  " Integral = ", result, " Error: ", result-ground_truth)
    print(result - ground_truth)

ground_truth = -1
print ("Integral of ln(x): ", ground_truth) 
for i in range (10, 210, 10):  
    result =  integral (math.log,0,1,True,i)
    #print ("i = " ,i, " Integral = ",result, " Error: " , result-ground_truth)
    print(result - ground_truth)

print ('Uneven integral')
print ("Integral of ln(x): ")
for width_factor in [1,1.2,1.3,1.4,1.45,1.5,1.6,1.7,1.8,1.9,2]:    
    result =  integral (math.log,0,1,False,width_factor)
    print ("Width factor = ", width_factor," Integral = ",result," Error: ",result - ground_truth)


