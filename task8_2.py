import matplotlib.pyplot as plt 
import numpy as np
import math
import random
from functools import cmp_to_key
import time
import memory_profiler

#jarvis algorithm
# ---------------------------

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y
 
def Left_index(points):
     
    '''
    Finding the left most point
    '''
    minn = 0
    for i in range(1,len(points)):
        if points[i].x < points[minn].x:
            minn = i
        elif points[i].x == points[minn].x:
            if points[i].y > points[minn].y:
                minn = i
    return minn
 
def orientation(p, q, r):
    '''
    To find orientation of ordered triplet (p, q, r).
    The function returns following values
    0 --> p, q and r are collinear
    1 --> Clockwise
    2 --> Counterclockwise
    '''
    val = (q.y - p.y) * (r.x - q.x) - \
          (q.x - p.x) * (r.y - q.y)
 
    if val == 0:
        return 0
    elif val > 0:
        return 1
    else:
        return 2
 
def jarvis_convexHull(points, n):
     
    # There must be at least 3 points
    if n < 3:
        return
 
    # Find the leftmost point
    l = Left_index(points)
 
    hull = []
     
    '''
    Start from leftmost point, keep moving counterclockwise
    until reach the start point again. This loop runs O(h)
    times where h is number of points in result or output.
    '''
    p = l
    q = 0
    while(True):
         
        # Add current point to result
        hull.append(p)
 
        '''
        Search for a point 'q' such that orientation(p, q,
        x) is counterclockwise for all points 'x'. The idea
        is to keep track of last visited most counterclock-
        wise point in q. If any point 'i' is more counterclock-
        wise than q, then update q.
        '''
        q = (p + 1) % n
 
        for i in range(n):
             
            # If i is more counterclockwise
            # than current q, then update q
            if(orientation(points[p],
                           points[i], points[q]) == 2):
                q = i
 
        '''
        Now q is the most counterclockwise with respect to p
        Set p as q for next iteration, so that q is added to
        result 'hull'
        '''
        p = q
 
        # While we don't come to first point
        if(p == l):
            break
    convex_hull = []
    # Print Result
    for each in hull:
        point=[]
        # print(points[each].x, points[each].y)
        point.append(points[each].x)
        point.append(points[each].y)
        convex_hull.append(point)
    return convex_hull


#graham algorithm 
# ---------------------------

 
# A class used to store the x and y coordinates of points
class Point:
    def __init__(self, x = None, y = None):
        self.x = x
        self.y = y
 
# A global point needed for sorting points with reference
# to the first point
p0 = Point(0, 0)
 
# A utility function to find next to top in a stack
def nextToTop(S):
    return S[-2]
 
# A utility function to return square of distance
# between p1 and p2
def distSq(p1, p2):
    return ((p1.x - p2.x) * (p1.x - p2.x) +
            (p1.y - p2.y) * (p1.y - p2.y))
 
# To find orientation of ordered triplet (p, q, r).
# The function returns following values
# 0 --> p, q and r are collinear
# 1 --> Clockwise
# 2 --> Counterclockwise
def orientation_(p, q, r):
    val = ((q.y - p.y) * (r.x - q.x) -
           (q.x - p.x) * (r.y - q.y))
    if val == 0:
        return 0  # collinear
    elif val > 0:
        return 1  # clock wise
    else:
        return 2  # counterclock wise
 
# A function used by cmp_to_key function to sort an array of
# points with respect to the first point
def compare(p1, p2):
   
    # Find orientation
    o = orientation_(p0, p1, p2)
    if o == 0:
        if distSq(p0, p2) >= distSq(p0, p1):
            return -1
        else:
            return 1
    else:
        if o == 2:
            return -1
        else:
            return 1
 
# Prints convex hull of a set of n points.

def get_slope(p1, p2):
    if p1[0] == p2[0]:
        return float('inf')
    else:
        return 1.0*(p1[1]-p2[1])/(p1[0]-p2[0])

def get_cross_product(p1,p2,p3):
    return ((p2[0] - p1[0])*(p3[1] - p1[1])) - ((p2[1] - p1[1])*(p3[0] - p1[0]))

def graham_convexHull(points):
    hull = []
    points.sort(key=lambda x:[x[0],x[1]])
    start = points.pop(0)
    
    points.sort(key=lambda p: (get_slope(p,start), -p[1],p[0]))
    hull = [start]

    for p in points:
        hull.append(p)
        while len(hull) > 2 and get_cross_product(hull[-3], hull[-2], hull[-1]) < 0:
            hull.pop(-2)
    
    while len(hull) > 2 and get_cross_product(hull[-3],hull[-2],hull[-1]) < 0:
        hull.pop(-2)
        
    return hull
        


# Monotone_chain algoritgm 
def Monotone_chain_convex_hull(points):


    # Sort the points lexicographically (tuples are compared lexicographically).
    # Remove duplicates to detect the case we have just one unique point.
    points = sorted(set(points))

    # Boring case: no points or a single point, possibly repeated multiple times.
    if len(points) <= 1:
        return points

    # 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
    # Returns a positive value, if OAB makes a counter-clockwise turn,
    # negative for clockwise turn, and zero if the points are collinear.
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    # Build lower hull 
    lower = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    # Build upper hull
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    # Concatenation of the lower and upper hulls gives the convex hull.
    # Last point of each list is omitted because it is repeated at the beginning of the other list. 
    return lower[:-1] + upper[:-1]


def main():

    points = []
    jarvis_x = []
    jarvis_y = []
    graham_x = []
    graham_y = []
    jarvis_time = []
    graham_time = []
    Monotone_time = []
    jarvis_space = []
    graham_space = []
    Monotone_space = []
    # points.append((0, 3))
    # points.append((2, 2))
    # points.append((1, 1))
    # points.append((2, 1))
    # points.append((3, 0))
    # points.append((0, 0))
    # points.append((3, 3))
    indexs = []
        
#----------------------------------------------------------------------- 

    for j in range(20,1000,20):
        indexs.append(j)
        points = [(random.randint(0,100),random.randint(0,100)) for i in range(j)]
        
        self_points = []
        for point in points:
            self_points.append(Point(point[0], point[1]))
        
        # for i in range(len(points)):
        #     plt.scatter(points[i][0],points[i][1])
        # print(type(points))
        # print(points)
        
        n = len(points)
        
        init = memory_profiler.memory_usage()
        start = time.time()
        jarvis = jarvis_convexHull(self_points,n)
        end = time.time()
        finish = memory_profiler.memory_usage()
        jarvis_time.append(end-start)
        jarvis_space.append(finish[0]-init[0])
        
                
        self_points = []
        for point in points:
            self_points.append(Point(point[0], point[1]))
        
        init = memory_profiler.memory_usage()
        start = time.time()
        graham = graham_convexHull(points)
        end = time.time()
        finish = memory_profiler.memory_usage()
        graham_time.append(end-start)
        graham_space.append(finish[0]-init[0])
        
        init = memory_profiler.memory_usage()
        start = time.time()
        Monotone_chain = Monotone_chain_convex_hull(points)
        end = time.time()
        finish = memory_profiler.memory_usage()
        Monotone_time.append(end-start)
        Monotone_space.append(finish[0]-init[0])
        
        
        if(j == 20):
            plt.subplot(3,1,1)
                
            for i in range(len(points)):
                plt.scatter(points[i][0],points[i][1])
            for i in range(len(jarvis)):
                jarvis_x.append(jarvis[i][0])
                jarvis_y.append(jarvis[i][1])
            jarvis_x.append(jarvis_x[0])
            jarvis_y.append(jarvis_y[0])
            plt.plot(jarvis_x,jarvis_y,label='Jarvis')
            plt.legend()
        
            print('this is ultimate answer by jarvis algorithm')
            print(jarvis) 
            print()
        

        
            print('this is ultimate answer by graham algorithm')
            print(graham)
            print()
            graham_x = []
            graham_y = []
            plt.subplot(3,1,2)
                
            for i in range(len(points)):
                plt.scatter(points[i][0],points[i][1])
            for i in range(len(graham)):
                graham_x.append(graham[i][0])
                graham_y.append(graham[i][1])
            graham_x.append(graham_x[0])
            graham_y.append(graham_y[0])
            plt.plot(graham_x,graham_y,label='Graham')
            plt.legend()
        
            print('this is ultimate answer by graham algorithm')
            print(Monotone_chain)
            print()
            
            Monotone_chain_x = []
            Monotone_chain_y = []
            plt.subplot(3,1,3)
                
            for i in range(len(points)):
                plt.scatter(points[i][0],points[i][1])
            for i in range(len(graham)):
                Monotone_chain_x.append(graham[i][0])
                Monotone_chain_y.append(graham[i][1])
            Monotone_chain_x.append(Monotone_chain_x[0])
            Monotone_chain_y.append(Monotone_chain_y[0])
            plt.plot(Monotone_chain_x,Monotone_chain_y,label='Monotone_chain')
            plt.legend()
    number = len(jarvis_time)
    index = range(0 , number)
    plt.figure(figsize=(20, 8))
    plt.subplot(2,1,1)

    plt.xticks(range(0,len(indexs)), indexs)
    plt.plot(index , jarvis_time , label = 'jarvis_time')
    plt.plot(index , graham_time , label = 'graham_time')
    plt.plot(index , Monotone_time , label = 'Monotone_time')
    plt.legend()
    
    
    plt.subplot(2,1,2)
    plt.xticks(range(0,len(indexs)), indexs)
    plt.plot(index , jarvis_space , label = 'jarvis_space')
    plt.plot(index , graham_space , label = 'graham_space')
    plt.plot(index , Monotone_space , label = 'Monotone_space')    

    plt.legend()  
    
if __name__ == '__main__':
    main()
    
