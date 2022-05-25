import numpy as np
import math
import matplotlib.pyplot as plt



#----UAV離脱；従来手法----#
def conv_break(uav_loc, uav_break):

    uav_break_fig = uav_loc[uav_break]
    uav_loc = np.delete(uav_loc, uav_break, 0)
    n = len(uav_loc)
    
    return n, uav_loc, uav_break_fig



#----UAV参加；従来手法----#
def conv_join(uav_loc, uav_join):

    uav_loc = np.vstack((uav_loc, uav_join))
    n = len(uav_loc)

    return n, uav_loc



#----通信遅延時間，ホップ数----#
def get_result():

    bnd = [(0.0, 0.0), (0.0, 1000.0), (1000.0, 0.0), (1000.0, 1000.0)]
    hop = (n - 1) / 2
    q = []
     
    for i in range(4):
        l = math.sqrt((bnd[i][0] - dest[0]) ** 2 + (bnd[i][1] - dest[1]) ** 2)
        q.append(l)
  
    t0 = 1000 * 1000 / s / n
    p = n / 100
    t1 = hop / p
    t2 = max(q) / 2
    time = t0 + t1 + t2

    return t0, t1, t2, time, hop



#----ボロノイ図を保存----#
def get_fig():
    
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.set_xlim(0, 1000)
    ax.set_ylim(0, 1000)
    ax.tick_params(direction = "in")
    ax.grid(linestyle=':')
    ax.set_axisbelow(True)
    
    ax.scatter(src[0], src[1], c = 'black', s = 120, marker = "X", clip_on = False, label = 'Source node')
    ax.scatter(dest[0], dest[1], c = 'firebrick', s = 120, marker = "X", clip_on = False, label = 'Destination node')
    ax.scatter(uav_loc[:,0], uav_loc[:,1], c = 'black', s = 100, marker = "+", clip_on = False, label = 'UAV')
    #ax.scatter(uav_join[0], uav_join[1], c = 'green', s = 100, marker = "+", clip_on = False, label = 'J_UAV')
    #ax.scatter(uav_break_loc[0], uav_break_loc[1], c = 'gray', s = 100, marker = "+", clip_on = False, label = 'B_UAV')
        
    ax.legend(loc='upper left')
    fig.savefig('rwp.png')
    
    return 0


v = 5 # Movement speed of UAV
r = 30 # Communication range of UAV
s = math.pi * r * r + 2 * r * v

n = 5 # Number of UAVs
src = 1000 * np.random.rand(2) # Location of source node
dest = 1000 * np.random.rand(2) # Location of destination node
uav_loc = 1000 * np.random.rand(n, 2) # Initial location


#----UAV離脱----#
uav_break = np.random.randint(0, n) # B_UAV
n, uav_loc, uav_break_fig = conv_break(uav_loc, uav_break)


#----UAV参加----#
uav_join = 1000 * np.random.rand(2) # Location of J_UAV
n, uav_loc = conv_join(uav_loc, uav_join)
    

t0, t1, t2, time, hop = get_result()


#get_fig()
#print(t0)
#print(t1)
#print(t2)
#print(time)
#print(hop)

