import numpy as np
import math
import heapq
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from scipy.spatial import Voronoi
from shapely.geometry import Polygon
from turfpy.measurement import boolean_point_in_polygon
from geojson import Point, Polygon as Poly, Feature



#----ボロノイ領域；頂点，頂点数----#
# ref：https://qiita.com/supbon2/items/30e0cb49c9338e721b8c
def get_vor(bnd, uav_loc):
    
    vor_poly = []
    all_loc = np.concatenate([uav_loc, np.array([[100000, 100000], [100000, -100000], [-100000, 0]])])
    vor = Voronoi(all_loc)
    bnd_poly = Polygon(bnd)
    
    for i in range(len(all_loc) - 3):
        p = [vor.vertices[v] for v in vor.regions[vor.point_region[i]]]
        i_cell = bnd_poly.intersection(Polygon(p))
        vor_poly.append(list(i_cell.exterior.coords[:-1]))
    
    vor_v_num = [len(n) for n in vor_poly]
    return vor_poly, vor_v_num



#----三角形；重心----#
def get_tricnt(p1, p2, p3):
    
    cx = (p1[0] + p2[0] + p3[0]) / 3
    cy = (p1[1] + p2[1] + p3[1]) / 3
    return (cx, cy)



#-----三角形；面積----#
def get_triarea(p1, p2, p3):

    a = ((p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0])) / 2
    return a



#----ボロノイ領域；重心----#
def get_cnt(vor_poly):
    
    vor_cnt = []
    
    for i in range(n):
        p1 = []; p2 = []
        q1 = []; q2 = []; q3 = []
        for j in range(vor_v_num[i] - 2):
            p1 = get_triarea(vor_poly[i][0], vor_poly[i][j+1], vor_poly[i][j+2])
            p2 = get_tricnt(vor_poly[i][0], vor_poly[i][j+1], vor_poly[i][j+2])
            q1.append(p1 * p2[0])
            q2.append(p1 * p2[1])
            q3.append(p1)    
        cx = np.sum(q1) / np.sum(q3)
        cy = np.sum(q2) / np.sum(q3)
        vor_cnt.append((cx, cy))
    
    return vor_cnt



#----ボロノイ領域；面積----#
def get_area(vor_poly):

    vor_area = []
    
    for i in range(n):
        p1 = []; p2 = []
        for j in range(vor_v_num[i] - 2):
            p1 = get_triarea(vor_poly[i][0], vor_poly[i][j+1], vor_poly[i][j+2])
            p2.append(abs(p1))
        vor_area.append(np.sum(p2))

    return vor_area



#----ボロノイ境界；長さ----#
def get_bndlen(vor_poly):
    
    vor_bndlen = []
    
    for i in range(n):
        p = []
        for j in range(vor_v_num[i]):
            if j < vor_v_num[i] - 1:
                l = math.sqrt((vor_poly[i][j+1][0] - vor_poly[i][j][0]) ** 2 + (vor_poly[i][j+1][1] - vor_poly[i][j][1]) ** 2)
            elif j == vor_v_num[i] - 1:
                l = math.sqrt((vor_poly[i][0][0] - vor_poly[i][j][0]) ** 2 + (vor_poly[i][0][1] - vor_poly[i][j][1]) ** 2)
            p.append(l)              
        vor_bndlen.append(p)
        
    return vor_bndlen



#----隣接するボロノイ領域(UAV)；番号，個数----#
def get_adj(vor_poly):

    vor_adj = []

    for i in range(n):
        p = []
        for j in range(vor_v_num[i]):
            for k in range(n):
                if i != k and vor_poly[i][j] in vor_poly[k] and k not in p:
                   p.append(k)
                   p.sort()
        vor_adj.append(p)
    
    vor_adj_num = [len(n) for n in vor_adj]
    
    return vor_adj, vor_adj_num



#----隣接する頂点；番号----#
def get_adjv(vor_poly):
    
    vor_adjv = []
    
    for i in range(n):
        p = []      
        for k in range(vor_adj_num[i]):
            q = []
            for j in range(vor_v_num[i]):
                if vor_poly[i][j] in vor_poly[vor_adj[i][k]]:
                    q.append(j)
            p.append(tuple(q))
        vor_adjv.append(p)
    
    return vor_adjv



#----2直線；交点----#
def get_cross(p1, p2, q1, q2):

    px = p2[0] - p1[0]; py = p2[1] - p1[1]
    qx = q2[0] - q1[0]; qy = q2[1] - q1[1]
    r1 = px * qy - qx * py
    
    if r1 == 0:
        return None
    
    r2 = ((q1[0] - p1[0]), (q1[1] - p1[1]))
    r3 = ((q2[1] - q1[1]) * r2[0] - (q2[0] - q1[0]) * r2[1]) / r1
    r4 = ((p2[0] - p1[0]) * r3, (p2[1] - p1[1]) * r3)

    return (p1[0] + r4[0]), (p1[1] + r4[1])



#----Waypoint(WP)，隣接するWP；番号，個数----#
def get_wp(vor_poly):

    vor_wp = []
    vor_wp_m = []
    
    r0 = []; r1 = []; r2 = []

    for i in range(n):
        p1 = []; p2 = []
        for j in range(vor_adj_num[i]):
            for k in range(len(vor_adjv[i][j])-1):
                cp = get_cross(vor_poly[i][vor_adjv[i][j][k]], vor_poly[i][vor_adjv[i][j][k+1]], vor_cnt[i], vor_cnt[vor_adj[i][j]])
                q1 = math.sqrt((vor_poly[i][vor_adjv[i][j][k+1]][0] - vor_poly[i][vor_adjv[i][j][k]][0]) ** 2 + (vor_poly[i][vor_adjv[i][j][k+1]][1] - vor_poly[i][vor_adjv[i][j][k]][1]) ** 2)
                q2 = math.sqrt((vor_poly[i][vor_adjv[i][j][k]][0] - cp[0]) ** 2 + (vor_poly[i][vor_adjv[i][j][k]][1] - cp[1]) ** 2)
                q3 = math.sqrt((vor_poly[i][vor_adjv[i][j][k+1]][0] - cp[0]) ** 2 + (vor_poly[i][vor_adjv[i][j][k+1]][1] - cp[1]) ** 2)
                if q1 > q2 and q1 > q3:
                    p1.append(cp)
                    p2.append(vor_adj[i][j])
        r1.append(p1)
        r2.append(p2)

    r1_num = [len(n) for n in r1]
    for i in range(n):
        for j in range(r1_num[i]):
            if (int(r1[i][j][0]), int(r1[i][j][1])) not in r0:
                r0.append((int(r1[i][j][0]), int(r1[i][j][1])))
                vor_wp.append(r1[i][j])
                vor_wp_m.append((i, r2[i][j]))
    
    vor_wp_num = len(vor_wp)

    return vor_wp, vor_wp_m, vor_wp_num



#----同一直線上；判定----#
def on_line(p1, p2, p3):

    if p1[0] == p3[0] == p2[0]:
        return True
    elif p1[1] == p3[1] == p2[1]:
        return True
    elif p1[0] - p3[0] != 0:
        q1 = (p1[1] - p3[1]) / (p1[0] - p3[0])
        q2 = p3[1] - (q1 * p1[0])
        if q1 * p2[0] + q2 == p2[1]:
            return True
        else:
            return False
    else:
        return False



#----UAV離脱；提案手法(1)----#
def prop_break1(uav_loc, vor_poly, uav_break):
 
    p = []
    vor_break = []
    uav_break_loc = uav_loc[uav_break]

    for j in range(vor_adj_num[uav_break]):
        p.append(vor_area[vor_adj[uav_break][j]])
    p1 = p.index(min(p))
    p2 = vor_adj[uav_break][p1] # 最小面積uav

    for j in range(vor_v_num[uav_break]):
        vor_break.append(vor_poly[uav_break][j])
    
    for j in range(vor_adj_num[p2]):
        if vor_adj[p2][j] == uav_break:
            q = j

    q1 = vor_break[vor_adjv[uav_break][p1][1]]
    for j in range(vor_v_num[uav_break]):
        if vor_break[0] != q1:
            vor_break.append(vor_break[0])
            del vor_break[0]     
        else:
            break              
    
    q2 = vor_adjv[p2][q][0]
    del vor_poly[p2][vor_adjv[p2][q][1]]
    del vor_poly[p2][vor_adjv[p2][q][0]]   
    for i in range(vor_v_num[uav_break]):
        vor_poly[p2].insert(q2, vor_break[i])
        q2 = q2 + 1

    vor_v_num1 = len(vor_poly[p2])
    r = []
    
    for j in range(vor_v_num1):
        if j < vor_v_num1-2:
            tf = on_line(vor_poly[p2][j], vor_poly[p2][j+1], vor_poly[p2][j+2])
            if tf == True:
                r.append(j+1)
        elif j == vor_v_num1-2:
            tf = on_line(vor_poly[p2][j], vor_poly[p2][j+1], vor_poly[p2][0])
            if tf == True:
                r.append(j+1)
        elif j == vor_v_num1-1: 
            tf = on_line(vor_poly[p2][j], vor_poly[p2][0], vor_poly[p2][1])
            if tf == True:
                r.append(0)
    r.sort(reverse=True)

    bnd1 = [(0.0, 0.0), (1000.0, 0.0), (1000.0, 1000.0), (0.0, 1000.0)]
    if len(r) > 0:
        for i in range(len(r)):
            if vor_poly[p2][r[i]] not in bnd1:
                del vor_poly[p2][r[i]]

    uav_loc = np.delete(uav_loc, uav_break, 0)
    n = len(uav_loc)
    del vor_poly[uav_break]
    vor_v_num2 = [len(n) for n in vor_poly]

    return n, uav_loc, uav_break_loc, vor_poly, vor_v_num2



#----UAV離脱；提案手法(2)----#
def prop_break2(uav_loc, vor_poly, uav_break):
    
    p = []
    vor_break = []
    uav_break_loc = uav_loc[uav_break]
    
    for j in range(vor_adj_num[uav_break]):
        l = math.sqrt((uav_loc[vor_adj[uav_break][j]][0] - uav_loc[uav_break][0]) ** 2 + (uav_loc[vor_adj[uav_break][j]][1] - uav_loc[uav_break][1]) ** 2)
        p.append(l)
    p1 = p.index(min(p))
    p2 = vor_adj[uav_break][p1]

    for j in range(vor_v_num[uav_break]):
        vor_break.append(vor_poly[uav_break][j])
    
    for j in range(vor_adj_num[p2]):
        if vor_adj[p2][j] == uav_break:
            q = j

    q1 = vor_break[vor_adjv[uav_break][p1][1]]
    for j in range(vor_v_num[uav_break]):
        if vor_break[0] != q1:
            vor_break.append(vor_break[0])
            del vor_break[0]     
        else:
            break              
    
    q2 = vor_adjv[p2][q][0]
    del vor_poly[p2][vor_adjv[p2][q][1]]
    del vor_poly[p2][vor_adjv[p2][q][0]]   
    for i in range(vor_v_num[uav_break]):
        vor_poly[p2].insert(q2, vor_break[i])
        q2 = q2 + 1

    vor_v_num1 = len(vor_poly[p2])
    r = []
    
    for j in range(vor_v_num1):
        if j < vor_v_num1-2:
            tf = on_line(vor_poly[p2][j], vor_poly[p2][j+1], vor_poly[p2][j+2])
            if tf == True:
                r.append(j+1)
        elif j == vor_v_num1-2:
            tf = on_line(vor_poly[p2][j], vor_poly[p2][j+1], vor_poly[p2][0])
            if tf == True:
                r.append(j+1)
        elif j == vor_v_num1-1: 
            tf = on_line(vor_poly[p2][j], vor_poly[p2][0], vor_poly[p2][1])
            if tf == True:
                r.append(0)
    r.sort(reverse=True)

    bnd1 = [(0.0, 0.0), (1000.0, 0.0), (1000.0, 1000.0), (0.0, 1000.0)]
    if len(r) > 0:
        for i in range(len(r)):
            if vor_poly[p2][r[i]] not in bnd1:
                del vor_poly[p2][r[i]]

    uav_loc = np.delete(uav_loc, uav_break, 0)
    n = len(uav_loc)
    del vor_poly[uav_break]
    vor_v_num2 = [len(n) for n in vor_poly]

    return n, uav_loc, uav_break_loc, vor_poly, vor_v_num2



#----UAV離脱；提案手法(3)----#
def prop_break3(uav_loc, vor_poly, uav_break):
    
    p = []
    vor_break = []
    uav_break_loc = uav_loc[uav_break]
    
    for j in range(vor_adj_num[uav_break]):
        p0 = []
        for k in range(len(vor_adjv[uav_break][j])-1):
            l = vor_bndlen[uav_break][vor_adjv[uav_break][j][k]]
            p0.append(l)   
        p.append(max(p0))
    p1 = p.index(max(p))
    p2 = vor_adj[uav_break][p1]

    for j in range(vor_v_num[uav_break]):
        vor_break.append(vor_poly[uav_break][j])
    
    for j in range(vor_adj_num[p2]):
        if vor_adj[p2][j] == uav_break:
            q = j

    q1 = vor_break[vor_adjv[uav_break][p1][1]]
    for j in range(vor_v_num[uav_break]):
        if vor_break[0] != q1:
            vor_break.append(vor_break[0])
            del vor_break[0]     
        else:
            break              
    
    q2 = vor_adjv[p2][q][0]
    del vor_poly[p2][vor_adjv[p2][q][1]]
    del vor_poly[p2][vor_adjv[p2][q][0]]   
    for i in range(vor_v_num[uav_break]):
        vor_poly[p2].insert(q2, vor_break[i])
        q2 = q2 + 1

    vor_v_num1 = len(vor_poly[p2])
    r = []
    
    for j in range(vor_v_num1):
        if j < vor_v_num1-2:
            tf = on_line(vor_poly[p2][j], vor_poly[p2][j+1], vor_poly[p2][j+2])
            if tf == True:
                r.append(j+1)
        elif j == vor_v_num1-2:
            tf = on_line(vor_poly[p2][j], vor_poly[p2][j+1], vor_poly[p2][0])
            if tf == True:
                r.append(j+1)
        elif j == vor_v_num1-1: 
            tf = on_line(vor_poly[p2][j], vor_poly[p2][0], vor_poly[p2][1])
            if tf == True:
                r.append(0)
    r.sort(reverse=True)

    bnd1 = [(0.0, 0.0), (1000.0, 0.0), (1000.0, 1000.0), (0.0, 1000.0)]
    if len(r) > 0:
        for i in range(len(r)):
            if vor_poly[p2][r[i]] not in bnd1:
                del vor_poly[p2][r[i]]

    uav_loc = np.delete(uav_loc, uav_break, 0)
    n = len(uav_loc)
    del vor_poly[uav_break]
    vor_v_num2 = [len(n) for n in vor_poly]

    return n, uav_loc, uav_break_loc, vor_poly, vor_v_num2



#----垂直二等分線----#
def get_perpbis(p1, p2):
    
    p = []
    q1 = (p2[1] - p1[1]) / (p2[0] - p1[0])
    q2 = -1 / q1
    
    cx = (p1[0] + p2[0]) / 2 # 中点；x座標
    cy = (p1[1] + p2[1]) / 2 # 中点；y座標
    r = cy - q2 * cx
    y = q2 * (cx + 1) + r

    p.append((cx, cy))
    p.append((cx + 1, y))

    return p



#----UAV参加；提案手法----#
def prop_join(n, uav_loc, vor_poly, uav_join):
  
    pnt = Feature(geometry = Point((uav_join[0], uav_join[1])))
    vor = []
    vor_join = []
    
    for i in range(n):
        vor.append(vor_poly[i])
        polygon = Poly(vor)
        if boolean_point_in_polygon(pnt, polygon) == True:
            minuav = i
        vor.clear()

    pb = get_perpbis(uav_loc[minuav], uav_join)
    for j in range(vor_v_num[minuav]):
        vor_join.append(vor_poly[minuav][j])
    
    cut_v = []
    p = []; q = []

    for j in range(vor_v_num[minuav]):
        if j < vor_v_num[minuav] - 1:
            cp = get_cross(vor_poly[minuav][j], vor_poly[minuav][j+1], pb[0], pb[1])
            p1 = math.sqrt((vor_poly[minuav][j+1][0] - vor_poly[minuav][j][0]) ** 2 + (vor_poly[minuav][j+1][1] - vor_poly[minuav][j][1]) ** 2)
            p2 = math.sqrt((vor_poly[minuav][j][0] - cp[0]) ** 2 + (vor_poly[minuav][j][1] - cp[1]) ** 2)
            p3 = math.sqrt((vor_poly[minuav][j+1][0] - cp[0]) ** 2 + (vor_poly[minuav][j+1][1] - cp[1]) ** 2)
        elif j == vor_v_num[minuav] - 1:
            cp = get_cross(vor_poly[minuav][j], vor_poly[minuav][0], pb[0], pb[1])
            p1 = math.sqrt((vor_poly[minuav][0][0] - vor_poly[minuav][j][0]) ** 2 + (vor_poly[minuav][0][1] - vor_poly[minuav][j][1]) ** 2)
            p2 = math.sqrt((vor_poly[minuav][j][0] - cp[0]) ** 2 + (vor_poly[minuav][j][1] - cp[1]) ** 2)
            p3 = math.sqrt((vor_poly[minuav][0][0] - cp[0]) ** 2 + (vor_poly[minuav][0][1] - cp[1]) ** 2)
        if p1 > p2 and p1 > p3:
            cut_v.append(cp)
            p.append(j)
    
    q1 = []; q2 = []; q3 = []; q4 = []; q5 = []

    if p[0] < vor_v_num[minuav] - 1:
        q.append((p[0], p[0]+1))
    elif p[0] == vor_v_num[minuav] - 1:
        q.append((0, p[0]+1))

    if p[1] < vor_v_num[minuav] - 1:
        q.append((p[1], p[1]+1))
    elif p[1] == vor_v_num[minuav] - 1:
        q.append((0, p[1]+1))

    for j in range(vor_adj_num[minuav]): 
        if q[0][0] in vor_adjv[minuav][j] and q[0][1] in vor_adjv[minuav][j]:  
            q1.append(vor_adj[minuav][j])
        if q[0][0] not in vor_adjv[minuav][j] or q[0][1] not in vor_adjv[minuav][j]:      
            q1.append(None)
        if q[1][0] in vor_adjv[minuav][j] and q[1][1] in vor_adjv[minuav][j]:  
            q2.append(vor_adj[minuav][j])
        if q[1][0] not in vor_adjv[minuav][j] or q[1][1] not in vor_adjv[minuav][j]:      
            q2.append(None)
    
    if all(v is None for v in q1):
        q3.append(None)
    else:
        for i in range(vor_adj_num[minuav]):
            if q1[i] != None:
                q3.append(q1[i])
    if all(v is None for v in q2):
        q3.append(None)
    else:
        for i in range(vor_adj_num[minuav]):
            if q2[i] != None:
                q3.append(q2[i])

    if q3[0] is not None:
        for i in range(vor_v_num[q3[0]]):
            if vor_poly[q3[0]][i] == vor_poly[minuav][q[0][0]]:
                q4.append(i)
            if vor_poly[q3[0]][i] == vor_poly[minuav][q[0][1]]:
                q4.append(i)
        q4.sort
        if q4[0] == 0 and q4[1] == vor_v_num[q3[0]]-1:
            vor_poly[q3[0]].append(cut_v[0])
        else:
            vor_poly[q3[0]].insert(q4[1], cut_v[0])

    if q3[1] is not None:
        for i in range(vor_v_num[q3[1]]):
            if vor_poly[q3[1]][i] == vor_poly[minuav][q[1][0]]:
                q5.append(i)
            if vor_poly[q3[1]][i] == vor_poly[minuav][q[1][1]]:
                q5.append(i)
        q5.sort
        if q5[0] == 0 and q5[1] == vor_v_num[q3[1]]-1:
            vor_poly[q3[1]].append(cut_v[1])
        else:
            vor_poly[q3[1]].insert(q5[1], cut_v[1])
    
    x1 = cut_v[1][0] - cut_v[0][0]
    y1 = cut_v[1][1] - cut_v[0][1]
    x2 = uav_loc[minuav][0] - cut_v[0][0]
    y2 = uav_loc[minuav][1] - cut_v[0][1]
    lr = x1 * y2 - y1 * x2
    
    if lr < 0:
        r1 = []
        r2 = []
        for j in range(vor_v_num[minuav]):
            if j <= p[0] or j > p[1]:
                r1.append(j)
            elif p[0] < j <= p[1]:
                r2.append(j)

        r1.sort(reverse=True)
        r2.sort(reverse=True)

        for i in range(len(r1)):
            del vor_join[r1[i]]
        for i in range(len(r2)):
            del vor_poly[minuav][r2[i]]
    
        vor_poly[minuav].insert(p[0]+1, cut_v[0])
        vor_poly[minuav].insert(p[0]+2, cut_v[1])
        vor_join.insert(0, cut_v[0])
        vor_join.append(cut_v[1])
    
    else:
        r1 = []
        r2 = []
        for j in range(vor_v_num[minuav]):
            if j <= p[0] or j > p[1]:
                r1.append(j)
            elif p[0] < j <= p[1]:
                r2.append(j)
        
        r1.sort(reverse=True)
        r2.sort(reverse=True)

        for i in range(len(r1)):
            del vor_poly[minuav][r1[i]]
        for i in range(len(r2)):
            del vor_join[r2[i]]
    
        vor_poly[minuav].insert(0, cut_v[0])
        vor_poly[minuav].append(cut_v[1])
        vor_join.insert(p[0]+1, cut_v[0])
        vor_join.insert(p[0]+2, cut_v[1])
    
    uav_loc = np.vstack((uav_loc, uav_join))
    n1 = len(uav_loc)
    vor_poly.append(vor_join)
    vor_v_num1 = [len(n) for n in vor_poly]
    
    return n1, uav_loc, vor_poly, vor_v_num1



#----UAV離脱；ボロノイ領域を再計算する手法----#
def cover_break(uav_loc, vor_poly, uav_break):

    uav_break_loc = uav_loc[uav_break]
    uav_loc = np.delete(uav_loc, uav_break, 0)
    n = len(uav_loc)
    vor_poly, vor_v_num = get_vor(bnd, uav_loc) 
    
    return n, uav_loc, uav_break_loc, vor_poly, vor_v_num



#----UAV参加；ボロノイ領域を再計算する手法----#
def cover_join(uav_loc, vor_poly, uav_join):

    uav_loc = np.vstack((uav_loc, uav_join))
    n = len(uav_loc)
    vor_poly, vor_v_num = get_vor(bnd, uav_loc) 
    
    return n, uav_loc, vor_poly, vor_v_num



#----送信元・宛先を含むボロノイ領域；番号----#
def get_be(vor_poly): 
  
    pnt1 = Feature(geometry = Point((src[0], src[1])))
    pnt2 = Feature(geometry = Point((dest[0], dest[1])))
    vor = []
  
    for i in range(n):
        vor.append(vor_poly[i])
        polygon = Poly(vor)
        if boolean_point_in_polygon(pnt1, polygon) == True:
            begin = i
        if boolean_point_in_polygon(pnt2, polygon) == True:
            end = i
        vor.clear() 
    
    return begin, end



#----目的地数，経路----#
def get_edge(vor_poly):

    edge = []
    p = []

    for i in range(n):
        if i == begin:
            for j in range(vor_wp_num):
                if i in vor_wp_m[j]:
                    l = math.sqrt((vor_cnt[begin][0] - vor_wp[j][0]) ** 2 + (vor_cnt[begin][1] - vor_wp[j][1]) ** 2)
                    p.append([0, j + 1, l])        
        elif i == end:
            for j in range(vor_wp_num):
                if i in vor_wp_m[j]:
                    l = math.sqrt((vor_cnt[end][0] - vor_wp[j][0]) ** 2 + (vor_cnt[end][1] - vor_wp[j][1]) ** 2)
                    p.append([j + 1, vor_wp_num + 1, l])
        else:
            for j in range(vor_wp_num):
                if i in vor_wp_m[j]:
                    for k in range(vor_wp_num - j - 1):
                        if k + j + 1 < vor_wp_num and i in vor_wp_m[k + j + 1]:
                            l = math.sqrt((vor_wp[j][0] - vor_wp[k + j + 1][0]) ** 2 + (vor_wp[j][1] - vor_wp[k + j + 1][1]) ** 2)
                            p.append([j+1, k+j+2, l])
                            p.append([k+j+2, j+1, l])

    for i in range(vor_wp_num + 1):
        q = []
        for j in range(len(p)):
            if p[j][0] == i:
                q.append([p[j][1], p[j][2]])
        edge.append(q)

    node_num = vor_wp_num + 2

    return node_num, edge



#----ダイクストラ法----#
# ref：https://qiita.com/Yuya-Shimizu/items/eefdc6f854534e90c988
def dijkstra(edge, node_num):
    
    node = [float('inf')] * node_num
    node[0] = 0

    node_name = []
    heapq.heappush(node_name, [0, [0]])

    while len(node_name) > 0:
        _, min_point = heapq.heappop(node_name)
        last = min_point[-1]
        
        if last == node_num - 1:
            return min_point, node

        for factor in edge[last]:
            goal = factor[0]
            cost  = factor[1]

            if node[last] + cost < node[goal]:
                node[goal] = node[last] + cost
                heapq.heappush(node_name, [node[last] + cost, min_point + [goal]])
    
    return []



#----通信遅延時間，ホップ数----#
def get_result():

    if begin == end:
        hop = 1
        t0 = vor_area[begin] / s / 2
        t1 = 0
        t2 = (math.sqrt((vor_cnt[end][0] - dest[0]) ** 2 + (vor_cnt[end][1] - dest[1]) ** 2)) / v
        time = t0 + t1 + t2
    
    else:
        node_num, edge = get_edge(vor_poly)
        root, cost = dijkstra(edge, node_num) 
        hop = len(root) + 1
        t0 = vor_area[begin] / s / 2
        t1 = cost[-1] / v
        t2 = (math.sqrt((vor_cnt[end][0] - dest[0]) ** 2 + (vor_cnt[end][1] - dest[1]) ** 2)) / v
        time = t0 + t1 + t2

    return t0, t1, t2, time, hop



#----ボロノイ図を保存----#
def get_fig():
    
    p1 = []; p2 = []; p3 = []; p4 = []
    
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
    
    vor = PolyCollection(vor_poly, edgecolor="black", facecolors="None", linewidth = 1.0)
    ax.add_collection(vor)
        
    for i in range(n):
        p1.append(vor_cnt[i][0])
        p2.append(vor_cnt[i][1])
    for i in range(vor_wp_num):
        p3.append(vor_wp[i][0])
        p4.append(vor_wp[i][1])   
    ax.scatter(p1, p2, c = 'firebrick', s = 100, marker = "p", clip_on = False, label = 'Voronoi centroid')
    ax.scatter(p3, p4, c = 'mediumblue', s = 50, marker = "D", clip_on = False, label = 'WP')
    
    ax.legend(loc='upper left')
    fig.savefig('voronoi.png')
    
    return 0



bnd = np.array([[0, 0], [1000, 0], [1000, 1000], [0, 1000]]) # Simulation area
v = 5 # Movement speed of UAV
r = 30 # Communication range of UAV
s = math.pi * r * r + 2 * r * v

n = 5 # Number of UAVs
src = 1000 * np.random.rand(2) # Location of source node
dest = 1000 * np.random.rand(2) # Location of destination node
uav_loc = 1000 * np.random.rand(n, 2) # Initial location

vor_poly, vor_v_num = get_vor(bnd, uav_loc) 
vor_cnt = get_cnt(vor_poly)
vor_area = get_area(vor_poly)
vor_bndlen = get_bndlen(vor_poly)
vor_adj, vor_adj_num = get_adj(vor_poly) 
vor_adjv = get_adjv(vor_poly)
vor_wp, vor_wp_m, vor_wp_num = get_wp(vor_poly)


#----提案手法；UAV離脱----#
uav_break = np.random.randint(0, n) # B_UAV
n, uav_loc, uav_break_loc, vor_poly, vor_v_num = prop_break1(uav_loc, vor_poly, uav_break) 
#n, uav_loc, uav_break_loc, vor_poly, vor_v_num = prop_break2(uav_loc, vor_poly, uav_break)
#n, uav_loc, uav_break_loc, vor_poly, vor_v_num = prop_break3(uav_loc, vor_poly, uav_break)
vor_cnt = get_cnt(vor_poly)
vor_area = get_area(vor_poly)
vor_bndlen = get_bndlen(vor_poly)
vor_adj, vor_adj_num = get_adj(vor_poly) 
vor_adjv = get_adjv(vor_poly)
vor_wp, vor_wp_m, vor_wp_num = get_wp(vor_poly)


#----提案手法；UAV参加----#
uav_join = 1000 * np.random.rand(2) # Location of J_UAV
n, uav_loc, vor_poly, vor_v_num = prop_join(n, uav_loc, vor_poly, uav_join)
vor_cnt = get_cnt(vor_poly)
vor_area = get_area(vor_poly)
vor_bndlen = get_bndlen(vor_poly)
vor_adj, vor_adj_num = get_adj(vor_poly) 
vor_adjv = get_adjv(vor_poly)
vor_wp, vor_wp_m, vor_wp_num = get_wp(vor_poly)


#----従来手法；UAV離脱----#
uav_break = np.random.randint(0, n) # B_UAV
n, uav_loc, uav_break_loc, vor_poly, vor_v_num = cover_break(uav_loc, vor_poly, uav_break)
vor_cnt = get_cnt(vor_poly)
vor_area = get_area(vor_poly)
vor_bndlen = get_bndlen(vor_poly)
vor_adj, vor_adj_num = get_adj(vor_poly) 
vor_adjv = get_adjv(vor_poly)
vor_wp, vor_wp_m, vor_wp_num = get_wp(vor_poly)


#----従来手法；UAV参加----#
uav_join = 1000 * np.random.rand(2) # Location of J_UAV
n, uav_loc, vor_poly, vor_v_num = cover_join(uav_loc, vor_poly, uav_join)
vor_cnt = get_cnt(vor_poly)
vor_area = get_area(vor_poly)
vor_bndlen = get_bndlen(vor_poly)
vor_adj, vor_adj_num = get_adj(vor_poly) 
vor_adjv = get_adjv(vor_poly)
vor_wp, vor_wp_m, vor_wp_num = get_wp(vor_poly)


begin, end = get_be(vor_poly)
t0, t1, t2, time, hop = get_result()


#get_fig()
#print(vor_poly)
#print(vor_v_num)
#print(vor_cnt)
#print(vor_area)
#print(vor_bndlen)
#print(vor_adj)
#print(vor_adj_num)
#print(vor_adjv)
#print(vor_wp)
#print(vor_wp_m)
#print(vor_wp_num)

#print(t0)
#print(t1)
#print(t2)
#print(time)
#print(hop)
#print(begin)
#print(end)