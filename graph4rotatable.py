import random
import math
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as patches
import numpy as np

# Input parameters
RoomWidth = 500.0     # room width in centi metres
RoomHeight = 500.0    # room height in centi metres
grid_size = 5     # centi-metres

MAX_COLUMN = math.ceil(RoomWidth/grid_size)
MAX_ROW = math.ceil(RoomWidth/grid_size)

allOBSvertex = [[(70,120),(70,300),(200,300),(200,220),(350,220),(350,120)],
                [(15,435),(70,490),(120,445),(65,390)],
                [(485,435),(430,490),(380,445),(435,390)]]

TxLocation = (250,500)
c1 = []
c2 = 0
cover2 = []
cov = []
T = []
Txcoverage = 0
for k in range(MAX_ROW):
    c1.append([0]*MAX_COLUMN)
    T.append([0]*MAX_COLUMN)
for x in range(MAX_COLUMN) :
  for y in range(MAX_ROW) :
       c1[x][y] = 0.0
       T[x][y] = 0.0
#mirror = (['R',250,10,150,0])  #wallSide,position,distance,mirror_width,angle

wall_axis = [[[0,0],[RoomWidth,0]],[[0,0],[0,RoomHeight]],[[0,RoomHeight],[RoomWidth,RoomHeight]],[[RoomWidth,RoomHeight],[RoomWidth,0]]]

class PixelPlane() :
    def __init__(self,MAX_ROW,MAX_COLUMN) :
        self.MAX_ROW = MAX_ROW
        self.MAX_COLUMN = MAX_COLUMN 
        self.pixel_intensity = []
        for k in range(self.MAX_ROW):
            self.pixel_intensity.append([0]*self.MAX_COLUMN)
        for x in range(self.MAX_COLUMN) :
            for y in range(self.MAX_ROW) :
                self.pixel_intensity[x][y] = 0.0


def intersect(x1,y1,x2,y2,x3,y3,x4,y4,debug=False):
    if debug:
        print("========= Start Intersect =========")
        print(x1,y1,x2,y2,x3,y3,x4,y4)
    denominator = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)
    if debug:
        print(denominator)
    if denominator!=0:
        Px = (x1*y2 - y1*x2)*(x3-x4) - (x1-x2)*(x3*y4-y3*x4)
        Px = Px/denominator
        Py = (x1*y2 - y1*x2)*(y3-y4) - (y1-y2)*(x3*y4-y3*x4)
        Py = Py/denominator
        if debug:
            print(Px,Py)
        onSegment1 = False
        onSegment2 = False

        alpha = 0.00001
        if ((Px>=x1-alpha and Px <= x2+alpha) or (Px>=x2-alpha and Px <= x1+alpha)) and ((Py>=y1-alpha and Py <= y2+alpha) or (Py>=y2-alpha and Py <= y1+alpha)):
                onSegment2 = True
        if ((Px>=x3-alpha and Px <= x4+alpha) or (Px>=x4-alpha and Px <= x3+alpha)) and ((Py>=y3-alpha and Py <= y4+alpha) or (Py>=y4-alpha and Py <= y3+alpha)):
                onSegment1 = True

        if debug:
            print(onSegment1,onSegment2)

        if onSegment1 == True and onSegment2 == True:
            return (True,Px,Py)
        else:
            return (False,Px,Py)
    else:
        return (False,'Parallel')

class Mirror():
    def __init__(self,wallSide,position,distance,mirror_width,angle):
        self.wallSide = wallSide
        self.position = position
        self.distance = distance
        self.mirror_width = mirror_width
        self.angle = angle
        mirror_loc = self.mirror_install(wallSide,position,distance,mirror_width,angle)
        self.mirror_loc = mirror_loc
        self.mirror_pixel_list = self.mirror_pixel_plane(mirror_loc[0],mirror_loc[1],mirror_loc[2],mirror_loc[3])
        self.project_pixel_list = self.pixel_mirror_projection()
        
        # ** TEST MULTIPLE RAYS ADD
        #self.project_pixel_list = self.shotgun(2) #Experiment (Multiple rays)
        
        self.mirror_reflect_points = self.mirror_reflect_points()
        self.tx_mirror_perp = self.perp_toMirror(TxLocation[0],TxLocation[1],mirror_loc[0],mirror_loc[1],mirror_loc[2],mirror_loc[3])
        ##print("self.tx_mirror_perp =",self.tx_mirror_perp)       
        all_mirror_reflect_points = self.mirror_reflect_points
        self.ray_to_object_list = self.ray_to_object()            
        self.all_ray_coverage_list = self.all_ray_coverage()
        #self.count_coverage()
        
    def mirror_install(self,wallSide,position,distance,mirror_width,angle): # (Wall,location,distance from wall,mirror length, angle)       
        if wallSide == 'Right' or wallSide == 'R':
            #print('R')
            #At angle = 0          
            x1 = RoomWidth - distance
            y1 = position - mirror_width/2
            x2 = RoomWidth - distance
            y2 = position + mirror_width/2
      
            #With angle
            x1 = x1 - np.cos(np.deg2rad(90+angle))*mirror_width/2
            y1 = y1 + ((mirror_width/2)-np.cos(np.deg2rad(angle))*mirror_width/2)
            x2 = x2 + np.cos(np.deg2rad(90+angle))*mirror_width/2
            y2 = y2 - ((mirror_width/2)-np.cos(np.deg2rad(angle))*mirror_width/2)
            
        elif wallSide == 'Left' or wallSide == 'L':
            #print('LEFT')
            x1 = distance
            y1 = position + mirror_width/2
            x2 = distance
            y2 = position - mirror_width/2

            #with angle
            x1 = x1 - np.cos(np.deg2rad(90-angle))*mirror_width/2
            y1 = y1 - ((mirror_width/2)-np.cos(np.deg2rad(angle))*mirror_width/2)
            x2 = x2 + np.cos(np.deg2rad(90-angle))*mirror_width/2
            y2 = y2 + ((mirror_width/2)-np.cos(np.deg2rad(angle))*mirror_width/2)

        elif wallSide == 'Top' or wallSide == 'T' :
            print('Top')
            x1 = position - mirror_width/2
            y1 = RoomHeight - distance
            x2 = position + mirror_width/2
            y2 = RoomHeight - distance
            
            #with angle
            x1 = x1 + ((mirror_width/2)-np.cos(np.deg2rad(angle))*mirror_width/2)
            y1 = y1 - np.cos(np.deg2rad(90-angle))*mirror_width/2
            x2 = x2 - ((mirror_width/2)-np.cos(np.deg2rad(angle))*mirror_width/2)
            y2 = y2 + np.cos(np.deg2rad(90-angle))*mirror_width/2

        elif wallSide == 'Bottom' or wallSide == 'B' : 
            print('Bottom')
            x1 = position - mirror_width/2
            y1 = distance
            x2 = position + mirror_width/2
            y2 = distance

            #with angle
            x1 = x1 + mirror_width/2 - (mirror_width/2*(np.cos(np.deg2rad(angle))))
            y1 = y1 - ((mirror_width/2)*np.cos(np.deg2rad(90-angle)))
            x2 = x2 - (mirror_width/2) + (mirror_width/2*(np.cos(np.deg2rad(angle))))
            y2 = y2 + (np.cos(np.deg2rad(90-angle)))*mirror_width/2

        #print(x1/grid_size,y1/grid_size,x2/grid_size,y2/grid_size)
        #print(x1,y1,x2,y2)
        return x1,y1,x2,y2
        
    def mirror_pixel_plane(self,x1,y1,x2,y2):
        mirror_plane = []
        for px in range(int(math.floor(min([x1,x2])/grid_size)),(int(math.ceil(max([x1,x2])/grid_size)+1))):
            for py in range(int(math.floor(min([y1,y2])/grid_size)),(int(math.ceil(max([y1,y2])/grid_size)+1))):
                pix_sides = [[px,py,px+1,py],[px,py,px,py+1],[px,py+1,px+1,py+1],[px+1,py,px+1,py+1]]
                ##print("pix sides =",pix_sides)
                pix_sides = np.multiply(pix_sides,grid_size)      
                for side in pix_sides:
                    mirror_on_grid = intersect(x1,y1,x2,y2,side[0],side[1],side[2],side[3])
                    if mirror_on_grid[0] == True:
                        if px< MAX_COLUMN and py < MAX_ROW:
                            mirror_plane.append([px,py])
                        break
        return mirror_plane

    def pixel_mirror_projection(self):
        project_pixel_list=[] #List for projection pixel position into exact mirror line
        #mid_pixel_list = [[x+0.5,y+0.5] for [x,y] in self.mirror_pixel_list]  #<-- If want all rays from mirror
        mid_pixel_list = [[x-0.5,y-0.5] for [x,y] in self.light_on_mirror(tx_id)] #<-- only light part can reflect
        #print("Mid pix =",mid_pixel_list)   
        for i in range(len(mid_pixel_list)):
            project_pixel_list.append(mid_pixel_list[i])
            if i+1 < len(mid_pixel_list):
                 project_pixel_list.append([np.mean([mid_pixel_list[i][0],mid_pixel_list[i+1][0]]),np.mean([mid_pixel_list[i][1],mid_pixel_list[i+1][1]])]) 
        #print('mirror_pixel_list',self.mirror_pixel_list)
        #print('project_pixel_list',project_pixel_list) 
        return project_pixel_list        
    
    def mirror_reflect_points(self):
        all_reflect_points = []
        pixel_to_loc = [[px*grid_size,py*grid_size] for [px,py] in self.project_pixel_list]
        for loc in pixel_to_loc:
            point = self.perp_toMirror(loc[0],loc[1],self.mirror_loc[0],self.mirror_loc[1],self.mirror_loc[2],self.mirror_loc[3])
            all_reflect_points.append(point)
        
        #print("all reflect points = ",all_reflect_points)
        return all_reflect_points
    
    def reference_point_mirror(self,x1,y1):
        reference_point_mirror = [x1-(self.tx_mirror_perp[0]-x1),y1-(self.tx_mirror_perp[1]-y1)]
        #print("ref point mirror = ",reference_point_mirror)        
        return reference_point_mirror
        
    def reference_point_reflect(self,x1,y1):
        reference_point_mirror = self.reference_point_mirror(x1,y1)
        reference_point_reflect = [reference_point_mirror[0]-(self.tx_mirror_perp[0]-TxLocation[0]),reference_point_mirror[1]+(TxLocation[1]-self.tx_mirror_perp[1])]   
        return reference_point_reflect
    
    def all_mirror_plane_to_wall(self):   #find intersect points between mirror plane and wall 
        x1 = self.mirror_loc[0]
        y1 = self.mirror_loc[1]
        x2 = self.mirror_loc[2]
        y2 = self.mirror_loc[3]
        all_mirror_plane_to_wall = []
        for wall in wall_axis:
            mirror_plane_to_wall = self.room_wall_intersect(np.array([x1,y1]),np.array([x2,y2]),np.array(wall[0]),np.array(wall[1]))
            if mirror_plane_to_wall[0]>=0 and mirror_plane_to_wall[0]<=RoomWidth and mirror_plane_to_wall[1]>=0 and mirror_plane_to_wall[1]<=RoomHeight:
               all_mirror_plane_to_wall.append(mirror_plane_to_wall)
        #print(all_mirror_plane_to_wall)
        return all_mirror_plane_to_wall
    
    def wall_position(self,x1,y1):       
        all_mirror_plane_to_wall = self.all_mirror_plane_to_wall()
        reference_point_reflect = self.reference_point_reflect(x1,y1)
        for wall in wall_axis:
            wall_position1 = self.room_wall_intersect(np.array([x1,y1]),np.array([reference_point_reflect[0],reference_point_reflect[1]]),np.array(wall[0]),np.array(wall[1]))
            #print("wall chk :",wall_position1)
            if wall_position1[0]>=0 and wall_position1[0]<=RoomWidth and wall_position1[1]>=0 and wall_position1[1]<=RoomHeight:
                mirror_intersect = intersect(all_mirror_plane_to_wall[0][0],all_mirror_plane_to_wall[0][1],all_mirror_plane_to_wall[1][0],all_mirror_plane_to_wall[1][1],TxLocation[0],TxLocation[1],wall_position1[0],wall_position1[1])
                refPoint_mir_wall_intersect = intersect(reference_point_reflect[0],reference_point_reflect[1],wall_position1[0],wall_position1[1],all_mirror_plane_to_wall[0][0],all_mirror_plane_to_wall[0][1],all_mirror_plane_to_wall[1][0],all_mirror_plane_to_wall[1][1])
                
                if mirror_intersect[0] == False:
                    if refPoint_mir_wall_intersect[0]==False:
                        break
                else:
                    pass
                    #print("but through mirror !! sorry let's seek another wall... ")
        #print(wall_position1)
        #plt.plot([wall_position1[0],x1],[wall_position1[1],y1],'r:',alpha=0.5)
        #print(wall)
        return wall_position1  #return to wall in all_ray_list
    
    def all_ray_list(self):
        ray_list =[]
        for p in self.mirror_reflect_points:      
            wall = self.wall_position(p[0],p[1])
            if wall[0]>=0 and wall[0]<=500 and wall[1]>=0 and wall[1]<=500 :
                ray_list.append([p[0],p[1],wall[0],wall[1]])
        return ray_list
     
    def light_on_mirror(self,tx_id): 
        light_on_mirror = []
        for px,py in self.mirror_pixel_list:
            cover = True
            for tx_plane in tx_id.all_tx_plane_id:
                if tx_plane.pixel_intensity[px][py]==0.0:
                    cover = False
                    break
            if cover==True:
               light_on_mirror.append([px,py])  
        return light_on_mirror

    def ray_to_object(self):
        ray_to_object_list = self.all_ray_list() #Initial rays from ray list
        i = 0
        for ray in ray_to_object_list:
            ray_distance = np.hypot(ray[2]-ray[0],ray[3]-ray[1])  
            for obs in all_obs_id:
                for k in range(len(obs.obsX)):
                    x3 = obs.obsX[k]
                    y3 = obs.obsY[k]
                    x4 = obs.obsX[(k+1)%len(obs.obsX)]
                    y4 = obs.obsY[(k+1)%len(obs.obsY)]
                    result = intersect(ray[0],ray[1],ray[2],ray[3],x3,y3,x4,y4)
                    if result[0] == True:
                        distance = np.hypot(ray[0]-result[1],ray[1]-result[2])  
                        if distance < ray_distance:
                            ray_distance = distance  
                            ray_to_object_list[i][2] = result[1]  #change wall[0] of ray_list[i] to be result[1]
                            ray_to_object_list[i][3] = result[2]
            i+=1
        
        #print("ray_to_object_list =",ray_to_object_list)                
        return ray_to_object_list  
    
    def all_ray_coverage(self):
        ray_coverage_list = []
        for ray in self.ray_to_object_list:              
            ray_coverage_list.append(self.mirror_pixel_plane(ray[0],ray[1],ray[2],ray[3]))  
        #print(ray_coverage_list)
        return ray_coverage_list
    
    def count_coverage(self):
        self.ray_pixel_coverage = PixelPlane(MAX_ROW,MAX_COLUMN)    
        for ray in self.all_ray_coverage_list:            #self.all_ray_coverage_list = self.all_ray_coverage()
            #print(ray)
            for gx,gy in ray:
                 if not math.isclose(c1[gx][gy],5.0) :   
                    if not math.isclose(c1[gx][gy],2.0) :
                        global c2 
                        c2 = c2 + 1
                        self.ray_pixel_coverage.pixel_intensity[gx][gy] = 2.0
                        c1[gx][gy] = 2.0
                ##print("set coverage at ",gx,gy)                                 
        
        #count = sum(self.ray_pixel_coverage.pixel_intensity)  # <-- some computer cannot run "sum"
        count = 0
        for gx in range(MAX_COLUMN):
            for gy in range(MAX_ROW):
                if self.ray_pixel_coverage.pixel_intensity[gx][gy] == 2.0:
                    count +=1
        
        
        count_percent = (count*100)/(MAX_ROW*MAX_COLUMN)
        
        # For compare with light coverage without mirror
        tx_cov_count = 0
        improve_coverage_from_light = count
        for px in range(MAX_COLUMN):
            for py in range(MAX_ROW):
                cover = True
                for tx_plane in tx_id.all_tx_plane_id:
                    if tx_plane.pixel_intensity[px][py]==0.0:
                        cover = False
                        break
                if cover==True:
                   T[px][py] = 2.0
                   tx_cov_count+=1
                   if self.ray_pixel_coverage.pixel_intensity[px][py] == 2.0:
                       improve_coverage_from_light -= 1
        global Txcoverage
        Txcoverage = tx_cov_count                    
        tx_cov_count_percent = (tx_cov_count*100)/(MAX_ROW*MAX_COLUMN)
        improve_percent = improve_coverage_from_light*100/(MAX_ROW*MAX_COLUMN)
        new_overall_covearge =  tx_cov_count +  improve_coverage_from_light
        new_overall_covearge_percent = (new_overall_covearge*100)/(MAX_ROW*MAX_COLUMN)
        cover2.append(new_overall_covearge)          
        #print("Mirror's Light Coverage Area (Pixels) =",count)
        #print("Mirror's Light Coverage Percent =",count_percent,"%")
        #print("Transmitter Coverage Area =",tx_cov_count)
        #print("Transmitter Coverage Percent =",tx_cov_count_percent,"%")
        #print("Mirror's Improved Area (reflect to shadow area) =", improve_coverage_from_light)
        #print("Mirror's Improved Percent =",improve_percent,"%")
        #print("New Coverage Area =",new_overall_covearge)
        #print("New Covearge Area Percent =",new_overall_covearge_percent,"%")
           
        return count,count_percent,tx_cov_count,improve_coverage_from_light,improve_percent,new_overall_covearge,new_overall_covearge       
     

    def perp_toMirror(self,tx,ty,rx1,ry1,rx2,ry2): #Find the point that perpendicular to mirror
        x1=rx1
        y1=ry1
        x2=rx2
        y2=ry2
        x3=tx
        y3=ty
        k = ((y2-y1) * (x3-x1) - (x2-x1) * (y3-y1)) / ((y2-y1)**2 + (x2-x1)**2)
        x4 = x3 - k * (y2-y1)
        y4 = y3 + k * (x2-x1)
        ##print(x4,y4) # point perpend to mirror
        return (x4,y4)
    
    def perp(self,a) : #vector
        b = np.empty_like(a)
        b[0] = -a[1]
        b[1] = a[0]
        return b

    def room_wall_intersect(self,p1,p2, q1,q2) :
        da = p2-p1
        db = q2-q1
        dp = p1-q1
        dap = self.perp(da)
        denom = np.dot( dap, db)
        num = np.dot( dap, dp )
        return (num / denom.astype(float))*db + q1
    
    def shotgun(self,multiple): # On experiment // build multiple rays between 2 rays
        new_mirror_reflect_point_list = []
        for k in range(len(self.project_pixel_list)-1):
            px = self.project_pixel_list[k][0]
            pxx = self.project_pixel_list[k+1][0]
            py = self.project_pixel_list[k][1]
            pyy = self.project_pixel_list[k+1][1]
            add_px=[]
            add_py=[]
            for m in range(multiple):
                x_step = (pxx-px)/multiple
                add_px.append(px+x_step)
                y_step = (pyy-py)/multiple
                add_py.append(py+y_step)

            #print("add_px =",add_px)
            #add_py = np.arange(py,pyy,(pyy-py)/multiple)
            for i in range(multiple):
                new_mirror_reflect_point_list.append([add_px[i],add_py[i]])

        new_mirror_reflect_point_list.append(self.project_pixel_list[k+1])
        #print("New multiple rays (shotgun) =",new_mirror_reflect_point_list )
        return new_mirror_reflect_point_list

class Obstacle():
    def __init__(self,obsX,obsY,grid_size):
        self.obsX = obsX    #obsX
        self.obsY = obsY    #obsY
        self.grid_size = grid_size
        self.obs_plane = PixelPlane(MAX_ROW,MAX_COLUMN)
        self.px_max = math.ceil(max(obsX)/grid_size)
        self.py_max = math.ceil(max(obsY)/grid_size)
        self.px_min = math.floor(min(obsX)/grid_size)
        self.py_min = math.floor(min(obsY)/grid_size)
        self.pixel_count = 0
        self.setPixel()
        self.countPixel()

    def setPixel(self):
        for py in range(self.py_min,self.py_max+1):
            for px in range(self.px_min,self.px_max+1):
                cnt_crosspoint_r = 0
                cnt_crosspoint_l = 0
                x = px*self.grid_size
                y = py*self.grid_size       
                is_edge = False
                for k in range(len(self.obsX)):  
                    x1 = x
                    y1 = y
                    x2 = self.px_max*self.grid_size
                    y2 = y
                    x3 = self.obsX[k]
                    y3 = self.obsY[k]
                    x4 = self.obsX[(k+1)%len(self.obsX)]
                    y4 = self.obsY[(k+1)%len(self.obsY)]
                    
                    result = intersect(x1,y1,x2,y2,x3,y3,x4,y4)
                    if result[0]==True:
                        Xx = result[1]   
                        Xy = result[2]   
                        if Xx == x3 and Xy == y3:   
                            if y4<y1:
                                cnt_crosspoint_r += 1
                        elif Xx == x4 and Xy == y4:
                            if y3<y1:
                                cnt_crosspoint_r += 1
                        else:
                            cnt_crosspoint_r += 1
                        
                    else:  # two lines are parallel
                        if y1==y3 and y1==y4:  # two lines are overlapped
                            if (x1 >= x3 and x1 <= x4) or (x1>=x4 and x1 <= x3):
                                is_edge = True
                                self.obs_plane.pixel_intensity[int(x1/self.grid_size)][int(y1/self.grid_size)] = 1.0
                                #print([int(x1/self.grid_size)][int(y1/self.grid_size)])
                                break

                    x2 = self.px_min*self.grid_size
                    result = intersect(x1,y1,x2,y2,x3,y3,x4,y4)
                    if result[0]==True:
                        Xx = result[1]
                        Xy = result[2]                               
                        if Xx == x3 and Xy == y3:
                            if y4<y1:
                                cnt_crosspoint_l += 1
                        elif Xx == x4 and Xy == y4:
                            if y3<y1:
                                cnt_crosspoint_l += 1
                        else:
                            cnt_crosspoint_l += 1

                    else:  # two lines are parallel
                        if y1==y3 and y1==y4:  # two lines are overlapped
                            if (x1 >= x3 and x1 <= x4) or (x1>=x4 and x1 <= x3):
                                is_edge = True
                                self.obs_plane.pixel_intensity[int(x1/self.grid_size)][int(y1/self.grid_size)] = 1.0
                                break
                    
                if is_edge==False:      
                    if cnt_crosspoint_r%2==1 or cnt_crosspoint_l%2==1:
                        self.obs_plane.pixel_intensity[int(x1/self.grid_size)][int(y1/self.grid_size)] = 1.0    

    def countPixel(self):
        px_cnt = 0
        for gx in range(MAX_COLUMN):
            for gy in range(MAX_ROW):
                if self.obs_plane.pixel_intensity[gx][gy] == 1.0:
                    px_cnt+=1
        #print("Pixel Count =",px_cnt)
                
        self.pixel_count = px_cnt
### create Pixel Class
class Display() :
    def __init__(self,MAX_ROW,MAX_COLUMN,grid_size) :
        self.MAX_ROW = MAX_ROW
        self.MAX_COLUMN = MAX_COLUMN
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111,aspect='equal')
        self.ax.set_xlim([-grid_size,MAX_COLUMN*grid_size])
        self.ax.set_ylim([-grid_size,MAX_ROW*grid_size])
        self.grid_size = grid_size
        self.pixel_id = []
        self.createpixel()
        plt.show(block=False)
        
    def createpixel(self):       
        for k in range (self.MAX_COLUMN):
            self.pixel_id.append([0]*self.MAX_ROW)
        for x in range (self.MAX_COLUMN) :
            for y in range (self.MAX_ROW) :
                self.pixel_id[x][y] = self.ax.add_patch(
                  patches.Rectangle(
                  ((x-0.5)*self.grid_size,(y-0.5)*self.grid_size),   # (x,y)
                    self.grid_size,          # width
                    self.grid_size,          # height
                    color = 'salmon',alpha=0.1,ec='lightgrey'))

    def draw_obstacle(self,obs_id,color):
        obstacle_plane = obs_id.obs_plane
        
        for x in range (self.MAX_COLUMN):
            for y in range (self.MAX_ROW):
                if obstacle_plane.pixel_intensity[x][y]>0.0:
                    c1[x][y] = 5.0
                    self.pixel_id[x][y].set_facecolor(color)
                    self.pixel_id[x][y].set_alpha(0.5)

        X = [k for k in obs_id.obsX]
        X.append(obs_id.obsX[0])
        Y = [k for k in obs_id.obsY]
        Y.append(obs_id.obsY[0])      
        self.ax.plot(X,Y,'-b')

      
        X =  [obs_id.px_min,obs_id.px_max,obs_id.px_max,obs_id.px_min,obs_id.px_min]
        Y =  [obs_id.py_min,obs_id.py_min,obs_id.py_max,obs_id.py_max,obs_id.py_min]
        X = [k*self.grid_size for k in X]
        Y = [k*self.grid_size for k in Y]
        self.ax.plot(X,Y,color='gray',alpha=0.3)
        self.fig.show()
        
    def draw_tx(self,tx_id):
        self.ax.plot(tx_id.xloc,tx_id.yloc,marker='o',color='red',alpha=1)
        for px in range (self.MAX_COLUMN):
            for py in range (self.MAX_ROW):
                cover = True
                for tx_plane in tx_id.all_tx_plane_id:
                    if tx_plane.pixel_intensity[px][py]==0.0:
                        cover = False
                        break
                if cover==True:
                    self.pixel_id[px][py].set_facecolor('g')
                    self.pixel_id[px][py].set_alpha(0.5)
        self.fig.show()
        
    def draw_mirror(self,mirror_id): 
        wallSide=mirror_id.wallSide
        position =mirror_id.position
        distance  =mirror_id.distance
        mirror_width =mirror_id.mirror_width
        angle  =mirror_id.angle
        
        mirro_loc = mirror_id.mirror_loc
        x1 = mirro_loc[0]
        y1 = mirro_loc[1]
        x2 = mirro_loc[2]
        y2 = mirro_loc[3]
        
        self.ax.plot([x1,x2],[y1,y2],marker='s',color='red') 
        if mirror_id.wallSide == 'R' :
            self.ax.plot([np.average([x1,x2]),RoomWidth],[position,position],marker='o',color='b')
        elif mirror_id.wallSide == 'L' :
            self.ax.plot([0,np.average([x1,x2])],[position,position],marker='o',color='b')
        elif mirror_id.wallSide == 'T' :
            self.ax.plot([position,position],[np.average([y1,y2]),RoomHeight],marker='o',color='b')
        elif mirror_id.wallSide == 'B' :
            self.ax.plot([position,position],[0,np.average([y1,y2])],marker='o',color='b')
        

        for gx,gy in mirror_id.mirror_pixel_list:
            self.pixel_id[gx][gy].set_facecolor('yellow')
            
        self.fig.show()   
        
    def draw_ray(self,mirror_id):
        ray_list = mirror_id.all_ray_list()  
        for ray in ray_list:
            x1=ray[0] 
            y1=ray[1]
            x2=ray[2]  
            y2=ray[3]
            self.ax.plot([x1,x2],[y1,y2],color='gold',linestyle=':',alpha=0.5)
            
    def draw_ray_to_object(self,mirror_id):
        ray_list = mirror_id.ray_to_object()  
        for ray in ray_list:
            x1=ray[0]
            y1=ray[1]
            x2=ray[2]
            y2=ray[3]
            self.ax.plot([x1,x2],[y1,y2],color='red',linestyle=':',alpha=0.2)
      
    def draw_ray_coverage(self,mirror_id):
        for ray_cov in mirror_id.all_ray_coverage_list:   #/grid
            for gx,gy in ray_cov:
                if not math.isclose(c1[gx][gy],5.0) :
                   self.pixel_id[gx][gy].set_facecolor('blue')
                   self.pixel_id[gx][gy].set_alpha(0.2)
        self.fig.show()
         
        
class Transmitter():
    def __init__(self,TxLocation,grid_size,all_obs_id):
        self.xloc = TxLocation[0] 
        self.yloc = TxLocation[1]
        self.grid_size = grid_size
        self.all_obs_id = all_obs_id
        
        self.all_tx_plane_id = []
        for obs in all_obs_id:
            self.set_pixel_coverage(obs)
        self.count_coverage()
    def set_pixel_coverage(self,obs):   
        tx_plane_id = PixelPlane(MAX_ROW,MAX_COLUMN)
        self.all_tx_plane_id.append(tx_plane_id) 
        #### Find coverage area
        x1 = self.xloc
        y1 = self.yloc   
        for px2 in range(MAX_COLUMN):
            for py2 in range(MAX_ROW):
                if obs.obs_plane.pixel_intensity[px2][py2]==0.0:
                    x2 = px2*grid_size
                    y2 = py2*grid_size
                    light_blocked = False
                    for k in range(len(obs.obsX)):
                        x3 = obs.obsX[k]
                        y3 = obs.obsY[k]
                        x4 = obs.obsX[(k+1)%len(obs.obsX)]
                        y4 = obs.obsY[(k+1)%len(obs.obsY)]
                        result = intersect(x1,y1,x2,y2,x3,y3,x4,y4)
                        if result[0]==True:
                            light_blocked = True       
                            break
                    if light_blocked == False:
                         tx_plane_id.pixel_intensity[int(x2/grid_size)][int(y2/grid_size)] = 1.0
    def count_coverage(self):
        pass


myDisplay = Display(MAX_ROW,MAX_COLUMN,grid_size) 


all_obs_id = []
for obs_vertex in allOBSvertex:
    obsX = [k[0] for k in obs_vertex]
    obsY = [k[1] for k in obs_vertex]
    all_obs_id.append(Obstacle(obsX,obsY,grid_size))
   
myDisplay.draw_obstacle(all_obs_id[0],'tan')
myDisplay.draw_obstacle(all_obs_id[1],'chocolate')
myDisplay.draw_obstacle(all_obs_id[2],'chocolate')

tx_id = Transmitter(TxLocation,grid_size,all_obs_id)
myDisplay.draw_tx(tx_id)
'''
num_of_mirror = 4
mirrorDisplays = []

for i in range(num_of_mirror):
    myDisplay2 = Display(MAX_ROW,MAX_COLUMN,grid_size)
    myDisplay2.draw_tx(tx_id)
    
    myDisplay2.draw_obstacle(all_obs_id[0],'tan')
    myDisplay2.draw_obstacle(all_obs_id[1],'chocolate')
    myDisplay2.draw_obstacle(all_obs_id[2],'chocolate')
    mirrorDisplays.append(myDisplay2)
'''   
obs_pixels =0
for obs in all_obs_id:
    #print("Pixel count =",obs.pixel_count)
    obs_pixels += obs.pixel_count
    
 
fig = plt.figure()

ax = fig.add_subplot(111)
#ax1 = fig.add_subplot(111)
colorlist = ['blue','green','magenta','orange','red','brown','pink','grey','black']

ax.legend()

for index,max_angle in enumerate(range(5,50,5)):
    size = []
    covlist = []
    all_myMirror = []
    for s in range(0,50,10):
        c1 = []
        c2 = 0
        cover2 = []
        cov = []
        T = []
        Txcoverage = 0
        
        for k in range(MAX_ROW):
            c1.append([0]*MAX_COLUMN)
            T.append([0]*MAX_COLUMN)
        for x in range(MAX_COLUMN) :
            for y in range(MAX_ROW) :
                c1[x][y] = 0.0
                T[x][y] = 0.0

        for angle in range(-max_angle,max_angle+1,1):
            mirror = (['R',125,30,s,angle],['L',125,30,s,angle],['L',375,30,s,angle],['L',375,30,s,angle])
            for i in range (len(mirror)):
                mir = mirror[i]
                myMirror = Mirror(mir[0],mir[1],mir[2],mir[3],mir[4]) 
                all_myMirror.append(myMirror)
                covv = myMirror.count_coverage()
            
        size.append(s)
        c22 = c2
        for x in range(MAX_COLUMN) :
            for y in range(MAX_ROW) :
                if math.isclose(T[x][y],2.0) :   
                    if math.isclose(c1[x][y],2.0) :
                        c22 -= 1
                        
        cov = Txcoverage+c22
        #print(obs_pixels)
        #print(Txcoverage)
        cov_percent = cov*100/((MAX_ROW*MAX_COLUMN)-obs_pixels)
        #print(cov_percent)
        covlist.append(cov_percent)
    #print(covlist)
    #print(covlist)
    ax.plot(size,covlist, color = colorlist[index], marker = 'D', linestyle = '-',label='R&L-max angle ={}'.format(max_angle))
    #line = mlines.Line2D(color=colorlist[index], marker = 'D', label='R&L-max angle ={}'.format(max_angle))


ax.set_xlabel('Reflector size (cm)')
ax.set_ylabel('Signal Coverage (%)')
plt.yticks(range(60,105,5))
plt.xticks(range(0,75,5))
            
plt.grid()
plt.legend()

plt.show()

#
