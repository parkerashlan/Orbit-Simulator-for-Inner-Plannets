"""
Created on Thu May 23 06:49:26 2019
    
@author: Ashlan Parker
"""
from numba import jit
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
import time as time

class Body():
    def __init__(self,distance,velocity,period):
        self.distance = distance
        self.velocity = np.sqrt((6.67e-11*1.989e30)/self.distance)
        self.period = np.sqrt(((4*np.pi**2)/(6.67e-11*1.989e30))*distance**3)

Mercury = Body(57.9e9,0,0)
Venus = Body(108.2e9,0,0)
Earth = Body(1.5e11,0,0)
Mars = Body(2.279e11,0,0)
Jupiter = Body(778.3e10,0,0)
Saturn = Body(1472.0e10,0,0)
Uranus = Body(2871.0e10,0,0)
Neptune = Body(4497e10,0,0)

planetdict ={'Mercury':Mercury,'Venus':Venus,'Earth':Earth,'Mars':Mars,'Jupiter':Jupiter,
             'Saturn':Saturn,'Uranus':Uranus,'Neptune':Neptune}

@jit(nopython=True,parallel=True)
def SunGravity(x0,y0,vx0,vy0,t,dt):
    
    G = 6.67e-11
    m = 1.989e30
    n = int(t/dt)
    
    xvals = np.zeros(n+1)
    yvals = np.zeros(n+1)
    
    ax = lambda x,y: -((G*m*x)/((x**2+y**2)**(3/2)))
    ay = lambda x,y: -((G*m*y)/((x**2+y**2)**(3/2)))
    
    x=x0
    y=y0
    
    vx=vx0
    vy=vy0
    
    xvals[0] = x0
    yvals[0] = y0
    
    for i in range(1,n+1):
        k1x = vx
        k1y = vy

        k1vx = ax(x,y)
        k1vy = ay(x,y)
        
        k2x = (vx+0.5*dt*k1vx)
        k2y = (vy+0.5*dt*k1vy)
        
        k2vx = ax(x+dt*0.5*k1x, y+dt*0.5*k1y)
        k2vy = ay(x+dt*0.5*k1x, y+dt*0.5*k1y)
        
        k3x = (vx + 0.5*k2vx)
        k3y = (vy +0.5*k2vy)
        
        k3vx = ax(x+dt*0.5*k2x, y+dt*0.5*k2y)
        k3vy = ay(x+dt*0.5*k2x, y+dt*0.5*k2y)
        
        k4x = vx+(dt*k3vx)
        k4y = vy+(dt*k3vy)
        
        k4vx = ax(x+dt*k3x,y+dt*k3y)
        k4vy = ay(x+dt*k3x,y+dt*k3y)
        
        
        x = x + (dt/6)*(k1x + 2*k2x + 2*k3x +k4x)
        y = y + (dt/6)*(k1y + 2*k2y + 2*k3y +k4y)
        vx = vx + (dt/6)*(k1vx + 2*k2vx +2*k3vx + k4vx)
        vy = vy + (dt/6)*(k1vy + 2*k2vy +2*k3vy + k4vy)
       
        xvals[i] = x
        yvals[i] = y
 
    return xvals,yvals

def ChooseBody(planet):
    
    if planetdict[planet].period >= 1e8:
        planetmotion = SunGravity(0,planetdict[planet].distance,planetdict[planet].velocity,0,planetdict[planet].period,10000)
    else:
        planetmotion = SunGravity(0,planetdict[planet].distance,planetdict[planet].velocity,0,planetdict[planet].period,100)

    return planetmotion

numberofplanets = input('Please input the number of planets:')
    
planets = []

for i in range(1,int(numberofplanets)+1):
    planetinput = input('Please input planet %s:' %i)
    planets.append(planetinput.capitalize())

planetlist = [ChooseBody(i) for i in planets] 
planetlistX,planetlistY = zip(*planetlist)
planetlistX = [i/1.496e11 for i in planetlistX]
planetlistY = [j/1.496e11 for j in planetlistY]

    
f1 = plt.figure()

if 'Jupiter' in planets or 'Saturn' in planets:
    ax1 = plt.axes(xlim = (-130,130),ylim = (-130,130))
elif 'Uranus' in planets or 'Neptune' in planets:
    ax1 = plt.axes(xlim = (-310,310),ylim = (-310,310))
else:
    ax1 = plt.axes(xlim = (-3,3),ylim = (-3,3))

line1, = ax1.plot([],[],lw=2,color='blue')
line2, = ax1.plot([],[],lw=2,color='green')
line3, = ax1.plot([],[],lw=2,color='orange')
line4, = ax1.plot([],[],lw=2,color='red')
line5, = ax1.plot([],[],lw=2,color='yellow')
line6, = ax1.plot([],[],lw=2,color='violet')
line7, = ax1.plot([],[],lw=2,color='cyan')
line8, = ax1.plot([],[],lw=2,color='black')

sun = patches.Circle((0,0),0.1,color='Yellow',fill=True)
ax1.add_patch(sun)

p1text = ax1.text(0.02,0.95,'', transform = ax1.transAxes)
p2text = ax1.text(0.02,0.90,'', transform = ax1.transAxes)
p3text = ax1.text(0.02,0.85,'', transform = ax1.transAxes)
p4text = ax1.text(0.02,0.80,'', transform = ax1.transAxes)
p5text = ax1.text(0.02,0.75,'', transform = ax1.transAxes)
p6text = ax1.text(0.02,0.70,'', transform = ax1.transAxes)
p7text = ax1.text(0.02,0.65,'', transform = ax1.transAxes)
p8text = ax1.text(0.02,0.60,'', transform = ax1.transAxes)

def init():
    line1.set_data([],[])
    line2.set_data([],[])
    line3.set_data([],[])
    line4.set_data([],[])
    line5.set_data([],[])
    line6.set_data([],[])
    line7.set_data([],[])
    line8.set_data([],[])
    p1text.set_text('')
    p2text.set_text('')
    p3text.set_text('')
    p4text.set_text('')
    p5text.set_text('')
    p6text.set_text('')
    p7text.set_text('')
    p8text.set_text('')
    return line1, line2, line3, line4, line5, line6, line7, line8, p1text, p2text, p3text, p4text, p5text, p6text, p7text, p8text

def updatePlanet1(i):
    orbitpx = planetlistX[0]
    orbitpy = planetlistY[0]
    line1.set_data(orbitpx[:1000*i],orbitpy[:1000*i])
    p1text.set_text('Blue Orbit: %s'.format(line1)%planets[0])
    return line1,p1text

aniplanet1 = animation.FuncAnimation(f1,updatePlanet1,init_func = init, 
                                     frames = int(np.size(planetlistX[0])/1000),interval = 1)

if int(numberofplanets) >= 2:
    def updatePlanet2(i):
        orbitp2x = planetlistX[1]
        orbitp2y = planetlistY[1]
        line2.set_data(orbitp2x[:1000*i],orbitp2y[:1000*i])
        p2text.set_text('Green Orbit: %s'.format(line2)%planets[1])
        return line2,p2text
    aniplanet2 = animation.FuncAnimation(f1,updatePlanet2,init_func = init, 
                                     frames = int(np.size(planetlistX[1])/1000),interval = 1)
if int(numberofplanets) >= 3:
    def updatePlanet3(i):
        orbitpx = planetlistX[2]
        orbitpy = planetlistY[2]
        line3.set_data(orbitpx[:1000*i],orbitpy[:1000*i])
        p3text.set_text('Orange Orbit: %s'.format(line3)%planets[2])
        return line3,p3text
    aniplanet3 = animation.FuncAnimation(f1,updatePlanet3,init_func = init, 
                                     frames = int(np.size(planetlistX[2])/1000),interval = 1)

if int(numberofplanets) >= 4:
    def updatePlanet4(i):
        orbitp2x = planetlistX[3]
        orbitp2y = planetlistY[3]
        line4.set_data(orbitp2x[:1000*i],orbitp2y[:1000*i])
        p4text.set_text('Red Orbit: %s'.format(line4)%planets[3])
        return line4,p4text
    aniplanet4 = animation.FuncAnimation(f1,updatePlanet4,init_func = init, 
                                     frames = int(np.size(planetlistX[3])/1000),interval = 1)
if int(numberofplanets) >= 5:
    def updatePlanet5(i):
        orbitp2x = planetlistX[4]
        orbitp2y = planetlistY[4]
        line5.set_data(orbitp2x[:1000*i],orbitp2y[:1000*i])
        p5text.set_text('Yellow Orbit: %s'.format(line5)%planets[4])
        return line5,p5text
    aniplanet5 = animation.FuncAnimation(f1,updatePlanet5,init_func = init, 
                                     frames = int(np.size(planetlistX[4])/1000),interval = 1)
if int(numberofplanets) >= 6:
    def updatePlanet6(i):
        orbitp2x = planetlistX[5]
        orbitp2y = planetlistY[5]
        line6.set_data(orbitp2x[:1000*i],orbitp2y[:1000*i])
        p6text.set_text('Violet Orbit: %s'.format(line6)%planets[5])
        return line6,p6text
    aniplanet6 = animation.FuncAnimation(f1,updatePlanet2,init_func = init, 
                                     frames = int(np.size(planetlistX[5])/1000),interval = 1)
if int(numberofplanets) >= 7:
    def updatePlanet7(i):
        orbitp2x = planetlistX[6]
        orbitp2y = planetlistY[6]
        line7.set_data(orbitp2x[:1000*i],orbitp2y[:1000*i])
        p7text.set_text('Cyan Orbit: %s'.format(line7)%planets[6])
        return line7,p7text
    aniplanet7 = animation.FuncAnimation(f1,updatePlanet7,init_func = init, 
                                     frames = int(np.size(planetlist[6])/1000),interval = 1)
if int(numberofplanets) == 8:
    def updatePlanet8(i):
        orbitp2x = planetlist[7]
        orbitp2y = planetlist[7]
        line8.set_data(orbitp2x[:1000*i],orbitp2y[:1000*i])
        p8text.set_text('Black Orbit: %s'.format(line8)%planets[7])
        return line8,p8text
    aniplanet8 = animation.FuncAnimation(f1,updatePlanet2,init_func = init, 
                                     frames = int(np.size(planetlist[7])/1000),interval = 1)

ax1.set_xlabel('Distance in X direction (AU)')
ax1.set_ylabel('Distance in Y direction (AU)')
ax1.set_title('Planetary Orbits')
plt.tight_layout
plt.show()
