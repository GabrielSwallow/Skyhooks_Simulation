# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 16:19:01 2020

@author: AlexanderSinclairTeixeira
"""

import scipy as sp
#import scipy.integrate as itg
import matplotlib
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import matplotlib.patches as pat #shapes

#%%
TIME_STEP = 10 #time resolution in seconds
FRAMES=4000 #num frames in animation
INTERVAL=5 #millis between each frame in animation
G = 6.673e-11 #grav const, N m^2 kg^-2
zHat = sp.array([0,0,1]) #unit vector lol

def mag(r1, r2): #fast global computation of size of vector
    return sp.sqrt(sp.dot(r2-r1, r2-r1))

def unitV(r1, r2): #unit vector from r1 to r2
    sep = r2 - r1
    return sep / mag(r1, r2)


class _RigidBody: 
    #Parent class for all objects. Do not instantiate directly, rather use 
    #one of its subclasses (either Large or Small)
    def __init__(self, mass, pos, vel, rad, theta, angV, pointFmt): #on start:
        self.mass = mass
        self.pos = sp.array(pos) #posn vector
        self.vel = sp.array(vel) #vel vector
        self.rad = rad #radius
        self.theta = theta #angle from global +ve x-axis in radians
        self.angV = angV #angular velocity, +ve is anticlockwise
        self.I = 2 / 5 * self.mass * self.rad ** 2 #sphere moment of inertia
        #should this be of a thin disc or a circle instead?!?
        self.patch = pat.Circle(xy=self.pos, radius=self.rad, alpha=0.3) 
        #circles are nice, default, hitbox
        self.point, = ax1.plot(self.pos[0], self.pos[1], pointFmt)

    def pointUpdate(self):
        self.point.set_data(self.pos[0], self.pos[1])
        return self.point

    def move(self): #add steps number?
        #for i in range(steps):
        self.pos = self.pos + self.vel * TIME_STEP
        self.theta = self.theta + self.angV * TIME_STEP
        return self.patchUpdate(), self.pointUpdate()

    def contact(self, other): #is centre of other inside radius?
        return True if mag(self.pos, other.pos) <= self.rad else False


class Large(_RigidBody): #has gravity, fixed circular path

    def __init__(self, mass, pos, vel, rad, theta, angV, pointFmt):
        super().__init__(mass, pos, vel, rad, theta, angV, pointFmt) #as above

    def patchUpdate(self): #for animation
        self.patch.set_center(self.pos)
        return self.patch

    def g(self, r): #gravitational field strength at vector r
        g = - G * self.mass / mag(self.pos, r)**2 #-GM/R^2 
        return g * unitV(self.pos, r) #with direction

    def circular(self, r, antiCW): 
        #velocity for a circular orbit around body at vector r
        #Anticlockwise is default, enter -1 for clockwise
        v = sp.sqrt(G * self.mass / mag(self.pos, r)) #mag of vel v=(GM/R)^0.5
        rHat = sp.append(unitV(self.pos, r), 0) #make it 3D
        vHat = antiCW * sp.cross(zHat, rHat)[:2] #return to 2D
        return v * vHat


class Small(_RigidBody): #variable v, no gravity

    def __init__(self, mass, pos, vel, rad, theta, angV, pointFmt):
        super().__init__(mass, pos, vel, rad, theta, angV, pointFmt)
        self.patch = pat.Arrow(self.pos[0], self.pos[1], 
                               1e4*self.vel[0], 1e4*self.vel[1], 
                               width=lim2/50) #position at tail of arrow
        self.a = sp.array([0,0]) #initialise acceleration

    def coalesce(self, other): #a completely useless function
        mass = self.mass + other.mass
        pos = self.pos + other.rad * unitV(self.pos, other.pos) #posn avg
        vel = (self.mass * self.vel
               + other.mass * other.vel) / mass #conservation of momentum
        rad = (self.rad ** 2 + other.rad ** 2) ** 0.5 #assuming equal densities
        L = self.I * self.angV + other.I * other.angV #sum of amgular momenta
        angV = L / (2 / 5 * mass * rad ** 2) #find new ang vel
        return _RigidBody(mass, pos, vel, rad, self.theta, angV) #returns class

    def patchUpdate(self): #objects too small to see so vel shown
        ax2.patches.remove(self.patch)
        self.patch = pat.Arrow(self.pos[0], self.pos[1],
                               50*self.vel[0], 50*self.vel[1],
                               width=lim2/50) #redraw arrow
        ax2.add_patch(self.patch) #add to axes
        return self.patch #for animation func


#%%
fig = plt.figure(figsize=[16,8])
lim1 = 2.5e11
ax1 = fig.add_subplot(121, xlim=(-lim1, lim1), ylim=(-lim1, lim1))
lim2 = 12e6
ax2 = fig.add_subplot(122, xlim=(-lim2, lim2), ylim=(-lim2, lim2))
fig.tight_layout()
t = 0.0

#%%
Sun = Large(mass=1.989e30, pos=[0,0], vel=[0,0], rad=6.963e8, 
            theta=0, angV=0, pointFmt=".y")
Earth = Large(mass=5.972e24, pos=[1.496e11,0], vel=[0,0], rad=6.37e6, 
              theta=0, angV=2*sp.pi/60/60/24, pointFmt=".b") #24 hour day
Mars = Large(mass=6.39e23, pos=[-2.289e11, 0], vel=[0,0], rad=3.389e6, 
             theta=0, angV=2*sp.pi/60/60/(24+2/3), pointFmt=".r") #24h40min day
Demo2 = Small(mass=1.5e4, pos=[1.496086e11, 0], vel=[-2e3,3.7e4], rad=6, 
              theta=0, angV=0, pointFmt="xg") #mid earth orbit: 2230km
LargeList = [Sun, Earth, Mars]
SmallList = [Demo2]

#%%
def focus(body):
#    xmin, xmax = ax2.get_xlim()
#    ymin, ymax = ax2.get_ylim()
#    if (body.pos[0] >= xmax or body.pos[0] <= xmin
#    or body.pos[1] >= ymax or body.pos[1] <= ymin):
    ax2.set_xlim(body.pos[0]-lim2/2, body.pos[0]+lim2/2)
    ax2.set_ylim(body.pos[1]-lim2/2, body.pos[1]+lim2/2)

def init():
    for i in ax2.patches:
        ax2.patches.remove(i)
    ax2.add_patch(Sun.patch)
    ax2.add_patch(Earth.patch)
    ax2.add_patch(Mars.patch)
    ax2.add_patch(Demo2.patch)
    focus(Earth)
    return (Sun.point, Earth.point, Mars.point, Demo2.point, *ax2.patches)

def animate(frame):
    Earth.vel = Sun.circular(Earth.pos, 1)
    Mars.vel = Sun.circular(Mars.pos, 1)
    Demo2.a = sum(body.g(Demo2.pos) for body in LargeList)
    Demo2.vel = Demo2.vel + Demo2.a * TIME_STEP
    focus(Demo2)
    return (Sun.point, *Earth.move(), *Mars.move(), *Demo2.move())

#%%
anim = FuncAnimation(fig, animate, frames=FRAMES, 
                     init_func=init, interval=INTERVAL, 
                     repeat=False, blit=True)

#init()
fig.show()

