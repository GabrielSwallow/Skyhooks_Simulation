# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 15:10:26 2020

@author: AlexanderSinclairTeixeira
"""

import scipy as sp
import scipy.integrate as spi
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import matplotlib.patches as pat #shapes

#%%
TIME_STEP = 20 #time resolution in seconds
FRAMES=5000 #num frames in animation
INTERVAL=10 #millis between each frame in animation
t = sp.linspace(0.0, TIME_STEP * (FRAMES - 1), num = FRAMES)
G = 6.673e-11 #grav const, N m^2 kg^-2
zHat = sp.array([0,0,1]) #unit vector lol

def mag(r1, r2): #fast global computation of size of vector
    return sp.sqrt(sp.dot(r2-r1, r2-r1))

def unitV(r1, r2): #unit vector !!!FROM r1 TO r2!!!
    sep = r2 - r1
    return sep / mag(r1, r2)

def g(Mass, rMass, r): #test grav acceleration at r from mass at rMass
    g = - G * Mass / mag(rMass, r)**2
    return g * unitV(rMass, r)


class _RigidBody:
    #Parent class for all objects. Do not instantiate directly, rather use
    #one of its subclasses (either Large or Small)
    def __init__(self, mass, pos, rad, pointFmt): #on start:
        self.mass = mass
        self.pos = sp.array(pos) #posn vector
        self.rad = rad #radius
        self.I = 2 / 5 * self.mass * self.rad ** 2 #sphere moment of inertia
        #should this be of a thin disc or a circle instead?!?
        self.point, = ax1.plot(self.pos[0], self.pos[1], pointFmt)

    def pointUpdate(self):
        self.point.set_data(self.pos[0], self.pos[1])
        return self.point

    def contact(self, other): #is centre of other inside radius?
        return True if mag(self.pos, other.pos) <= self.rad else False

    def update(self): #package artists into list to extend updated list
        return [self.patchUpdate(), self.pointUpdate()]

    def insideEvent(self, t, y): #convenience func
        # return zero if point is inside the body for integrator to stop
        return 0.0 if mag(self.orbitPos(t), y[:2]) < self.rad else 1.0
    insideEvent.terminal = True



class Large(_RigidBody): #has gravity, fixed circular/static path

    def __init__(self, mass, rad, orbitR, theta, pointFmt):
        self.orbitR = orbitR #orbital radius
        if orbitR != 0.0: #from Kep's 3rd law
            self.omega = sp.sqrt( G * Sun.mass / self.orbitR**3 )
        else: #if its the sun, or at origin
            self.omega = 0.0
        self.theta = theta #angle from global +ve x-axis, rads
        pos = self.orbitPos(0.0) #now we can define init pos
        super().__init__(mass, pos, rad, pointFmt) #continue as above
        self.patch = pat.Circle(xy=self.pos, radius=self.rad, alpha=0.4)
        #circles are nice, default, hitbox

    def move(self):
        pass #why did i do that

    def patchUpdate(self): #for animation
        self.patch.set_center(self.pos)
        return self.patch

    def gravFS(self, r): #grav field strength at r from THIS BODY
        g = - G * self.mass / mag(self.pos, r)**2 #-GM/R^2
        return g * unitV(self.pos, r) #with direction
        #repeated function.... take this one out?

    def orbitVel(self, r, antiCW):  #legacy
        #velocity for a circular orbit around body at vector r
        #Anticlockwise is default, enter -1 for clockwise
        v = sp.sqrt(G * self.mass / mag(self.pos, r)) #mag of vel v=(GM/R)^0.5
        rHat = sp.append(unitV(self.pos, r), 0) #make it 3D
        vHat = antiCW * sp.cross(zHat, rHat)[:2] #return to 2D
        return v * vHat

    def orbitPos(self, t): #exact orbit position as f(t)
        x = self.orbitR * sp.cos(self.omega * t + self.theta)
        y = self.orbitR * sp.sin(self.omega * t + self.theta)
        return sp.array([x,y])


class Small(_RigidBody): #variable v, no gravity

    def __init__(self, mass, pos, vel, rad, angle, angV, pointFmt):
        super().__init__(mass, pos, rad, pointFmt)
        self.vel = sp.array(vel) #vel vector
        self.angle = angle #angle from local +ve x-axis, radians
        self.angV = angV #angular velocity, radians per second
        self.patch = pat.Arrow(self.pos[0], self.pos[1],
                               1e4*self.vel[0], 1e4*self.vel[1],
                               width=lim2/50) #position at tail of arrow
        self.a = sp.array([0,0]) #initialise acceleration

    def move(self): #add steps number?
        #for i in range(steps):
        self.pos = self.pos + self.vel * TIME_STEP
        self.angle = self.angle + self.angV * TIME_STEP
        return self.pos, self.angle

    def coalesce(self, other): #useless function to define new arguments
        mass = self.mass + other.mass
        pos = self.pos + other.rad * unitV(self.pos, other.pos) #posn avg
        vel = (self.mass * self.vel
               + other.mass * other.vel) / mass #conservation of momentum
        rad = (self.rad ** 2 + other.rad ** 2) ** 0.5 #assuming equal densities
        L = self.I * self.angV + other.I * other.angV #sum of amgular momenta
        angV = L / (2 / 5 * mass * rad ** 2) #find new ang vel
        return Small(mass, pos, vel, rad, self.angle, angV, pointFmt="ok")

    def patchUpdate(self): #objects too small to see so vel shown
        ax2.patches.remove(self.patch)
        self.patch = pat.Arrow(self.pos[0], self.pos[1],
                               50*self.vel[0], 50*self.vel[1],
                               width=lim2/50) #redraw arrow
        ax2.add_patch(self.patch) #add to axes
        return self.patch #for animation func


#%% put into function?
fig = plt.figure(figsize=[18,8])
lim1 = 2.5e11
ax1 = fig.add_subplot(121, xlim=(-lim1, lim1), ylim=(-lim1, lim1))
lim2 = 4e6
ax2 = fig.add_subplot(122, xlim=(-lim2, lim2), ylim=(-lim2, lim2))
fig.tight_layout()

#%% put into function?
Sun = Large(mass=1.989e30, rad=6.963e8, orbitR=0.0,
            theta=0, pointFmt=".y")
Earth = Large(mass=5.972e24, rad=6.37e6, orbitR=1.496e11,
              theta=0, pointFmt=".b") #24 hour day
Mars = Large(mass=6.39e23, rad=3.389e6, orbitR=2.289e11,
             theta=sp.pi/4, pointFmt=".r") #24h40min day
Demo2 = Small(mass=1.5e4, pos=[1.496086e11, 10], vel=[-2e3,3.6431e4], rad=20,
              angle=0, angV=0, pointFmt="xg") #mid earth orbit: 2230km
LargeList = [Sun, Earth, Mars] #sun MUST be first
SmallList = [Demo2]
AllList = LargeList + SmallList

#%% solving
def gSum(t, r):
    allForces = sum( g(body.mass, body.orbitPos(t), r) for body in LargeList )
    #print(t, allForces)
    return allForces

def dU_dt(t, U): #vector U = (x, y, xdot, ydot)
    # returns Udot = (xdot, ydot, xdotdot, ydotdot)
    return (*U[2:], *gSum(t, U[:2]))

path = spi.solve_ivp(dU_dt, t_span=(t[0], t[-1]), y0=(*Demo2.pos, *Demo2.vel),
                     atol=(1e-1,1e-1, 1e-3, 1e-3), dense_output=True,
                     events=list(map(lambda x: x.insideEvent, LargeList)),
                     max_step=3600 ) #, t_eval=t)

#%%
def focus(body): #centre axis on body (or reset box)
#    xmin, xmax = ax2.get_xlim()
#    ymin, ymax = ax2.get_ylim()
#    if (body.pos[0] >= xmax or body.pos[0] <= xmin
#    or body.pos[1] >= ymax or body.pos[1] <= ymin):
    ax2.set_xlim(body.pos[0]-lim2/2, body.pos[0]+lim2/2)
    ax2.set_ylim(body.pos[1]-lim2/2, body.pos[1]+lim2/2)
    return ax2.get_xlim(), ax2.get_ylim()

def init():
    for i in ax2.patches:
        ax2.patches.remove(i)
    for body in AllList:
        ax2.add_patch(body.patch)
    focus(Demo2)
    return (Sun.point, Earth.point, Mars.point, Demo2.point, *ax2.patches)

def animate(frame):
    updated = [Sun.point, Sun.patch] #sun doesnt move
    for planet in LargeList[1:]: #large bodies after sun
        planet.pos = planet.orbitPos(t[frame])
        updated.extend(planet.update())
    Demo2.pos, Demo2.vel = sp.split(path.sol(t[frame]), 2)#x, y, v_x, v_y
    updated.extend(Demo2.update())
    focus(Demo2)
    return tuple(updated)

#%%
if path.status == 1: #temination events
    tmax = min(sp.concatenate(path.t_events)) #stop at first event
    tStopIndex = sp.searchsorted(t, tmax)
else:
    tStopIndex = len(t)
if path.success:
    anim = FuncAnimation(fig, animate, frames=len(t[:tStopIndex]),
                         init_func=init, interval=INTERVAL,
                         repeat=False, blit=True)
    fig.show()
