# -*- coding: utf-8 -*-
"""
Created on Wed May 13 16:28:18 2020

@author: AlexanderSinclairTeixeira
"""

import scipy as sp
#import scipy.integrate as itg
import matplotlib
import matplotlib.animation
import matplotlib.pyplot as plt
import matplotlib.patches as pat #shapes

#%%
TIME_STEP = 0.1 #time resolution in seconds
G = 6.673e-11 #grav const, N m^2 kg^-2
class PointMass(object): #Parent class for all objects. Do not instantiate directly
                         #rather, use one of its subclasses (either Large or Small)
    def __init__(self, mass, r, v): #on startup:
        self.mass = mass
        self.r = sp.array(r) #posn vector
        self.v = sp.array(v) #vel vector
    def move(self): #all objects move
        self.v = self.v + self.a * TIME_STEP #incrementing v
        self.r = self.r + self.v * TIME_STEP #incrementing r
        return self.patchUpdate() #redraw patch for use in animation function

class Large(PointMass): #fixed v, has gravity
    def __init__(self, mass, r, v, R):
        super().__init__(mass, r, v) #as above
        self.R = R #radius
        self.patch = pat.Circle(xy=self.r, radius=self.R, alpha=0.3) #planets are circles
        self.a = sp.array([0,0]) #no accel for now, change to centripetal off orbit radius param
    def patchUpdate(self): #for animation
        self.patch.set_center(self.r)
        return self.patch
    def g(self, s): #gravitational field strength at vector s
        sep = s - self.r #sep vector, r
        mag = sp.sqrt(sep.dot(sep)) #magnitude, ¦sep¦
        return - sep * G * self.mass / (mag ** 3) # g=-GM rHat/¦r¦^2, with rHat = r/¦r¦

class Small(PointMass): #variable v, no gravity
    def __init__(self, mass, r, v):
        super().__init__(mass, r, v)
        self.patch = pat.Arrow(self.r[0], self.r[1], self.v[0], self.v[1], width=3e5) #position at tail of arrow
        self.a = sp.array([0,0]) #initialise
    def patchUpdate(self): #ditto, small objects are too small to see so v vector shown
        ax1.patches.remove(self.patch)
        self.patch = pat.Arrow(self.r[0], self.r[1], 100*self.v[0], 100*self.v[1], width=3e5) #redraw arrow
        ax1.add_patch(self.patch) #add to axes
        return self.patch #for animation func

#%%
fig = plt.figure(figsize=[8,8]) #set up figure
ax1 = plt.axes(xlim=(-8e6, 8e6), ylim=(-8e6, 8e6) ) #set up first set of axes
A = Large(mass=5e24, r=[0,0], v=[0, 0], R=6.38e6) #massive bod
B = Small(mass=250, r=[0, 6.5e6], v=[7164, 0]) #petite lil satellite, moving at orbital v
#plt.grid()

def init():
    for i in ax1.patches: #using cla() or clear() doesnt work so -_-
        ax1.patches.remove(i)
    A.r = sp.array([0,0])#reinitilise positions
    B.r = sp.array([0, 6.8e6])
    ax1.add_patch(A.patch) #draw initial states
    ax1.add_patch(B.patch)
    return A.patch, B.patch #patches returned as iterable for blitting alg

def animate(frame):
    B.a = A.g(B.r)
    return A.move(), B.move() #for blitting alg

anim = matplotlib.animation.FuncAnimation(fig, animate, frames=7000,
                                          init_func=init, interval=1*TIME_STEP,
                                          blit=True, repeat = False) #interval in millis
#INCLUDE LINE BELOW TO SAVE
#anim.save("orbit.mp4", fps=100)
plt.show()
