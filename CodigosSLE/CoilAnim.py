import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
axis = plt.axes(xlim=(-50,50), ylim=(-50,50))
line, = axis.plot([],[], lw = 2)

def init():
    line.set_data([],[])
    return line,

xdata, ydata = [], []

def animate(i):
    t = 0.1 * i
    x = t * np.sin(t)
    y = t * np.cos(t)

    xdata.append(x)
    ydata.append(y)
    line.set_data(xdata, ydata)
    return line,

anim = animation.FuncAnimation(fig, animate, init_func = init,
                               frames = 500, interval = 20, blit = True)
anim.save("/mnt/c/Users/Esteban/Desktop/LearningPython/growingCoil.mp4", writer = 'ffmpeg', fps = 30)
