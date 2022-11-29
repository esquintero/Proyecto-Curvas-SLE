import matplotlib.pyplot as plt
import numpy as np
import math
import cmath as cm

ReAxis = np.linspace(0, 10, 10)
ImAxis = 0 + ReAxis * 1j
dt = 1
tmax = 1
N = 1000
# time = np.arange(dt, tmax, dt)
time = np.arange(0, N, 1)  # ahora en funcion del numero de puntos


mallacompleja = np.zeros((np.size(ReAxis), np.size(ImAxis)), dtype=complex)
for ix in range(np.size(ReAxis)):
    for iy in range(np.size(ImAxis)):
        mallacompleja[ix, iy] = ReAxis[ix] + ImAxis[iy]

# ++++++++++++++ DE LA FUNCION DIRECTORA A LA TRAZA +++++++++++++++++++++++++
# Metodo numerico del paper: 'Exact Solutions for Loewner Evolutions'
def drivingfunction(t):
    return t


def backpropagation_method(t, dt):
    z = 0j
    Rez = 0
    Imz = 0

    A = 1
    B = drivingfunction(t) + drivingfunction(t - dt)
    C = 2 * dt + drivingfunction(t) * drivingfunction(t - dt)
    discriminant = -4 * A * C + B**2
    if discriminant < 0:
        Rez = (1 / 2 * A) * (-B)
        Imz = (1 / 2 * A) * math.sqrt(-discriminant)
        z = complex(Rez, Imz)
        return z
    else:
        Rez = (1 / 2 * A) * (-B)
        Imz = (1 / 2 * A) * math.sqrt(discriminant)
        z = complex(Rez + Imz, 0)
        return z

    # while t > 0:
    #     A = 1
    #     B = drivingfunction(t) + drivingfunction(t - dt)
    #     C = 2 * dt + drivingfunction(t) * drivingfunction(t - dt)
    #     discriminant = -4 * A * C + B**2
    #     if discriminant < 0:
    #         Rez = (1 / 2 * A) * (-B)
    #         Imz = (1 / 2 * A) * math.sqrt(-discriminant)
    #         z = complex(Rez, Imz)
    #     else:
    #         Rez = (1 / 2 * A) * (-B)
    #         Imz = (1 / 2 * A) * math.sqrt(discriminant)
    #         z = complex(Rez + Imz, 0)
    #     t -= dt
    # return z


# Metodo numerico vertical slits
def verticalslits(rez, imz, delta, t, dt):
    z = complex(rez, imz)
    # hz = 1j * cm.sqrt(-4 * dt - (z - delta) ** 2)
    hz = cm.sqrt(z * z - 4 * dt) + delta
    if t < dt:
        t = 0
        return hz
    else:
        return verticalslits(np.real(hz), np.imag(hz), delta, t - dt, dt)


# Metodo numerico tilted slits
def tiltedslits(rez, imz, delta, t, dt):
    z = complex(rez, imz)
    v = delta * delta / dt
    # v = 8 / 3
    alpha = 0
    if delta > 0:
        alpha = 0.5 - 0.5 * math.sqrt(v / (16 + v))
    else:
        alpha = 0.5 + 0.5 * math.sqrt(v / (16 + v))
    A = (z + 2 * cm.sqrt(dt * (1 - alpha) / alpha)) ** (1 - alpha)
    B = (z - 2 * cm.sqrt(dt * alpha / (1 - alpha))) ** alpha
    hz = A * B
    # if t == 0:
    #     # t = 0
    #     return hz
    # else:
    #     return tiltedslits(np.real(hz), np.imag(hz), delta, t - 1, dt)
    return hz


# +++++++++++++++++++++++ FUNCION DIRECTORA +++++++++++++++++++++++++++++++++
# Obtencion del delta_k por random walk
def dk_randomwalk(kappa, deltat):
    delta_k = 0
    sign = np.random.uniform(0, 1)
    # sign = np.random.rand(1)
    if sign > 0.5:
        delta_k = math.sqrt(kappa * deltat)
    else:
        delta_k = -math.sqrt(kappa * deltat)
    return delta_k


# Obtenciond el delta_k por distribucion normal
def dk_normal(kappa, deltat):
    delta_k = 0
    mu = 0  # media de la gaussiana
    sigma = math.sqrt(kappa * deltat)  # desviacion estandar de la gaussiana
    return np.random.normal(mu, sigma)


# +++++++++++++++++++++++ MAIN +++++++++++++++++++++++++++++++++

tracex = np.zeros(np.size(time))
tracey = np.zeros(np.size(time))
aux = 0j
compvector = np.zeros(
    np.size(time), dtype=complex
)  # Almacena la traza en complejos, para analizar resultados

# Intento con el metodo de backpropagation, status: failed
# Comentarios: La funcion directora termina siendo siempre una constantes, falta de entendimiento del metodo numerico empleado.

# for i in range(np.size(time)):
#     aux = backpropagation_method(time[i], dt)
#     tracex[i] = np.real(aux)
#     tracey[i] = np.imag(aux)
#     compvector[i] = aux

# Intento con el metodo de vertical slits, status: failed
# Comentarios: Se comprobo que el Delta_k (nuestro dt en la simulacion) es 1/N, con N el numero de particiones del tiempo. La curva generada no se queda en el semiplano superior y describe un triangulo. Al menos ya se tiene un metodo recursivo que funciona. Se agrega una condicion a la fuerza para que nunca pase al semiplano inferior. Con esto la traza nunca pasa mucho mas alla del eje real negativo. No se ha encontrado la explicacion.
# for j in range(np.size(time)):
#     # dk = dk_randomwalk(8 / 3, dt)
#     dk = dk_normal(8 / 3, dt)
#     # print(dk)
#     aux = verticalslits(0, 0, dk, time[j], dt)  # verticalslits(rez, imz, delta, t, dt)
#     # condicion de rebote en el eje real
#     if np.imag(aux) < 0:
#         aux = aux.conjugate()
#     tracex[j] = np.real(aux)
#     tracey[j] = np.imag(aux)
#     compvector[j] = aux

# plt.scatter(tracex, tracey)
# plt.plot(tracex, tracey)
# numero_iter = 0
# delta_k = dk_randomwalk(8 / 3, 5)
# verticalslits(0, 0, delta_k, 1, dt, numero_iter)

# Intento con el metodo tilted slits, status: Ongoing
# Comentarios: Las graficas no se ven igual que en matlab. Existe un problema en que la curva nunca se dobla sobre si misma.
for k in range(np.size(time)):
    if k > 0:
        dk = dk_randomwalk(8 / 3, dt)
        # dk = dk_normal(8 / 3, dt)
        # print(dk)
        aux = tiltedslits(
            tracex[k - 1], tracey[k - 1], dk, time[k], dt
        )  # verticalslits(rez, imz, delta, t, dt)
        # condicion de rebote en el eje real
        # if np.imag(aux) < 0:
        #     aux = aux.conjugate()
        tracex[k] = np.real(aux)
        tracey[k] = np.imag(aux)
        compvector[k] = aux

# plt.scatter(tracex, tracey)
# plt.plot(tracex, tracey, "-", linewidth=1.0)
# # plt.axis("equal")
# plt.set_aspect("equal", "box")

fig, ax = plt.subplots()
ax.plot(tracex, tracey, "-", linewidth=1.0)
ax.axis("equal")
fig.tight_layout()
plt.show()
