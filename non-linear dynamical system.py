import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

class NL_DS:

    def __init__(self, m=0.425 , g = 9.81, d = 0.42 , delta = 0.65, r = 0.125, R = 53, L0 = 120, L1 = 25,t_final=2, alpha = 1.2, c = 6815, k = 1880, b = 10.4, phi = 42, V=8, x1=0., x2=0., x3=0.):
        self.__m = m
        self.__g = g
        self.__d = d
        self.__r = r
        self.__alpha = alpha
        self.__R = R
        self.__L0 = L0
        self.__L1 = L1
        self.__c = c
        self.__k = k
        self.__delta = delta
        self.__b= b
        self.__phi = phi
        self.__V= V
        self.__x1= x1
        self.__x2= x2
        self.__x3= x3
        self.__t_final = t_final

        def NL_Dynamical_System(t, z):
            theta = z[2]
            t_final= 2
            return [self.__x2,
                   (7*m)/5 * (c * ((x3)^2))/(delta-x1) + (m*g*np.sin(phi)) - k*(x1-d) - b*(x2),
                   (V - R*x3)/(L0 +L1*exp^(-alpha*(delta-x1)))]

        z_initial = [self.__x1, self.__x2, self.__x3]
        solution = solve_ivp(NL_Dynamical_System,
                             [0, t_final],
                             z_initial,
                             t_eval=np.linspace(0, t_final, 1000))
        self.__x1final = solution.y[0][-1]
        self.__x2final = solution.y[1][-1]
        self.__x3final = solution.y[2][-1]
        times= solution.t

        return solution

    def x1(self):
        return self.__x1final

    def x2(self):
        return self.__x2final

    def x3(self):
        return self.__x3final

plt.plot(x1final, times,)
plt.ylabel("x1 trajectory (m)")
plt.xlabel("time (s)")
plt.title("The Trajectory of the x1 over a Period of Time")
plt.show()

plt.plot(times, x2final)
plt.xlabel('Time (s)')
plt.ylabel('x2 (m) ')
plt.title("The Trajectory of the x2 over a Period of Time")
plt.show()

plt.plot(times, x3final)
plt.xlabel('Time (s)')
plt.ylabel('x3 (m) ')
plt.title("The Trajectory of the xx3 over a Period of Time")
plt.show()



