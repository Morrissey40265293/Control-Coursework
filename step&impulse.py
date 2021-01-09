import numpy as np
import sympy as sym
from scipy import integrate
import matplotlib.pyplot as plt

# Input constants
m =  0.425 # mass (kg)
g = 9.81 # acceleration (m/s^2)
d = 0.42  #
delta = 0.65 #m
r = 0.125 #m
R = 53 #Ω
L0 = 120 #mH
L1 = 25 #mH
alpha = 1.2 #m−1
c = 6815 #(g·m3)/(A^2 .s^2)
k = 1880 #N/m
b = 10.4 #Ns/m
phi = 42
V= (34.45-53*(47.86))*(1880*(47.86)-792.4)

x1, x2, x3, V = sym.symbols('x1, x2, x3,V')
phi = (7*m)/5 * (c * ((x3)^2))/(delta-x1) + (m*g*np.sin(phi)) - k*(x1-d) - b*(x2)
psi= (V - R*x3)
psi/= L0 +L1*2.718281828^(-alpha*(delta-x1))
# Determine the partial derivatives of φ, w.r.t to x1, x2, x3
phi_deriv_x1 = phi.diff(x1)
phi_deriv_x2 = phi.diff(x2)
phi_deriv_x3 = phi.diff(x3)


psi_deriv_x1 = psi.diff(x1)
psi_deriv_x3 = psi.diff(x3)
psi_deriv_V= psi.diff(V)

x1e = 0
x2e = 0
x3e = 0
xe = 0
Ve = V

# Substitute symbols with equilibrium points
phi_deriv_x1_at_equlibrium = phi_deriv_x1.subs([(x1, x1e), (x2, x2e), (x3, x3e)])
phi_deriv_x2_at_equlibrium = phi_deriv_x2.subs([(x1, x1e), (x2, x2e), (x3, x3e)])
phi_deriv_x3_at_equlibrium = phi_deriv_x3.subs([(x1, x1e), (x2, x2e), (x3, x3e)])

psi_deriv_x1_at_equlibrium = psi_deriv_x1.subs([(x1, x1e), (x3, x3e), (V, Ve)])
psi_deriv_x3_at_equlibrium = psi_deriv_x3.subs([(x1, x1e), (x3, x3e), (V, Ve)])
psi_deriv_V_at_equlibrium = psi_deriv_V.subs([(x1, x1e), (x3, x3e), (V, Ve)])

# x2' = a1x1 + a2x2 + a3x3
a1= phi_deriv_x1_at_equlibrium
a2= phi_deriv_x2_at_equlibrium
a3= phi_deriv_x3_at_equlibrium

#x3' = b1x1 + b2x3 + b3V
b1= psi_deriv_x1_at_equlibrium
b2= psi_deriv_x3_at_equlibrium
b3= psi_deriv_V_at_equlibrium

s, t = sym.symbols('s, t')
a1, a2, a3, b1, b2, b3 = sym.symbols('a1, a2, a3, b1, b2, b3', real=True, positive=True)

# Define G_theta
G_theta= 1/b1 * (((-b1*x1)-(b3*V))/V)

F_s_impulse = 1
F_s_step = 1 / s


X1_s_impulse_G_theta = G_theta * F_s_impulse
x1_t_impulse_G_theta = sym.inverse_laplace_transform(X1_s_impulse_G_theta, s, t)
X1_s_step_G_theta = G_theta * F_s_step
x1_t_step_G_theta = sym.inverse_laplace_transform(X1_s_step_G_theta, s, t)

n_points = 500
t_final = 50
t_span = np . linspace (0 , t_final , n_points )


# Print the symbolic expressions
sym.pprint(x1_t_impulse_G_theta.simplify())
sym.pprint(x1_t_step_G_theta.simplify())

plt.plot(x1_t_impulse_G_theta, t_span )
plt.plot(x1_t_step_G_theta, t_span )
plt.show()
