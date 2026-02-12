## simple cut and paste of functions from henry's code, to be used in the notebook.
import math
import numpy as np

#### Universal constants

Mp = 1.220890e19
G = 1 / Mp**2
M_sun = 2e30 * 5.61e26

M_ann = 3.1e-04  # in solar masses
mua_ann = 2e-16  # in GeV
a_ann = 0.999999
alpha_ann = G*M_ann*mua_ann

M_tr = 1e-06     # in solar masses
mua_tr = 1e-13   # in GeV
a_tr = 0.999999
alpha_tr = G*M_tr*mua_tr



### Universal functions



def Omega_plus(a, r):
    return a / (2 * r)

def r_plus(M, a):
    return M * (1 + np.sqrt(1 - a**2)) * G

def gamma_a(l, alpha, M):
    if l == 1:
        p = 17
    else:
        p = 4 * l + 1
    rg = G * M
    return G * 1e-10 / rg**3 * (((alpha / l) * 0.5)**p + ((alpha / l) * 0.5)**(p + 1))

def gamma_t(ng,ne, mu, alpha, M):
    rg = G*M
    rc = ne**2 / alpha**2 * rg
    omega = 0.5* mu * alpha**2*(1/ng**2-1/ne**2)
    return 2 * G * omega**5 / 5 * mu**2 * rc**4


def glm(a, m, l, r, omega):
    g = 1
    if l == 1:
        g *= l**2 * (1 - a**2) + (a * m - 2 * r * omega)**2
    if l != 1:
        for k in range(1, l):
            g *= k**2 * (1 - a**2) + (a * m - 2 * r * omega)**2
    return g

def super_gamma(n, l, m, omega, mu, M, r, a):
    C = (
        2**(4 * l - 1)
        * math.factorial(n + l)
        / (n**(2 * l + 4) * math.factorial(n - l - 1))
        * (math.factorial(l) / (math.factorial(2 * l) * math.factorial(2 * l + 1)))**2
    )
    g = glm(a, m, l, r, omega)
    omega_plus = Omega_plus(a, r)
    gamma_nlm = (
        2 * r / M * C * g * (m * omega_plus - omega) * (G * mu * M)**(4 * l + 5) / G 
    )
    return gamma_nlm

# Annihilation Functions

def gamma_a(l, alpha, M):
    if l == 1:
        p = 17
    else:
        p = 4 * l + 1
    rg = G * M
    return G * 1e-10 / rg**3 * (((alpha / l) * 0.5)**p + ((alpha / l) * 0.5)**(p + 1))

def omega_ann(mua, alpha, n):
    return 2 * mua * (1 - alpha**2 / (2 * n**2))


def N_max(M):
    # Arvanitaki Eq 8, approximation
    return 10**(76)*(M/10)**2





##### Level Transition Frequency Domain Strain


def N_g_omega (omega, gamma_g,gamma_t, N_e,T):
    return 1/(gamma_g+gamma_t*N_e- 1j*omega)/(2*np.pi)*(1-np.exp(gamma_g+gamma_t*N_e- 1j*omega/(2*np.pi))*T)

def N_e_omega (omega, gamma_e,gamma_t, N_g,T):
    return 1/(gamma_e+gamma_t*N_g- 1j*omega/(2*np.pi))*(1-np.exp(gamma_e+gamma_t*N_g- 1j*omega/(2*np.pi))*T)

def h_tr(r,omega_tr,gamma_e,gamma_g,gamma_tr, T, omega,N_e,N_g): 
    # may need to fix expression due to complex factor in exponential
   
    return 1/np.sqrt(2*np.pi)*np.sqrt(4*G/(r**2*omega_tr) *gamma_tr *Ng*Ne) / ((gamma_g+gamma_e+gamma_tr*N_e-gamma_tr*N_g)- 1j*omega/(2*np.pi))*(1-np.exp(  (gamma_g+gamma_e+gamma_tr*N_e-gamma_tr*N_g)*T - 1j* omega/(2*np.pi) ))


####