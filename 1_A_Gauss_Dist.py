import numpy as np
import matplotlib.pyplot as plt

class Gaussian_Distribution:
    def __init__(self, mu=0, k=1.3806505e-23, T=3000, N=1e5, m=1.67355767e-27):
        self.mu = mu                 # mu is the middle value
        self.k = k                   # Boltzmanns constant
        self.T = T                   # temperature in Kelvin
        self.N = N                   # Number of molecules
        self.m = 2*m                   # mass
        self.sigma = np.sqrt(k*T/(2*m))  # sigma is the standard deviation

    def f(self, x):
        '''Normal probability distribution function'''
        return 1 / (np.sqrt(2*np.pi) * self.sigma) * np.exp(-0.5*((x - self.mu) / self.sigma)**2)

    def set_vx(self, vx):
        '''Gives instance vx list'''
        self.vx = vx
        return self.vx

    def Pvx(self, vx):
        k, T, N, m = self.k, self.T, self.N, self.m
        return np.sqrt(m/(2*np.pi*k*T)) * np.exp(-0.5 * ((m*vx**2)/(k*T)))


    def P(self, no, N=10001):
        '''Evaluates integral of f(x)'''
        self.a = self.mu - no*self.sigma    # lower limit
        self.b = self.mu + no*self.sigma    # upper limit
        ans = 0
        dx = (self.b - self.a)/N
        for i in range(N):
            ans += self.f(self.a + i*dx) * dx
        return ans

    def integral(self, a, b, N=10001):
        ans = 0
        dvx = (b - a)/N
        for i in range(N):
            ans += self.Pvx(a + i*dvx) * dvx
        return ans

    def abs_vel(self, v):
        m, k, T = self.m, self.k, self.T
        return (m/(2*np.pi*k*T))**(3/2) * np.exp(-0.5 * ((m*v**2)/(k*T))) * 4*np.pi*v**2

    def plot_vx(self):
        '''Method plots Pvx'''
        Pvx = self.Pvx(self.vx)
        plt.plot(vx, Pvx)
        plt.title('$v_x$ velocity distribution')
        plt.xlabel('$v_x$ m/s')
        plt.ylabel('P($v_x$)')
        plt.savefig('vx_velocity_distribution.png')

    def plot_abs_vel(self, v):
        abs = self.abs_vel(v)
        plt.plot(v, abs)
        plt.ticklabel_format(style='plain', useOffset=None)


ex_input = Gaussian_Distribution()

for i in range(1,4):
    P = ex_input.P(i)
    print(f'P(x in [-{i}sigma, {i}sigma]) = {P:.4f} = {P*100:.2f}%')

vx = np.linspace(-2.5e4, 2.5e4, 10001)    # 10^4 m/s

ex_input.set_vx(vx)
ex_input.plot_vx()
probability = ex_input.integral(5e3, 30e3)
# plt.show()
print(probability)
print(probability*ex_input.N)
v = np.linspace(0, 3e4, 10001)
abs_vel = ex_input.plot_abs_vel(v)
plt.show()
