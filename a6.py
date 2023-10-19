import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as l
import control
from scipy.integrate import odeint

c1 = 1/60
c2 = 0.2
c3 = 0.1
ybar = 110

def ex3():
    c1 = 1/60
    c2 = 0.2
    c3 = 0.1

    A = np.array([[0,1],[0,0]])
    B = np.array([[0],[1]])

    # print(Be)
    # Ee = l.eigvals(Ae-0.05*Be @ np.array([[0,0,0,0,1]]))
    # print(Ee)
    Q = np.diag([1,1])
    R = 1

    L,_,_ = control.lqr(A,B,Q,R)
    print(L)

def ex4():
    def dx(x, u):
        # u is a 1-vector
        # u = u[0]

        x1,x2,x3,x4 = x[0],x[1],x[2],x[3]
        
        dx1 = -u*x1*x3 + c1*x4
        dx2 = u*x1*x3 - c2*x2
        dx3 = -c3*x3 + c2*x2
        dx4 = c3*x3 - c1*x4
        return [dx1,dx2,dx3,dx4]

    def simulate(u_fun, x0, t):
        def odefun(x,_):
            return dx(x, u_fun(x))
        return odeint(odefun, x0, t)

    x0 = np.array([60,15,20,5]).T
    t = np.linspace(0,20,1000)
    def alpha(x):
        return -((-c2**2-c2*c3)*x[1]+c3**2*x[2])/(c2*x[0]*x[2])
    def beta(x):
        return 1/(c2*x[0]*x[2])
    def v_fun(x):
        dx3 = -c3*x[2]+c2*x[1]
        e = x[2]-ybar
        return -(e+np.sqrt(3)*dx3)
    def u_fun(x):
        return alpha(x) + beta(x)*v_fun(x)
    solution = simulate(u_fun, x0, t)
    
    fig, ax = plt.subplots()
    ax.plot(t,solution, label=["Susceptible (x1)","Exposed (x2)","Infected (x3)","Recovered (x4)"])
    ax.legend()
    ax.grid()
    ax.set_title(f"ref = {ybar}")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Percent (%)")
    plt.show()

ex4()
