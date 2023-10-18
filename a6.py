import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as l
import control
from scipy.integrate import odeint

def ex3():
    c1 = 1/60
    c2 = 0.2
    c3 = 0.1
    A = np.array([[0,0,0,c1],[0,-c2,0,0],[0,c2,-c3,0],[0,0,c3,-c1]])
    E = l.eigvals(A)
    # print(E)

    # Ae = np.array([[0,0,0,c1,0],[0,-c2,0,0,0],[0,c2,-c3,0,0],[0,0,c3,-c1,0],[0,0,-1,0,0]])
    al = 1.1
    Ae = np.array([[-al,0,0,c1],[al,-c2,0,0],[0,c2,-c3,0],[0,0,c3,-c1]])
    # print(Ae)
    # Be = np.array([[-1],[1],[0],[0],[0]])
    Be = np.array([[-1],[1],[0],[0]])
    # print(Be)
    # Ee = l.eigvals(Ae-0.05*Be @ np.array([[0,0,0,0,1]]))
    # print(Ee)
    Q = np.diag([0,0,1,0])
    R = 1

    L = control.lqr(Ae,Be,Q,R)
    print(L)

def ex4():
    def dx(x, u):
        # u is a 1-vector
        u = u[0]

        c1 = 1/60
        c2 = 0.2
        c3 = 0.1
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

    x0 = np.array([0.1,0.2,0.3,0.4]).T
    t = np.linspace(0,10,100)
    K = np.array([[1,2,3,4]])
    u_fun = lambda x : -K @ x
    solution = simulate(u_fun, x0, t)
    
    fig, ax = plt.subplots()
    ax.plot(t,solution, label=["x1","x2","x3","x4"])
    ax.legend()
    ax.grid()
    ax.set_xlabel("Time [s]")
    plt.show()

ex4()
