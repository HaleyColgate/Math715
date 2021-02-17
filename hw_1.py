# functions to help you
import numpy as np
import scipy.sparse as spsp
from scipy.sparse.linalg import spsolve
import scipy.integrate as integrate

import matplotlib
from matplotlib import pyplot as plt
import pandas as pd
from pandas import Series, DataFrame

class Mesh:
  def __init__(self, points):
    # self.n_p  number of node points               type : int
    # self.p    array with the node points (sorted) type : np.array dim: (n_p)
    # self.n_s  number of segments                  type : int
    # self.s    array with indices of points per    type : np.array dim: (n_s, 2) 
    #           segment  
    # self.bc.  array with the indices of boundary  type : np.array dim: (2)
    #           points

    self.n_p = np.shape(points)[0]
    self.p   = np.sort(points)

    self.n_s = np.shape(points)[0]-1
    self.s   = np.transpose(np.array([self.p.tolist()[:-1],self.p.tolist()[1:]]))

    self.bc  = np.array([self.p[0],self.p[-1]])

class V_h:
  def __init__(self, mesh):
    # self.mesh Mesh object containg geometric info type: Mesh
    # self.dim  dimension of the space              type: int

    self.mesh = mesh
    self.dim  = self.mesh.n_s+1

  def eval(self, xi, x):
    """ evaluation of the piece wise local polynomial given by
       the coefficients xi, at the point x 
    """

    # compute the index of the interval in which x is contained
    for i in range(self.mesh.n_s):
        if x >= self.mesh.p[i]:
            xIndex = i + 1
    
    # compute the size of the interval
    intervalSize = self.mesh.s[xIndex-1][1] - self.mesh.s[xIndex-1][0] 

    # here return the value of the function 
    return xi[xIndex-1]*(self.mesh.s[xIndex-1][1]-x)/intervalSize + xi[xIndex]*(x-self.mesh.s[xIndex-1][0])/intervalSize


class Function:
  def __init__(self, xi, v_h):
    self.xi  = xi
    self.v_h = v_h

  def __call__(self,x):
    # wrapper for calling eval in V_h
    
    # use the function defined in v_h
    return self.v_h.eval(self.xi, x)


def mass_matrix(v_h):

    # sparse matrix easy to change sparsity pattern
    # this initializes an empty sparse matrix of 
    # size v_h.dim x v_h.dim
    M = spsp.lil_matrix((v_h.dim,v_h.dim))

    #for loop
    for i in range(v_h.mesh.n_s):
        # extract the indices
        x1 = v_h.mesh.s[i][0]
        x2 = v_h.mesh.s[i][1]

        # compute the lengh of the segment
        h = x2-x1

        # add the values to the matrix
        subM = np.array([[h/3, h/6],[h/6, h/3]])
        M[i:i+2, [i, i+1]] += subM
    return M


def stiffness_matrix(v_h, sigma):
    # matrix easy to change sparsity pattern
    S = spsp.lil_matrix((v_h.dim,v_h.dim))

    # for loop
    for i in range(v_h.mesh.n_s):
        # extract the indices
        x1 = v_h.mesh.s[i][0]
        x2 = v_h.mesh.s[i][1]
        
        # compute the length of the segments
        h = x2-x1

        # sample sigma
        a = sigma((x1+x2)/2)

        # update the stiffness matrix
        subS = np.array([[a/h, -a/h],[-a/h, a/h]])
        S[i:i+2, [i, i+1]] += subS

    return S


# show differences between Trapezoidal rule and Simpson rule
def load_vector(v_h, f):
    # allocate the vector
    b = np.zeros(v_h.dim)
    # for loop over the segments
    for i in range(v_h.mesh.n_s):
        # extracting the indices and
        x1 = v_h.mesh.s[i][0]
        x2 = v_h.mesh.s[i][1]
        m = (x1+x2)/2
            
        # computing the length of the interval 
        h = x2-x1

        # update b
        b[i:(i+2)] += np.array([h/6*f(x1)+h/3*f(m),h/6*f(x2)+h/3*f(m)])
    return b




def source_assembler(v_h, f, sigma, u_dirichlet):
  # computing the load vector (use the function above)
    b = load_vector(v_h, f)

  # extract the interval index for left boundary
    x1L = v_h.mesh.s[0][0]
    x2L = v_h.mesh.s[0][1]
    mL = (x1L+x2L)/2

  # compute the lenght of the interval
    hL = x2L - x1L

  # sample sigma at the middle point
    lsigma = sigma(mL)

  # update the source_vector
    b[1] += u_dirichlet[0] * lsigma / hL


  # extract the interval index for the right boudanry
    x1R = v_h.mesh.s[-1][0]
    x2R = v_h.mesh.s[-1][1]
    mR = (x1R + x2R)/2

  # compute the length of the interval
    hR = x2R - x1R

  # sample sigma at the middle point
    rsigma = sigma(hR)

  # update the source_vector
    b[-2] += u_dirichlet[0] * rsigma / hR

  # return only the interior nodes
    return b[1:-1]




def solve_poisson_dirichelet(v_h, f, sigma, u_dirichlet=np.zeros((2)) ):
    """ function to solbe the Poisson equation with 
    Dirichlet boundary conditions
    input:  v_h         function space
            f           load (python function)
            sigma       conductivity
            u_dirichlet boundary conditions
    output: u           approximation (Function class)
    """  
    # we compute the stiffness matrix, we only use the  
    # the interior dof, and we need to convert it to 
    # a csc_matrix
    S = spsp.csc_matrix(stiffness_matrix(v_h, sigma))[1:-1,1:-1]
    
    # we build the source
    b = source_assembler(v_h, f, sigma, u_dirichlet)

    # solve for the interior degrees of freedom
    u_interior = spsolve(S,b)

    # concatenate the solution to add the boundary 
    # conditions
    xi_u = np.concatenate([u_dirichlet[:1], 
                           u_interior, 
                           u_dirichlet[1:]])

    # return the function
    return Function(xi_u, v_h)




def pi_h(v_h, f):
    """interpolation function
        input:  v_h   function space
                f     function to project
        output: pih_f function that is the interpolation 
                    of f into v_h
    """
    #pi_h_f = 
    pass

#  return pi_h_f




def p_h(v_h, f):
    """projection function
        input:  v_h   function space
                f     function to project
        output: ph_f  function that is the projection 
                    of f into v_h
    """
    # compute load vector
    b = load_vector(v_h, f)

    # compute Mass matrix and convert it to csc type
    M = spsp.csc_matrix(mass_matrix(v_h))

    # solve the system
    xi = spsolve(M,b)

    # create the new function (this needs to be an instance)
    # of a Function class
    ph_f = Function(xi, v_h)

    return ph_f

#Problem c part (b)
def compResiduals(mesh, f):
    eta = []
    maxEta = 0
    for i in range(mesh.n_s):
        x0 = mesh.s[i][0]
        x1 = mesh.s[i][1]
        h = x1 - x0
        m = (x0+x1)/2
        quad = (h/6)*((f(x0))**2+4*(f(m))**2+(f(x1))**2)
        eta.append(h*np.sqrt(quad))
    return eta

f = lambda x: np.exp(-50*(x-1/2)**2)

def refine_grid(tolerance, f, low, high, alpha):
    #fine grid solve:
    fineMesh = Mesh(np.linspace(low,high,1001))
    v_hFine = V_h(fineMesh)
    sigma = lambda x : 1 + 0*x  

    fine_sol = solve_poisson_dirichelet(v_hFine, f, sigma)

    x = Mesh([low, high])
    error = 2*tolerance
    while error > tolerance:
        xNew = x.p.tolist()
        resids = compResiduals(x, f)
        error = 0
        refineBound = alpha*max(resids)
        for i in range(len(resids)):
            error += resids[i]**2
            if resids[i] > refineBound:
                x0 = x.s[i][0]
                x1 = x.s[i][1]
                m = (x0+x1)/2
                xNew.append(m)
        x = Mesh(xNew)
        v_h = V_h(x)
        u = solve_poisson_dirichelet(v_h, f, sigma, [0,0])
        errF = lambda x: np.square(fine_sol(x) - u(x))
        error = np.sqrt(integrate.quad(errF, 0.0,1., limit = 100)[0])

    return x, error
    
alpha = .9
newMesh, error = refine_grid(.001, f, 0, 1, alpha)
print(newMesh.p, error, len(newMesh.p), 'alpha:', alpha)

alpha = 0
newMesh, error = refine_grid(.001, f, 0, 1, alpha)
print(newMesh.p, error, len(newMesh.p), 'alpha:', alpha)



if __name__ == "__main__":

  """ This is the main function, which will run 
  if you try to run this script, this code was provided 
  to you to help you debug your code. 
  """ 

  x = np.linspace(0,1,11)

  mesh = Mesh(x)
  v_h  = V_h(mesh)

  
  f_load = lambda x: 2+0*x
  xi = f_load(x) # linear function

  u = Function(xi, v_h) 

  assert np.abs(u(x[5]) - f_load(x[5])) < 1.e-6

  # check if this is projection
  ph_f = p_h(v_h, f_load)
  ph_f2 = p_h(v_h, ph_f)

  # 
  assert np.max(ph_f.xi - ph_f2.xi) < 1.e-6

  # using analytical solution
  u = lambda x : -np.sin(np.pi*x)/(np.pi)**2 + 1
  # building the correct source file
  f = lambda x : np.sin(x)
  # conductivity is constant
  sigma = lambda x : 1 + 0*x  

  u_sol = solve_poisson_dirichelet(v_h, f, sigma, [1,1])

  err = lambda x: np.square(u_sol(x) - u(x))
  # we use an fearly accurate quadrature 
  l2_err = np.sqrt(integrate.quad(err, 0.0,1., limit = 100)[0])

  print("L^2 error using %d points is %.6f"% (v_h.dim, l2_err))
  # this should be quite large

  # define a finer partition 
  x = np.linspace(0,1,21)
  # init mesh and fucntion space
  mesh = Mesh(x)
  v_h  = V_h(mesh)

  u_sol = solve_poisson_dirichelet(v_h, f, sigma, [1,1])

  err = lambda x: np.square(u_sol(x) - u(x))
  # we use an fearly accurate quadrature 
  l2_err = np.sqrt(integrate.quad(err, 0.0,1., limit = 500)[0])

  # print the error
  print("L^2 error using %d points is %.6f"% (v_h.dim, l2_err))

#B part 1
def showHconvergence(numSteps):
    hvals = []
    errs = []
    for step in numSteps:
        hvals.append(np.log(1/(step-1)))
        x = np.linspace(0,1,step)
        mesh = Mesh(x)
        v_h  = V_h(mesh)
        f = lambda x : (4*np.pi)**2*np.sin(4*np.pi*x)
        ph_f = p_h(v_h, f)

        err = lambda x: np.square(f(x) - ph_f(x))
        # we use an fearly accurate quadrature 
        errs.append(np.log(integrate.quad(err, 0.0,1., limit = 500)[0]))
    errSeries = Series(errs, hvals)
    ax = errSeries.plot()
    ax.set_title("Mesh Size vs L2 Error for Projection")
    ax.set_ylabel("Log of Squared L2 Error of Projection")
    ax.set_xlabel("Log of Mesh Size")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.savefig('ProblemB_hGraph.png', bbox_inches='tight')
    
#numSteps = range(6, 101, 1)
#showHconvergence(numSteps)

#B part 3
def counterHconvergence(numSteps):
    hvals = []
    errs = []
    for step in numSteps:
        hvals.append(np.log(1/(step-1)))
        x = np.linspace(0,1,step)
        mesh = Mesh(x)
        v_h  = V_h(mesh)
        f = lambda x : 2*x
        ph_f = p_h(v_h, f)

        err = lambda x: np.square(f(x) - ph_f(x))
        # we use an fearly accurate quadrature 
        errs.append(np.log(integrate.quad(err, 0.0,1., limit = 500)[0]))
    errSeries = Series(errs, hvals)
    ax = errSeries.plot()
    ax.set_title("Mesh Size vs L2 Error for Projection")
    ax.set_ylabel("Log of Squared L2 Error of Projection")
    ax.set_xlabel("Log of Mesh Size")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.savefig('ProblemB_hCounter_Graph.png', bbox_inches='tight')
#counterHconvergence(numSteps)

#B part 2
def showDoubleDerivative_convergence(kVals):
    doubleDeriv = []
    errs = []
    step = 51
    for k in kVals:
        k = k/100
        l2norm = 8*k**4-4*k**3*np.sin(2*k)*np.cos(2*k)
        doubleDeriv.append(np.sqrt(l2norm))
        x = np.linspace(0,1,step)
        mesh = Mesh(x)
        v_h  = V_h(mesh)
        f = lambda x : np.sin(2*k*x)
        ph_f = p_h(v_h, f)
        err = lambda x: np.square(f(x) - ph_f(x))
        # we use an fearly accurate quadrature 
        errs.append(np.sqrt(integrate.quad(err, 0.0,1., limit = 500)[0]))
    errSeries = Series(errs, doubleDeriv)
    ax = errSeries.plot()
    ax.set_title("L2 Norm of the Second Derivative vs L2 Error")
    ax.set_ylabel("L2 Error of Projection")
    ax.set_xlabel("L2 Norm of Second Derivative")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.savefig('ProblemB_fGraph.png')


#kVals = range(1, 100)
#showDoubleDerivative_convergence(kVals)


#c(a)
def derivConvergence(numSteps):
    hvals = []
    errs = []
    sigma = lambda x : 1 + 0*x  
    u_dirichlet = np.zeros((2))
    for step in numSteps:
        hvals.append(np.log(np.pi/(step-1)))
        x = np.linspace(0, np.pi, step)
        mesh = Mesh(x)
        v_h  = V_h(mesh)
        # using analytical solution
        uPrime = lambda x : np.cos(x)
        # building the correct source file
        f = lambda x : np.sin(x)
        u_sol = solve_poisson_dirichelet(v_h, f, sigma)
        err = 0
        for i in range(v_h.mesh.n_s):
            lower = v_h.mesh.s[i][0] 
            upper = v_h.mesh.s[i][1]
            h = upper - lower
            slope = (u_sol(upper) - u_sol(lower))/h
            derivFunction = lambda x : np.square(uPrime(x) - slope)
            # we use an fearly accurate quadrature 
            err += integrate.quad(derivFunction, lower, upper, limit = 100)[0]
        errs.append(np.log(err))
    errSeries = Series(errs, hvals)
    ax = errSeries.plot()
    ax.set_title("Mesh Size vs L2 Error of (u_h - u)'")
    ax.set_ylabel("Log of Squared L2 Error")
    ax.set_xlabel("Log of Mesh Size")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.savefig('ProblemC_hConvergence', bbox_inches='tight')
    
numSteps = range(6, 10, 1)
derivConvergence(numSteps)
