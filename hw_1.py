# functions to help you
import numpy as np
import scipy.sparse as spsp
from scipy.sparse.linalg import spsolve
import scipy.integrate as integrate

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


myMesh = Mesh(np.array([0,.55,.75, .25, .8]))

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

functionSpace = V_h(myMesh)
#print(functionSpace.eval([1,2,3,4,5], 0))

class Function:
  def __init__(self, xi, v_h):
    self.xi  = xi
    self.v_h = v_h

  def __call__(self,x):
    # wrapper for calling eval in V_h
    
    # use the function defined in v_h
    return self.v_h.eval(self.xi, x)

function1 = Function([1,2,3,4,5], functionSpace)


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

#print(mass_matrix(functionSpace).todense())

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
        subS = np.array([[a/h, -a/h],[-a/h, h]])
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

def f(x):
    return x

#print(load_vector(functionSpace, f))



def source_assembler(v_h, f, u_dirichlet):
  # computing the load vector (use the function above)
    pass

  # extract the interval index for left boundary


  # compute the lenght of the interval


  # sample sigma at the middle point


  # update the source_vector



  # extract the interval index for the right boudanry



  # compute the length of the interval



  # sample sigma at the middle point



  # update the source_vector


  # return only the interior nodes
  #return b[1:-1]




#def solve_poisson_dirichelet(v_h, f, sigma, u_dirichlet=np.zeros((2)) ):
#    pass
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
    #S = 
    
    # we build the source
    #b = 

    # solve for the interior degrees of freedom
    #u_interior = spsolve(S,b)

    # concatenate the solution to add the boundary 
    # conditions
    #xi_u = np.concatenate([u_dirichlet[:1], 
    #                       u_interior, 
    #                       u_dirichlet[1:]])

    # return the function
    #return Function(xi_u, v_h)




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

print(p_h(functionSpace,f)(0))

if __name__ == "__main__":

  """ This is the main function, which will run 
  if you try to run this script, this code was provided 
  to you to help you debug your code. 
  """ 

  x = np.linspace(0,1,11)

  #mesh = Mesh(x)
  #v_h  = V_h(mesh)

  
  #f_load = lambda x: 2+0*x
  #xi = f_load(x) # linear function

  #u = Function(xi, v_h) 

  #assert np.abs(u(x[5]) - f_load(x[5])) < 1.e-6

  # check if this is projection
  #ph_f = p_h(v_h, f_load)
  #ph_f2 = p_h(v_h, ph_f)

  # 
  #assert np.max(ph_f.xi - ph_f2.xi) < 1.e-6

  # using analytical solution
  #u = lambda x : np.sin(4*np.pi*x)
  # building the correct source file
  #f = lambda x : (4*np.pi)**2*np.sin(4*np.pi*x)
  # conductivity is constant
  #sigma = lambda x : 1 + 0*x  

  #u_sol = solve_poisson_dirichelet(v_h, f, sigma)

  #err = lambda x: np.square(u_sol(x) - u(x))
  # we use an fearly accurate quadrature 
  #l2_err = np.sqrt(integrate.quad(err, 0.0,1.)[0])

  #print("L^2 error using %d points is %.6f"% (v_h.dim, l2_err))
  # this should be quite large

  # define a finer partition 
  #x = np.linspace(0,1,21)
  # init mesh and fucntion space
  #mesh = Mesh(x)
  #v_h  = V_h(mesh)

  #u_sol = solve_poisson_dirichelet(v_h, f, sigma)

  #err = lambda x: np.square(u_sol(x) - u(x))
  # we use an fearly accurate quadrature 
  #l2_err = np.sqrt(integrate.quad(err, 0.0,1.)[0])

  # print the error
  #print("L^2 error using %d points is %.6f"% (v_h.dim, l2_err))






