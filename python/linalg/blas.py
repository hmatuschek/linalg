import linalg._blas as _blas
import numpy as np


def dot(x,y):
	"""	Implements the some dot-product of two vectors x and y. """
	assert (x.dtype == y.dtype)
	if (np.float64 == x.dtype):
		return _blas.ddot(x,y)
	else:
		raise TypeError("Can not call dot() with vector of type {0}".format(x.dtype));
		
		
def nrm2(x):
	"""	Implements the 2-norm of a vector. """
	if (np.float64 == x.dtype):
		return _blas.nrm2(x)
	else:
		raise TypeError("Can not call nrm2() with vector of type {0}".format(x.dtype));
		
		
def gemv(alpha, A, x, beta, y):
	"""	Implements general matrix-vector product so y' = alphaA*x + beta*y. """
	assert (A.dtype == x.dtype)
	assert (A.dtype == x.dtype)
	
	if (np.float64 == x.dtype):
		return _blas.nrm2(x)
	else:
		raise TypeError("Can not call gemv() with arrays of type {0}".format(x.dtype));