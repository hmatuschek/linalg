import unittest
from numpy import *
import linalg._blas



class GEMVTest (unittest.TestCase):
    def setUp(self):
        self.Arow = array( [[ 1, 2, 3],
                            [ 4, 5, 6],
                            [ 7, 8, 9]], dtype=float, order='C' )

        self.Acol = array( self.Arow, dtype=float, order='F' )
        
        self.x = array( [10, 11, 12], dtype=float )


    def testRowOrder(self):
        y = empty( (3,) )

        linalg._blas.dgemv(1., self.Arow, self.x, 0., y);
        self.assertTrue((y == inner(self.Arow, self.x)).all())

        linalg._blas.dgemv(1., self.Arow.T, self.x, 0., y);
        self.assertTrue((y == inner(self.Arow.T, self.x)).all())

        linalg._blas.dgemv(1., self.Arow[:,:2], self.x[:2], 0., y);
        self.assertTrue((y == inner(self.Arow[:,:2], self.x[:2])).all())
        
        linalg._blas.dgemv(1., self.Arow[:,1:], self.x[1:], 0., y);
        self.assertTrue((y == inner(self.Arow[:,1:], self.x[1:])).all())

        linalg._blas.dgemv(1., self.Arow[:2,:].T, self.x[:2], 0., y);
        self.assertTrue((y == inner(self.Arow[:2,:].T, self.x[:2])).all())
        
        linalg._blas.dgemv(1., self.Arow[1:,:].T, self.x[1:], 0., y);
        self.assertTrue((y == inner(self.Arow[1:,:].T, self.x[1:])).all())


    def testColOrder(self):
        y = empty( (3,) )
        linalg._blas.dgemv(1., self.Acol, self.x, 0., y);
        self.assertTrue((y == inner(self.Acol, self.x)).all())

        linalg._blas.dgemv(1., self.Acol.T, self.x, 0., y);
        self.assertTrue((y == inner(self.Acol.T, self.x)).all())

        linalg._blas.dgemv(1., self.Acol[:,:2], self.x[:2], 0., y);
        self.assertTrue((y == inner(self.Acol[:,:2], self.x[:2])).all())
        
        linalg._blas.dgemv(1., self.Acol[:,1:], self.x[1:], 0., y);
        self.assertTrue((y == inner(self.Acol[:,1:], self.x[1:])).all())

        linalg._blas.dgemv(1., self.Acol[:2,:].T, self.x[:2], 0., y);
        self.assertTrue((y == inner(self.Acol[:2,:].T, self.x[:2])).all())
        
        linalg._blas.dgemv(1., self.Acol[1:,:].T, self.x[1:], 0., y);
        self.assertTrue((y == inner(self.Acol[1:,:].T, self.x[1:])).all())
    
