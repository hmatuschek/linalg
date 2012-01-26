import unittest
from numpy import *
import linalg


class DOTTest (unittest.TestCase):
    """ Tests the python-interface to Linalg::Blas::dot(). """
    
    def test_simple(self):
        x = array( [1,2,3,4], dtype=float)
        y = array( [5,6,7,8], dtype=float)
        self.assertEqual(linalg.ddot(x,y), dot(x,y))

    def test_array_view(self):
        x = array( [1,2,3,4], dtype=float)
        y = array( [5,6,7,8], dtype=float)
        self.assertEqual(linalg.ddot(x[:3],y[:3]), dot(x[:3], y[:3]))
        self.assertEqual(linalg.ddot(x[1:],y[1:]), dot(x[1:], y[1:]))
        self.assertEqual(linalg.ddot(x[1:3],y[1:3]), dot(x[1:3], y[1:3]))

    
