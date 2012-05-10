import unittest
from numpy import *
import linalg._blas
import sys


class DOTTest (unittest.TestCase):
	""" Tests the python-interface to Linalg::Blas::dot(). """
    
	def test_simple(self):
		x = array( [1,2,3,4], dtype=float)
		y = array( [5,6,7,8], dtype=float)
		self.assertEqual(linalg._blas.ddot(x,y), dot(x,y))

	def test_array_view(self):
		x = array( [1,2,3,4], dtype=float)
		y = array( [5,6,7,8], dtype=float)
		self.assertEqual(linalg._blas.ddot(x[:3],y[:3]), dot(x[:3], y[:3]))
		self.assertEqual(linalg._blas.ddot(x[1:],y[1:]), dot(x[1:], y[1:]))
		self.assertEqual(linalg._blas.ddot(x[1:3],y[1:3]), dot(x[1:3], y[1:3]))

	def test_refcount(self):
		x = array( [1,2,3,4], dtype=float)
		y = array( [5,6,7,8], dtype=float)
		x_count, y_count = sys.getrefcount(x), sys.getrefcount(y);
		linalg._blas.ddot(x,y)
		self.assertEqual(x_count, sys.getrefcount(x))
		self.assertEqual(y_count, sys.getrefcount(y))