import unittest
import numpy as np
from numpy import array, sqrt, sum
import linalg

class NRM2Test (unittest.TestCase):

    def test_double(self):
        x = array([1,2,3,4,5], dtype=float)

        self.assertEqual(linalg.dnrm2(x), sqrt(sum(x**2)))
        self.assertEqual(linalg.dnrm2(x[:4]), sqrt(sum(x[:4]**2)))
        self.assertEqual(linalg.dnrm2(x[1:]), sqrt(sum(x[1:]**2)))
        self.assertEqual(linalg.dnrm2(x[1:4]), sqrt(sum(x[1:4]**2)))
