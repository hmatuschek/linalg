#include "lapack.hh"
#include "exception.hh"
#include <iostream>



/*
 * C interface to FORTRAN functions (direct access).
 */
extern "C" {
extern void dgetc2_(const int *n, double *a, const int *lda,
                    int *ipiv, int *jpiv, int *info);
}


using namespace Fluc;
using namespace Fluc::Linalg;


void
Linalg::getc2(Matrix<double> &a, Vector<int> &ipiv, Vector<int> &jpiv, size_t &rank)
{
  // Assert A being a N-by-N matrix:
  if (a.rows() != a.cols())
  {
    Linalg::IndexError err;
    err << "a must be a quadratic matrix!";
    throw err;
  }

  int n   = a.rows();
  int lda = a.leading_dimension();
  int info = 0;

  // Calling fortran:
  dgetc2_(&n, *a, &lda, *ipiv, *jpiv, &info);

  // Check for errors:
  if (0 > info)
  {
    Linalg::LapackError err("Lapack: ");
    err << "LAPACK returned error: " << -info;
    throw err;
  }

  rank = n-info;
}



LUDec::LUDec(const Matrix<double> &a)
  : lu(std::max(a.rows(), a.cols()), std::max(a.rows(), a.cols())),
    ipiv(std::max(a.rows(), a.cols())), jpiv(std::max(a.rows(), a.cols()))
{
  // Determine N as max(rows or columns of A):
  size_t N = std::max(a.rows(), a.cols());

  // Allocate (zeros) N-by-N matrix lu:
  this->lu.set(0,0, this->lu.rows(), this->lu.cols(), Matrix<double>::zeros(N,N));

  // And store A in upper-left "corner"
  lu.set(0,0, a.rows(), a.cols(), a);

  // perform LU decomposition
  getc2(this->lu, this->ipiv, this->jpiv, this->rank);

  //std::cerr << "LU rank:" << this->rank << " from " << N << std::endl;
}


Matrix<double>
LUDec::getL()
{
  Matrix<double> L(this->lu.copy());

  // Set upper triangular to 0.0 and diagonal to 1.0
  for (size_t i=0; i<this->lu.rows(); i++)
  {
    for (size_t j=i; j<this->lu.cols(); j++)
    {
      L(i,j) = 0.0;
    }

    L(i,i) = 1.0;
  }

  return L;
}


Matrix<double>
LUDec::getU()
{
  Matrix<double> U(this->lu.copy());

  // Set lower part 0:
  for (size_t i=1; i<this->lu.rows(); i++)
  {
    for (size_t j=0; j<i; j++)
    {
      U(i,j) = 0.0;
    }
  }

  return U;
}


Vector<int>
LUDec::getRowPerm()
{
  return this->ipiv;
}


Vector<int>
LUDec::getColPerm()
{
  return this->jpiv;
}


Matrix<double>
LUDec::getP()
{
  Matrix<double> P = Matrix<double>::zeros(this->lu.cols(), this->lu.rows());

  for (size_t i=0; i<this->lu.cols(); i++)
  {
    P(i, this->jpiv(i));
  }

  return P;
}


Matrix<double>
LUDec::getQ()
{
  Matrix<double> Q = Matrix<double>::zeros(this->lu.rows(), this->lu.cols());

  for (size_t i=0; i<this->lu.rows(); i++)
  {
    Q(this->ipiv(i), i);
  }

  return Q;
}
