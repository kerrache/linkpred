/*
 * This file is part of LinkPred.
 *
 * LinkPred: A high performance library for link prediction in complex networks.
 * Copyright (C) 2017  by Said Kerrache.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * \file
 * @ingroup Numerical
 * @brief Contains the implementation of a full matrix.
 */

#ifndef GFMATRIX_HPP_
#define GFMATRIX_HPP_

#include <cstdlib>
#include <vector>
#include "linkpred/numerical/linear/vec.hpp"

namespace LinkPred {

class Vec;

/**
 * @brief Generalized full matrix. The storage scheme used is column-major.
 */
class GFMatrix {

	/**
	 * GFMatrix-vector multiplication.
	 * @param mat The matrix.
	 * @param vec The vector that will be multiplied by the matrix.
	 * @return The resulting vector.
	 */
	friend Vec operator*(const GFMatrix &mat, const Vec &vec);

	/**
	 * GFMatrix * GFMatrix.
	 * @param mat1 The first matrix.
	 * @param mat2 The second matrix.
	 * @return The resulting matrix.
	 */
	friend GFMatrix operator*(const GFMatrix & mat1, const GFMatrix & mat2);

	/**
	 * GFMatrix + GFMatrix.
	 * @param mat1 The first matrix.
	 * @param mat2 The second matrix.
	 * @return The resulting matrix.
	 */
	friend GFMatrix operator+(const GFMatrix & mat1, const GFMatrix & mat2);

	/**
	 * GFMatrix - GFMatrix.
	 * @param mat1 The first matrix.
	 * @param mat2 The second matrix.
	 * @return The resulting matrix.
	 */
	friend GFMatrix operator-(const GFMatrix & mat1, const GFMatrix & mat2);

	/**
	 * GFMatrix / GFMatrix (element-wise division).
	 * @param mat1 The first matrix.
	 * @param mat2 The second matrix.
	 * @return The resulting matrix.
	 */
	friend GFMatrix operator/(const GFMatrix & mat1, const GFMatrix & mat2);

	/**
	 * GFMatrix + Vec.
	 * @param mat1 The first matrix.
	 * @param mat2 The second matrix (diagonal).
	 * @return The resulting matrix.
	 */
	friend GFMatrix operator+(const GFMatrix & mat1, const Vec & mat2);

	/**
	 * Vec + GFMatrix.
	 * @param mat1 The first matrix (diagonal).
	 * @param mat2 The second matrix.
	 * @return The resulting matrix.
	 */
	friend GFMatrix operator+(const Vec & mat1, const GFMatrix & mat2);

	/**
	 * GFMatrix - Vec.
	 * @param mat1 The first matrix.
	 * @param mat2 The second matrix (diagonal).
	 * @return The resulting matrix.
	 */
	friend GFMatrix operator-(const GFMatrix & mat1, const Vec & mat2);

	/**
	 * Vec - GFMatrix.
	 * @param mat1 The first matrix (diagonal).
	 * @param mat2 The second matrix.
	 * @return The resulting matrix.
	 */
	friend GFMatrix operator-(const Vec & mat1, const GFMatrix & mat2);

	/**
	 * a * GFMatrix.
	 * @param a scalar.
	 * @param mat The second matrix.
	 * @return The resulting matrix.
	 */
	friend GFMatrix operator*(double a, const GFMatrix & mat);

	/**
	 * GFMatrix * a.
	 * @param mat The second matrix.
	 * @param a scalar.
	 * @return The resulting matrix.
	 */
	friend GFMatrix operator*(const GFMatrix & mat, double a);

	/**
	 * a + GFMatrix.
	 * @param a scalar.
	 * @param mat The second matrix.
	 * @return The resulting matrix.
	 */
	friend GFMatrix operator+(double a, const GFMatrix & mat);

	/**
	 * GFMatrix + a.
	 * @param mat The second matrix.
	 * @param a A scalar.
	 * @return The resulting matrix.
	 */
	friend GFMatrix operator+(const GFMatrix & mat, double a);

	/**
	 * a - GFMatrix.
	 * @param a scalar.
	 * @param mat The second matrix.
	 * @return The resulting matrix.
	 */
	friend GFMatrix operator-(double a, const GFMatrix & mat);

	/**
	 * GFMatrix - a.
	 * @param mat The second matrix.
	 * @param a scalar.
	 * @return The resulting matrix.
	 */
	friend GFMatrix operator-(const GFMatrix & mat, double a);

protected:
	std::string name = ""; /**< Name of the matrix. */
	int m = 0; /**< Number of rows. */
	int n = 0; /**< Number of columns. */
	double * values = nullptr; /**< The values. */

	/**
	 * Matrix multiplication. Both matrices are not transposed.
	 * @param values1 Values of matrix 1.
	 * @param values2 Values of matrix 2.
	 * @param values3 Values of the result matrix.
	 * @param m1 Number of rows of matrix 1.
	 * @param n1 Number of columns of matrix 1.
	 * @param m2 Number of rows of matrix 2.
	 * @param n2 Number of columns of matrix 2.
	 */
	static void matMultFF(double const *values1, double const *values2,
			double * values3, int m1, int n1, int m2, int n2);

	/**
	 * Matrix multiplication. Second matrix only is transposed.
	 * @param values1 Values of matrix 1.
	 * @param values2 Values of matrix 2.
	 * @param values3 Values of the result matrix.
	 * @param m1 Number of rows of matrix 1.
	 * @param n1 Number of columns of matrix 1.
	 * @param m2 Number of rows of matrix 2.
	 * @param n2 Number of columns of matrix 2.
	 */
	static void matMultFT(double const *values1, double const *values2,
			double * values3, int m1, int n1, int m2, int n2);

	/**
	 * Matrix multiplication. First matrix only is transposed.
	 * @param values1 Values of matrix 1.
	 * @param values2 Values of matrix 2.
	 * @param values3 Values of the result matrix.
	 * @param m1 Number of rows of matrix 1.
	 * @param n1 Number of columns of matrix 1.
	 * @param m2 Number of rows of matrix 2.
	 * @param n2 Number of columns of matrix 2.
	 */
	static void matMultTF(double const *values1, double const *values2,
			double * values3, int m1, int n1, int m2, int n2);

	/**
	 * Matrix multiplication. Both matrices are transposed.
	 * @param values1 Values of matrix 1.
	 * @param values2 Values of matrix 2.
	 * @param values3 Values of the result matrix.
	 * @param m1 Number of rows of matrix 1.
	 * @param n1 Number of columns of matrix 1.
	 * @param m2 Number of rows of matrix 2.
	 * @param n2 Number of columns of matrix 2.
	 */
	static void matMultTT(double const *values1, double const *values2,
			double * values3, int m1, int n1, int m2, int n2);

public:

	/**
	 * Default constructor.
	 */
	GFMatrix() = delete;

	/**
	 * Constructor.
	 * @param m The number of rows.
	 * @param n The number of columns.
	 * @param initZero If set to true, the matrix is initialized to zero.
	 */
	GFMatrix(int m, int n, bool initZero = false);

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	GFMatrix(GFMatrix const & that);

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	GFMatrix & operator =(GFMatrix const & that);

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	GFMatrix(GFMatrix && that);

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	GFMatrix & operator =(GFMatrix && that);

	/**
	 * @param that The other matrix.
	 * @return True if the two matrices are equal.
	 */
	bool operator==(const GFMatrix &that) const;

	/**
	 * @param that The other matrix.
	 * @return True if the two matrices are not equal.
	 */
	bool operator!=(const GFMatrix &that) const;

	/**
	 * Set a matrix entry.
	 * @param i Row index.
	 * @param j Column index.
	 * @param v The new value.
	 */
	inline void set(int i, int j, double v) {
		values[i + j * m] = v;
	}

	/**
	 * Get a matrix entry.
	 * @param i Row index.
	 * @param j Column index.
	 */
	inline double get(int i, int j) const {
		return values[i + j * m];
	}

	/**
	 * Set row vector.
	 * @param i The row index.
	 * @param vec The row vector to set.
	 */
	void setRow(int i, Vec const & vec);

	/**
	 * @param i The row index.
	 * @return The row vector i.
	 */
	Vec getRow(int i) const;

	/**
	 * Set column vector.
	 * @param j The column index.
	 * @param vec The column vector to set.
	 */
	void setCol(int j, Vec const & vec);

	/**
	 * @param j The column index.
	 * @return The column vector j.
	 */
	Vec getCol(int j) const;

	/**
	 * @param rowInd The row indexes.
	 * @return The matrix composed from the specified rows.
	 */
	GFMatrix getRows(std::vector<int> const & rowInd) const;

	/**
	 * @param colInd The column indexes.
	 * @return The matrix composed from the specified columns.
	 */
	GFMatrix getCols(std::vector<int> const & colInd) const;

	/**
	 * @return The sum of the columns of the matrix.
	 */
	Vec sumCols() const;

	/**
	 * @return The sum of the rows of the matrix.
	 */
	Vec sumRows() const;

	/**
	 * @return The sum of all elements.
	 */
	double sum();

	/**
	 * Print matrix to standard output.
	 */
	void print() const;

	/**
	 * Print matrix to standard output adding the name.
	 * @param name Name added to the output.
	 */
	void print(std::string name) const;

	/**
	 * @return The diagonal of the matrix.
	 */
	Vec diag() const;

	/**
	 * Set diagonal.
	 * @param v The value to be put in the diagonal.
	 */
	void setDiag(double v);

	/**
	 * Matrix-matrix multiplication with possibility of transposing.
	 * @param mat1 First matrix.
	 * @param mat2 Second matrix.
	 * @param trans1 Transposing mat1.
	 * @param trans2 Transposing mat2.
	 * @return The resulting matrix.
	 */
	static GFMatrix mult(GFMatrix const & mat1, GFMatrix const & mat2,
			bool trans1 = false, bool trans2 = false);

	/**
	 * Element-wise matrix multiplication.
	 * @param mat1 First matrix.
	 * @param mat2 Second matrix.
	 * @return The resulting matrix.
	 */
	static GFMatrix elemMult(GFMatrix const & mat1, GFMatrix const & mat2);

	/**
	 * Vector'-Vector multiplication.
	 * @param v1 First vector.
	 * @param v2 Second vector.
	 * @return The resulting matrix.
	 */
	static GFMatrix mult(Vec const & v1, Vec const & v2);

	/**
	 * Remove NaN entries.
	 */
	void removeNaN();

	/**
	 * @return The number of rows.
	 */
	int getM() const {
		return m;
	}

	/**
	 * @return The number of columns.
	 */
	int getN() const {
		return n;
	}

	/**
	 * @return The name of the matrix.
	 */
	std::string getName() const {
		return name;
	}

	/**
	 * Set the name of the matrix.
	 * @param name The new matrix name.
	 */
	void setName(const std::string& name) {
		this->name = name;
	}

	/**
	 * Destructor.
	 */
	virtual ~GFMatrix();
};

} /* namespace LinkPred */

#endif /* GFMATRIX_HPP_ */
