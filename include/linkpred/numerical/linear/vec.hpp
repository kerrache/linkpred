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
 * @brief Contains the implementation of a vector (in the sense of linear algebra).
 */

#ifndef INCLUDE_VEC_HPP_
#define INCLUDE_VEC_HPP_

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include <cstdlib>

namespace LinkPred {

/**
 * Forward declarations.
 */
class SMatrix;
class Vec;
class GFMatrix;

/**
 * Vector addition.
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @return v1 + v2.
 */
Vec operator +(Vec const & v1, Vec const & v2);

/**
 * Vector subtraction.
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @return v1 - v2.
 */
Vec operator -(Vec const & v1, Vec const & v2);

/**
 * Vector element-wise multiplication.
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @return v1 .* v2.
 */
Vec operator *(Vec const & v1, Vec const & v2);

/**
 * Vector element-wise division.
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @return v1 ./ v2.
 */
Vec operator /(Vec const & v1, Vec const & v2);

/**
 * Scalar-vector addition.
 * @param a The scalar.
 * @param v The vector.
 * @return a + v.
 */
Vec operator +(double a, Vec const & v);

/**
 * Scalar-vector subtraction.
 * @param a The scalar.
 * @param v The vector.
 * @return a - v.
 */
Vec operator -(double a, Vec const & v);

/**
 * Scalar-vector multiplication.
 * @param a The scalar.
 * @param v The vector.
 * @return a * v.
 */
Vec operator *(double a, Vec const & v);

/**
 * Scalar-vector division.
 * @param a The scalar.
 * @param v The vector.
 * @return a / v.
 */
Vec operator /(double a, Vec const & v);

/**
 * Vector-scalar addition.
 * @param v The vector.
 * @param a The scalar.
 * @return v + a.
 */
Vec operator +(Vec const & v, double a);

/**
 * Vector-scalar subtraction.
 * @param v The vector.
 * @param a The scalar.
 * @return v - a.
 */
Vec operator -(Vec const & v, double a);

/**
 * Vector-scalar multiplication.
 * @param v The vector.
 * @param a The scalar.
 * @return v * a.
 */
Vec operator *(Vec const & v, double a);

/**
 * Vector-scalar division.
 * @param v The vector.
 * @param a The scalar.
 * @return v / a.
 */
Vec operator /(Vec const & v, double a);

/**
 * Dot product.
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @return v1' * v2.
 */
double operator ^(Vec const & v1, Vec const & v2);

/**
 * @brief This class represent a vector (in the sense of linear algebra).
 */
class Vec {
	friend Vec operator +(Vec const & v1, Vec const & v2);
	friend Vec operator -(Vec const & v1, Vec const & v2);
	friend Vec operator *(double a, Vec const & v);
	friend Vec operator *(Vec const & v1, Vec const & v2);
	friend Vec operator /(Vec const & v1, Vec const & v2);
	friend Vec operator +(double a, Vec const & v);
	friend Vec operator -(double a, Vec const & v);
	friend Vec operator *(double a, Vec const & v);
	friend Vec operator /(double a, Vec const & v);
	friend Vec operator +(Vec const & v, double a);
	friend Vec operator -(Vec const & v, double a);
	friend Vec operator *(Vec const & v, double a);
	friend Vec operator /(Vec const & v, double a);
	friend double operator ^(Vec const & v1, Vec const & v2);

	friend GFMatrix operator+(const GFMatrix & mat1, const Vec & mat2);
	friend GFMatrix operator+(const Vec & mat1, const GFMatrix & mat2);
	friend GFMatrix operator-(const GFMatrix & mat1, const Vec & mat2);
	friend GFMatrix operator-(const Vec & mat1, const GFMatrix & mat2);

	friend class GFMatrix;

	/**
	 * GFMatrix-vector multiplication.
	 * @param mat The matrix.
	 * @param vec The vector that will be multiplied by the matrix.
	 * @return The resulting vector.
	 */
	friend Vec operator*(const GFMatrix &mat, const Vec &vec);

protected:
	int n; /**< The dimension. */
	double *values; /**< The values. */
	bool constant; /**< True if the vector is constant. */
	// TODO use constant in all methods.
public:

	/**
	 * Default constructor.
	 */
	Vec();

	/**
	 * Constructor.
	 * @param _n The dimension.
	 */
	Vec(int _n);

	/**
	 * Constructor.
	 * @param n The dimension.
	 * @param vals The values.
	 */
	Vec(int n, double * vals);

	/**
	 * Create a new vector from an existing one by copying specified elements.
	 * @param v The vector from which the data is taken.
	 * @param ind Index of the elements to be copied.
	 */
	Vec(Vec const & v, std::shared_ptr<std::vector<int>> ind);

	/**
	 * Create a new vector from two existing ones by concatenation.
	 * @param v1 The first vector.
	 * @param v2 The second vector.
	 */
	Vec(Vec const & v1, Vec const & v2);

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	Vec(Vec const & that);

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	Vec & operator =(Vec const & that);

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	Vec(Vec && that);

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	Vec & operator =(Vec && that);

	/**
	 * @return The size of the vector.
	 */
	int size() const {
		return n;
	}

	/**
	 * @return Array of values.
	 */
	double* getValues() const {
		return values;
	}

	/**
	 * @return True if all entries of the vector are equal.
	 */
	bool isConstant() const {
		return constant;
	}

	/**
	 * @param constant New value for constant.
	 */
	void setConstant(bool constant) {
		this->constant = constant;
	}

	/**
	 * Check if the vector is constant (all entries are equal) up to a specific tolerance.
	 * @param tol The tolerance.
	 * @return True if the vector is constant.
	 */
	bool makeConstant(double tol);

	/**
	 * Resize the vector.
	 * @param newSize The new size. If this is smaller than the current size, nothing happens.
	 */
	void resize(int newSize);

	/**
	 * Resize the vector.
	 * @param newSize The new size. If this is smaller than the current size, nothing happens.
	 * @param val Initialization value for new entries.
	 */
	void resize(int newSize, double val);

	/**
	 * Overloads operator ==.
	 * 	@param other The other vector.
	 * 	@return True if the the two vectors are equal.
	 */
	bool operator==(const Vec &other) const;

	/**
	 * Overloads operator !=.
	 * 	@param other The other vector.
	 * 	@return True if the the two vectors are not equal.
	 */
	bool operator!=(const Vec &other) const;

	/**
	 * Overloaded subscript.
	 * 	@param k index.
	 * 	@return reference to the k'th entry.
	 */
	inline double& operator[](const int k) {
		return values[k];
	}

	/**
	 * Overloaded subscript.
	 * 	@param k index.
	 * 	@return the k'th entry.
	 */
	inline double operator[](const int k) const {
		return values[k];
	}

	/**
	 * @return Indices of non-zero elements.
	 */
	std::vector<int> findNz() const;

	/**
	 * @param val The avlue we are lookign for.
	 * @return Indices of elements equal to val.
	 */
	std::vector<int> find(double val) const;

	/**
	 * Reset vector to 0.
	 */
	void reset();

	/**
	 * Reset vector to a specified value.
	 * @param val All vector elements are set to val.
	 */
	void reset(double val);

	/**
	 * Read the values from a text file.
	 * 	@param in Input file.
	 */
	void read(std::ifstream *in);

	/**
	 * Write values to file.
	 * 	@param out Output file.
	 */
	void write(std::ofstream *out) const;

	/**
	 * Read the values from a text file.
	 * 	@param fileName Input file name.
	 */
	void read(std::string fileName);

	/**
	 * Write values to file.
	 * 	@param fileName Output file name.
	 */
	void write(std::string fileName) const;

	/**
	 * Computes the norm of the vector.
	 * 	@return The norm of the vector.
	 */
	double norm() const;

	/**
	 * Computes the mean of the vector.
	 * 	@return The mean of the vector.
	 */
	double mean() const;

	/**
	 * @return The minimum absolute value in the vector.
	 */
	double amin() const;

	/**
	 * @return The maximum absolute value in the vector.
	 */
	double amax() const;

	/**
	 * @return The minimum value in the vector.
	 */
	double min() const;

	/**
	 * @return The maximum value in the vector.
	 */
	double max() const;

	/**
	 * @return The index of the minimum value in the vector.
	 */
	int imin() const;

	/**
	 * @return The index of the maximum value in the vector.
	 */
	int imax() const;

	/**
	 * @return The index of the first minimum value in the vector.
	 */
	int ifmin() const;

	/**
	 * @return The index of the first maximum value in the vector.
	 */
	int ifmax() const;

	/**
	 * @return The index of the last minimum value in the vector.
	 */
	int ilmin() const;

	/**
	 * @return The index of the last maximum value in the vector.
	 */
	int ilmax() const;

	/**
	 * Center vector.
	 */
	void center();

	/**
	 * Computes the norm of the vector.
	 * @param lim The limit.
	 * 	@return The norm of the vector.
	 */
	double norm(int lim) const;

	/**
	 * Computes the mean of the vector.
	 * 	@return The mean of the vector.
	 */
	double mean(int lim) const;

	/**
	 * @return The minimum absolute value in the vector up to lim number of elements.
	 * @param lim The limit.
	 */
	double amin(int lim) const;

	/**
	 * @return The maximum absolute value in the vector up to lim number of elements.
	 * @param lim The limit.
	 */
	double amax(int lim) const;

	/**
	 * @return The minimum value in the vector up to lim number of elements.
	 * @param lim The limit.
	 */
	double min(int lim) const;

	/**
	 * @return The maximum value in the vector up to lim number of elements.
	 * @param lim The limit.
	 */
	double max(int lim) const;

	/**
	 * Center vector up to lim.
	 * @param lim The limit.
	 */
	void center(int lim);

	/**
	 * Computes y = a * x + y, where y is the current object (this) and x is a full vector.
	 * @param a The value of a.
	 * @param x The other vector.
	 */
	void axpy(double a, Vec const & x);

	/**
	 * Scale the vector by a scalar.
	 * @param a The scale.
	 */
	void scale(double a);

	/**
	 * Scale part of the vector.
	 * @param start The starting index.
	 * @param nb Number of elements to be scaled.
	 * @param a The scale factor.
	 */
	void partScale(int start, int nb, double a);

	/**
	 * Scatter a vector into the current vector starting from given position.
	 * @param start The starting position.
	 * @param v The vector to be scattered.
	 * @param ind The index. Only entries in this array are affected by this operation.
	 */
	void scatter(int start, Vec const & v, int const * ind);

	/**
	 * @return The sum of the elements of the vector.
	 */
	double sum() const;

	/**
	 * Print vector to standard output.
	 */
	void print() const;

	/**
	 * Print vector to standard output.
	 * @param header Header.
	 */
	void print(std::string header) const;

	/**
	 * Print vector to output stream.
	 * @param out The output stream.
	 */
	void print(std::ostream & out) const;

	/**
	 * Print vector to output stream.
	 * @param out The output stream.
	 * @param header Header.
	 */
	void print(std::ostream & out, std::string header) const;

	/**
	 * Destructor.
	 */
	virtual ~Vec();
};

} /* namespace LinkPred */

#endif /* INCLUDE_VEC_HPP_ */
