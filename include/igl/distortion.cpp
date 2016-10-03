// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Gabriel Cordes <gabriel.cordes@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "distortion.h"

#include <cassert>
#include <vector>

template <
	typename DerivedV,
	typename DerivedF>
	IGL_INLINE double igl::compute_angle_distortion(
		const Eigen::PlainObjectBase<DerivedV> & V_uv,
		const Eigen::PlainObjectBase<DerivedV> & V,
		const Eigen::PlainObjectBase<DerivedF> & F)
{
	assert(V_uv.cols() == 2 && "V_vu should have 2 cols");
	assert(V.cols() == 3 && "V should have 3 cols");
	assert(F.cols() == 3 && "F should have 3 cols (i.e. every face should be a triangle)");

	using namespace std;
	using namespace Eigen;

	double total_area = 0;
	//double total_area_uv = 0;
	double distortion = 0;

	const unsigned int rowsF = F.rows();

	for (unsigned int i = 0; i < rowsF; i++) {
		double area = igl::compute_area_face(V, F, i);
		double area_uv = igl::compute_area_face(V_uv, F, i);

		total_area += area;
		//total_area_uv += area_uv;

		distortion += igl::compute_angle_distortion_triangle_squared(V_uv, V, F, i) * area;// *area_uv;
	}

	distortion /= total_area;// *total_area_uv;

	return sqrt(distortion);
}

template <
	typename DerivedV,
	typename DerivedF>
	IGL_INLINE double igl::compute_angle_distortion_triangle_squared(
		const Eigen::PlainObjectBase<DerivedV> & V_uv,
		const Eigen::PlainObjectBase<DerivedV> & V,
		const Eigen::PlainObjectBase<DerivedF> & F,
		unsigned int face_index)
{
	assert(V_uv.cols() == 2 && "V_vu should have 2 cols");
	assert(V.cols() == 3 && "V should have 3 cols");
	assert(F.cols() == 3 && "F should have 3 cols (i.e. every face should be a triangle)");
	assert((face_index >= 0 && face_index < F.rows()) && "face index out of bound");

	using namespace std;

	double distortion = 0;
	
	const unsigned int colsF = F.cols();
	const unsigned int colsV = V.cols();
	const unsigned int colsV_uv = V_uv.cols();

	vector<vector<double>> vec_v_uv(colsF, vector<double>(colsV_uv));
	vector<vector<double>> vec_v(colsF, vector<double>(colsV));
	vector<double> alphas(colsF);
	vector<double> betas(colsF);
	
	for (unsigned int v_index = 0; v_index < colsF; v_index++) {
		for (unsigned int v_coord = 0; v_coord < colsV; v_coord++) {
			vec_v[v_index][v_coord] = V(F(face_index, v_index), v_coord);
		}

		for (unsigned int v_coord = 0; v_coord < colsV_uv; v_coord++) {
			vec_v_uv[v_index][v_coord] = V_uv(F(face_index, v_index), v_coord);
		}
	}

	for (unsigned int v_index = 0; v_index < colsF; v_index++) {
		alphas[v_index] = p_vec_angle(vec_v_uv[v_index], vec_v_uv[(v_index + 1) % colsF], vec_v_uv[(v_index + 2) % colsF]);
		betas[v_index] = p_vec_angle(vec_v[v_index], vec_v[(v_index + 1) % colsF], vec_v[(v_index + 2) % colsF]);
	}

	for (unsigned int l = 0; l < colsF; l++) {
		double diff = alphas[l] - betas[l];

		distortion += diff * diff / (betas[l] * betas[l]);
	}

	return distortion;
}

template <
	typename DerivedV,
	typename DerivedF>
	IGL_INLINE double igl::compute_stretch_distortion(
		const Eigen::PlainObjectBase<DerivedV> & V_uv,
		const Eigen::PlainObjectBase<DerivedV> & V,
		const Eigen::PlainObjectBase<DerivedF> & F)
{
	assert(V_uv.cols() == 2 && "V_vu should have 2 cols");
	assert(V.cols() == 3 && "V should have 3 cols");
	assert(F.cols() == 3 && "F should have 3 cols (i.e. every face should be a triangle)");

	using namespace std;
	using namespace Eigen;

	double total_area = 0;
	double total_area_uv = 0;
	double distortion = 0;

	const unsigned int rowsF = F.rows();

	for (unsigned int i = 0; i < rowsF; i++) {
		double area = igl::compute_area_face(V, F, i);
		double area_uv = igl::compute_area_face(V_uv, F, i);

		total_area += area;
		total_area_uv += area_uv;

		distortion += igl::compute_stretch_distortion_triangle_squared(V_uv, V, F, i) * area;
	}
	
	distortion *= total_area_uv / (total_area * total_area);

	return sqrt(distortion);
}

template <
	typename DerivedV,
	typename DerivedF>
	IGL_INLINE double igl::compute_stretch_distortion_triangle_squared(
		const Eigen::PlainObjectBase<DerivedV> & V_uv,
		const Eigen::PlainObjectBase<DerivedV> & V,
		const Eigen::PlainObjectBase<DerivedF> & F,
		unsigned int face_index)
{
	assert(V_uv.cols() == 2 && "V_vu should have 2 cols");
	assert(V.cols() == 3 && "V should have 3 cols");
	assert(F.cols() == 3 && "F should have 3 cols (i.e. every face should be a triangle)");
	assert((face_index >= 0 && face_index < F.rows()) && "face index out of bound");

	double partial_s[3];
	double partial_t[3];

	compute_partial(V_uv, V, F, face_index, partial_s, partial_t);
		
	return (dot_v3v3(partial_s, partial_s) + dot_v3v3(partial_t, partial_t)) * 0.5;
}

template <
	typename DerivedV,
	typename DerivedF>
	IGL_INLINE double igl::compute_shear_distortion(
		const Eigen::PlainObjectBase<DerivedV> & V_uv,
		const Eigen::PlainObjectBase<DerivedV> & V,
		const Eigen::PlainObjectBase<DerivedF> & F)
{
	assert(V_uv.cols() == 2 && "V_vu should have 2 cols");
	assert(V.cols() == 3 && "V should have 3 cols");
	assert(F.cols() == 3 && "F should have 3 cols (i.e. every face should be a triangle)");

	using namespace std;
	using namespace Eigen;

	double total_area = 0;
	double distortion = 0;

	const unsigned int rowsF = F.rows();

	for (unsigned int i = 0; i < rowsF; i++) {
		double area = igl::compute_area_face(V, F, i);

		total_area += area;

		distortion += igl::compute_shear_distortion_triangle_squared(V_uv, V, F, i) * area;
	}

	distortion /= total_area;

	return sqrt(distortion);
}

template <
	typename DerivedV,
	typename DerivedF>
	IGL_INLINE double igl::compute_shear_distortion_triangle_squared(
		const Eigen::PlainObjectBase<DerivedV> & V_uv,
		const Eigen::PlainObjectBase<DerivedV> & V,
		const Eigen::PlainObjectBase<DerivedF> & F,
		unsigned int face_index)
{
	assert(V_uv.cols() == 2 && "V_vu should have 2 cols");
	assert(V.cols() == 3 && "V should have 3 cols");
	assert(F.cols() == 3 && "F should have 3 cols (i.e. every face should be a triangle)");
	assert((face_index >= 0 && face_index < F.rows()) && "face index out of bound");

	double partial_s[3];
	double partial_t[3];

	compute_partial(V_uv, V, F, face_index, partial_s, partial_t);

	double shear = dot_v3v3(partial_s, partial_t) / (norm_v3(partial_s) * norm_v3(partial_t));
	
	return shear * shear;
}

template <
	typename DerivedV,
	typename DerivedF>
	IGL_INLINE double igl::compute_area_distortion(
		const Eigen::PlainObjectBase<DerivedV> & V_uv,
		const Eigen::PlainObjectBase<DerivedV> & V,
		const Eigen::PlainObjectBase<DerivedF> & F)
{
	assert(V_uv.cols() == 2 && "V_vu should have 2 cols");
	assert(V.cols() == 3 && "V should have 3 cols");
	assert(F.cols() == 3 && "F should have 3 cols (i.e. every face should be a triangle)");

	using namespace std;
	using namespace Eigen;

	double total_area = 0;
	double total_area_uv = 0;
	double distortion = 0;

	const unsigned int rowsF = F.rows();

	for (unsigned int i = 0; i < rowsF; i++) {
		double area_uv = igl::compute_area_face(V_uv, F, i);
		double area = igl::compute_area_face(V, F, i);

		double area_diff = area_uv - area;

		total_area += area;
		total_area_uv += area_uv;
	}
	
	for (unsigned int i = 0; i < rowsF; i++) {
		// duplicate calculation
		double area_uv = igl::compute_area_face(V_uv, F, i);
		double area = igl::compute_area_face(V, F, i);

		double area_diff = area_uv / total_area_uv - area / total_area;
		
		distortion += area_diff * area_diff * area;
	}

	distortion /= total_area;

	return sqrt(distortion);
}


template <
	typename DerivedV,
	typename DerivedF>
	IGL_INLINE double igl::compute_area_face(
		const Eigen::PlainObjectBase<DerivedV> & V,
		const Eigen::PlainObjectBase<DerivedF> & F,
		unsigned int face_index)
{
	assert((V.cols() >= 2 && V.cols() <= 3) && "V should have two or three cols");
	assert(F.cols() == 3 && "F should have 3 cols (i.e. every face should be a triangle)");
	assert((face_index >= 0 && face_index < F.rows()) && "face index out of bound");

	using namespace std;

	const unsigned int colsF = F.cols();
	const unsigned int colsV = V.cols();

	vector<vector<double>> vec_v(colsF, vector<double>(colsV));

	for (unsigned int v_index = 0; v_index < colsF; v_index++) {
		for (unsigned int v_coord = 0; v_coord < colsV; v_coord++) {
			vec_v[v_index][v_coord] = V(F(face_index, v_index), v_coord);
		}
	}

	double area = 0;

	for (unsigned int v_index = 1; v_index < colsF - 1; v_index++) {
		area += p_area(vec_v[0], vec_v[v_index], vec_v[v_index + 1]);
	}

	return area;
}

template <
	typename DerivedV,
	typename DerivedF>
	IGL_INLINE double igl::compute_area(
		const Eigen::PlainObjectBase<DerivedV> & V,
		const Eigen::PlainObjectBase<DerivedF> & F)
{
	const unsigned int rowsF = F.rows();

	assert((V.cols() >= 2 && V.cols() <= 3) && "V should have two or three cols");
	assert(F.cols() == 3 && "F should have 3 cols (i.e. every face should be a triangle)");

	using namespace std;

	double total_area = 0;

	for (unsigned int i = 0; i < rowsF; i++) {
		double area = igl::compute_area_face(V, F, i);

		total_area += area;
	}

	return total_area;
}

/* Geometry */

inline void copy_v3_v3(double r[3], const double a[3])
{
	r[0] = a[0];
	r[1] = a[1];
	r[2] = a[2];
}

inline double dot_v2v2(const double a[2], const double b[2])
{
	return a[0] * b[0] + a[1] * b[1];
}

inline double dot_v3v3(const double a[3], const double b[3])
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline void mul_v2_fl(double a[2], const double f)
{
	a[0] *= f;
	a[1] *= f;
}

inline void mul_v3_fl(double a[3], const double f)
{
	a[0] *= f;
	a[1] *= f;
	a[2] *= f;
}

inline void mul_v3_v3fl(const double a[3], const double f, double r[3])
{
	r[0] = a[0] * f;
	r[1] = a[1] * f;
	r[2] = a[2] * f;
}

inline void add_v3_v3(double a[3], const double b[3])
{
	a[0] += b[0];
	a[1] += b[1];
	a[2] += b[2];
}

inline double norm_v3(const double a[3]) {
	return sqrt(dot_v3v3(a, a));
}

inline double norm_v2(const double a[2]) {
	return sqrt(dot_v2v2(a, a));
}

inline void normalize_v2(double a[2]) {
	double norm = norm_v2(a);

	mul_v2_fl(a, 1.0 / norm);
}

inline void normalize_v3(double a[3]) {
	double norm = norm_v3(a);

	mul_v3_fl(a, 1.0 / norm);
}

inline void cross_product3(const double a[3], const double b[3], double r[3])
{
	r[0] = a[1] * b[2] - a[2] * b[1];
	r[1] = a[2] * b[0] - a[0] * b[2];
	r[2] = a[0] * b[1] - a[1] * b[0];
}

inline double cross_product2(const double a[2], const double b[2])
{
	return a[0] * b[1] - a[1] * b[0];
}

inline double p_vec_angle_cos(std::vector<double> v1, std::vector<double> v2, std::vector<double> v3)
{
	double d1[3], d2[3];

	d1[0] = v1[0] - v2[0];
	d1[1] = v1[1] - v2[1];
	
	d2[0] = v3[0] - v2[0];
	d2[1] = v3[1] - v2[1];

	if (v1.size() > 2) {
		d1[2] = v1[2] - v2[2];
		d2[2] = v3[2] - v2[2];
	}
	else {
		d1[2] = 0;
		d2[2] = 0;
	}

	normalize_v3(d1);
	normalize_v3(d2);

	return dot_v3v3(d1, d2);
}

inline double p_vec_angle(std::vector<double> v1, std::vector<double> v2, std::vector<double> v3)
{

	double dot = p_vec_angle_cos(v1, v2, v3);

	if (dot <= -1.0f)
		return M_PI;
	else if (dot >= 1.0f)
		return 0.0f;
	else
		return acos(dot);
}

inline double p_area(std::vector<double> v1, std::vector<double> v2, std::vector<double> v3)
{
	double d1[3], d2[3];

	d1[0] = v1[0] - v2[0];
	d1[1] = v1[1] - v2[1];

	d2[0] = v3[0] - v2[0];
	d2[1] = v3[1] - v2[1];

	if (v1.size() > 2) {
		d1[2] = v1[2] - v2[2];
		d2[2] = v3[2] - v2[2];
	}
	else {
		d1[2] = 0;
		d2[2] = 0;
	}

	double cross[3];
	cross_product3(d1, d2, cross);

	return 0.5 * norm_v3(cross);
}

template <
	typename DerivedV,
	typename DerivedF>
	IGL_INLINE void compute_partial(
		const Eigen::PlainObjectBase<DerivedV> & V_uv,
		const Eigen::PlainObjectBase<DerivedV> & V,
		const Eigen::PlainObjectBase<DerivedF> & F,
		unsigned int face_index,
		double partial_s[3],
		double partial_t[3]) {

	const unsigned int colsV_uv = V_uv.cols();
	const unsigned int colsV = V.cols();
	const unsigned int colsF = F.cols();

	assert(colsV_uv == 2 && "V_vu should have 2 cols");
	assert(colsV == 3 && "V should have 3 cols");
	assert(colsF == 3 && "F should have 3 cols (i.e. every face should be a triangle)");
	assert((face_index >= 0 && face_index < F.rows()) && "face index out of bound");

	using namespace std;

	vector<vector<double>> vec_v_uv(colsF, vector<double>(colsV_uv));
	vector<vector<double>> vec_v(colsF, vector<double>(colsV));

	for (unsigned int v_index = 0; v_index < colsF; v_index++) {
		for (unsigned int v_coord = 0; v_coord < colsV; v_coord++) {
			vec_v[v_index][v_coord] = V(F(face_index, v_index), v_coord);
		}

		for (unsigned int v_coord = 0; v_coord < colsV_uv; v_coord++) {
			vec_v_uv[v_index][v_coord] = V_uv(F(face_index, v_index), v_coord);
		}
	}

	double temp[3] = {0,0,0};

	// init to 0
	copy_v3_v3(partial_s, temp);
	copy_v3_v3(partial_t, temp);

	copy_v3_v3(temp, &vec_v[0][0]);
	mul_v3_fl(temp, vec_v_uv[1][1] - vec_v_uv[2][1]);
	add_v3_v3(partial_s, temp);

	copy_v3_v3(temp, &vec_v[1][0]);
	mul_v3_fl(temp, vec_v_uv[2][1] - vec_v_uv[0][1]);
	add_v3_v3(partial_s, temp);

	copy_v3_v3(temp, &vec_v[2][0]);
	mul_v3_fl(temp, vec_v_uv[0][1] - vec_v_uv[1][1]);
	add_v3_v3(partial_s, temp);


	copy_v3_v3(temp, &vec_v[0][0]);
	mul_v3_fl(temp, vec_v_uv[2][0] - vec_v_uv[1][0]);
	add_v3_v3(partial_t, temp);

	copy_v3_v3(temp, &vec_v[1][0]);
	mul_v3_fl(temp, vec_v_uv[0][0] - vec_v_uv[2][0]);
	add_v3_v3(partial_t, temp);

	copy_v3_v3(temp, &vec_v[2][0]);
	mul_v3_fl(temp, vec_v_uv[1][0] - vec_v_uv[0][0]);
	add_v3_v3(partial_t, temp);

	double area_parallelogram = (vec_v_uv[1][0] - vec_v_uv[0][0]) * (vec_v_uv[2][1] - vec_v_uv[0][1]) - (vec_v_uv[2][0] - vec_v_uv[0][0]) * (vec_v_uv[1][1] - vec_v_uv[0][1]);

	mul_v3_fl(partial_s, 1.0 / area_parallelogram);
	mul_v3_fl(partial_t, 1.0 / area_parallelogram);
}


#ifdef IGL_STATIC_LIBRARY
template double igl::compute_area_distortion<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
template double igl::compute_area<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
template double igl::compute_area_face<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, unsigned int);
template double igl::compute_angle_distortion<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
template double igl::compute_angle_distortion_triangle_squared<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, unsigned int);
template double igl::compute_stretch_distortion<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
template double igl::compute_stretch_distortion_triangle_squared<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, unsigned int);
template double igl::compute_shear_distortion<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
template double igl::compute_shear_distortion_triangle_squared<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, unsigned int);
#endif