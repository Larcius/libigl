// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_DISTORTION_ANGLE_H
#define IGL_DISTORTION_ANGLE_H
#include "igl_inline.h"
#include <Eigen/Geometry>
namespace igl
{
	template <
		typename DerivedV,
		typename DerivedF>
		IGL_INLINE double compute_angle_distortion(
			const Eigen::PlainObjectBase<DerivedV> & V_uv,
			const Eigen::PlainObjectBase<DerivedV> & V,
			const Eigen::PlainObjectBase<DerivedF> & F);

	template <
		typename DerivedV,
		typename DerivedF>
		IGL_INLINE double compute_angle_distortion_triangle_squared(
			const Eigen::PlainObjectBase<DerivedV> & V_uv,
			const Eigen::PlainObjectBase<DerivedV> & V,
			const Eigen::PlainObjectBase<DerivedF> & F,
			unsigned int face_index);

	template <
		typename DerivedV,
		typename DerivedF>
		IGL_INLINE double compute_stretch_distortion(
			const Eigen::PlainObjectBase<DerivedV> & V_uv,
			const Eigen::PlainObjectBase<DerivedV> & V,
			const Eigen::PlainObjectBase<DerivedF> & F);

	template <
		typename DerivedV,
		typename DerivedF>
		IGL_INLINE double compute_stretch_distortion_triangle_squared(
			const Eigen::PlainObjectBase<DerivedV> & V_uv,
			const Eigen::PlainObjectBase<DerivedV> & V,
			const Eigen::PlainObjectBase<DerivedF> & F,
			unsigned int face_index);

	template <
		typename DerivedV,
		typename DerivedF>
		IGL_INLINE double compute_shear_distortion(
			const Eigen::PlainObjectBase<DerivedV> & V_uv,
			const Eigen::PlainObjectBase<DerivedV> & V,
			const Eigen::PlainObjectBase<DerivedF> & F);

	template <
		typename DerivedV,
		typename DerivedF>
		IGL_INLINE double compute_shear_distortion_triangle_squared(
			const Eigen::PlainObjectBase<DerivedV> & V_uv,
			const Eigen::PlainObjectBase<DerivedV> & V,
			const Eigen::PlainObjectBase<DerivedF> & F,
			unsigned int face_index);

	template <
		typename DerivedV,
		typename DerivedF>
		IGL_INLINE double compute_area_distortion(
			const Eigen::PlainObjectBase<DerivedV> & V_uv,
			const Eigen::PlainObjectBase<DerivedV> & V,
			const Eigen::PlainObjectBase<DerivedF> & F);

	template <
		typename DerivedV,
		typename DerivedF>
		IGL_INLINE double compute_area(
			const Eigen::PlainObjectBase<DerivedV> & V,
			const Eigen::PlainObjectBase<DerivedF> & F);

	template <
		typename DerivedV,
		typename DerivedF>
		IGL_INLINE double compute_area_face(
			const Eigen::PlainObjectBase<DerivedV> & V,
			const Eigen::PlainObjectBase<DerivedF> & F,
			unsigned int face_index);
}

#ifndef IGL_STATIC_LIBRARY
#include "distortion.cpp"
#endif

#endif
