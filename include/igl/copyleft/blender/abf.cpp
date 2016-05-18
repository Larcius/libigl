// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Gabriel Cordes <gabriel.cordes@googlemail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "abf.h"
#include "source\blender\editors\uvedit\uvedit_parametrizer.h"
#include "BLI_alloca.h"
#include <cassert>
#include <iostream>

template <
	typename DerivedV,
	typename DerivedF>
	IGL_INLINE bool igl::copyleft::abf_solve(
		const Eigen::PlainObjectBase<DerivedV> & V,
		const Eigen::PlainObjectBase<DerivedF> & F,
		Eigen::PlainObjectBase<DerivedV> & V_uv,
		bool preFillHoles)
{
	using namespace std;
	using namespace Eigen;

	ParamHandle *handle;

	handle = param_construct_begin();

	Matrix<double, Dynamic, Dynamic, RowMajor> uv_all;

	add_faces<DerivedV, DerivedF>(handle, V, F, uv_all);

	param_construct_end(handle, preFillHoles ? PARAM_TRUE : PARAM_FALSE, PARAM_FALSE);

	if(true) {
		param_abf_begin(handle);
		param_abf_solve(handle);
		param_abf_end(handle);
	} else {
		param_lscm_begin(handle, PARAM_FALSE, PARAM_TRUE);
		param_lscm_solve(handle);
		param_lscm_end(handle);
	}

	//param_average(handle);
	//param_pack(handle, 0, false);

	param_flush(handle);

	param_delete(handle);

	const unsigned int rows = uv_all.rows();
	const unsigned int cols = uv_all.cols();
	V_uv.resize(rows, cols);
	
	for (unsigned int i = 0; i < rows; i++) {
		for (unsigned int j = 0; j < cols; j++) {
			V_uv(i, j) = uv_all(i, j);
		}
	}

	return true;
}

template <
	typename DerivedV,
	typename DerivedF>
	IGL_INLINE void add_faces(
		ParamHandle* handle,
		const Eigen::PlainObjectBase<DerivedV> & V,
		const Eigen::PlainObjectBase<DerivedF> & F,
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & V_uv)
{
	using namespace std;
	using namespace Eigen;
	
	ParamBool tri_pin[3] = { PARAM_FALSE, PARAM_FALSE, PARAM_FALSE };
	ParamBool tri_select[3] = { PARAM_FALSE, PARAM_FALSE, PARAM_FALSE };

	const size_t numVertices = V.rows();
	const size_t uvCols = 2;
	V_uv.resize(numVertices, uvCols);
	double* uv_data = V_uv.data();
	
	for (size_t i = 0, numVerticesPerFace = F.cols(), numFaces = F.rows(); i < numFaces; i++)
	{
		/*
		ParamKey key = i;

		/
		vector<ParamKey> vkeys;

		vector<double*> cos;// [numVerticesPerFace][4];
		vector<double*> uvs;// [numVerticesPerFace][4];
		/

		ParamKey *vkeys = (ParamKey *)BLI_array_alloca(vkeys, numVerticesPerFace);
		ParamBool *pin = (ParamBool*)BLI_array_alloca(pin, numVerticesPerFace);
		ParamBool *select = (ParamBool*)BLI_array_alloca(select, numVerticesPerFace);
		double **co = (double **)BLI_array_alloca(co, numVerticesPerFace);
		double **uv = (double **)BLI_array_alloca(uv, numVerticesPerFace);

		for (size_t j = 0; j < numVerticesPerFace; j++)
		{
			unsigned int vertex = F(i, j);

			vkeys[j] = (ParamKey)vertex;

			double temp_co[4] = {
				(double)V(vertex, 0),
				(double)V(vertex, 1),
				(double)V(vertex, 2),
				(double)0
			};

			co[j] = temp_co;

			double temp_uv[4] = {
				(double)V(vertex, 0),
				(double)V(vertex, 1),
				(double)V(vertex, 2),
				(double)0
			};

			uv[j] = temp_uv;

			pin[j] = PARAM_FALSE;
			select[j] = PARAM_FALSE;
		}

		//param_face_add(handle, key, numVerticesPerFace, &vkeys[0], &cos[0], &uvs[0], &pin, &select, NULL);
		param_face_add(handle, key, numVerticesPerFace, vkeys, co, uv, pin, select, NULL);
		*/

		ParamKey tri_vkeys[3] = { F(i, 0), F(i, 1), F(i, 2) };

		double tri_co0[3] = { V(tri_vkeys[0], 0), V(tri_vkeys[0], 1), V(tri_vkeys[0], 2) };
		double tri_co1[3] = { V(tri_vkeys[1], 0), V(tri_vkeys[1], 1), V(tri_vkeys[1], 2) };
		double tri_co2[3] = { V(tri_vkeys[2], 0), V(tri_vkeys[2], 1), V(tri_vkeys[2], 2) };
		double *tri_co[3] = { tri_co0, tri_co1, tri_co2 };

		double *tri_uv[3] = { &uv_data[tri_vkeys[0] * uvCols], &uv_data[tri_vkeys[1] * uvCols], &uv_data[tri_vkeys[2] * uvCols] };

		param_face_add(handle, i, 3, tri_vkeys, tri_co, tri_uv, tri_pin, tri_select, NULL);
	}
}

#ifdef IGL_STATIC_LIBRARY
//template bool igl::copyleft::abf_solve<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&);
template bool igl::copyleft::abf_solve<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, bool);
#endif
