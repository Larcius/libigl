// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Gabriel Cordes <gabriel.cordes@googlemail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_ABF_H
#define IGL_ABF_H
#include "igl/igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace igl
{
	namespace copyleft
	{
		// Compute necessary information to start using an ARAP deformation
		//
		// Inputs:
		//   V  #V by dim list of mesh positions
		//   F  #F by simplex-size list of triangle|tet indices into V
		//   dim  dimension being used at solve time. For deformation usually dim =
		//     V.cols(), for surface parameterization V.cols() = 3 and dim = 2
		//   b  #b list of "boundary" fixed vertex indices into V
		// Outputs:
		//   data  struct containing necessary precomputation
		/*
		template <
			typename DerivedV,
			typename DerivedF>
			IGL_INLINE ParamHandle abf_precomputation(
				const Eigen::PlainObjectBase<DerivedV> & V,
				const Eigen::PlainObjectBase<DerivedF> & F);
		*/
		// Inputs:
		//   bc  #b by dim list of boundary conditions
		//   data  struct containing necessary precomputation and parameters
		//   U  #V by dim initial guess
		/*
		template <
			typename Derivedbc,
			typename DerivedU>
			IGL_INLINE bool abf_solve(
				const Eigen::PlainObjectBase<Derivedbc> & bc,
				ParamHandle & data,
				Eigen::PlainObjectBase<DerivedU> & U);
		*/

		template <
			typename DerivedV,
			typename DerivedF >
			IGL_INLINE bool abf_solve(
				const Eigen::PlainObjectBase<DerivedV> & V,
				const Eigen::PlainObjectBase<DerivedF> & F,
				Eigen::PlainObjectBase<DerivedV> & V_uv,
				bool preFillHoles);
	};
};

#ifndef IGL_STATIC_LIBRARY
#include "abf.cpp"
#endif

#endif
