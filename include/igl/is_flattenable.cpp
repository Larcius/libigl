// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Gabriel Cordes <gabriel.cordes@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "is_flattenable.h"

#include "is_edge_manifold.h"
#include "is_vertex_manifold.h"
#include "boundary_loop.h"
#include "edges.h"
#include "facet_components.h"

template <typename DerivedV, typename DerivedF>
IGL_INLINE short igl::is_flattenable(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedF>& F)
{
	typedef typename DerivedF::Scalar Index;
	typedef typename DerivedF::Index FIndex;
	Eigen::Matrix<Index, -1, 1, 0, -1, 1> C;
	Eigen::Matrix<Index, -1, 2, 0, -1, 2> E;
	Eigen::Matrix<Index, -1, 1, 0, -1, 1> B;

	short ret = 0;

	if (!is_edge_manifold(V, F)) {
		ret |= flattenable::NO_EDGE_MANIFOLD;
	}

	if (!is_vertex_manifold(F, B)) {
		ret |= flattenable::NO_VERTEX_MANIFOLD;
	}

	facet_components(F, C);
	if (C.minCoeff() != C.maxCoeff()) {
		ret |= flattenable::NOT_CONNECTED;
	}
	
	std::vector<std::vector<Index> > L;
	boundary_loop(F, L);
	if (L.size() == 0) {
		ret |= flattenable::EMPTY_BOUNDARY;
	}
	
	edges(F, E);
	if (V.rows() - E.rows() + F.rows() != 2 - L.size()) {
		ret |= flattenable::NO_EULER_CHARACTERISTIC;
	}

	return ret;
}

#ifdef IGL_STATIC_LIBRARY
template short igl::is_flattenable<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
#endif
