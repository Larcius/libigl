// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Gabriel Cordes <gabriel.cordes@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_IS_FLATTENABLE_H
#define IGL_IS_FLATTENABLE_H
#include "igl_inline.h"

#include <Eigen/Core>
#include <Eigen/Dense>

namespace igl
{
  // check if the mesh is a flattenable surface
  //
  // Definition (surface):
  // A connected 2-manifold is called a surface.
  //
  // Definition (flat surface):
  // A surface that is isometric to a surface embedded in R^2 is called a flat surface
  //
  // Definition (flattenable surface):
  // A surface that is homeomprhic to a flat surface is called a flattenable surface.
  //
  //
  // Inputs:
  //   V  #V by dim list of mesh vertex positions
  //   F  #F by 3 list of triangle indices
  // Output:
  //   Flags that indicate which properties are violated in order to be a flattenable surface
  //   or 0 if it is a flattenable surface. See enum igl::flattenable.
  //
  // See also: is_edge_manifold
  //           is_vertex_manifold
  //           boundary_loop
  template <typename DerivedV, typename DerivedF>
  IGL_INLINE short is_flattenable(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedF>& F);

  enum flattenable : short
  {
	NO_EDGE_MANIFOLD = 1,
	NO_VERTEX_MANIFOLD = 2,
	EMPTY_BOUNDARY = 4,
	NO_EULER_CHARACTERISTIC = 8,
	NOT_CONNECTED = 16
  };
}

#ifndef IGL_STATIC_LIBRARY
#  include "is_flattenable.cpp"
#endif

#endif
