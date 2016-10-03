// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "writeSVG.h"
#include <cstdio>
#include <fstream>
#include <Eigen/Sparse>

// write mesh to an ascii svg file
template <typename DerivedV, typename DerivedF>
IGL_INLINE bool igl::writeSVG(
	const std::string fname,
	const Eigen::PlainObjectBase<DerivedV>& V,
	const Eigen::PlainObjectBase<DerivedF>& F)
{
	using namespace std;
	using namespace Eigen;
	assert(V.cols() >= 2 && "V should have at least 2 columns");
	ofstream s(fname);
	if (!s.is_open())
	{
		fprintf(stderr, "IOError: writeSVG() could not open %s\n", fname.c_str());
		return false;
	}

	int strokeWidthBorder = 5;
	int strokeWidthInterior = 4;
	int strokeWidth = max(strokeWidthBorder, strokeWidthInterior);

	double minX = V.col(0).minCoeff();
	double minY = -V.col(1).maxCoeff();

	double maxX = V.col(0).maxCoeff();
	double maxY = -V.col(1).minCoeff();

	s << "<svg xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" xmlns:cc=\"http://creativecommons.org/ns#\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\" viewBox=\"" << (minX - strokeWidth) << " " << (minY - strokeWidth) << " " << (maxX - minX + 2 * strokeWidth) << " " << (maxY - minY + 2 * strokeWidth) << "\">" << endl;

	unsigned int numFaces = F.rows();
	unsigned int verticesPerFace = F.cols();
	unsigned int numVertices = V.rows();
	
	typedef Eigen::Triplet<short> T;
	std::vector<T> tripletList;
	unsigned int estimation_of_entries = numFaces + numVertices; // - 2
	tripletList.reserve(estimation_of_entries);

	for (unsigned int i = 0; i < numFaces; i++) {
		unsigned int lastVertexIndex;

		for (unsigned int j = 0; j < verticesPerFace + 1; j++) {
			unsigned int vertexIndex = F(i, j % verticesPerFace);
			
			if (j > 0) {
				tripletList.push_back(T(min(vertexIndex, lastVertexIndex), max(vertexIndex, lastVertexIndex), (short)1));
			}

			lastVertexIndex = vertexIndex;
		}
	}

	SparseMatrix<short> interiorEdges(numVertices, numVertices);
	interiorEdges.setFromTriplets(tripletList.begin(), tripletList.end());

	s << "\t<g stroke-linejoin=\"round\" stroke-linecap=\"round\" fill=\"none\">" << endl;

	for (int type = 2; type > 0; type--) {
		s << "\t\t<g stroke-width=\"" << (type == 1 ? strokeWidthBorder : strokeWidthInterior) << "\" stroke=\"#" << (type == 1 ? "0f0" : "000") << "\">" << endl;

		for (int k = 0; k < interiorEdges.outerSize(); ++k) {
			for (Eigen::SparseMatrix<short>::InnerIterator it(interiorEdges, k); it; ++it) {
				short count = it.value();

				assert(count > 0 && "this iterator should only list non-zero values");

				if (it.value() != type) {
					continue;
				}

				VectorXd vertex1 = V.row(it.row());
				VectorXd vertex2 = V.row(it.col());

				s << "\t\t\t<path d=\"M " << vertex1[0] << " " << (-vertex1[1]) << " L " << vertex2[0] << " " << (-vertex2[1]) << "\" />" << endl;
			}
		}

		s << "\t\t</g>" <<endl;
	}

	s << "\t</g>" << endl << "</svg>";

	return true;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
// generated by autoexplicit.sh
template bool igl::writeSVG<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> >, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
/*
template bool igl::writeSVG<Eigen::Matrix<double, -1, 2, 1, -1, 2>, Eigen::Matrix<unsigned int, -1, -1, 1, -1, -1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> >, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 2, 1, -1, 2> > const&, Eigen::PlainObjectBase<Eigen::Matrix<unsigned int, -1, -1, 1, -1, -1> > const&);
template bool igl::writeSVG<Eigen::Matrix<float, -1, 2, 1, -1, 2>, Eigen::Matrix<unsigned int, -1, -1, 1, -1, -1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> >, Eigen::PlainObjectBase<Eigen::Matrix<float, -1, 2, 1, -1, 2> > const&, Eigen::PlainObjectBase<Eigen::Matrix<unsigned int, -1, -1, 1, -1, -1> > const&);
template bool igl::writeSVG<Eigen::Matrix<double, -1, 2, 0, -1, 2>, Eigen::Matrix<int, -1, 2, 0, -1, 2> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> >, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 2, 0, -1, 2> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 2, 0, -1, 2> > const&);
template bool igl::writeSVG<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 2, 0, -1, 2> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> >, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 2, 0, -1, 2> > const&);
*/
#endif