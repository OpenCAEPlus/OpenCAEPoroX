/*! \file    GmshGrid.cpp
 *  \brief   GmshGrid class declaration
 *  \author  Shizhe Li
 *  \date    Sep/10/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */



#ifdef WITH_GMSH


#include "GmshGrid.hpp"

Polygon::Polygon(const vector<OCP_USI>& pIndex, const OCP_USI& tag_in, const string& phyinfo) 
{
	p        = pIndex;
	tag      = tag_in;
	physical = phyinfo;
}


void Polygon::CalCenter(const vector<OCP_DBL>& points) 
{
	const USI np = p.size();
	center.Reset();
	for (USI i = 0; i < np; i++) {
		center.x += points[3 * p[i] + 0];
		center.y += points[3 * p[i] + 1];
		center.z += points[3 * p[i] + 2];
	}
	center *= 1.0 / np;
}


void Polygon::CalArea(const vector<OCP_DBL>& points) 
{
	if (p.size() == 3) {
		const OCP_DBL* p0 = &points[3 * p[0]];
		const OCP_DBL* p1 = &points[3 * p[1]];
		const OCP_DBL* p2 = &points[3 * p[2]];
		const USI      x = 0;
		const USI      y = 1;
		area = 0.5 * fabs((p2[x] - p1[x]) * (p0[y] - p1[y]) - (p2[y] - p1[y]) * (p0[x] - p1[x]));
	}
	else {
		const OCP_DBL* p0 = &points[3 * p[0]];
		const OCP_DBL* p1 = &points[3 * p[1]];
		const OCP_DBL* p2 = &points[3 * p[2]];
		const OCP_DBL* p3 = &points[3 * p[3]];
		const USI      x = 0;
		const USI      y = 1;
		area = 0.5 * (fabs((p2[x] - p1[x]) * (p0[y] - p1[y]) - (p2[y] - p1[y]) * (p0[x] - p1[x]))
			+ fabs((p0[x] - p3[x]) * (p2[y] - p3[y]) - (p0[y] - p3[y]) * (p2[x] - p3[x])));
	}
}


OCP_BOOL Polygon::IfPointInElement(const Point3D& objP, const vector<OCP_DBL>& points)
{
	OCP_DBL   tmpArea = 0;
	const USI numP    = p.size();
	for (USI i = 0; i < numP; i++) {
		const Point3D tmpPoint = CrossProduct(objP - Point3D(&points[3 * p[i % numP]]), 
			                                  objP - Point3D(&points[3 * p[(i + 1) % numP]]));
		tmpArea += sqrt(tmpPoint * tmpPoint);
	}
	tmpArea *= 0.5;

	if (fabs(area - tmpArea) < TINY) {
		return OCP_TRUE;
	}
	else {
		return OCP_FALSE;
	}
}


void GMSHGrid::InputGrid(const string& file)
{

	gmsh::initialize();
	gmsh::open(file);

	// Print the model name and dimension:
	std::string name;
	gmsh::model::getCurrent(name);
	std::cout << "Model " << name << " (" << gmsh::model::getDimension() << "D)\n";
	dimen = gmsh::model::getDimension();

	if (dimen == 2) {
		InputGrid2D(file);
	}
	gmsh::finalize();

	Setup();
}

void GMSHGrid::InputGrid2D(const string& file)
{

	// Get all mesh nodes
	std::vector<std::size_t> nodeTags;
	std::vector<double> nodeParams;
	gmsh::model::mesh::getNodes(nodeTags, points, nodeParams, -1, -1);

	// Get all the elementary entities in the model, as a vector of (dimension, tag) pairs:
	std::vector<std::pair<int, int>> entities;
	gmsh::model::getEntities(entities);

	OCP_USI lineTag   = 1;
	OCP_USI faceIndex = 0;

	for (const auto& e : entities) {
		// Dimension and tag of the entity:
		const int dim = e.first, tag = e.second;


		std::vector<int> physicalTags;
		gmsh::model::getPhysicalGroupsForEntity(dim, tag, physicalTags);

		if (physicalTags.empty()) {
			continue;
		}

		string physicalName;
		gmsh::model::getPhysicalName(dim, physicalTags[0], physicalName);

		// Get the mesh elements for the entity (dim, tag):
		std::vector<int> elemTypes;
		std::vector<std::vector<std::size_t> > elemTags, elemNodeTags;
		gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim, tag);

		if (dim == 1) {
			// for boundary lines
			for (USI t = 0; t < elemTypes.size(); t++) {
				for (OCP_USI l = 0; l < elemTags[t].size(); l++) {
					const OCP_USI bId = elemNodeTags[t][2 * l];
					const OCP_USI eId = elemNodeTags[t][2 * l + 1];
					edges.insert(Edge(bId - 1, eId - 1, elemTags[t][l], physicalName));
					lineTag++;
				}
			}			
		}
		else if (dim == 2) {

			INT faciesIndex = -1;
			for (INT f = 0; f < facies.size(); f++) {
				if (facies[f].name == physicalName) {
					faciesIndex = f;
					break;
				}				
			}
			if (faciesIndex == -1) {
				facies.push_back(Facies(physicalName));
				faciesIndex = facies.size() - 1;
			}

			// for triangle and quadrangle
			USI np = 0;
			for (USI t = 0; t < elemTypes.size(); t++) {
				if (elemTypes[t] == 2)  np = 3;  // for triangle
				else                    np = 4;  // for quadrangle

				vector<OCP_USI> indexFace(np);
				for (OCP_USI l = 0; l < elemTags[t].size(); l++) {
					indexFace.clear();					

					for (USI i = 0; i < np; i++) {
						const OCP_USI bId = elemNodeTags[t][np * l + (i % np)];
						const OCP_USI eId = elemNodeTags[t][np * l + (i + 1) % np];

						auto iter = edges.find(Edge(bId - 1, eId - 1));
						if (iter == edges.end()) {
							edges.insert(Edge(bId - 1, eId - 1, lineTag++, faceIndex, i));
						}
						else {
							iter->faceIndex.push_back(faceIndex);
							iter->faceIndex.push_back(i);
						}

						indexFace.push_back(elemNodeTags[t][np * l + i] - 1);
					}
					faceIndex++;
					elements.push_back(Polygon(indexFace, elemTags[t][l], physicalName));
					faciesNum.push_back(faciesIndex);
				}
			}
		}
	}
}


void GMSHGrid::Setup()
{
	if (dimen == 2) {
		CalAreaCenter2D();
		SetupConnAreaAndBoundary2D();
	}
}


void GMSHGrid::CalAreaCenter2D()
{
	for (auto& e : elements) {
		e.CalCenter(points);
		e.CalArea(points);
	}
}


void GMSHGrid::SetupConnAreaAndBoundary2D()
{
	for (const auto& e : edges) {

		if (e.faceIndex.size() == 2) {
			// boundary
			Polygon& element = elements[e.faceIndex[0]];
			element.location += (e.physical) + " & ";
			// Calculate effective area
			const OCP_USI   bId       = element.p[e.faceIndex[1]];
			const OCP_USI   eId       = element.p[(e.faceIndex[1] + 1) % (element.p.size())];
			const Point3D&& edgeNode0 = Point3D(&points[3 * bId]);
			const Point3D&& edgeNode1 = Point3D(&points[3 * eId]);
			const Point3D&& edgeNormal{ -(edgeNode0 - edgeNode1).y, (edgeNode0 - edgeNode1).x, 0 };
			const Point3D&& center2edge = 0.5 * (edgeNode0 + edgeNode1) - element.center;
			e.area.push_back(fabs((center2edge * edgeNormal) / sqrt(center2edge * center2edge)));
		}
		else {
			// internal edge
			for (USI i = 0; i < 2; i++) {
				const Polygon& element = elements[e.faceIndex[2 * i]];
				// Calculate effective area
				const OCP_USI   bId       = element.p[e.faceIndex[2 * i + 1]];
				const OCP_USI   eId       = element.p[(e.faceIndex[2 * i + 1] + 1) % (element.p.size())];
				const Point3D&& edgeNode0 = Point3D(&points[3 * bId]);
				const Point3D&& edgeNode1 = Point3D(&points[3 * eId]);
				const Point3D&& edgeNormal{ -(edgeNode0 - edgeNode1).y, (edgeNode0 - edgeNode1).x, 0 };
				const Point3D&& center2edge = 0.5 * (edgeNode0 + edgeNode1) - element.center;
				e.area.push_back(fabs((center2edge * edgeNormal) / sqrt(center2edge * center2edge)));
			}
		}
	}
}


void GMSHGrid::InputProperty(ifstream& ifs)
{
	if (elements.empty()) {
		OCP_ABORT("INPUT KEYWORD GMSH FIRST!");
	}

	while (!ifs.eof()) {

		string fbuf;
		getline(ifs, fbuf);

		if (fbuf == "GMSHPROEND") break;

		vector<string> vbuf;
		for (USI f = 0; f < facies.size(); f++) {
			if (fbuf == facies[f].name) {
				mapF2F.push_back(f);
				// Input for current facies
				while (true) {
					ReadLine(ifs, vbuf);
					if (vbuf[0] == "END") break;

					if (vbuf[0] == "*PORO")      facies[f].poro = stod(vbuf[1]);
					else if (vbuf[0] == "*PERM") {
						facies[f].kx = stod(vbuf[1]);
						facies[f].ky = stod(vbuf[1]);
						facies[f].kz = stod(vbuf[1]);
					}
				}
			}
		}
	}
}


#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Sep/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/
