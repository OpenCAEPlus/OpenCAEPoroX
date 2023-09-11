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


void GMSHGrid::Input(const string& file)
{

	gmsh::initialize();
	gmsh::open(file);

	// Print the model name and dimension:
	std::string name;
	gmsh::model::getCurrent(name);
	std::cout << "Model " << name << " (" << gmsh::model::getDimension() << "D)\n";
	dimen = gmsh::model::getDimension();

	if (dimen == 2) {
		Input2D(file);
	}
	gmsh::finalize();

	Setup();
}

void GMSHGrid::Input2D(const string& file)
{
	// Get all mesh nodes
	std::vector<std::size_t> nodeTags;
	std::vector<double> nodeParams;
	gmsh::model::mesh::getNodes(nodeTags, points, nodeParams, -1, -1);

	// Get all the elementary entities in the model, as a vector of (dimension, tag) pairs:
	std::vector<std::pair<int, int> > entities;
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
					edges.insert(Edge(bId, eId, elemTags[t][l], physicalName));
					lineTag++;
				}
			}			
		}
		else if (dim == 2) {
			// for triangle and quadrangle
			USI np = 0;
			for (USI t = 0; t < elemTypes.size(); t++) {
				if (elemTypes[t] == 2) {
					np = 3;  // for triangle
				}
				else {
					np = 4;  // for quadrangle
				}

				vector<OCP_USI> indexFace(np);
				for (OCP_USI l = 0; l < elemTags[t].size(); l++) {
					indexFace.clear();					

					for (USI i = 0; i < np; i++) {
						const OCP_USI bId = elemNodeTags[t][np * l + (i % np)];
						const OCP_USI eId = elemNodeTags[t][np * l + (i + 1) % np];

						auto iter = edges.find(Edge(bId, eId));
						if (iter == edges.end()) {
							edges.insert(Edge(bId, eId, lineTag++, faceIndex, i));
						}
						else {
							iter->faceIndex.push_back(faceIndex);
							iter->faceIndex.push_back(i);
						}

						indexFace.push_back(elemNodeTags[t][np * l + i]);
					}
					faceIndex++;
					elements.push_back(Polygon(indexFace, elemTags[t][l], physicalName));
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
			const Point3D&& center2edge = 0.5 * (edgeNode0 + edgeNode1) - Point3D(element.center.data());
			e.area.push_back(abs((center2edge * edgeNormal) / sqrt(center2edge * center2edge)));
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
				const Point3D&& center2edge = 0.5 * (edgeNode0 + edgeNode1) - Point3D(element.center.data());
				e.area.push_back(abs((center2edge * edgeNormal) / sqrt(center2edge * center2edge)));
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
