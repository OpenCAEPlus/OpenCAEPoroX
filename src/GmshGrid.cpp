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
	std::vector<double> nodeCoords, nodeParams;
	gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, -1, -1);

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
					edges.insert(Edge(elemNodeTags[t][2 * l], elemNodeTags[t][2 * l + 1], elemTags[t][l], physicalName));
					lineTag++;
				}
			}			
		}
		else if (dim == 2) {
			// for triangle and quadrangle
			for (USI t = 0; t < elemTypes.size(); t++) {
				if (elemTypes[t] == 2) {
					for (OCP_USI l = 0; l < elemTags[t].size(); l++) {
						auto iter = edges.find(Edge(elemNodeTags[t][3 * l + 0], elemNodeTags[t][3 * l + 1]));
						if (iter == edges.end()) {
							edges.insert(Edge(elemNodeTags[t][3 * l + 0], elemNodeTags[t][3 * l + 1], lineTag++, faceIndex, 0));
						}
						else {
							iter->faceIndex.push_back(faceIndex);
							iter->faceIndex.push_back(0);
						}

						iter = edges.find(Edge(elemNodeTags[t][3 * l + 1], elemNodeTags[t][3 * l + 2]));
						if (iter == edges.end()) {
							edges.insert(Edge(elemNodeTags[t][3 * l + 1], elemNodeTags[t][3 * l + 2], lineTag++, faceIndex, 1));
						}
						else {
							iter->faceIndex.push_back(faceIndex);
							iter->faceIndex.push_back(1);
						}

						iter = edges.find(Edge(elemNodeTags[t][3 * l + 2], elemNodeTags[t][3 * l + 0]));
						if (iter == edges.end()) {
							edges.insert((Edge(elemNodeTags[t][3 * l + 2], elemNodeTags[t][3 * l + 0], lineTag++, faceIndex, 2)));
						}
						else {
							iter->faceIndex.push_back(faceIndex);
							iter->faceIndex.push_back(2);
						}
						faceIndex++;
						elements.push_back(Polygon(&nodeCoords[3 * (elemNodeTags[t][3 * l + 0] - 1)],
												   &nodeCoords[3 * (elemNodeTags[t][3 * l + 1] - 1)],
												   &nodeCoords[3 * (elemNodeTags[t][3 * l + 2] - 1)], elemTags[t][l], physicalName));
					}
				}
				else {
					// for quadrangle
					for (OCP_USI l = 0; l < elemTags[t].size(); l++) {
						auto iter = edges.find(Edge(elemNodeTags[t][4 * l + 0], elemNodeTags[t][4 * l + 1]));
						if (iter == edges.end()) {
							edges.insert(Edge(elemNodeTags[t][4 * l + 0], elemNodeTags[t][4 * l + 1], lineTag++, faceIndex, 0));
						}
						else {
							iter->faceIndex.push_back(faceIndex);
							iter->faceIndex.push_back(0);
						}

						iter = edges.find(Edge(elemNodeTags[t][4 * l + 1], elemNodeTags[t][4 * l + 2]));
						if (iter == edges.end()) {
							edges.insert(Edge(elemNodeTags[t][4 * l + 1], elemNodeTags[t][4 * l + 2], lineTag++, faceIndex, 1));
						}
						else {
							iter->faceIndex.push_back(faceIndex);
							iter->faceIndex.push_back(1);
						}

						iter = edges.find(Edge(elemNodeTags[t][4 * l + 2], elemNodeTags[t][4 * l + 3]));
						if (iter == edges.end()) {
							edges.insert((Edge(elemNodeTags[t][4 * l + 2], elemNodeTags[t][4 * l + 3], lineTag++, faceIndex, 2)));
						}
						else {
							iter->faceIndex.push_back(faceIndex);
							iter->faceIndex.push_back(2);
						}

						iter = edges.find(Edge(elemNodeTags[t][4 * l + 3], elemNodeTags[t][4 * l + 0]));
						if (iter == edges.end()) {
							edges.insert((Edge(elemNodeTags[t][4 * l + 3], elemNodeTags[t][4 * l + 0], lineTag++, faceIndex, 3)));
						}
						else {
							iter->faceIndex.push_back(faceIndex);
							iter->faceIndex.push_back(3);
						}

						faceIndex++;
						elements.push_back(Polygon(&nodeCoords[3 * (elemNodeTags[t][4 * l + 0] - 1)],
							&nodeCoords[3 * (elemNodeTags[t][4 * l + 1] - 1)],
							&nodeCoords[3 * (elemNodeTags[t][4 * l + 2] - 1)],
							&nodeCoords[3 * (elemNodeTags[t][4 * l + 3] - 1)], elemTags[t][l], physicalName));
					}
				}
				
			}
		}
	}
}



void GMSHGrid::Setup()
{
	if (dimen == 2) {
		SetupConnAreaAndBoundary2D();
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
			const Point3D& edgeNode0 = element.p[e.faceIndex[1]];
			const Point3D& edgeNode1 = element.p[(e.faceIndex[1] + 1) % (element.p.size())];
			Point3D        edgeNormal{ -(edgeNode0 - edgeNode1).y, (edgeNode0 - edgeNode1).x, 0 };
			Point3D        center2edge = 0.5 * (edgeNode0 + edgeNode1) - element.center;
			e.area.push_back(abs((center2edge * edgeNormal) / sqrt(center2edge * center2edge)));
		}
		else {
			// internal edge
			for (USI i = 0; i < 2; i++) {
				const Polygon& element = elements[e.faceIndex[2 * i]];
				// Calculate effective area
				const Point3D& edgeNode0 = element.p[e.faceIndex[2 * i + 1]];
				const Point3D& edgeNode1 = element.p[(e.faceIndex[2 * i + 1] + 1) % (element.p.size())];
				Point3D        edgeNormal{ -(edgeNode0 - edgeNode1).y, (edgeNode0 - edgeNode1).x, 0 };
				Point3D        center2edge = 0.5 * (edgeNode0 + edgeNode1) - element.center;
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
