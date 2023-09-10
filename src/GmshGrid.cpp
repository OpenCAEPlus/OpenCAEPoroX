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


	// Get all mesh nodes
	std::vector<std::size_t> nodeTags;
	std::vector<double> nodeCoords, nodeParams;
	gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, -1, -1);

	// Get all the elementary entities in the model, as a vector of (dimension, tag) pairs:
	std::vector<std::pair<int, int> > entities;
	gmsh::model::getEntities(entities);


	OCP_USI lineTag = 1;


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
							edges.insert(Edge(elemNodeTags[t][3 * l + 0], elemNodeTags[t][3 * l + 1], lineTag++, elemTags[t][l]));
						}
						else {
							iter->faceTag.push_back(elemTags[t][l]);
						}

						iter = edges.find(Edge(elemNodeTags[t][3 * l + 1], elemNodeTags[t][3 * l + 2]));
						if (iter == edges.end()) {
							edges.insert(Edge(elemNodeTags[t][3 * l + 1], elemNodeTags[t][3 * l + 2], lineTag++, elemTags[t][l]));
						}
						else {
							iter->faceTag.push_back(elemTags[t][l]);
						}

						iter = edges.find(Edge(elemNodeTags[t][3 * l + 2], elemNodeTags[t][3 * l + 0]));
						if (iter == edges.end()) {
							edges.insert((Edge(elemNodeTags[t][3 * l + 2], elemNodeTags[t][3 * l + 0], lineTag++, elemTags[t][l])));
						}
						else {
							iter->faceTag.push_back(elemTags[t][l]);
						}

						elements.push_back(Polygon(&nodeCoords[3 * (elemNodeTags[t][3 * l + 0] - 1)],
												   &nodeCoords[3 * (elemNodeTags[t][3 * l + 1] - 1)],
												   &nodeCoords[3 * (elemNodeTags[t][3 * l + 2] - 1)], elemTags[t][l]));
					}
				}
				else {
					// for quadrangle
					for (OCP_USI l = 0; l < elemTags[t].size(); l++) {
						auto iter = edges.find(Edge(elemNodeTags[t][4 * l + 0], elemNodeTags[t][4 * l + 1]));
						if (iter == edges.end()) {
							edges.insert(Edge(elemNodeTags[t][4 * l + 0], elemNodeTags[t][4 * l + 1], lineTag++, elemTags[t][l]));
						}
						else {
							iter->faceTag.push_back(elemTags[t][l]);
						}

						iter = edges.find(Edge(elemNodeTags[t][4 * l + 1], elemNodeTags[t][4 * l + 2]));
						if (iter == edges.end()) {
							edges.insert(Edge(elemNodeTags[t][4 * l + 1], elemNodeTags[t][4 * l + 2], lineTag++, elemTags[t][l]));
						}
						else {
							iter->faceTag.push_back(elemTags[t][l]);
						}

						iter = edges.find(Edge(elemNodeTags[t][4 * l + 2], elemNodeTags[t][4 * l + 3]));
						if (iter == edges.end()) {
							edges.insert((Edge(elemNodeTags[t][4 * l + 2], elemNodeTags[t][4 * l + 3], lineTag++, elemTags[t][l])));
						}
						else {
							iter->faceTag.push_back(elemTags[t][l]);
						}

						iter = edges.find(Edge(elemNodeTags[t][4 * l + 3], elemNodeTags[t][4 * l + 0]));
						if (iter == edges.end()) {
							edges.insert((Edge(elemNodeTags[t][4 * l + 3], elemNodeTags[t][4 * l + 0], lineTag++, elemTags[t][l])));
						}
						else {
							iter->faceTag.push_back(elemTags[t][l]);
						}

						elements.push_back(Polygon(&nodeCoords[3 * (elemNodeTags[t][4 * l + 0] - 1)],
							&nodeCoords[3 * (elemNodeTags[t][4 * l + 1] - 1)],
							&nodeCoords[3 * (elemNodeTags[t][4 * l + 2] - 1)],
							&nodeCoords[3 * (elemNodeTags[t][4 * l + 3] - 1)], elemTags[t][l]));
					}
				}
				
			}
		}
	}

	gmsh::finalize();
}


#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Sep/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/
