/*! \file    BulkInitializer.cpp
 *  \brief   BulkInitializer class definition
 *  \author  Shizhe Li
 *  \date    Aug/26/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "BulkInitializer.hpp"


void BulkInitializer::Setup(const ParamReservoir& rs_param, const OCPMixtureType& mixType)
{
	initType = rs_param.initType;
	rsTemp   = rs_param.rsTemp;
	// for hydrostatic equilibrium
	if (initType == InitType::EQUIL) {

		Equil tmpEquil;

		if (rs_param.PBVD_T.data.size() > 0) tmpEquil.PBVD.Setup(rs_param.PBVD_T.data[0]);
		switch (mixType)
		{
		case OCPMixtureType::BO_OW:
		case OCPMixtureType::THERMALK_OW:
			tmpEquil.Dref = rs_param.EQUIL[0].data[0];
			tmpEquil.Pref = rs_param.EQUIL[0].data[1];
			tmpEquil.DOWC = rs_param.EQUIL[0].data[2];
			tmpEquil.PcOW = rs_param.EQUIL[0].data[3];
			break;
		case OCPMixtureType::BO_OGW:
		case OCPMixtureType::COMP:
			tmpEquil.Dref = rs_param.EQUIL[0].data[0];
			tmpEquil.Pref = rs_param.EQUIL[0].data[1];
			tmpEquil.DOWC = rs_param.EQUIL[0].data[2];
			tmpEquil.PcOW = rs_param.EQUIL[0].data[3];
			tmpEquil.DGOC = rs_param.EQUIL[0].data[4];
			tmpEquil.PcGO = rs_param.EQUIL[0].data[5];
			break;
		default:
			OCP_ABORT("Wrong Initialization Type!");
			break;
		}

		EQUIL.push_back(tmpEquil);

		// Zi distribution
		for (auto& z : rs_param.ZMFVD_T.data) {
			initZi_Tab.push_back(OCPTable(z));
		}
		if (mixType == OCPMixtureType::COMP) {
			if (initZi_Tab.empty()) {
				OCP_ABORT("ZMFVD is MISSING in COMPOSITIONAL MODEL!");
			}
		}

		// Temperature distribution
		for (auto& t : rs_param.TEMPVD_T.data) {
			initT_Tab.push_back(OCPTable(t));
		}
		if (mixType == OCPMixtureType::COMP ||
			mixType == OCPMixtureType::THERMALK_OW) {
			if (initT_Tab.empty()) {
				// Use RTEMP
				initT_Tab.push_back(OCPTable(vector<vector<OCP_DBL>> {
					vector<OCP_DBL>{0}, vector<OCP_DBL>{rs_param.rsTemp}
				}));
			}
		}
	}
}


void BulkInitializer::Initialize(BulkVarSet& bvs, const PVTModule& pvtm, const SATModule& satm, const BulkOptionalModules& optMs, const Domain& domain)
{
	bvs.initType = initType;
	// for hydrostatic equilibrium
	if (initType == InitType::EQUIL) {
		InitHydroEquil(bvs, pvtm, satm, domain);
	}
	else if (initType == InitType::PTN) {
		InitPTNi(bvs);
	}
	else {
		InitHydroEquilW(bvs, pvtm, satm, optMs, domain);
	}
}


/// initialize reservoir with given P,T,Ni
void BulkInitializer::InitPTNi(BulkVarSet& bvs)
{
	bvs.P = P;

	if (T.empty()) {
		fill(bvs.T.begin(), bvs.T.end(), rsTemp);
	}
	else {
		bvs.T = T;
	}

	if (Ni.size() != bvs.nc) {
		OCP_ABORT("Ni is not given correctly!");
	}
	for (OCP_USI n = 0; n < bvs.nb; n++) {
		for (USI i = 0; i < bvs.nc; i++) {
			bvs.Ni[n * bvs.nc + i] = Ni[i][n];
		}
	}

	if (Pj.size() != bvs.np) {
		OCP_ABORT("Pj is not given correctly!");
	}
	for (OCP_USI n = 0; n < bvs.nb; n++) {
		for (USI j = 0; j < bvs.np; j++) {
			bvs.Pj[n * bvs.np + j] = Pj[j][n];
		}
	}

	// free memory
	vector<OCP_DBL>().swap(P);
	vector<OCP_DBL>().swap(T);
	for (USI i = 0; i < bvs.nc; i++) {
		vector<OCP_DBL>().swap(Ni[i]);
	}
	for (USI j = 0; j < bvs.np; j++) {
		vector<OCP_DBL>().swap(Pj[j]);
	}

	fill(bvs.phaseExist.begin(), bvs.phaseExist.end(), OCP_FALSE);
}


// for spe11 now
void BulkInitializer::InitHydroEquilW(BulkVarSet& bvs, const PVTModule& PVTm, const SATModule& SATm, const BulkOptionalModules& optMs, const Domain& domain)
{
	OCP_FUNCNAME;

	// all water phase and all water component
	const OCP_DBL         Dref = 0;
	const OCP_DBL         Pref = PRESSURE_STD;
	const vector<OCP_DBL> tmpInitZi{ 1, 0 };

	auto PVT = PVTm.GetPVT(0);

	/// calculate initial pressure for each bulk
	numNodes = 100;
	const OCP_DBL tabdz = (1.2 * 1E2 - 0) / (numNodes - 1);
	vector<OCP_DBL> Ztmp(numNodes, Dref);
	vector<OCP_DBL> Potmp(numNodes, 0);
	vector<OCP_DBL> Pgtmp(numNodes, 0);
	vector<OCP_DBL> Pwtmp(numNodes, 0);
	for (USI i = 1; i < numNodes; i++) {
		Ztmp[i] = Ztmp[i - 1] + tabdz;
	}

	const OCP_DBL Ttmp = 20;
	Pwtmp[0] = Pref;

	OCP_DBL gammaWtmp;
	// find the water pressure
	for (USI id = 0; id < numNodes - 1; id++) {
		gammaWtmp = GRAVITY_FACTOR * PVT->RhoPhase(Pwtmp[id], 0, Ttmp, tmpInitZi, PhaseType::wat);
		Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
	}

	Potmp = Pwtmp;
	Pgtmp = Pwtmp;
	OCPTable DepthP(vector<vector<OCP_DBL>>{Ztmp, Potmp, Pgtmp, Pwtmp});
	if (CURRENT_RANK == MASTER_PROCESS)
		DepthP.Display();

	std::vector<OCP_DBL> data(4, 0);
	for (OCP_USI n = 0; n < bvs.nb; n++) {

		for (USI i = 0; i < bvs.nc; i++) {
			bvs.Ni[n * bvs.nc + i] = tmpInitZi[i];
		}
		
		DepthP.Eval_All(0, bvs.depth[n], data);
		bvs.P[n] = data[3];
		bvs.Pj[n * bvs.np + bvs.g] = data[2];
		bvs.Pj[n * bvs.np + bvs.w] = data[3];

		bvs.S[n * bvs.np + bvs.g] = 0;
		bvs.S[n * bvs.np + bvs.w] = 1;
	}
}


void BulkInitializer::InitHydroEquil(BulkVarSet& bvs, const PVTModule& PVTm, const SATModule& SATm, const Domain& domain)
{
	OCP_FUNCNAME;

	OCP_DBL Dref = EQUIL[0].Dref;
	OCP_DBL Pref = EQUIL[0].Pref;
	OCP_DBL DOWC = EQUIL[0].DOWC;
	OCP_DBL PcOW = EQUIL[0].PcOW;
	OCP_DBL DOGC = EQUIL[0].DGOC;
	OCP_DBL PcGO = EQUIL[0].PcGO;
	OCP_DBL zRange[2];
	OCP_DBL zRangeTmp[2] = { 1E8,0 };  // min , max

	for (OCP_USI n = 0; n < bvs.nb; n++) {
		OCP_DBL temp1 = bvs.depth[n] - bvs.dz[n] / 2;
		OCP_DBL temp2 = bvs.depth[n] + bvs.dz[n] / 2;
		zRangeTmp[0] = zRangeTmp[0] < temp1 ? zRangeTmp[0] : temp1;
		zRangeTmp[1] = zRangeTmp[1] > temp2 ? zRangeTmp[1] : temp2;
	}
	zRangeTmp[1] *= -1;

	MPI_Allreduce(&zRangeTmp, &zRange, 2, OCPMPI_DBL, MPI_MIN, domain.global_comm);
	const OCP_DBL Zmin = zRange[0];
	const OCP_DBL Zmax = -zRange[1];
	OCP_DBL tabdz = (Zmax - Zmin) / (numNodes - 1);

	vector<OCP_DBL> Ztmp(numNodes, 0);
	vector<OCP_DBL> Potmp(numNodes, 0);
	vector<OCP_DBL> Pgtmp(numNodes, 0);
	vector<OCP_DBL> Pwtmp(numNodes, 0);

	vector<OCP_DBL> tmpInitZi(bvs.nc, 0);

	// cal Tab_Ztmp
	Ztmp[0] = Zmin;
	for (USI i = 1; i < numNodes; i++) {
		Ztmp[i] = Ztmp[i - 1] + tabdz;
	}

	OCP_DBL myTemp;

	// find the RefId
	USI beginId = 0;
	if (Dref <= Ztmp[0]) {
		beginId = 0;
	}
	else if (Dref >= Ztmp[numNodes - 1]) {
		beginId = numNodes - 1;
	}
	else {
		beginId =
			distance(Ztmp.begin(), find_if(Ztmp.begin(), Ztmp.end(),
				[s = Dref](auto& t) { return t > s; }));
		beginId--;
	}

	// begin calculating oil pressure:
	OCP_DBL Pbb = Pref;
	OCP_DBL gammaOtmp, gammaWtmp, gammaGtmp;
	OCP_DBL Ptmp;
	USI     mynum = 10;
	OCP_DBL mydz = 0;
	OCP_DBL Poref, Pgref, Pwref;
	OCP_DBL Pbegin = 0;

	const auto initZi_flag = initZi_Tab.size() > 0 ? OCP_TRUE : OCP_FALSE;
	const auto initT_flag = initT_Tab.size() > 0 ? OCP_TRUE : OCP_FALSE;
	const auto PBVD_flag = EQUIL[0].PBVD.IsEmpty() ? OCP_FALSE : OCP_TRUE;

	auto PVT = PVTm.GetPVT(0);

	if (Dref < DOGC && bvs.g >= 0) {
		// reference pressure is gas pressure
		Pgref = Pref;
		if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
		if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);
		if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, Dref, 1);


		gammaGtmp = GRAVITY_FACTOR *
			PVT->RhoPhase(Pgref, Pbb, myTemp, tmpInitZi, PhaseType::gas);
		Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
		Pgtmp[beginId] = Pbegin;

		// find the gas pressure
		for (USI id = beginId; id > 0; id--) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
			if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, Ztmp[id], 1);

			gammaGtmp = GRAVITY_FACTOR * PVT->RhoPhase(Pgtmp[id], Pbb, myTemp,
				tmpInitZi, PhaseType::gas);
			Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
		}

		for (USI id = beginId; id < numNodes - 1; id++) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
			if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, Ztmp[id], 1);

			gammaGtmp = GRAVITY_FACTOR * PVT->RhoPhase(Pgtmp[id], Pbb, myTemp,
				tmpInitZi, PhaseType::gas);
			Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
		}

		// find the oil pressure in Dref by Pgref
		Poref = 0;
		Ptmp = Pgref;
		mydz = (DOGC - Dref) / mynum;
		OCP_DBL myz = Dref;

		for (USI i = 0; i < mynum; i++) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
			if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, myz, 1);

			gammaGtmp = GRAVITY_FACTOR *
				PVT->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi, PhaseType::gas);
			Ptmp += gammaGtmp * mydz;
			myz += mydz;
		}
		Ptmp -= PcGO;
		for (USI i = 0; i < mynum; i++) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
			if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, myz, 1);

			gammaOtmp = GRAVITY_FACTOR *
				PVT->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi, PhaseType::oil);
			Ptmp -= gammaOtmp * mydz;
			myz -= mydz;
		}
		Poref = Ptmp;

		// find the oil pressure in tab
		if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
		if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);
		if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, Dref, 1);

		gammaOtmp = GRAVITY_FACTOR *
			PVT->RhoPhase(Poref, Pbb, myTemp, tmpInitZi, PhaseType::oil);
		Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
		Potmp[beginId] = Pbegin;

		for (USI id = beginId; id > 0; id--) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
			if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, Ztmp[id], 1);

			gammaOtmp = GRAVITY_FACTOR * PVT->RhoPhase(Potmp[id], Pbb, myTemp,
				tmpInitZi, PhaseType::oil);
			Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
		}

		for (USI id = beginId; id < numNodes - 1; id++) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
			if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, Ztmp[id], 1);

			gammaOtmp = GRAVITY_FACTOR * PVT->RhoPhase(Potmp[id], Pbb, myTemp,
				tmpInitZi, PhaseType::oil);
			Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
		}

		// find the water pressure in Dref by Poref
		Pwref = 0;
		Ptmp = Poref;
		mydz = (DOWC - Dref) / mynum;
		myz = Dref;

		for (USI i = 0; i < mynum; i++) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
			if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, myz, 1);

			gammaOtmp = GRAVITY_FACTOR * PVT->RhoPhase(Poref, Pbb, myTemp,
				tmpInitZi, PhaseType::oil);
			Ptmp += gammaOtmp * mydz;
			myz += mydz;
		}
		Ptmp -= PcOW;
		for (USI i = 0; i < mynum; i++) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);

			gammaWtmp = GRAVITY_FACTOR * PVT->RhoPhase(Ptmp, Pbb, myTemp,
				tmpInitZi, PhaseType::wat);
			Ptmp -= gammaWtmp * mydz;
			myz -= mydz;
		}
		Pwref = Ptmp;

		// find the water pressure in tab
		if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
		if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);

		gammaWtmp = GRAVITY_FACTOR *
			PVT->RhoPhase(Pwref, Pbb, myTemp, tmpInitZi, PhaseType::wat);
		Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
		Pwtmp[beginId] = Pbegin;

		for (USI id = beginId; id > 0; id--) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

			gammaWtmp = GRAVITY_FACTOR * PVT->RhoPhase(Pwtmp[id], Pbb, myTemp,
				tmpInitZi, PhaseType::wat);
			Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
		}

		for (USI id = beginId; id < numNodes - 1; id++) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

			gammaWtmp = GRAVITY_FACTOR * PVT->RhoPhase(Pwtmp[id], Pbb, myTemp,
				tmpInitZi, PhaseType::wat);
			Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
		}
	}
	else if (Dref > DOWC) {
		OCP_DBL myz;
		// reference pressure is water pressure
		if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
		if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);

		Pwref = Pref;
		gammaWtmp = GRAVITY_FACTOR *
			PVT->RhoPhase(Pwref, Pbb, myTemp, tmpInitZi, PhaseType::wat);
		Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
		Pwtmp[beginId] = Pbegin;

		// find the water pressure
		for (USI id = beginId; id > 0; id--) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

			gammaWtmp = GRAVITY_FACTOR * PVT->RhoPhase(Pwtmp[id], Pbb, myTemp,
				tmpInitZi, PhaseType::wat);
			Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		for (USI id = beginId; id < numNodes - 1; id++) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

			gammaWtmp = GRAVITY_FACTOR * PVT->RhoPhase(Pwtmp[id], Pbb, myTemp,
				tmpInitZi, PhaseType::wat);
			Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
		}

		// find the oil pressure in Dref by Pwref
		Poref = 0;
		Ptmp = Pwref;
		mydz = (DOWC - Dref) / mynum;
		myz = Dref;

		for (USI i = 0; i < mynum; i++) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);

			gammaWtmp = GRAVITY_FACTOR * PVT->RhoPhase(Ptmp, Pbb, myTemp,
				tmpInitZi, PhaseType::wat);
			Ptmp += gammaWtmp * mydz;
			myz += mydz;
		}
		Ptmp += PcOW;

		for (USI i = 0; i < mynum; i++) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
			if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, myz, 1);

			gammaOtmp = GRAVITY_FACTOR *
				PVT->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi, PhaseType::oil);
			Ptmp -= gammaOtmp * mydz;
			myz -= mydz;
		}
		Poref = Ptmp;

		// find the oil pressure in tab
		if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
		if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);
		if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, Dref, 1);

		gammaOtmp = GRAVITY_FACTOR *
			PVT->RhoPhase(Poref, Pbb, myTemp, tmpInitZi, PhaseType::oil);
		Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
		Potmp[beginId] = Pbegin;

		for (USI id = beginId; id > 0; id--) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
			if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, Ztmp[id], 1);

			gammaOtmp = GRAVITY_FACTOR * PVT->RhoPhase(Potmp[id], Pbb, myTemp,
				tmpInitZi, PhaseType::oil);
			Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
		}

		for (USI id = beginId; id < numNodes - 1; id++) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
			if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, Ztmp[id], 1);

			gammaOtmp = GRAVITY_FACTOR * PVT->RhoPhase(Potmp[id], Pbb, myTemp,
				tmpInitZi, PhaseType::oil);
			Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
		}

		if (bvs.g >= 0) {
			// find the gas pressure in Dref by Poref
			Pgref = 0;
			Ptmp = Poref;
			mydz = (DOGC - Dref) / mynum;
			myz = Dref;

			for (USI i = 0; i < mynum; i++) {
				if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
				if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
				if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, myz, 1);

				gammaOtmp =
					GRAVITY_FACTOR *
					PVT->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi, PhaseType::oil);
				Ptmp += gammaOtmp * mydz;
				myz += mydz;
			}
			Ptmp += PcGO;
			for (USI i = 0; i < mynum; i++) {
				if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
				if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
				if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, myz, 1);

				gammaGtmp =
					GRAVITY_FACTOR *
					PVT->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi, PhaseType::gas);
				Ptmp -= gammaGtmp * mydz;
				myz -= mydz;
			}
			Pgref = Ptmp;

			// find the gas pressure in tab
			if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);
			if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, Dref, 1);

			gammaGtmp = GRAVITY_FACTOR * PVT->RhoPhase(Pgref, Pbb, myTemp,
				tmpInitZi, PhaseType::gas);
			Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
			Pgtmp[beginId] = Pbegin;

			for (USI id = beginId; id > 0; id--) {
				if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
				if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
				if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, Ztmp[id], 1);

				gammaGtmp =
					GRAVITY_FACTOR * PVT->RhoPhase(Pgtmp[id], Pbb, myTemp,
						tmpInitZi, PhaseType::gas);
				Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
			}
			for (USI id = beginId; id < numNodes - 1; id++) {
				if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
				if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
				if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, Ztmp[id], 1);

				gammaGtmp =
					GRAVITY_FACTOR * PVT->RhoPhase(Pgtmp[id], Pbb, myTemp,
						tmpInitZi, PhaseType::gas);
				Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
			}
		}

	}
	else {
		OCP_DBL myz;
		// reference pressure is oil pressure
		Poref = Pref;
		if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
		if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);
		if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, Dref, 1);

		gammaOtmp = GRAVITY_FACTOR *
			PVT->RhoPhase(Poref, Pbb, myTemp, tmpInitZi, PhaseType::oil);
		Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
		Potmp[beginId] = Pbegin;

		// find the oil pressure
		for (USI id = beginId; id > 0; id--) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
			if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, Ztmp[id], 1);

			gammaOtmp = GRAVITY_FACTOR * PVT->RhoPhase(Potmp[id], Pbb, myTemp,
				tmpInitZi, PhaseType::oil);
			Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		for (USI id = beginId; id < numNodes - 1; id++) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
			if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, Ztmp[id], 1);

			gammaOtmp = GRAVITY_FACTOR * PVT->RhoPhase(Potmp[id], Pbb, myTemp,
				tmpInitZi, PhaseType::oil);
			Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
		}

		if (bvs.g >= 0) {
			// find the gas pressure in Dref by Poref
			Pgref = 0;
			Ptmp = Poref;
			mydz = (DOGC - Dref) / mynum;
			myz = Dref;

			for (USI i = 0; i < mynum; i++) {
				if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
				if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
				if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, myz, 1);

				gammaOtmp =
					GRAVITY_FACTOR *
					PVT->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi, PhaseType::oil);
				Ptmp += gammaOtmp * mydz;
				myz += mydz;
			}
			Ptmp += PcGO;
			for (USI i = 0; i < mynum; i++) {
				if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
				if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
				if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, myz, 1);

				gammaGtmp =
					GRAVITY_FACTOR *
					PVT->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi, PhaseType::gas);
				Ptmp -= gammaGtmp * mydz;
				myz -= mydz;
			}
			Pgref = Ptmp;

			// find the gas pressure in tab
			if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);
			if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, Dref, 1);

			gammaGtmp = GRAVITY_FACTOR * PVT->RhoPhase(Pgref, Pbb, myTemp,
				tmpInitZi, PhaseType::gas);
			Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
			Pgtmp[beginId] = Pbegin;

			for (USI id = beginId; id > 0; id--) {
				if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
				if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
				if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, Ztmp[id], 1);

				gammaGtmp =
					GRAVITY_FACTOR * PVT->RhoPhase(Pgtmp[id], Pbb, myTemp,
						tmpInitZi, PhaseType::gas);
				Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
			}

			for (USI id = beginId; id < numNodes - 1; id++) {
				if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
				if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
				if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, Ztmp[id], 1);

				gammaGtmp =
					GRAVITY_FACTOR * PVT->RhoPhase(Pgtmp[id], Pbb, myTemp,
						tmpInitZi, PhaseType::gas);
				Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
			}
		}

		// find the water pressure in Dref by Poref
		Pwref = 0;
		Ptmp = Poref;
		mydz = (DOWC - Dref) / mynum;
		myz = Dref;

		for (USI i = 0; i < mynum; i++) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
			if (PBVD_flag) Pbb = EQUIL[0].PBVD.Eval(0, myz, 1);

			gammaOtmp = GRAVITY_FACTOR *
				PVT->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi, PhaseType::oil);
			Ptmp += gammaOtmp * mydz;
			myz += mydz;
		}
		Ptmp -= PcOW;
		for (USI i = 0; i < mynum; i++) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);

			gammaWtmp = GRAVITY_FACTOR * PVT->RhoPhase(Ptmp, Pbb, myTemp,
				tmpInitZi, PhaseType::wat);
			Ptmp -= gammaWtmp * mydz;
			myz -= mydz;
		}
		Pwref = Ptmp;

		// find the water pressure in tab
		if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
		if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);

		gammaWtmp = GRAVITY_FACTOR *
			PVT->RhoPhase(Pwref, Pbb, myTemp, tmpInitZi, PhaseType::wat);
		Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
		Pwtmp[beginId] = Pbegin;

		for (USI id = beginId; id > 0; id--) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

			gammaWtmp = GRAVITY_FACTOR * PVT->RhoPhase(Pwtmp[id], Pbb, myTemp,
				tmpInitZi, PhaseType::wat);
			Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
		}

		for (USI id = beginId; id < numNodes - 1; id++) {
			if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
			if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

			gammaWtmp = GRAVITY_FACTOR * PVT->RhoPhase(Pwtmp[id], Pbb, myTemp,
				tmpInitZi, PhaseType::wat);
			Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
		}
	}

	OCPTable DepthP(vector<vector<OCP_DBL>>{Ztmp, Potmp, Pgtmp, Pwtmp});

	if (CURRENT_RANK == MASTER_PROCESS)
		DepthP.Display();

	// calculate Pc from DepthP to calculate Sj
	std::vector<OCP_DBL> data(4, 0);

	for (OCP_USI n = 0; n < bvs.nb; n++) {
		if (initZi_flag) {
			initZi_Tab[0].Eval_All0(bvs.depth[n], tmpInitZi);
			for (USI i = 0; i < bvs.nc; i++) {
				bvs.Ni[n * bvs.nc + i] = tmpInitZi[i];
			}
		}
		if (initT_flag) {
			myTemp = initT_Tab[0].Eval(0, bvs.depth[n], 1);
			bvs.initT[n] = myTemp;
			bvs.T[n] = myTemp;
		}


		DepthP.Eval_All(0, bvs.depth[n], data);
		const auto SAT = SATm.GetSAT(n);
		OCP_DBL Po = data[1];
		OCP_DBL Pg = data[2];
		OCP_DBL Pw = data[3];
		OCP_DBL Pcgo = Pg - Po;
		OCP_DBL Pcow = Po - Pw;
		OCP_DBL Sw = SAT->CalSwByPcow(Pcow);
		OCP_DBL Sg = 0;
		if (bvs.g >= 0) {
			Sg = SAT->CalSgByPcgo(Pcgo);
		}
		if (Sw + Sg > 1) {
			// should be modified
			OCP_DBL Pcgw = Pcow + Pcgo;
			Sw = SAT->CalSwByPcgw(Pcgw);
			Sg = 1 - Sw;
		}

		if (1 - Sw < TINY) {
			// all water
			Po = Pw + SAT->CalPcowBySw(1.0);
		}
		else if (1 - Sg < TINY) {
			// all gas
			Po = Pg - SAT->CalPcgoBySg(1.0);
		}
		else if (1 - Sw - Sg < TINY) {
			// water and gas
			Po = Pg - SAT->CalPcgoBySg(Sg);
		}
		bvs.P[n] = Po;

		if (bvs.o >= 0) {
			bvs.Pj[n * bvs.np + bvs.o] = Po;
		}
		if (bvs.g >= 0) {
			bvs.Pj[n * bvs.np + bvs.g] = Pg;
		}
		if (bvs.w >= 0) {
			bvs.Pj[n * bvs.np + bvs.w] = Pw;
		}

		if (bvs.depth[n] < DOGC) {
			Pbb = Po;
		}
		else if (PBVD_flag) {
			Pbb = EQUIL[0].PBVD.Eval(0, bvs.depth[n], 1);
		}
		bvs.Pb[n] = Pbb;

		// cal Sw
		const OCP_DBL swco = SAT->GetSwco();
		if (fabs(SAT->CalPcowBySw(0.0 - TINY)) < TINY &&
			fabs(SAT->CalPcowBySw(1.0 + TINY) < TINY)) {
			bvs.S[n * bvs.np + bvs.w] = swco;
			if (bvs.g >= 0) {
				bvs.S[n * bvs.np + bvs.g] = Sg;
			}
			if (bvs.o >= 0) {
				bvs.S[n * bvs.np + bvs.o] = 1 - Sg - swco;
			}
			continue;
		}

		Sw = 0;
		Sg = 0;
		const USI ncut    = 10;
		OCP_DBL   avePcow = 0;

		for (USI k = 0; k < ncut; k++) {
			OCP_DBL tmpSw = 0;
			OCP_DBL tmpSg = 0;
			OCP_DBL dep = bvs.depth[n] + bvs.dz[n] / ncut * (k - (ncut - 1) / 2.0);
			DepthP.Eval_All(0, dep, data);
			Po = data[1];
			Pg = data[2];
			Pw = data[3];
			Pcow = Po - Pw;
			Pcgo = Pg - Po;
			avePcow += Pcow;
			tmpSw = SAT->CalSwByPcow(Pcow);
			if (bvs.g >= 0) {
				tmpSg = SAT->CalSgByPcgo(Pcgo);
			}
			if (tmpSw + tmpSg > 1) {
				// should be modified
				OCP_DBL Pcgw = Pcow + Pcgo;
				tmpSw = SAT->CalSwByPcgw(Pcgw);
				tmpSg = 1 - tmpSw;
			}
			Sw += tmpSw;
			Sg += tmpSg;
		}
		Sw /= ncut;
		Sg /= ncut;
		avePcow /= ncut;

		SAT->SetupScale(n, Sw, avePcow);
		bvs.S[n * bvs.np + bvs.w] = Sw;
		if (bvs.g >= 0) {
			bvs.S[n * bvs.np + bvs.g] = Sg;
		}
		if (bvs.o >= 0) {
			bvs.S[n * bvs.np + bvs.o] = 1 - Sg - Sw;
		}
	}
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/26/2023      Create file                          */
/*----------------------------------------------------------------------------*/