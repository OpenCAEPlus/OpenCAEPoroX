/*! \file    Reservoir.cpp
 *  \brief   Reservoir class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "Reservoir.hpp"

/////////////////////////////////////////////////////////////////////
// General
/////////////////////////////////////////////////////////////////////


void Reservoir::InputParam(PreProcess& prepro, ParamRead& param)
{
    SetupDomain(prepro.domain);
    InputDistParamGrid(param.paramRs, prepro.preParamGridWell);
    InputDistParamOthers(param);
}


OCP_INT Reservoir::GetSendVarInfo(const VarInfoGrid& varInfo, VarInfoGrid& send_var)
{
    OCP_INT send_var_info = 0;
    USI count = 0;
   
    for (const auto& v : varInfo.var_dbl) {
        if (v.src_ptr->size() > 0) {
            send_var.var_dbl.push_back(v);
            send_var_info += pow(2, count);
        }
        count++;
    }
    
    for (const auto& v : varInfo.var_usi) {
        if (v.src_ptr->size() > 0) {
            send_var.var_usi.push_back(v);
            send_var_info += pow(2, count);
        }
        count++;
    }

    // the last bit represents the type of grid
    if (count >= sizeof(OCP_INT) * 8 - 1) {
        OCP_ABORT("Too many send vars!");
    }

    send_var.CalNumByte();

    return send_var_info;
}


void Reservoir::InputDistParamGrid(ParamReservoir& rsparam, PreParamGridWell& mygrid)
{

    const PreParamGridWell& grid = mygrid;
    const ParamReservoir& rs     = rsparam;

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////

    const VarInfoGrid varInfo(vector<VarInfo<vector<OCP_DBL>>> {
            VarInfo<vector<OCP_DBL>>{ "DX", &grid.dx, &bulk.dx },
            VarInfo<vector<OCP_DBL>>{ "DY", &grid.dy, &bulk.dy },
            VarInfo<vector<OCP_DBL>>{ "DZ", &grid.dz, &bulk.dz },
            VarInfo<vector<OCP_DBL>>{ "V", &grid.v, &bulk.v },
            VarInfo<vector<OCP_DBL>>{ "DEPTH", &grid.depth, &bulk.depth },
            VarInfo<vector<OCP_DBL>>{ "PORO", &grid.poro, &bulk.poroInit },
            VarInfo<vector<OCP_DBL>>{ "NTG", &grid.ntg, &bulk.ntg },
            VarInfo<vector<OCP_DBL>>{ "PERMX", &rs.permX, &bulk.rockKx },
            VarInfo<vector<OCP_DBL>>{ "PERMY", &rs.permY, &bulk.rockKy },
            VarInfo<vector<OCP_DBL>>{ "PERMZ", &rs.permZ, &bulk.rockKz },
            VarInfo<vector<OCP_DBL>>{ "SWAT", &rs.Swat, &optFeatures.scalePcow.swatInit },
            VarInfo<vector<OCP_DBL>>{ "THCONR", &rs.thconr, &bulk.thconr },
    }, vector<VarInfo<vector<OCP_USI>>>{
        VarInfo<vector<OCP_USI>>{ "SATNUM", &rs.SATNUM, &bulk.SATNUM },
        VarInfo<vector<OCP_USI>>{ "PVTNUM", &rs.PVTNUM, &bulk.PVTNUM },
        VarInfo<vector<OCP_USI>>{ "ROCKNUM", &rs.ROCKNUM, &bulk.ROCKNUM },
    });

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////

    // Get grid-based vars from params
    optFeatures.numBulk  = domain.numGridLocal;
    bulk.numBulk         = domain.numGridLocal;   // Interior + ghost
    bulk.numBulkInterior = domain.numGridInterior;

    MPI_Comm      myComm  = domain.myComm;
    const OCP_INT numproc = domain.numproc;
    const OCP_INT myrank  = domain.myrank;

    // numGridLocal, numGridInterior, num edge
    OCP_INT numIgridEdge[3] = { 0,0,0 };
    // coded information
    OCP_INT send_var_info   = 0;

    /////////////////////////////////////////////////////////////////////
    // Grid-Based && Conn-Based
    /////////////////////////////////////////////////////////////////////
    
    if (myrank == MASTER_PROCESS) {
        
        // Calculate the number of kind of grid data needs to be sent to other process
        VarInfoGrid send_var;
        send_var_info = GetSendVarInfo(varInfo, send_var);
        if (grid.activeGridNum < grid.numGrid) {
            send_var_info += pow(2, sizeof(OCP_INT) * 8 - 2);
            domain.allActive = OCP_FALSE;
        }
                
        MPI_Bcast(&send_var_info, 1, MPI_INT, MASTER_PROCESS, myComm);

        // Master proc is excluded
        OCP_INT* numIgridEdgeproc = new OCP_INT[3 * numproc]();
        MPI_Gather(numIgridEdge, 3, MPI_INT, numIgridEdgeproc, 3, MPI_INT, MASTER_PROCESS, myComm);        
      
        OCP_INT* dislps               = new OCP_INT[numproc + 1]();
        OCP_USI* numGridInterior_proc = new OCP_USI[numproc]();
        OCP_INT maxNumElement = 0;
        OCP_INT maxNumEdge    = 0;
        OCP_INT* numIgridproc = numIgridEdgeproc;
        for (USI p = 1; p < numproc; p++) {
            numGridInterior_proc[p] = numIgridEdgeproc[3 * p + 1];
            maxNumElement = max(maxNumElement, numIgridEdgeproc[3 * p]);
            maxNumEdge    = max(maxNumEdge, numIgridEdgeproc[3 * p + 2]);

            numIgridproc[p] = numIgridEdgeproc[3 * p];
            dislps[p + 1]   = numIgridproc[p] + dislps[p];
        }

        OCP_INT* iGridproc = new OCP_INT[dislps[numproc]];
        MPI_Gatherv(0, 0, MPI_INT, iGridproc, numIgridproc, dislps, MPI_INT, MASTER_PROCESS, myComm);
        delete[] numIgridEdgeproc;

        // Send vars and conns to other process      
        vector<vector<OCP_CHAR>> send_buffer(1);
        // send_var_value, send_edge(direction, areaB, areaE)
        OCP_USI maxSendSize = maxNumElement * send_var.numByte_total + maxNumEdge * 3 * sizeof(OCP_DBL);
        if (!domain.allActive) {
            maxSendSize += maxNumElement * sizeof(OCP_USI);
        }
        send_buffer[0].resize(maxSendSize);
        vector<OCP_BOOL>        work_state{ OCP_FALSE };
        OCP_CHAR*               work_buffer;
        vector<OCP_USI>         proc_buffer_l;
        vector<OCP_USI>         proc_buffer_c; // which process use which buffer

        MPI_Request* request = new MPI_Request[numproc - 1];
        int          send_flag = 1;

        for (OCP_USI p = 1; p < numproc; p++) {
            // get non-working buffer
            OCP_USI ws;
            for (ws = 0; ws < send_buffer.size(); ws++) {
                if (!work_state[ws]) {
                    work_buffer    = send_buffer[ws].data();
                    work_state[ws] = OCP_TRUE;
                    break;
                }
            }
            if (ws == send_buffer.size()) {
                send_buffer.push_back(send_buffer[0]);
                work_buffer = send_buffer.back().data();
                work_state.push_back(OCP_TRUE);
            }
            proc_buffer_c.push_back(p);
            proc_buffer_c.push_back(ws);

            // grid
            const OCP_USI  numGrid   = dislps[p + 1] - dislps[p];
            const OCP_INT* gridIndex = &iGridproc[dislps[p]];
            vector<OCP_USI> gridIndexA(numGrid);
            for (OCP_USI n = 0; n < numGrid; n++) {
                gridIndexA[n] = grid.map_Act2All[gridIndex[n]];
            }

            // dbl
            OCP_DBL* dbl_ptr = (OCP_DBL*)work_buffer;
            for (USI s = 0; s < send_var.numvar_dbl; s++) {
                const OCP_DBL* src = send_var.var_dbl[s].src_ptr->data();
                OCP_DBL* dst = dbl_ptr + s * numGrid;
                for (OCP_USI n = 0; n < numGrid; n++) {
                    dst[n] = src[gridIndexA[n]];
                }
            }
            dbl_ptr += send_var.numvar_dbl * numGrid;
            // usi
            OCP_USI* usi_ptr = (OCP_USI*)dbl_ptr;
            for (USI s = 0; s < send_var.numvar_usi; s++) {
                const OCP_USI* src = send_var.var_usi[s].src_ptr->data();
                OCP_USI* dst = usi_ptr + s * numGrid;
                for (OCP_USI n = 0; n < numGrid; n++) {
                    dst[n] = src[gridIndexA[n]];
                }
            }
            usi_ptr += send_var.numvar_usi * numGrid;
            // grid-conn
            OCP_DBL* conn_ptr = (OCP_DBL*)usi_ptr;
            OCP_USI send_size = numGrid * send_var.numByte_total;
            OCP_USI conn_size = 0;
            for (OCP_USI n = 0; n < numGridInterior_proc[p]; n++) {
                for (const auto& gn : grid.gNeighbor[gridIndex[n]]) {
                    conn_ptr[conn_size++] = gn.direction;
                    conn_ptr[conn_size++] = gn.areaB;
                    conn_ptr[conn_size++] = gn.areaE;
                }
            }        
            send_size += conn_size * sizeof(OCP_DBL);

            // global index(include inactive grid)
            if (!domain.allActive) {
                OCP_USI* agi = (OCP_USI*)(conn_ptr + conn_size);
                copy(gridIndexA.data(), gridIndexA.data() + numGrid, agi);
                send_size += numGrid * sizeof(gridIndexA[0]);
            }
           
            // send
            MPI_Isend((void*)work_buffer, send_size, MPI_BYTE, p, 0, myComm, &request[p - 1]);
            // MPI_Send((void*)work_buffer, send_size, MPI_BYTE, p, 0, myComm);
            cout << "Third stage : 0 sends " << send_size << "b to " << p << endl;

            // update work_state(request)
            for (USI i = 0; i < proc_buffer_c.size(); i+=2) {

                MPI_Test(&request[proc_buffer_c[i] - 1], &send_flag, MPI_STATUS_IGNORE);
                if (send_flag) {
                    work_state[proc_buffer_c[i + 1]] = OCP_FALSE;
                }else{
                    proc_buffer_l.push_back(proc_buffer_c[i]);
                    proc_buffer_l.push_back(proc_buffer_c[i + 1]);
                }
            }
            proc_buffer_c.swap(proc_buffer_l);
            proc_buffer_l.clear();
        }

        delete[] iGridproc;       
        delete[] dislps;
        delete[] numGridInterior_proc;

        vector<OCP_BOOL>().swap(work_state);
        vector<OCP_USI>().swap(proc_buffer_c);
        vector<OCP_USI>().swap(proc_buffer_l);

        // set self
        // dbl
        for (auto& s : send_var.var_dbl) {
            s.dst_ptr->resize(bulk.numBulk);
            const OCP_DBL* src = s.src_ptr->data();           
            OCP_DBL*       dst = s.dst_ptr->data();
            for (OCP_USI n = 0; n < bulk.numBulk; n++) {
                dst[n] = src[grid.map_Act2All[domain.grid[n]]];
            }
        }
        // usi
        for (auto& s : send_var.var_usi) {
            s.dst_ptr->resize(bulk.numBulk);
            const OCP_USI* src = s.src_ptr->data();
            OCP_USI*       dst = s.dst_ptr->data();
            for (OCP_USI n = 0; n < bulk.numBulk; n++) {
                dst[n] = src[grid.map_Act2All[domain.grid[n]]];
            }
        }

        // Get conn-based vars from grid
        vector<BulkPair>* dst = &conn.iteratorConn;
        domain.neighborNum.reserve(domain.numGridInterior);
        const OCP_USI     global_well_start = domain.numElementTotal - domain.numWellTotal;
        const map<OCP_USI, OCP_USI>& init2local = domain.init_global_to_local;
        OCP_USI           bId, eId;

        for (OCP_USI n = 0; n < bulk.numBulkInterior; n++) {
            bId = n;
            domain.neighborNum.push_back(grid.gNeighbor[domain.grid[bId]].size() + 1);
            for (const auto& gn : grid.gNeighbor[domain.grid[bId]]) {
                if (gn.id >= global_well_start)
                    continue; // well is exculude

                eId = init2local.at(gn.id);
                if (eId > bId)
                    dst->push_back(BulkPair(bId, eId, gn.direction, gn.areaB, gn.areaE));
            }
        }
        conn.numConn = conn.iteratorConn.size();

        // global index(include inactive grid)
        if (!domain.allActive) {
            domain.gridAllIndex.resize(bulk.numBulk);
            for (OCP_USI n = 0; n < bulk.numBulk; n++) {
                domain.gridAllIndex[n] = grid.map_Act2All[domain.grid[n]];
            }
        }

        if (myrank == 1) {
            //for (OCP_USI n = 0; n < bulk.numBulk; n++) {
            //    cout << domain.grid[n] << "  " << bulk.poroInit[n] - 1 << endl;
            //}
            for (const auto& c : conn.iteratorConn) {
                cout << domain.grid[c.BId()] << "   " << domain.grid[c.EId()] << "   " << c.Direction() << "   "
                    << c.AreaB() << "   " << c.AreaE() << endl;
            }
            cout << conn.numConn << endl;
        }

        // Free grid Memory
        mygrid.FreeMemory();
        rsparam.FreeGridMemory();

        // wait   
        MPI_Waitall(numproc - 1, request, MPI_STATUS_IGNORE);
        delete[] request;
    }
    else {
        // Get num of kind of vars
        MPI_Bcast(&send_var_info, 1, MPI_INT, MASTER_PROCESS, myComm);
        
        numIgridEdge[0] = domain.numGridLocal;
        numIgridEdge[1] = domain.numGridInterior;
        for (const auto& e : domain.elementCSR) {
            numIgridEdge[2] += e[2];    // well may be included
        }

        MPI_Gather(numIgridEdge, 3, MPI_INT, 0, 3, MPI_INT, MASTER_PROCESS, myComm);
            
        MPI_Gatherv(domain.grid.data(), numIgridEdge[0], MPI_INT, 0, 0, 0, MPI_INT, MASTER_PROCESS, myComm);

        // recv vars to from master process
        // decode first
        vector<OCP_INT> recv_var_info;
        recv_var_info.reserve(sizeof(OCP_INT) * 8);
        for (OCP_INT i = sizeof(OCP_INT) * 8 - 1; i >= 0; i--) {
            recv_var_info.push_back(((send_var_info >> i) & 1));
        }
        if (recv_var_info[1] == 1)  domain.allActive = OCP_FALSE;
        reverse(recv_var_info.begin(), recv_var_info.end());

        VarInfoGrid recv_var;
        USI count = 0;
        for (const auto& v : varInfo.var_dbl) {
            if (recv_var_info[count++]) {
                recv_var.var_dbl.push_back(v);
            }
        }
        for (const auto& v : varInfo.var_usi) {
            if (recv_var_info[count++]) {
                recv_var.var_usi.push_back(v);
            }
        }
        recv_var.CalNumByte();

        MPI_Status status;
        OCP_USI   recv_size   = numIgridEdge[0] * recv_var.numByte_total + numIgridEdge[2] * 3 * sizeof(OCP_DBL);
        if (!domain.allActive) {
            // recv global index(include inactive grid)
            recv_size += numIgridEdge[0] * sizeof(OCP_USI);
        }
        OCP_CHAR* recv_buffer = new OCP_CHAR[recv_size]();

        MPI_Recv((void*)recv_buffer, recv_size, MPI_BYTE, MASTER_PROCESS, 0, myComm, &status);       
        cout << "Third stage : " << myrank <<  " receives " << recv_size << "b from 0" << endl;
             
        // dbl
        OCP_DBL* dbl_ptr = (OCP_DBL*)recv_buffer;
        for (USI r = 0; r < recv_var.numvar_dbl; r++) {
            recv_var.var_dbl[r].dst_ptr->resize(bulk.numBulk);
            const OCP_DBL* src = dbl_ptr + r * bulk.numBulk;
            OCP_DBL*       dst = recv_var.var_dbl[r].dst_ptr->data();
            copy(src, src + bulk.numBulk, dst);
        }      
        dbl_ptr += recv_var.numvar_dbl * bulk.numBulk;
        // usi
        OCP_USI* usi_ptr = (OCP_USI*)dbl_ptr;
        for (USI r = 0; r < recv_var.numvar_usi; r++) {
            recv_var.var_usi[r].dst_ptr->resize(bulk.numBulk);
            const OCP_USI* src = usi_ptr + r * bulk.numBulk;
            OCP_USI*       dst = recv_var.var_usi[r].dst_ptr->data();
            copy(src, src + bulk.numBulk, dst);
        }
        usi_ptr += recv_var.numvar_usi * bulk.numBulk;

        // Get Conn from recv_buffer, only Interior grids' neighbors are passed
        OCP_DBL*     conn_ptr = (OCP_DBL*)usi_ptr;
        vector<BulkPair>* dst = &conn.iteratorConn;       
        domain.neighborNum.reserve(domain.numGridInterior);
        const OCP_USI     global_well_start = domain.numElementTotal - domain.numWellTotal;
        const map<OCP_USI, OCP_USI>& init2local = domain.init_global_to_local;
        OCP_USI           bId, eId;
        // Traverse elementCSR in ascending order of process number
        set<OCP_USI> recv_proc;
        for (const auto& e : domain.elementCSR) {
            recv_proc.insert(e[0]);
        }
        for (const auto& s : recv_proc) {
            for (const auto& e : domain.elementCSR) {
                if (e[0] == s) {
                    const idx_t* my_vtx  = &e[3];
                    const idx_t* my_xadj = &e[3 + e[1]];
                    const idx_t* my_edge = &e[3 + e[1] + e[1] + 1];
                    for (USI i = 0; i < e[1]; i++) {
                        if (my_vtx[i] >= global_well_start) {                   
                            continue;  // well is excluded
                        }
                        domain.neighborNum.push_back(my_xadj[i + 1] - my_xadj[i] + 1);
                        bId = init2local.at(my_vtx[i]);
                        for (USI j = my_xadj[i]; j < my_xadj[i + 1]; j++) {
                            const auto& iter = init2local.find(my_edge[j]);
                            if (iter != init2local.end()) {
                                // well is exculude
                                eId = iter->second;
                                if (eId > bId) {
                                    dst->push_back(BulkPair(bId, eId, conn_ptr[0], conn_ptr[1], conn_ptr[2]));
                                }
                            }
                            conn_ptr += 3;
                        }
                    }
                }
            }
        }
        conn.numConn = conn.iteratorConn.size();
        set<OCP_USI>().swap(recv_proc);

        // record global index of grid if inactive grid exist
        if (!domain.allActive) {
            domain.gridAllIndex.resize(bulk.numBulk);
            const OCP_USI* src = (OCP_USI*)conn_ptr;
            OCP_USI*       dst = domain.gridAllIndex.data();
            copy(src, src + bulk.numBulk, dst);
        }

        delete[] recv_buffer;

        if (myrank == 0) {
            //cout << bulk.numBulk << endl;
            //for (OCP_USI n = 0; n < bulk.numBulk; n++) {
            //    cout << domain.grid[n] << "  " << bulk.dx[n] << endl;
            //}
            //for (const auto& c : conn.iteratorConn) {
            //    cout << c.BId() << "   " << c.EId() << "   " << c.Direction() << "   "
            //        << c.AreaB() << "   " << c.AreaE() << endl;
            //    
            //}
            //cout << conn.numConn << endl;
        }
    }

    // Free memory
    vector<vector<idx_t>>().swap(domain.elementCSR);
}


void Reservoir::InputDistParamOthers(const ParamRead& param)
{
    domain.SetGirdDimens(param.paramRs.dimens.nx, param.paramRs.dimens.ny, param.paramRs.dimens.nz);
    bulk.InputParam(param.paramRs);
    allWells.InputParam(param.paramWell, domain);
    optFeatures.InputParam(param.paramRs);
}


void Reservoir::SetupIsoT()
{
    OCP_FUNCNAME;

    bulk.SetupIsoT(domain);
    conn.SetupIsoT(bulk);
    allWells.Setup(bulk);
    bulk.SetupOptionalFeatures(optFeatures);
}

void Reservoir::SetupT()
{
    bulk.SetupT(domain);
    conn.SetupT(bulk);
    allWells.Setup(bulk);
    bulk.SetupOptionalFeatures(optFeatures);
}

void Reservoir::ApplyControl(const USI& i)
{
    OCP_FUNCNAME;

    allWells.ApplyControl(i);
    allWells.SetupWellGroup(bulk);
}

void Reservoir::CalMaxChange()
{
    OCP_FUNCNAME;

    bulk.CalMaxChange();
    allWells.CalMaxBHPChange();
}

void Reservoir::CalIPRT(const OCP_DBL& dt)
{
    OCP_FUNCNAME;
    // Calculate injection / production rate for current step
    allWells.CalIPRT(bulk, dt);
    // Calculate Reinjection fluid for next step
    allWells.CalReInjFluid(bulk);
}

OCP_DBL Reservoir::CalCFL(const OCP_DBL& dt, const OCP_BOOL& ifComm) const
{
    fill(bulk.cfl.begin(), bulk.cfl.end(), 0.0);
    const USI np = bulk.numPhase;

    for (OCP_USI c = 0; c < conn.numConn; c++) {
        for (USI j = 0; j < np; j++) {
            const OCP_USI uId = conn.bcval.upblock[c * np + j];

            if (bulk.phaseExist[uId * np + j]) {
                bulk.cfl[uId * np + j] += fabs(conn.bcval.velocity[c * np + j]) * dt;
            }
        }
    }

    for (const auto& wl : allWells.wells) {
        if (wl.IsOpen() && wl.WellType() == PROD) {
            for (USI p = 0; p < wl.PerfNum(); p++) {
                if (wl.PerfState(p) == OPEN) {
                    const OCP_USI k = wl.PerfLocation(p);

                    for (USI j = 0; j < np; j++) {
                        bulk.cfl[k * np + j] += fabs(wl.PerfProdQj_ft3(p, j)) * dt;
                    }
                }
            }
        }
    }

    bulk.maxCFL_loc       = 0;
    const OCP_USI len = bulk.numBulk * np;
    for (OCP_USI n = 0; n < len; n++) {
        if (bulk.phaseExist[n] && bulk.vj[n] > TINY) {
            bulk.cfl[n] /= bulk.vj[n];
#ifdef DEBUG
            if (!isfinite(bulk.cfl[n])) {
                OCP_ABORT("cfl is nan!");
            }
#endif // DEBUG
            if (bulk.maxCFL_loc < bulk.cfl[n]) bulk.maxCFL_loc = bulk.cfl[n];
        }
    }
    if (ifComm) {

        GetWallTime timer;
        timer.Start();

        MPI_Allreduce(&bulk.maxCFL_loc, &bulk.maxCFL, 1, MPI_DOUBLE, MPI_MAX, domain.myComm);

        OCPTIME_COMM_COLLECTIVE += timer.Stop() / 1000;
    }
    else {
        bulk.maxCFL = bulk.maxCFL_loc;
    }
    return bulk.maxCFL;
}

void Reservoir::PrintSolFIM(const string& outfile) const
{
    ofstream outu(outfile);
    if (!outu.is_open()) cout << "Can not open " << outfile << endl;
    const OCP_USI nb = bulk.numBulk;
    const OCP_USI nc = bulk.numCom;

    for (OCP_USI n = 0; n < nb; n++) {
        // Pressure
        outu << bulk.P[n] << "\n";
        // Ni
        for (USI i = 0; i < nc; i++) {
            outu << bulk.Ni[n * nc + i] << "\n";
        }
    }
    // Well Pressure
    for (USI w = 0; w < allWells.numWell; w++) {
        outu << allWells.GetWBHP(w) << "\n";
    }
    outu.close();
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/