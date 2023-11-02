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
    InputDistParamGrid(prepro.preParamGridWell);
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


void Reservoir::InputDistParamGrid(PreParamGridWell& mygrid)
{

    if (CURRENT_RANK == MASTER_PROCESS) {
        OCP_INFO("Master Process Distribute Grid Params -- begin");
    }

    const PreParamGridWell& grid = mygrid;

    /////////////////////////////////////////////////////////////////////////
    // Distribute Grid-based data (lengths of all data are numGrid)
    /////////////////////////////////////////////////////////////////////////

    const VarInfoGrid varInfo(vector<VarInfo<vector<OCP_DBL>>> {
            VarInfo<vector<OCP_DBL>>{ "DX", &grid.dx, &bulk.vs.dx },
            VarInfo<vector<OCP_DBL>>{ "DY", &grid.dy, &bulk.vs.dy },
            VarInfo<vector<OCP_DBL>>{ "DZ", &grid.dz, &bulk.vs.dz },
            VarInfo<vector<OCP_DBL>>{ "V", &grid.v, &bulk.vs.v },
            VarInfo<vector<OCP_DBL>>{ "DEPTH", &grid.depth, &bulk.vs.depth },
            VarInfo<vector<OCP_DBL>>{ "PORO", &grid.poro, &bulk.vs.poroInit },
            VarInfo<vector<OCP_DBL>>{ "NTG", &grid.ntg, &bulk.vs.ntg },
            VarInfo<vector<OCP_DBL>>{ "PERMX", &grid.kx, &bulk.vs.rockKx },
            VarInfo<vector<OCP_DBL>>{ "PERMY", &grid.ky, &bulk.vs.rockKy },
            VarInfo<vector<OCP_DBL>>{ "PERMZ", &grid.kz, &bulk.vs.rockKz },
            VarInfo<vector<OCP_DBL>>{ "SIGMAV", &grid.sigma, &bulk.vs.sigma },
                VarInfo<vector<OCP_DBL>>{ "DZMTRXV", &grid.dzMtrx, &bulk.vs.dzMtrx },
                VarInfo<vector<OCP_DBL>>{ "SWAT", & grid.initR.swat, & bulk.INITm.GetSwat()},
            VarInfo<vector<OCP_DBL>>{ "SWATINIT", &grid.initR.swatInit, &bulk.optMs.scalePcow.GetSwatInit() },
            VarInfo<vector<OCP_DBL>>{ "BOUNDARYAREA", &grid.boundArea, &bulk.optMs.boundary.GetBoundaryArea() },
    }, vector<VarInfo<vector<OCP_USI>>>{
        VarInfo<vector<OCP_USI>>{"BOUNDARY", &grid.boundIndex, &bulk.optMs.boundary.GetBoundaryIndex() },
        VarInfo<vector<OCP_USI>>{ "SATNUM", &grid.SATNUM, &bulk.SATm.GetSATNUM() },
        VarInfo<vector<OCP_USI>>{ "PVTNUM", &grid.PVTNUM, &bulk.PVTm.GetPVTNUM() },
        VarInfo<vector<OCP_USI>>{ "ROCKNUM", &grid.ROCKNUM, &bulk.ROCKm.GetROCKNUM() },
    });


    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////

    // Get grid-based vars from params
    bulk.optMs.nb = domain.numGridLocal;
    bulk.vs.nb    = domain.numGridLocal;   // Interior + ghost
    bulk.vs.nbI   = domain.numGridInterior;

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
                
        MPI_Bcast(&send_var_info, 1, OCPMPI_INT, MASTER_PROCESS, myComm);

        // Master proc is excluded
        OCP_INT* numIgridEdgeproc = new OCP_INT[3 * numproc]();
        MPI_Gather(numIgridEdge, 3, OCPMPI_INT, numIgridEdgeproc, 3, OCPMPI_INT, MASTER_PROCESS, myComm);
      

        OCP_INT  maxNumElement = 0;
        OCP_INT  maxNumEdge    = 0;
        for (USI p = 1; p < numproc; p++) {       
            maxNumElement = max(maxNumElement, numIgridEdgeproc[3 * p]);
            maxNumEdge    = max(maxNumEdge, numIgridEdgeproc[3 * p + 2]);
        }
        OCP_ULL* iGridproc = new OCP_ULL[maxNumElement];

        // Send vars and conns to other process      
        vector<vector<OCP_CHAR>> send_buffer(1);
        // send_var_value, send_edge(direction, areaB, areaE, transmult)
        const USI varNumEdge = 4;
        OCP_USI maxSendSize = maxNumElement * send_var.numByte_total + maxNumEdge * varNumEdge * sizeof(OCP_DBL);

        send_buffer[0].resize(maxSendSize);
        vector<OCP_BOOL>        work_state{ OCP_FALSE };
        OCP_CHAR*               work_buffer;
        vector<OCP_USI>         proc_buffer_l;
        vector<OCP_USI>         proc_buffer_c; // which process use which buffer

        MPI_Request* request = new MPI_Request[numproc - 1];
        int          send_flag = 1;

        for (OCP_USI p = 1; p < numproc; p++) {

            MPI_Status status;
            MPI_Recv(iGridproc, numIgridEdgeproc[3 * p], OCPMPI_ULL, p, 0, myComm, &status);

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
            const OCP_USI  numGrid   = numIgridEdgeproc[3 * p];
            const OCP_ULL* gridIndex = iGridproc;
            vector<OCP_ULL> gridIndexA(numGrid);
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
            for (OCP_USI n = 0; n < numIgridEdgeproc[3 * p + 1]; n++) {
                for (const auto& gn : grid.gNeighbor[gridIndex[n]]) {
                    conn_ptr[conn_size++] = static_cast<OCP_DBL>(gn.Direct());
                    conn_ptr[conn_size++] = gn.AreaB();
                    conn_ptr[conn_size++] = gn.AreaE();
                    conn_ptr[conn_size++] = gn.TransMult();
                }
            }        
            send_size += conn_size * sizeof(OCP_DBL);        
            // send
            MPI_Isend((void*)work_buffer, send_size, OCPMPI_BYTE, p, 0, myComm, &request[p - 1]);
            // MPI_Send((void*)work_buffer, send_size, OCPMPI_BYTE, p, 0, myComm);
            // cout << "Third stage : 0 sends " << send_size << "b to " << p << endl;

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

        delete[] numIgridEdgeproc;
        delete[] iGridproc;

        vector<OCP_BOOL>().swap(work_state);
        vector<OCP_USI>().swap(proc_buffer_c);
        vector<OCP_USI>().swap(proc_buffer_l);

        // set self
        // dbl
        for (auto& s : send_var.var_dbl) {
            s.dst_ptr->resize(bulk.vs.nb);
            const OCP_DBL* src = s.src_ptr->data();           
            OCP_DBL*       dst = s.dst_ptr->data();
            for (OCP_USI n = 0; n < bulk.vs.nb; n++) {
                dst[n] = src[grid.map_Act2All[domain.grid[n]]];
            }
        }
        // usi
        for (auto& s : send_var.var_usi) {
            s.dst_ptr->resize(bulk.vs.nb);
            const OCP_USI* src = s.src_ptr->data();
            OCP_USI*       dst = s.dst_ptr->data();
            for (OCP_USI n = 0; n < bulk.vs.nb; n++) {
                dst[n] = src[grid.map_Act2All[domain.grid[n]]];
            }
        }

        // Get conn-based vars from grid
        vector<BulkConnPair>* dst = &conn.iteratorConn;
        domain.neighborNum.reserve(domain.numGridInterior);
        const auto  global_well_start = domain.numElementTotal - domain.numWellTotal;
        const auto& init2local = domain.init_global_to_local;
        OCP_USI           bId, eId;

        for (OCP_USI n = 0; n < bulk.vs.nbI; n++) {
            bId = n;
            domain.neighborNum.push_back(grid.gNeighbor[domain.grid[bId]].size() + 1);
            for (const auto& gn : grid.gNeighbor[domain.grid[bId]]) {
                if (gn.ID() >= global_well_start) {
                    // well connection
                    const USI wIndex = gn.ID() - global_well_start;
                    const USI len    = domain.wellWPB.size();
                    USI w = 0;
                    for (w = 0; w < len; w++) {
                        if (wIndex == domain.wellWPB[w][0]) {
                            domain.wellWPB[w].push_back(static_cast<OCP_USI>(gn.AreaB()));
                            domain.wellWPB[w].push_back(bId);
                            break;
                        }
                    }
                    if (w == len) {
                        domain.wellWPB.push_back(vector<OCP_USI>{
                            wIndex,
                            static_cast<OCP_USI>(gn.AreaB()), bId});
                    }
                    continue; 
                }
                eId = init2local.at(gn.ID());
                if (eId > bId)
                    dst->push_back(BulkConnPair(bId, eId, gn.Direct(), gn.AreaB(), gn.AreaE(), gn.TransMult()));
            }
        }
        conn.numConn = conn.iteratorConn.size();

        if (myrank == 1) {
            //for (OCP_USI n = 0; n < bulk.numBulk; n++) {
            //    cout << domain.grid[n] << "  " << bulk.poroInit[n] - 1 << endl;
            //}
            for (const auto& c : conn.iteratorConn) {
                cout << domain.grid[c.BId()] << "   " << domain.grid[c.EId()] << "   " << static_cast<USI>(c.Direction()) << "   "
                    << c.AreaB() << "   " << c.AreaE() << endl;
            }
            cout << conn.numConn << endl;
        }

        // Free grid Memory
        mygrid.FreeMemory();

        // wait   
        MPI_Waitall(numproc - 1, request, MPI_STATUS_IGNORE);
        delete[] request;
    }
    else {
        // Get num of kind of vars
        MPI_Bcast(&send_var_info, 1, OCPMPI_INT, MASTER_PROCESS, myComm);
        
        numIgridEdge[0] = domain.numGridLocal;
        numIgridEdge[1] = domain.numGridInterior;
        for (const auto& e : domain.elementCSR) {
            numIgridEdge[2] += e[2];    // well may be included
        }

        MPI_Gather(numIgridEdge, 3, OCPMPI_INT, 0, 3, OCPMPI_INT, MASTER_PROCESS, myComm);
            
        MPI_Send(domain.grid.data(), numIgridEdge[0], OCPMPI_ULL, MASTER_PROCESS, 0, myComm);

        // recv vars to from master process
        // decode first
        vector<OCP_INT> recv_var_info;
        recv_var_info.reserve(sizeof(OCP_INT) * 8);
        for (OCP_INT i = sizeof(OCP_INT) * 8 - 1; i >= 0; i--) {
            recv_var_info.push_back(((send_var_info >> i) & 1));
        }
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

        const USI varNumEdge = 4;
        MPI_Status status;
        OCP_USI   recv_size   = numIgridEdge[0] * recv_var.numByte_total + numIgridEdge[2] * varNumEdge * sizeof(OCP_DBL);
        OCP_CHAR* recv_buffer = new OCP_CHAR[recv_size]();

        MPI_Recv((void*)recv_buffer, recv_size, OCPMPI_BYTE, MASTER_PROCESS, 0, myComm, &status);
        // cout << "Third stage : " << myrank <<  " receives " << recv_size << "b from 0" << endl;
             
        // dbl
        OCP_DBL* dbl_ptr = (OCP_DBL*)recv_buffer;
        for (USI r = 0; r < recv_var.numvar_dbl; r++) {
            recv_var.var_dbl[r].dst_ptr->resize(bulk.vs.nb);
            const OCP_DBL* src = dbl_ptr + r * bulk.vs.nb;
            OCP_DBL*       dst = recv_var.var_dbl[r].dst_ptr->data();
            copy(src, src + bulk.vs.nb, dst);
        }      
        dbl_ptr += recv_var.numvar_dbl * bulk.vs.nb;
        // usi
        OCP_USI* usi_ptr = (OCP_USI*)dbl_ptr;
        for (USI r = 0; r < recv_var.numvar_usi; r++) {
            recv_var.var_usi[r].dst_ptr->resize(bulk.vs.nb);
            const OCP_USI* src = usi_ptr + r * bulk.vs.nb;
            OCP_USI*       dst = recv_var.var_usi[r].dst_ptr->data();
            copy(src, src + bulk.vs.nb, dst);
        }
        usi_ptr += recv_var.numvar_usi * bulk.vs.nb;

        // Get Conn from recv_buffer, only Interior grids' neighbors are passed
        OCP_DBL*     conn_ptr = (OCP_DBL*)usi_ptr;
        vector<BulkConnPair>* dst = &conn.iteratorConn;       
        domain.neighborNum.reserve(domain.numGridInterior);
        const auto  global_well_start = domain.numElementTotal - domain.numWellTotal;
        const auto& init2local = domain.init_global_to_local;
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
                                // bulk connection
                                eId = iter->second;
                                if (eId > bId) {
                                    dst->push_back(BulkConnPair(bId, eId, static_cast<ConnDirect>(conn_ptr[0]), conn_ptr[1], conn_ptr[2], conn_ptr[3]));
                                }
                            }
                            else {
                                // well connection
                                const USI wIndex = my_edge[j] - global_well_start;
                                const USI len    = domain.wellWPB.size();
                                USI w = 0;
                                for (w = 0; w < len; w++) {
                                    if (wIndex == domain.wellWPB[w][0]) {
                                        domain.wellWPB[w].push_back(static_cast<OCP_USI>(conn_ptr[1]));
                                        domain.wellWPB[w].push_back(bId);
                                        break;
                                    }
                                }
                                if (w == len) {
                                    domain.wellWPB.push_back(vector<OCP_USI>{
                                        wIndex,
                                        static_cast<OCP_USI>(conn_ptr[1]), bId});
                                }
                            }
                            conn_ptr += varNumEdge;
                        }
                    }
                }
            }
        }
        conn.numConn = conn.iteratorConn.size();
        set<OCP_USI>().swap(recv_proc);

        delete[] recv_buffer;
    }

    // Free memory
    vector<vector<idx_t>>().swap(domain.elementCSR);



#ifdef WITH_GMSH
    /////////////////////////////////////////////////////////////////////////
    // Distribute Boundary Name (used in Gmsh grid now)
    /////////////////////////////////////////////////////////////////////////


    if (myrank == MASTER_PROCESS) {
        vector<string> names;       
        INT numName = grid.gmshGrid.GetBoundaryName(names);
        MPI_Bcast(&numName, 1, OCPMPI_INT, MASTER_PROCESS, myComm);
        if (numName > 0) {
            vector<INT> lenName(numName);
            string      nameset;
            for (USI i = 0; i < numName; i++) {
                lenName[i] = names[i].size();
                nameset    += names[i];
            }
            MPI_Bcast(&lenName[0], numName, OCPMPI_INT, MASTER_PROCESS, myComm);
            INT len = 0;
            for (USI i = 0; i < numName; i++) {
                len += lenName[i];
            }                               
            MPI_Bcast(&nameset[0], len, OCPMPI_CHAR, MASTER_PROCESS, myComm);

            /// Set boundary name
            auto& boundaryName = bulk.optMs.boundary.GetBoundName();
            boundaryName       = names;
        }       
    }
    else {
        INT numName;
        MPI_Bcast(&numName, 1, OCPMPI_INT, MASTER_PROCESS, myComm);
        if (numName > 0) {
            vector<INT> lenName(numName);
            MPI_Bcast(&lenName[0], numName, OCPMPI_INT, MASTER_PROCESS, myComm);
            INT len = 0;
            for (USI i = 0; i < numName; i++) {
                len += lenName[i];
            }
            vector<OCP_CHAR> nameset(len);
            MPI_Bcast(&nameset[0], len, OCPMPI_CHAR, MASTER_PROCESS, myComm);

            /// Set boundary name
            auto& boundaryName = bulk.optMs.boundary.GetBoundName();
            USI bId = 0;
            for (USI i = 0; i < numName; i++) {
                boundaryName.push_back(string(&nameset[bId], lenName[i]));
                bId += lenName[i];
            }
        }
    }
#endif
  
    MPI_Barrier(myComm);


    if (CURRENT_RANK == MASTER_PROCESS) {
        OCP_INFO("Master Process Distribute Grid Params -- end");
    }
}


void Reservoir::InputDistParamOthers(const ParamRead& param)
{
    if (CURRENT_RANK == MASTER_PROCESS) {
        OCP_INFO("Input Reservoir Params -- begin");
    }

    bulk.InputParam(param.paramRs);
    conn.InputParam(param.paramRs, bulk);
    allWells.InputParam(param.paramWell, domain);

    if (CURRENT_RANK == MASTER_PROCESS) {
        OCP_INFO("Input Reservoir Params -- end");
    }
}


void Reservoir::Setup()
{
    OCP_FUNCNAME;

    bulk.Setup();
    allWells.Setup(bulk);
}

void Reservoir::ApplyControl(const USI& i)
{
    OCP_FUNCNAME;

    allWells.ApplyControl(i);
    allWells.SetupWellGroup(bulk);
}



void Reservoir::CalIPRT(const OCP_DBL& dt)
{
    OCP_FUNCNAME;
    // Calculate injection / production rate for current step
    allWells.CalIPRT(bulk, dt);
}


void Reservoir::PrintSolFIM(const string& outfile) const
{
    ofstream outu(outfile);
    if (!outu.is_open()) cout << "Can not open " << outfile << endl;
    const OCP_USI nb = bulk.vs.nb;
    const OCP_USI nc = bulk.vs.nc;

    for (OCP_USI n = 0; n < nb; n++) {
        // Pressure
        outu << bulk.vs.P[n] << "\n";
        // Ni
        for (USI i = 0; i < nc; i++) {
            outu << bulk.vs.Ni[n * nc + i] << "\n";
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