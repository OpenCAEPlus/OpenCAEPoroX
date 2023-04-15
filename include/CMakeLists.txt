target_sources(
  OpenCAEPoroX
  PUBLIC AcceleratePVT.hpp
		 AllWells.hpp
		 BulkConn.hpp
		 Bulk.hpp
		 CornerGrid.hpp
		 DenseMat.hpp
		 Domain.hpp
		 FaspSolver.hpp
		 FlowUnit.hpp
		 IsothermalSolver.hpp
		 LinearSolver.hpp
		 LinearSystem.hpp
		 MixtureBO.hpp
		 MixtureComp.hpp
		 Mixture.hpp
		 MixtureThermal.hpp
		 OCPConst.hpp
		 OCPControl.hpp
		 OCPFluidMethod.hpp
		 OCP.hpp
		 OCPOutput.hpp
		 OCPStructure.hpp
		 OCPTable.hpp
		 OCPTimeRecord.hpp
		 OptionalFeatures.hpp
		 Output4Vtk.hpp
		 ParamControl.hpp
		 ParamOutput.hpp
		 ParamRead.hpp
		 ParamReservoir.hpp
		 ParamWell.hpp
		 PardisoSolver.hpp
		 Partition.hpp
		 PhasePermeability.hpp
		 PreParamGridWell.hpp
		 PreProcess.hpp
		 Reservoir.hpp
		 Rock.hpp
		 SamgSolver.hpp
		 Solver.hpp
		 ThermalMethod.hpp
		 ThermalSolver.hpp
		 UtilError.hpp
		 UtilInput.hpp
		 UtilOutput.hpp
		 UtilTiming.hpp
		 Well.hpp
		 WellOpt.hpp
		 WellPerf.hpp)
         
         
target_include_directories(OpenCAEPoroX PUBLIC ${CMAKE_CURRENT_LIST_DIR})

