/*!
 * turbomachinery.cpp
* \brief Headers of the turbomachinery subroutines of the SU2 solvers.
 * \author S. Vitale, O. Pfeifle
 * \version 4.1.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */


#include "../include/turbomachinery.hpp"

CTurboMachinery::CTurboMachinery(void) {


	/*--- Turbo-Performance array initialization ---*/

	  TotalStaticEfficiency         = NULL;
	  TotalTotalEfficiency 					= NULL;
	  KineticEnergyLoss							= NULL;
	  TotalPressureLoss							= NULL;
	  MassFlowIn										= NULL;
	  MassFlowOut										= NULL;
	  FlowAngleIn										= NULL;
	  FlowAngleIn_BC 								= NULL;
	  FlowAngleOut									= NULL;
	  EulerianWork									= NULL;
	  TotalEnthalpyIn								= NULL;
	  TotalEnthalpyIn_BC 						= NULL;
	  EntropyIn 										= NULL;
	  EntropyIn_BC 									= NULL;
	  PressureRatio									= NULL;
	  TotalPresureIn								= NULL;
		TotalTemperatureIn						= NULL;
	  EnthalpyOut										= NULL;
	  MachIn												= NULL;
	  MachOut												= NULL;
	  VelocityOutIs									= NULL;
		DensityIn											= NULL;
		PressureIn										= NULL;
		TurboVelocityIn								= NULL;
		DensityOut										= NULL;
		PressureOut										= NULL;
		TurboVelocityOut							= NULL;
		EnthalpyOutIs									= NULL;
		EntropyGen										= NULL;
		AbsFlowAngleIn								= NULL;
		TotalEnthalpyOut							= NULL;
		TotalEnthalpyOutIs						= NULL;
		TotalRothalpyIn								= NULL;
		TotalRothalpyOut							= NULL;
		AbsFlowAngleOut								= NULL;
		PressureOut_BC								= NULL;



}

CTurboMachinery::CTurboMachinery( CConfig *config, CGeometry *geometry) : CTurboMachinery() {


	unsigned short nMarkerTurboPerf = config->GetnMarker_TurboPerformance();
	unsigned short iMarker, iDim;

	nDim = geometry->GetnDim();

	/*--- Initializate quantities for turboperformace ---*/

	  TotalStaticEfficiency         = new su2double[nMarkerTurboPerf];
	  TotalTotalEfficiency 					= new su2double[nMarkerTurboPerf];
	  KineticEnergyLoss							= new su2double[nMarkerTurboPerf];
	  TotalPressureLoss							= new su2double[nMarkerTurboPerf];
	  MassFlowIn										= new su2double[nMarkerTurboPerf];
	  MassFlowOut										= new su2double[nMarkerTurboPerf];
	  FlowAngleIn										= new su2double[nMarkerTurboPerf];
	  FlowAngleIn_BC 								= new su2double[nMarkerTurboPerf];
	  FlowAngleOut									= new su2double[nMarkerTurboPerf];
	  EulerianWork									= new su2double[nMarkerTurboPerf];
	  TotalEnthalpyIn								= new su2double[nMarkerTurboPerf];
	  TotalEnthalpyIn_BC 						= new su2double[nMarkerTurboPerf];
	  EntropyIn 										= new su2double[nMarkerTurboPerf];
	  EntropyIn_BC 									= new su2double[nMarkerTurboPerf];
	  PressureRatio									= new su2double[nMarkerTurboPerf];
	  TotalPresureIn								= new su2double[nMarkerTurboPerf];
		TotalTemperatureIn						= new su2double[nMarkerTurboPerf];
	  EnthalpyOut										= new su2double[nMarkerTurboPerf];
	  MachIn												= new su2double*[nMarkerTurboPerf];
	  MachOut												= new su2double*[nMarkerTurboPerf];
	  VelocityOutIs									= new su2double[nMarkerTurboPerf];
		DensityIn											= new su2double[nMarkerTurboPerf];
		PressureIn										= new su2double[nMarkerTurboPerf];
		TurboVelocityIn								= new su2double*[nMarkerTurboPerf];
		DensityOut										= new su2double[nMarkerTurboPerf];
		PressureOut										= new su2double[nMarkerTurboPerf];
		TurboVelocityOut							= new su2double*[nMarkerTurboPerf];
		EnthalpyOutIs									= new su2double[nMarkerTurboPerf];
		EntropyGen										= new su2double[nMarkerTurboPerf];
		AbsFlowAngleIn								= new su2double[nMarkerTurboPerf];
		TotalEnthalpyOut							= new su2double[nMarkerTurboPerf];
		TotalEnthalpyOutIs						= new su2double[nMarkerTurboPerf];
		TotalRothalpyIn								= new su2double[nMarkerTurboPerf];
		TotalRothalpyOut							= new su2double[nMarkerTurboPerf];
		AbsFlowAngleOut								= new su2double[nMarkerTurboPerf];
		PressureOut_BC								= new su2double[nMarkerTurboPerf];

		for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++){
			MachIn[iMarker]													= new su2double[nDim];
			MachOut[iMarker]												= new su2double[nDim];
			TurboVelocityIn[iMarker]								= new su2double[nDim];
			TurboVelocityOut[iMarker]								= new su2double[nDim];
			for (iDim = 0; iDim < nDim; iDim++){
				MachIn[iMarker][iDim]										= 0.0;
				MachOut[iMarker][iDim]									= 0.0;
				TurboVelocityIn[iMarker][iDim]					= 0.0;
				TurboVelocityOut[iMarker][iDim]					= 0.0;
			}
		}

	  for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++){
	    TotalStaticEfficiency[iMarker]					= 0.0;
	    TotalTotalEfficiency[iMarker]						= 0.0;
	    KineticEnergyLoss[iMarker]							= 0.0;
	    TotalPressureLoss[iMarker]							= 0.0;
	    MassFlowIn[iMarker]											= 0.0;
	    MassFlowOut[iMarker]										= 0.0;
	    FlowAngleIn[iMarker]										= 0.0;
	    FlowAngleIn_BC[iMarker]									= 0.0;
	    FlowAngleOut[iMarker]										= 0.0;
	    EulerianWork[iMarker]										= 0.0;
	    TotalEnthalpyIn[iMarker]								= 0.0;
	    TotalEnthalpyIn_BC[iMarker]							= 0.0;
	    EntropyIn [iMarker]											= 0.0;
	    EntropyIn_BC [iMarker]									= 0.0;
	    PressureRatio[iMarker]									= 0.0;
	    PressureOut[iMarker]										= 0.0;
			TotalPresureIn[iMarker]									= 0.0;
			TotalTemperatureIn[iMarker]							= 0.0;
	    EnthalpyOut[iMarker]										= 0.0;
	    VelocityOutIs[iMarker]									= 0.0;
			DensityIn[iMarker]											= 0.0;
			PressureIn[iMarker]											= 0.0;
			DensityOut[iMarker]											= 0.0;
			PressureOut[iMarker]										= 0.0;
			EnthalpyOutIs[iMarker]									= 0.0;
			EntropyGen[iMarker]											= 0.0;
			AbsFlowAngleIn[iMarker]									= 0.0;
			TotalEnthalpyOut[iMarker]								= 0.0;
			TotalEnthalpyOutIs[iMarker]							= 0.0;
			TotalRothalpyIn[iMarker]								= 0.0;
			TotalRothalpyOut[iMarker]								= 0.0;
			AbsFlowAngleOut[iMarker]								= 0.0;
			PressureOut_BC[iMarker]									= 0.0;
	  }




}

CTurboMachinery::~CTurboMachinery(void) {
  unsigned short iVar, iMarker;

  /*--- Array deallocation ---*/

  if (TotalStaticEfficiency != NULL) delete []TotalStaticEfficiency;
  if (TotalTotalEfficiency  != NULL) delete []TotalTotalEfficiency;
  if (KineticEnergyLoss	    != NULL) delete []KineticEnergyLoss;
  if (TotalPressureLoss			!= NULL) delete []TotalPressureLoss;
  if (MassFlowIn						!= NULL) delete []MassFlowIn;
  if (MassFlowOut						!= NULL) delete []MassFlowOut;
  if (FlowAngleIn						!= NULL) delete []FlowAngleIn;
  if (FlowAngleIn_BC 				!= NULL) delete []FlowAngleIn_BC;
  if (FlowAngleOut					!= NULL) delete []FlowAngleOut;
  if (EulerianWork					!= NULL) delete []EulerianWork;
  if (TotalEnthalpyIn				!= NULL) delete []TotalEnthalpyIn;
  if (TotalEnthalpyIn_BC 		!= NULL) delete []TotalEnthalpyIn_BC;
  if (EntropyIn 						!= NULL) delete []EntropyIn;
  if (EntropyIn_BC 					!= NULL) delete []EntropyIn_BC;
  if (PressureRatio					!= NULL) delete []PressureRatio;
  if (TotalPresureIn				!= NULL) delete []TotalPresureIn;
  if (TotalTemperatureIn		!= NULL) delete []TotalTemperatureIn;
  if (EnthalpyOut						!= NULL) delete []EnthalpyOut;
  if (MachIn								!= NULL) delete []MachIn;
  if (MachOut								!= NULL) delete []MachOut;
  if (VelocityOutIs					!= NULL) delete []VelocityOutIs;
  if (DensityIn							!= NULL) delete []DensityIn;
  if (PressureIn						!= NULL) delete []PressureIn;
  if (TurboVelocityIn				!= NULL) delete []TurboVelocityIn;
  if (DensityOut						!= NULL) delete []DensityOut;
  if (PressureOut						!= NULL) delete []PressureOut;
  if (TurboVelocityOut			!= NULL) delete []TurboVelocityOut;
  if (EnthalpyOutIs					!= NULL) delete []EnthalpyOutIs;
  if (EntropyGen						!= NULL) delete []EntropyGen;
  if (AbsFlowAngleIn				!= NULL) delete []AbsFlowAngleIn;
  if (TotalEnthalpyOut			!= NULL) delete []TotalEnthalpyOut;
  if (TotalEnthalpyOutIs		!= NULL) delete []TotalEnthalpyOutIs;
  if (TotalRothalpyIn				!= NULL) delete []TotalRothalpyIn;
  if (TotalRothalpyOut			!= NULL) delete []TotalRothalpyOut;
  if (AbsFlowAngleOut				!= NULL) delete []AbsFlowAngleOut;
  if (PressureOut_BC				!= NULL) delete []PressureOut_BC;

}

//
//void CTurboMachinery::TurboPerformance(CSolver *solver, CConfig *config, CGeometry *geometry){
//
//	unsigned short iMarker, iMarkerTP;
//	su2double  avgVel2In, avgVel2Out,avgVelRel2In, avgVelRel2Out, avgGridVel2In, avgGridVel2Out, avgTotalEnthalpyIn= 0.0,avgTotalRothalpyIn,
//	avgTotalEnthalpyOut, avgTotalRothalpyOut, avgTotalEnthalpyOutIs, avgEnthalpyOut, avgEnthalpyOutIs,
//	avgPressureOut, avgTotalRelPressureIn, avgTotalRelPressureOut, avgEntropyIn, avgEntropyOut, flowAngleIn, massFlowIn, tangMachIn, normalMachIn, 	flowAngleOut,
//	massFlowOut, tangMachOut, normalMachOut, avgTotTempIn, avgTotPresIn, P_Total, T_Total, *FlowDir, alphaIn_BC, entropyIn_BC, totalEnthalpyIn_BC, densityIn_Mix,
//	pressureIn_Mix, normalVelocityIn_Mix, tangVelocityIn_Mix, densityOut_Mix, pressureOut_Mix, normalVelocityOut_Mix, tangVelocityOut_Mix, absFlowAngleIn,
//	absFlowAngleOut, pressureOut_BC;
//	su2double pitch;
//	unsigned short iDim, i, n1, n2, n1t,n2t;
//
//	int rank = MASTER_NODE;
//	int size = SINGLE_NODE;
//	int markerTP;
//	string Marker_Tag;
//
//#ifdef HAVE_MPI
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//  MPI_Comm_size(MPI_COMM_WORLD, &size);
//  su2double *TurbPerfIn= NULL,*TurbPerfOut= NULL;
//  su2double *TotTurbPerfIn = NULL,*TotTurbPerfOut = NULL;
//  int *TotMarkerTP;
//  n1  = 18;
//  n2  = 18;
//  n1t = n1*size;
//  n2t = n2*size;
//  TurbPerfIn = new su2double[n1];
//  TurbPerfOut = new su2double[n2];
//
//  for (i=0;i<n1;i++)
//  	TurbPerfIn[i]= -1.0;
//  for (i=0;i<n2;i++)
//		TurbPerfOut[i]= -1.0;
//#endif
//
//  avgTotalRothalpyIn     = -1.0;
//  avgTotalEnthalpyIn		 = -1.0;
//  avgEntropyIn           = -1.0;
//  avgTotalRelPressureIn  = -1.0;
//	flowAngleIn						 = -1.0;
//	massFlowIn						 = -1.0;
//	tangMachIn						 = -1.0;
//	normalMachIn					 = -1.0;
//	avgTotTempIn					 = -1.0;
//	avgTotPresIn           = -1.0;
//	alphaIn_BC             = -1.0;
//	entropyIn_BC           = -1.0;
//	totalEnthalpyIn_BC     = -1.0;
//  avgTotalRothalpyOut    = -1.0;
//	avgTotalEnthalpyOut    = -1.0;
//	avgTotalRelPressureOut = -1.0;
//	avgPressureOut				 = -1.0;
//	avgEnthalpyOut				 = -1.0;
//	avgGridVel2Out				 = -1.0;
//	flowAngleOut					 = -1.0;
//	massFlowOut						 = -1.0;
//	tangMachOut						 = -1.0;
//	normalMachOut					 = -1.0;
//  densityIn_Mix					 = -1.0;
//	pressureIn_Mix				 = -1.0;
//	normalVelocityIn_Mix	 = -1.0;
//	tangVelocityIn_Mix		 = -1.0;
//	densityOut_Mix				 = -1.0;
//	pressureOut_Mix				 = -1.0;
//	normalVelocityOut_Mix  = -1.0;
//	tangVelocityOut_Mix		 = -1.0;
//	absFlowAngleIn         = -1.0;
//	absFlowAngleOut        = -1.0;
//	pressureOut_BC         = -1.0;
//	markerTP				       = -1;
//
//
//	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
//		for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
//			if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
//				Marker_Tag         = config->GetMarker_All_TagBound(iMarker);
//				pitch = config->GetPeriodicRotAngles(iMarkerTP -1);
//				// to avoid nan with 2D axial case.
//				if(pitch <= EPS) pitch=2*PI_NUMBER;
//				/*--- compute or retrieve inlet information ---*/
//				if (config->GetMarker_All_TurbomachineryFlag(iMarker) == INFLOW){
//					markerTP = iMarkerTP;
////					cout << markerTP -1 << " "<< config->GetMarker_All_TagBound(iMarker)<<endl;
//					avgVelRel2In= 0.0;
//					avgGridVel2In= 0.0;
//					avgVel2In= 0.0;
//					for (iDim = 0; iDim < nDim; iDim++){
//						if(AverageTurboVelocity[iMarker][0][1] >= 0.0){
//							avgVelRel2In +=( AverageVelocity[iMarker][0][iDim] - geometry->GetAverageGridVel(iMarker, 0)[iDim])*( AverageVelocity[iMarker][0][iDim] - geometry->GetAverageGridVel(iMarker, 0)[iDim]);
//						}
//						else{
//							avgVelRel2In +=( AverageVelocity[iMarker][0][iDim] + geometry->GetAverageGridVel(iMarker, 0)[iDim])*( AverageVelocity[iMarker][0][iDim] + geometry->GetAverageGridVel(iMarker, 0)[iDim]);
//						}
//						avgGridVel2In += geometry->GetAverageGridVel(iMarker, 0)[iDim]*geometry->GetAverageGridVel(iMarker, 0)[iDim];
//						avgVel2In += AverageVelocity[iMarker][0][iDim]*AverageVelocity[iMarker][0][iDim];
//					}
//
//					avgTotalRothalpyIn 			= AverageEnthalpy[iMarker][0] + 0.5*avgVelRel2In - 0.5*avgGridVel2In;
//					avgTotalEnthalpyIn 			= AverageEnthalpy[iMarker][0] + 0.5*avgVel2In;
//					avgEntropyIn 						= AverageEntropy[iMarker][0];
//					FluidModel->SetTDState_hs(avgTotalRothalpyIn, avgEntropyIn);
//					avgTotalRelPressureIn   = FluidModel->GetPressure();
//					flowAngleIn							= SpanFlowAngle[iMarker][0];
//					massFlowIn							= SpanMassFlow[iMarker][0]*2*PI_NUMBER/abs(pitch);
//					tangMachIn							= AverageTurboMach[iMarker][0][1];
//					normalMachIn						= AverageTurboMach[iMarker][0][0];
//					avgTotTempIn						= AverageTotTemperature[iMarker][0];
//					avgTotPresIn						= AverageTotPressure[iMarker][0];
//				  densityIn_Mix					  = AverageDensity[iMarker][0];
//					pressureIn_Mix				  = AveragePressure[iMarker][0];
//					normalVelocityIn_Mix	  = AverageTurboVelocity[iMarker][0][0];
//					tangVelocityIn_Mix		  = AverageTurboVelocity[iMarker][0][1];
//					absFlowAngleIn          = atan(AverageTurboVelocity[iMarker][0][1]/AverageTurboVelocity[iMarker][0][0]);
//
////TODO(turbo) better location has to be found for this computation, perhaps in the outputstructure file.
//					if(config->GetBoolNRBC() || config->GetBoolRiemann()){
//
//						if(config->GetBoolRiemann()){
//							P_Total  = config->GetRiemann_Var1(Marker_Tag);
//							T_Total  = config->GetRiemann_Var2(Marker_Tag);
//							FlowDir = config->GetRiemann_FlowDir(Marker_Tag);
//							alphaIn_BC = atan(FlowDir[1]/FlowDir[0]);
//							P_Total /= config->GetPressure_Ref();
//							T_Total /= config->GetTemperature_Ref();
//
//						}else{
//							if(config->GetKind_Data_NRBC(Marker_Tag) == TOTAL_CONDITIONS_PT){
//								P_Total  = config->GetNRBC_Var1(Marker_Tag);
//								T_Total  = config->GetNRBC_Var2(Marker_Tag);
//								FlowDir = config->GetNRBC_FlowDir(Marker_Tag);
//								alphaIn_BC = atan(FlowDir[1]/FlowDir[0]);
//								P_Total /= config->GetPressure_Ref();
//								T_Total /= config->GetTemperature_Ref();
//
//							}
//							else{
//								P_Total  = ExtAverageTotPressure[iMarker][0];
//								T_Total  = ExtAverageTotTemperature[iMarker][0];
//								alphaIn_BC = atan(ExtAverageTurboVelocity[iMarker][0][1]/ExtAverageTurboVelocity[iMarker][0][0]);
//							}
//						}
//
//
//
//						/* --- Computes the total state --- */
//						FluidModel->SetTDState_PT(P_Total, T_Total);
//						totalEnthalpyIn_BC= FluidModel->GetStaticEnergy()+ FluidModel->GetPressure()/FluidModel->GetDensity();
//						entropyIn_BC= FluidModel->GetEntropy();
//					}else{
//						cout << " Inlet BC convergence can't be checked "<<endl;
//						entropyIn_BC = 0.0;
//						totalEnthalpyIn_BC= 0.0;
//						alphaIn_BC =0.0;
//					}
//
//
//
//
//
//
//
//#ifdef HAVE_MPI
//					TurbPerfIn[0]  = avgTotalRothalpyIn;
//					TurbPerfIn[1]  = avgTotalEnthalpyIn;
//					TurbPerfIn[2]  = avgEntropyIn;
//					TurbPerfIn[3]  = avgTotalRelPressureIn;
//					TurbPerfIn[4]  = flowAngleIn;
//					TurbPerfIn[5]  = massFlowIn;
//					TurbPerfIn[6]  = tangMachIn;
//					TurbPerfIn[7]  = normalMachIn;
//					TurbPerfIn[8]  = avgTotTempIn;
//					TurbPerfIn[9]  = avgTotPresIn;
//					TurbPerfIn[10] = alphaIn_BC;
//					TurbPerfIn[11] = entropyIn_BC;
//					TurbPerfIn[12] = totalEnthalpyIn_BC;
//					TurbPerfIn[13] = densityIn_Mix;
//				  TurbPerfIn[14] = pressureIn_Mix;
//					TurbPerfIn[15] = normalVelocityIn_Mix;
//					TurbPerfIn[16] = tangVelocityIn_Mix;
//					TurbPerfIn[17] = absFlowAngleIn;
//#endif
//				}
//
//				/*--- compute or retrieve outlet information ---*/
//				if (config->GetMarker_All_TurbomachineryFlag(iMarker) == OUTFLOW){
//					avgVelRel2Out = 0.0;
//					avgGridVel2Out = 0.0;
//					avgVel2Out = 0.0;
//					for (iDim = 0; iDim < nDim; iDim++){
//						if(AverageTurboVelocity[iMarker][0][1]  >= 0.0){
//							avgVelRel2Out +=( AverageVelocity[iMarker][0][iDim] + geometry->GetAverageGridVel(iMarker, 0)[iDim])*( AverageVelocity[iMarker][0][iDim] + geometry->GetAverageGridVel(iMarker, 0)[iDim]);
//						}
//						else{
//							avgVelRel2Out +=( AverageVelocity[iMarker][0][iDim] - geometry->GetAverageGridVel(iMarker, 0)[iDim])*( AverageVelocity[iMarker][0][iDim] - geometry->GetAverageGridVel(iMarker, 0)[iDim]);
//						}
//						avgGridVel2Out += geometry->GetAverageGridVel(iMarker, 0)[iDim]*geometry->GetAverageGridVel(iMarker, 0)[iDim];
//						avgVel2Out += AverageVelocity[iMarker][0][iDim]*AverageVelocity[iMarker][0][iDim];
//					}
//					avgTotalRothalpyOut       = AverageEnthalpy[iMarker][0] + 0.5*avgVelRel2Out - 0.5*avgGridVel2Out;
//					avgTotalEnthalpyOut       = AverageEnthalpy[iMarker][0] + 0.5*avgVel2Out;
//					avgEntropyOut             = AverageEntropy[iMarker][0];
//					avgEnthalpyOut            = AverageEnthalpy[iMarker][0];
//					FluidModel->SetTDState_hs(avgTotalRothalpyOut, avgEntropyOut);
//					avgTotalRelPressureOut    =  FluidModel->GetPressure();
//					avgPressureOut						= AveragePressure[iMarker][0];
//					flowAngleOut							= SpanFlowAngle[iMarker][0];
//					massFlowOut							  = SpanMassFlow[iMarker][0]*2*PI_NUMBER/abs(pitch);
//					tangMachOut								= AverageTurboMach[iMarker][0][1];
//					normalMachOut						  = AverageTurboMach[iMarker][0][0];
//					densityOut_Mix					  = AverageDensity[iMarker][0];
//					pressureOut_Mix				  = AveragePressure[iMarker][0];
//					normalVelocityOut_Mix	  = AverageTurboVelocity[iMarker][0][0];
//					tangVelocityOut_Mix		  = AverageTurboVelocity[iMarker][0][1];
//					absFlowAngleOut          = atan(AverageTurboVelocity[iMarker][0][1]/AverageTurboVelocity[iMarker][0][0]);
//
//
//					if(config->GetBoolNRBC() || config->GetBoolRiemann()){
//
//						if(config->GetBoolRiemann()){
//							pressureOut_BC  = config->GetRiemann_Var1(Marker_Tag);
//							pressureOut_BC /= config->GetPressure_Ref();
//						}
//						else{
//							pressureOut_BC  = config->GetNRBC_Var1(Marker_Tag);
//							pressureOut_BC /= config->GetPressure_Ref();
//
//						}
//					}
//					else{
//						cout << " OUTLET BC convergence can't be checked "<<endl;
//						pressureOut_BC = 0.0;
//					}
//
//
//
//
//
//#ifdef HAVE_MPI
//					TurbPerfOut[0]  = avgTotalRothalpyOut;
//					TurbPerfOut[1]  = avgTotalEnthalpyOut;
//					TurbPerfOut[2]  = avgTotalRelPressureOut;
//					TurbPerfOut[3]  = avgPressureOut;
//					TurbPerfOut[4]  = avgEnthalpyOut;
//					TurbPerfOut[5]  = avgGridVel2Out;
//					TurbPerfOut[6]  = flowAngleOut;
//					TurbPerfOut[7]  = massFlowOut;
//					TurbPerfOut[8]  = tangMachOut;
//					TurbPerfOut[9]  = normalMachOut;
//					TurbPerfOut[10] = avgEntropyOut;
//					TurbPerfOut[11] = densityOut_Mix;
//					TurbPerfOut[12] = pressureOut_Mix;
//					TurbPerfOut[13] = normalVelocityOut_Mix;
//					TurbPerfOut[14] = tangVelocityOut_Mix;
//					TurbPerfOut[15] = absFlowAngleOut;
//					TurbPerfOut[16] = pressureOut_BC;
//					TurbPerfOut[17] = avgVel2Out;
//
//#endif
//
//				}
//			}
//		}
//	}
//
//#ifdef HAVE_MPI
//if (rank == MASTER_NODE){
//  TotTurbPerfIn = new su2double[n1t];
//  TotTurbPerfOut = new su2double[n2t];
//  for (i=0;i<n1t;i++)
//  	TotTurbPerfIn[i]= -1.0;
//  for (i=0;i<n2t;i++)
//		TotTurbPerfOut[i]= -1.0;
//	TotMarkerTP = new int[size];
//	for(i=0; i<size; i++){
//		TotMarkerTP[i] = -1;
//	}
//}
//	SU2_MPI::Gather(TurbPerfIn, n1, MPI_DOUBLE, TotTurbPerfIn, n1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
//	SU2_MPI::Gather(TurbPerfOut, n2, MPI_DOUBLE,TotTurbPerfOut, n2, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
//	SU2_MPI::Gather(&markerTP, 1, MPI_INT,TotMarkerTP, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
//	delete [] TurbPerfIn, delete [] TurbPerfOut;
//
//	if (rank == MASTER_NODE){
////		for (i=0;i<size;i++){
////			cout << TotMarkerTP[i] << endl;
////		}
//		for (i=0;i<size;i++){
//			if(TotTurbPerfIn[n1*i] > 0.0){
//				avgTotalRothalpyIn 		 = 0.0;
//				avgTotalRothalpyIn 		 = TotTurbPerfIn[n1*i];
//				avgTotalEnthalpyIn 		 = 0.0;
//				avgTotalEnthalpyIn 		 = TotTurbPerfIn[n1*i+1];
//				avgEntropyIn 					 = 0.0;
//				avgEntropyIn 			     = TotTurbPerfIn[n1*i+2];
//				avgTotalRelPressureIn  = 0.0;
//				avgTotalRelPressureIn  = TotTurbPerfIn[n1*i+3];
//				flowAngleIn						 = 0.0;
//				flowAngleIn						 = TotTurbPerfIn[n1*i+4];
//				massFlowIn						 = 0.0;
//				massFlowIn						 = TotTurbPerfIn[n1*i+5];
//				tangMachIn						 = 0.0;
//				tangMachIn						 = TotTurbPerfIn[n1*i+6];
//				normalMachIn					 = 0.0;
//				normalMachIn					 = TotTurbPerfIn[n1*i+7];
//				avgTotTempIn           = 0.0;
//				avgTotTempIn           = TotTurbPerfIn[n1*i+8];
//				avgTotPresIn           = 0.0;
//				avgTotPresIn           = TotTurbPerfIn[n1*i+9];
//				alphaIn_BC						 = 0.0;
//				alphaIn_BC 						 = TotTurbPerfIn[n1*i+10];
//				entropyIn_BC           = 0.0;
//				entropyIn_BC					 = TotTurbPerfIn[n1*i+11];
//				totalEnthalpyIn_BC		 = 0.0;
//				totalEnthalpyIn_BC     = TotTurbPerfIn[n1*i+12];
//				densityIn_Mix  				 = 0.0;
//				densityIn_Mix  				 = TotTurbPerfIn[n1*i+13];
//				pressureIn_Mix         = 0.0;
//				pressureIn_Mix         = TotTurbPerfIn[n1*i+14];
//				normalVelocityIn_Mix   = 0.0;
//				normalVelocityIn_Mix   = TotTurbPerfIn[n1*i+15];
//				tangVelocityIn_Mix     = 0.0;
//				tangVelocityIn_Mix     = TotTurbPerfIn[n1*i+16];
//				absFlowAngleIn         = 0.0;
//				absFlowAngleIn         = TotTurbPerfIn[n1*i+17];
//				markerTP               = -1;
//				markerTP               = TotMarkerTP[i];
////				cout << " I am assigning this value "<< TotMarkerTP[i] << endl;
//			}
//
//			if(TotTurbPerfOut[n2*i] > 0.0){
//				avgTotalRothalpyOut    = 0.0;
//				avgTotalRothalpyOut    = TotTurbPerfOut[n2*i];
//				avgTotalEnthalpyOut    = 0.0;
//				avgTotalEnthalpyOut    = TotTurbPerfOut[n2*i+1];
//				avgTotalRelPressureOut = 0.0;
//				avgTotalRelPressureOut = TotTurbPerfOut[n2*i+2];
//				avgPressureOut				 = 0.0;
//				avgPressureOut				 = TotTurbPerfOut[n2*i+3];
//				avgEnthalpyOut				 = 0.0;
//				avgEnthalpyOut				 = TotTurbPerfOut[n2*i+4];
//				avgGridVel2Out				 = 0.0;
//				avgGridVel2Out				 = TotTurbPerfOut[n2*i+5];
//				flowAngleOut					 = 0.0;
//				flowAngleOut					 = TotTurbPerfOut[n2*i+6];
//				massFlowOut						 = 0.0;
//				massFlowOut 					 = TotTurbPerfOut[n2*i+7];
//				tangMachOut						 = 0.0;
//				tangMachOut						 = TotTurbPerfOut[n2*i+8];
//				normalMachOut					 = 0.0;
//				normalMachOut					 = TotTurbPerfOut[n2*i+9];
//				avgEntropyOut					 = 0.0;
//				avgEntropyOut					 = TotTurbPerfOut[n2*i+10];
//				densityOut_Mix				 = 0.0;
//				densityOut_Mix         = TotTurbPerfOut[n2*i+11];
//				pressureOut_Mix				 = 0.0;
//				pressureOut_Mix        = TotTurbPerfOut[n2*i+12];
//				normalVelocityOut_Mix  = 0.0;
//				normalVelocityOut_Mix  = TotTurbPerfOut[n2*i+13];
//				tangVelocityOut_Mix    = 0.0;
//				tangVelocityOut_Mix    = TotTurbPerfOut[n2*i+14];
//				absFlowAngleOut        = 0.0;
//				absFlowAngleOut        = TotTurbPerfOut[n2*i+15];
//				pressureOut_BC         = 0.0;
//				pressureOut_BC         = TotTurbPerfOut[n2*i+16];
//				avgVel2Out						 = 0.0;
//				avgVel2Out						 = TotTurbPerfOut[n2*i+17];
//			}
//		}
//
//		delete [] TotTurbPerfIn, delete [] TotTurbPerfOut; delete [] TotMarkerTP;
//}
//
//#endif
//
//	if (rank == MASTER_NODE){
//
//
//		//IMPORTANT this approach of multi-zone performances rely upon the fact that turbomachinery markers follow the natural (stator-rotor) development of the real machine.
//
//		/*--- compute outlet isoentropic conditions ---*/
//		FluidModel->SetTDState_Ps(avgPressureOut, avgEntropyIn);
//		avgEnthalpyOutIs = FluidModel->GetStaticEnergy() + avgPressureOut/FluidModel->GetDensity();
//		avgTotalEnthalpyOutIs = avgEnthalpyOutIs + 0.5*avgVel2Out;
//
//		/*--- store turboperformance informations ---*/
//		PressureRatio[markerTP -1] = avgTotalRelPressureIn/avgPressureOut;
//
//
//		/*----Quantities needed for computing the turbomachinery performance -----*/
//		TotalPressureLoss[markerTP -1] 		= (avgTotalRelPressureIn - avgTotalRelPressureOut)/(avgTotalRelPressureOut - avgPressureOut) ;
//		KineticEnergyLoss[markerTP -1] 		= (avgEnthalpyOut - avgEnthalpyOutIs)/(avgTotalRothalpyIn - avgEnthalpyOut + 0.5*avgGridVel2Out);
//		EulerianWork[markerTP -1]      		= avgTotalEnthalpyIn - avgTotalEnthalpyOut;
//		TotalEnthalpyIn[markerTP -1]   		= avgTotalEnthalpyIn;
//		TotalEnthalpyOut[markerTP -1]   	= avgTotalEnthalpyOut;
//		TotalRothalpyIn[markerTP -1]   		= avgTotalRothalpyIn;
//		TotalRothalpyOut[markerTP -1]   	= avgTotalRothalpyOut;
//		TotalEnthalpyOutIs[markerTP -1]		=	avgTotalEnthalpyOutIs;
//		EntropyIn[markerTP -1]				 		= avgEntropyIn;
//		EntropyGen[markerTP -1]           = (avgEntropyOut - avgEntropyIn)/avgEntropyIn;
//		AbsFlowAngleIn[markerTP -1]       = absFlowAngleIn;
//		AbsFlowAngleOut[markerTP -1]      = absFlowAngleOut;
//		FlowAngleIn[markerTP -1]       		= flowAngleIn;
//		FlowAngleOut[markerTP -1]      		= flowAngleOut;
//		MassFlowIn[markerTP -1]        		= massFlowIn;
//		MassFlowOut[markerTP -1]       		= massFlowOut;
//		MachIn[markerTP -1][0]            = normalMachIn;
//		MachOut[markerTP -1][0]           = normalMachOut;
//		MachIn[markerTP -1][1]     				= tangMachIn;
//		MachOut[markerTP -1][1]    				= tangMachOut;
//		EnthalpyOut[markerTP -1]       		= avgEnthalpyOut;
//		EnthalpyOutIs[markerTP -1]        = avgEnthalpyOutIs;
//		VelocityOutIs[markerTP -1]    		= sqrt(2.0*(avgTotalRothalpyIn - avgEnthalpyOut + 0.5*avgGridVel2Out));
//
//		/*----Quantities needed for BC convergence test -----*/
//
//		TotalPresureIn[markerTP -1]       = avgTotPresIn;
//		TotalTemperatureIn[markerTP -1]   = avgTotTempIn;
//		FlowAngleIn_BC[markerTP -1] 			= alphaIn_BC;
//		EntropyIn_BC[markerTP -1]					= entropyIn_BC;
//		TotalEnthalpyIn_BC[markerTP -1]   = totalEnthalpyIn_BC;
//		DensityIn[markerTP -1]						= densityIn_Mix;
//		PressureIn[markerTP -1]					  = pressureIn_Mix;
//		TurboVelocityIn[markerTP -1][0] 	= normalVelocityIn_Mix;
//		TurboVelocityIn[markerTP -1][1]		= tangVelocityIn_Mix;
//		DensityOut[markerTP -1]						= densityOut_Mix;
//		PressureOut[markerTP -1]					= pressureOut_Mix;
//		TurboVelocityOut[markerTP -1][0] 	= normalVelocityOut_Mix;
//		TurboVelocityOut[markerTP -1][1]	= tangVelocityOut_Mix;
//		PressureOut_BC[markerTP -1]       = pressureOut_BC;
//
//
//
//
//	}
//}
//
//void CTurboMachinery::MultiStageTurboPerformance(CSolver *solver, CConfig *config){
//
//	unsigned short nBladesRow, nStages;
//	unsigned short iStage, iDim;
//	int rank = MASTER_NODE;
//	int size = SINGLE_NODE;
//#ifdef HAVE_MPI
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//  MPI_Comm_size(MPI_COMM_WORLD, &size);
//#endif
//
//
//	//IMPORTANT this approach of multi-zone performances rely upon the fact that turbomachinery markers follow the natural (stator-rotor) development of the real machine.
//
//  nBladesRow = config->GetnMarker_Turbomachinery();
//  nStages    = int(nBladesRow/2);
//  su2double  vel2out;
//  if (rank == MASTER_NODE){
//  	EulerianWork[nBladesRow + nStages]        = 0.0;
//  	/*---Comnpute performance for each stage---*/
//  	for(iStage = 0; iStage < nStages; iStage++ ){
//  		FluidModel->SetTDState_Ps(PressureOut[iStage*2 +1], EntropyIn[iStage*2]);
//  		EnthalpyOutIs[nBladesRow + iStage] 				 = FluidModel->GetStaticEnergy() + PressureOut[iStage*2 +1]/FluidModel->GetDensity();
//  		FluidModel->SetTDState_Prho(PressureOut[iStage*2 +1], DensityOut[iStage*2 +1]);
//  		vel2out = 0.0;
//  		for (iDim = 0; iDim<nDim;iDim++) vel2out += MachOut[iStage*2 +1][iDim]*MachOut[iStage*2 +1][iDim];
//  		vel2out /= FluidModel->GetSoundSpeed2();
//  		TotalEnthalpyOutIs[nBladesRow + iStage] 	 = EnthalpyOutIs[nBladesRow + iStage] + 0.5*vel2out;
//
//  		TotalTotalEfficiency[nBladesRow + iStage]  = (TotalEnthalpyIn[iStage*2] - TotalEnthalpyOut[iStage*2 + 1])/(TotalEnthalpyIn[iStage*2] - TotalEnthalpyOutIs[nBladesRow + iStage]);
//  		TotalStaticEfficiency[nBladesRow + iStage] = (TotalEnthalpyIn[iStage*2] - TotalEnthalpyOut[iStage*2 + 1])/(TotalEnthalpyIn[iStage*2] - EnthalpyOutIs[nBladesRow + iStage]);
//  		EntropyGen[nBladesRow + iStage]            = ((EntropyIn[iStage*2 + 1]*EntropyGen[iStage*2 + 1] + EntropyIn[iStage*2 + 1]) - EntropyIn[iStage*2])/abs(EntropyIn[iStage*2]);
//  		PressureRatio[nBladesRow + iStage]         = (PressureRatio[iStage*2]*PressureOut[iStage*2]/PressureOut[iStage*2 + 1]);
//  		MassFlowIn[nBladesRow + iStage]         	 = MassFlowIn[iStage*2];
//  		MassFlowOut[nBladesRow + iStage]         	 = MassFlowIn[iStage*2 + 1];
//
//  	}
//
//  	/*---Comnpute performance for full machine---*/
//  	FluidModel->SetTDState_Ps(PressureOut[nBladesRow-1], EntropyIn[0]);
//		EnthalpyOutIs[nBladesRow + nStages] 				 = FluidModel->GetStaticEnergy() + PressureOut[nBladesRow-1]/FluidModel->GetDensity();
//		FluidModel->SetTDState_Prho(PressureOut[nBladesRow-1], DensityOut[nBladesRow-1]);
//		vel2out = 0.0;
//		for (iDim = 0; iDim<nDim;iDim++) vel2out += MachOut[nBladesRow-1][iDim]*MachOut[nBladesRow-1][iDim];
//		vel2out /= FluidModel->GetSoundSpeed2();
//		TotalEnthalpyOutIs[nBladesRow + nStages] 	 = EnthalpyOutIs[nBladesRow + nStages] + 0.5*vel2out;
//
//		TotalTotalEfficiency[nBladesRow + nStages] = (TotalEnthalpyIn[0] - TotalEnthalpyOut[nBladesRow-1])/(TotalEnthalpyIn[0] - TotalEnthalpyOutIs[nBladesRow + nStages]);
//    TotalStaticEfficiency[nBladesRow +nStages] = (TotalEnthalpyIn[0] - TotalEnthalpyOut[nBladesRow-1])/(TotalEnthalpyIn[0] - EnthalpyOutIs[nBladesRow + nStages]);
//  	EntropyGen[nBladesRow + iStage]            = ((EntropyIn[nBladesRow-1]*EntropyGen[nBladesRow-1] + EntropyIn[nBladesRow-1]) - EntropyIn[0])/abs(EntropyIn[0]);
//    PressureRatio[nBladesRow + nStages]        = PressureRatio[0]*PressureOut[0]/PressureOut[nBladesRow-1];
//		MassFlowIn[nBladesRow + nStages]         	 = MassFlowIn[0];
//  	MassFlowOut[nBladesRow + nStages]          = MassFlowIn[nBladesRow-1];
//  }
//}
