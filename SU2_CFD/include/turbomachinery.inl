/*!
 * \file turbomachinery.inl
 * \brief In-Line subroutines of the <i>turbomachinery.hpp</i> file.
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

#pragma once


inline su2double CTurboMachinery::GetTotalPressureLoss(unsigned short inMarkerTP){return TotalPressureLoss[inMarkerTP];}

inline su2double CTurboMachinery::GetKineticEnergyLoss(unsigned short inMarkerTP){return KineticEnergyLoss[inMarkerTP];}

inline su2double CTurboMachinery::GetTotalStaticEfficiency(unsigned short inMarkerTP){return TotalStaticEfficiency[inMarkerTP];}

inline su2double CTurboMachinery::GetTotalTotalEfficiency(unsigned short inMarkerTP){return TotalTotalEfficiency[inMarkerTP];}

inline su2double CTurboMachinery::GetEulerianWork(unsigned short inMarkerTP){return EulerianWork[inMarkerTP];}

inline su2double CTurboMachinery::GetTotalEnthalpyIn(unsigned short inMarkerTP){return TotalEnthalpyIn[inMarkerTP];}

inline su2double CTurboMachinery::GetFlowAngleIn(unsigned short inMarkerTP){return FlowAngleIn[inMarkerTP];}

inline su2double CTurboMachinery::GetFlowAngleOut(unsigned short inMarkerTP){return FlowAngleOut[inMarkerTP];}

inline su2double CTurboMachinery::GetMassFlowIn(unsigned short inMarkerTP){return MassFlowIn[inMarkerTP];}

inline su2double CTurboMachinery::GetMassFlowOut(unsigned short inMarkerTP){return MassFlowOut[inMarkerTP];}

inline su2double* CTurboMachinery::GetMachIn(unsigned short inMarkerTP){return MachIn[inMarkerTP];}

inline su2double* CTurboMachinery::GetMachOut(unsigned short inMarkerTP){return MachOut[inMarkerTP];}

inline su2double CTurboMachinery::GetEnthalpyOut(unsigned short inMarkerTP){return EnthalpyOut[inMarkerTP];}

inline su2double CTurboMachinery::GetVelocityOutIs(unsigned short inMarkerTP){return VelocityOutIs[inMarkerTP];}

inline su2double CTurboMachinery::GetTotalPresureIn(unsigned short inMarkerTP){return TotalPresureIn[inMarkerTP];}

inline su2double CTurboMachinery::GetTotalTemperatureIn(unsigned short inMarkerTP){return TotalTemperatureIn[inMarkerTP];}

inline su2double CTurboMachinery::GetFlowAngleIn_BC(unsigned short inMarkerTP){return FlowAngleIn_BC[inMarkerTP];}

inline su2double CTurboMachinery::GetEntropyIn(unsigned short inMarkerTP){return EntropyIn[inMarkerTP];}

inline su2double CTurboMachinery::GetEntropyIn_BC(unsigned short inMarkerTP){return EntropyIn_BC[inMarkerTP];}

inline su2double CTurboMachinery::GetTotalEnthalpyIn_BC(unsigned short inMarkerTP){return TotalEnthalpyIn_BC[inMarkerTP];}

inline su2double CTurboMachinery::GetPressureRatio(unsigned short inMarkerTP){return PressureRatio[inMarkerTP];}

inline su2double CTurboMachinery::GetDensityIn(unsigned short inMarkerTP){return DensityIn[inMarkerTP];}

inline su2double CTurboMachinery::GetPressureIn(unsigned short inMarkerTP){return PressureIn[inMarkerTP];}

inline su2double* CTurboMachinery::GetTurboVelocityIn(unsigned short inMarkerTP){return TurboVelocityIn[inMarkerTP];}

inline su2double CTurboMachinery::GetDensityOut(unsigned short inMarkerTP){return DensityOut[inMarkerTP];}

inline su2double CTurboMachinery::GetPressureOut(unsigned short inMarkerTP){return PressureOut[inMarkerTP];}

inline su2double* CTurboMachinery::GetTurboVelocityOut(unsigned short inMarkerTP){return TurboVelocityOut[inMarkerTP];}

inline su2double CTurboMachinery::GetEnthalpyOutIs(unsigned short inMarkerTP){return EnthalpyOutIs[inMarkerTP];}

inline su2double CTurboMachinery::GetEntropyGen(unsigned short inMarkerTP){return EntropyGen[inMarkerTP];}

inline su2double CTurboMachinery::GetAbsFlowAngleIn(unsigned short inMarkerTP){return AbsFlowAngleIn[inMarkerTP];}

inline su2double CTurboMachinery::GetTotalEnthalpyOut(unsigned short inMarkerTP){return TotalEnthalpyOut[inMarkerTP];}

inline su2double CTurboMachinery::GetTotalEnthalpyOutIs(unsigned short inMarkerTP){return TotalEnthalpyOutIs[inMarkerTP];}

inline su2double CTurboMachinery::GetTotalRothalpyIn(unsigned short inMarkerTP){return TotalRothalpyIn[inMarkerTP];}

inline su2double CTurboMachinery::GetTotalRothalpyOut(unsigned short inMarkerTP){return TotalRothalpyOut[inMarkerTP];}

inline su2double CTurboMachinery::GetAbsFlowAngleOut(unsigned short inMarkerTP){return AbsFlowAngleOut[inMarkerTP];}

inline su2double CTurboMachinery::GetPressureOut_BC(unsigned short inMarkerTP){return PressureOut_BC[inMarkerTP];}

inline void CTurboMachinery::SetTotalPressureLoss(su2double value, unsigned short inMarkerTP){ TotalPressureLoss[inMarkerTP]=value;}

inline void CTurboMachinery::SetKineticEnergyLoss(su2double value, unsigned short inMarkerTP){ KineticEnergyLoss[inMarkerTP]=value;}

inline void CTurboMachinery::SetTotalStaticEfficiency(su2double value, unsigned short inMarkerTP){ TotalStaticEfficiency[inMarkerTP]=value;}

inline void CTurboMachinery::SetTotalTotalEfficiency(su2double value, unsigned short inMarkerTP){ TotalTotalEfficiency[inMarkerTP]=value;}

inline void CTurboMachinery::SetEulerianWork(su2double value, unsigned short inMarkerTP){ EulerianWork[inMarkerTP]=value;}

inline void CTurboMachinery::SetTotalEnthalpyIn(su2double value, unsigned short inMarkerTP){ TotalEnthalpyIn[inMarkerTP]=value;}

inline void CTurboMachinery::SetFlowAngleIn(su2double value, unsigned short inMarkerTP){ FlowAngleIn[inMarkerTP]=value;}

inline void CTurboMachinery::SetFlowAngleOut(su2double value, unsigned short inMarkerTP){ FlowAngleOut[inMarkerTP]=value;}

inline void CTurboMachinery::SetMassFlowIn(su2double value, unsigned short inMarkerTP){ MassFlowIn[inMarkerTP]=value;}

inline void CTurboMachinery::SetMassFlowOut(su2double value, unsigned short inMarkerTP){ MassFlowOut[inMarkerTP]=value;}
