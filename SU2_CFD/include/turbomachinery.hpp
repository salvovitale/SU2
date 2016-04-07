/*!
 * \file turbomachinery.hpp
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

#pragma once
#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../include/solver_structure.hpp"
/*!
 * \class CTurboMachinery
 * \brief Main class for Turbomachinery capabilities of the tool.
 * a child class for each particular tubomachinery case
 * \author: S.Vitale, O. Pfeifle
 * \version 4.1.0 "Cardinal"
 */

class CTurboMachinery {
protected:

	unsigned short nDim;

	su2double *TotalStaticEfficiency,
	  	  	  	*TotalTotalEfficiency,
				*KineticEnergyLoss,
				*TotalPressureLoss,
	  	  *MassFlowIn,
				*MassFlowOut,
				*FlowAngleIn,
				*FlowAngleIn_BC,
				*FlowAngleOut,
				*EulerianWork,
				*TotalEnthalpyIn,
				*TotalEnthalpyIn_BC,
				*EntropyIn,
				*EntropyIn_BC,
				*PressureRatio,
				*TotalPresureIn,
				*TotalTemperatureIn,
				*EnthalpyOut,
				**MachIn,
				**MachOut,
				*VelocityOutIs,
				*DensityIn,
				*PressureIn,
				**TurboVelocityIn,
				*DensityOut,
				*PressureOut,
				**TurboVelocityOut,
				*EnthalpyOutIs,
				*EntropyGen,
				*AbsFlowAngleIn,
				*TotalEnthalpyOut,
				*TotalRothalpyIn,
				*TotalRothalpyOut,
				*TotalEnthalpyOutIs,
				*AbsFlowAngleOut,
				*PressureOut_BC;


public:



	 /*!
	 	 * \brief Constructor of the class.
	 	 */
	 	CTurboMachinery(void);

	 	/*!
	 	 * \overload
	 	 * \param[in] geometry - Geometrical definition of the problem.
	 	 * \param[in] config - Definition of the particular problem.
	 	 */
	 	CTurboMachinery(CConfig *config, CGeometry *geometry);

	 	/*!
	 	 * \brief Destructor of the class.
	 	 */
	 	virtual ~CTurboMachinery(void);

	 	/*!
		 * \brief Compute turboperformance.
		 * \param[in] solver 		- solver container.
		 * \param[in] config 		- config file information.
		 * \param[in] geometry  - geometry information.
		 */
	 	void TurboPerformance(CSolver *solver, CConfig *config, CGeometry *geometry);
	 	/*!
		 * \brief Compute turboperformance.
		 * \param[in] solver 		- solver container.
		 * \param[in] config 		- config file information.
		 */

	 	void MultiStageTurboPerformance(CSolver *solver, CConfig *config);
		/*!
		 * \brief Provide Total Pressure Losses (turbomachinery performance).
		 * \param[in] inMarkerTP - turboperformance marker.
		 * \return Value of Total Pressure Losses.
		 */
		su2double GetTotalPressureLoss(unsigned short inMarkerTP);

		/*!
		 * \brief Provide Kinetic Energy Losses (turbomachinery performance).
		 * \param[in] inMarkerTP - turboperformance marker.
		 * \return Value of the Kinetic Energy Losses.
		 */
		su2double GetKineticEnergyLoss(unsigned short inMarkerTP);

		/*!
		 * \brief Provide Total-Total Efficiency (turbomachinery performance).
		 * \param[in] inMarkerTP - turboperformance marker.
		 * \return Value of the Total-Total Efficiency.
		 */
		su2double GetTotalTotalEfficiency(unsigned short inMarkerTP);

		/*!
		 * \brief Provide Total-Static Efficiency (turbomachinery performance).
		 * \param[in] inMarkerTP - turboperformance marker.
		 * \return Value of the Total-Static Efficiency.
		 */
		su2double GetTotalStaticEfficiency(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the Eulerian Work (turbomachinery performance).
		 * \param[in] inMarkerTP - turboperformance marker.
		 * \return Value of the Eulerian Work.
		 */
		su2double GetEulerianWork(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the Inlet Total Enthalpy (turbomachinery performance).
		 * \param[in] inMarkerTP - turboperformance marker.
		 * \return Value of the Inlet Total Enthalpy.
		 */
		su2double GetTotalEnthalpyIn(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the Inlet Flow Angle (turbomachinery performance).
		 * \param[in] inMarkerTP - turboperformance marker.
		 * \return Value of the Inlet Flow Angle.
		 */
		su2double GetFlowAngleIn(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the Outlet Flow Angle (turbomachinery performance).
		 * \param[in] inMarkerTP - turboperformance marker.
		 * \return Value of the Outlet FLow Angle.
		 */
		su2double GetFlowAngleOut(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the Inlet Mass Flow (turbomachinery performance).
		 * \param[in] inMarkerTP - turboperformance marker.
		 * \return Value of the Inlet Mass Flow.
		 */
		su2double GetMassFlowIn(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the Outlet Mass Flow (turbomachinery performance).
		 * \param[in] inMarkerTP - turboperformance marker.
		 * \return Value of the Outlet Mass Flow.
		 */
		su2double GetMassFlowOut(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the Inlet Mach number (turbomachinery performance).
		 * \param[in] inMarkerTP - turboperformance marker.
		 * \return Value of the Inlet Mach number.
		 */
		su2double* GetMachIn(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the Outlet Mach number (turbomachinery performance).
		 * \param[in] inMarkerTP - turboperformance marker.
		 * \return Value of the Outlet Mach number.
		 */
		su2double* GetMachOut(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the Outlet Enthalpy (turbomachinery performance).
		 * \param[in] inMarkerTP - turboperformance marker.
		 * \return Value of the Outlet Enthalpy.
		 */
		su2double GetEnthalpyOut(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the Isentropic Outlet Velocity (turbomachinery performance).
		 * \param[in] inMarkerTP - turboperformance marker.
		 * \return Value of the Isentropic Outlet Velocity.
		 */
		su2double GetVelocityOutIs(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the total pressure at the inlet for convergence monitoring.
		 * \param[in] val_marker - bound marker.
		 * \return Value of the Outlet Inlet total pressure.
		 */
		su2double GetTotalPresureIn(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the total temperature at the inlet for convergence monitoring.
		 * \param[in] val_marker - bound marker.
		 * \return Value of the inlet total temperature.
		 */
		su2double GetTotalTemperatureIn(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the inlet flow angle from BC for convergence monitoring.
		 * \param[in] val_marker - bound marker.
		 * \return Value of the inlet flow angle from BC.
		 */
		su2double GetFlowAngleIn_BC(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the inlet entropy for convergence monitoring.
		 * \param[in] val_marker - bound marker.
		 * \return Value of the inlet entropy.
		 */
		su2double GetEntropyIn(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the inlet entropy from BC for convergence monitoring.
		 * \param[in] val_marker - bound marker.
		 * \return Value of the inlet entropy from BC.
		 */
		su2double GetEntropyIn_BC(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the total inlet enthalpy from BC for convergence monitoring.
		 * \param[in] val_marker - bound marker.
		 * \return Value of the total enthalpy from BC.
		 */
		su2double GetTotalEnthalpyIn_BC(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the Pressure Ratio (tubomachinery performance).
		 * \param[in] inMarkerTP - turboperformance marker.
		 * \return Value of the Pressure Ratio.
		 */
		su2double GetPressureRatio(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the inlet density to check convergence of conservative mixing-plane.
		 * \param[in] inMarkerTP - bound marker.
		 * \return Value of the inlet density.
		 */
		 su2double GetDensityIn(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the inlet pressure to check convergence of conservative mixing-plane.
		 * \param[in] inMarkerTP - bound marker.
		 * \return Value of inlet pressure.
		 */
		su2double GetPressureIn(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the inlet normal velocity to check convergence of conservative mixing-plane.
		 * \param[in] inMarkerTP - bound marker.
		 * \return Value of the inlet normal velocity.
		 */
		su2double* GetTurboVelocityIn(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the outlet density to check convergence of conservative mixing-plane.
		 * \param[in] inMarkerTP - bound marker.
		 * \return Value of the outlet density.
		 */
		su2double GetDensityOut(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the outlet pressure to check convergence of conservative mixing-plane.
		 * \param[in] inMarkerTP - bound marker.
		 * \return Value of the outlet pressure.
		 */
		su2double GetPressureOut(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the outlet normal velocity to check convergence of conservative mixing-plane.
		 * \param[in] inMarkerTP - bound marker.
		 * \return Value of the outlet normal velocity.
		 */
		su2double* GetTurboVelocityOut(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the isentropic outlet Enthalpy for turbomachinery performance.
		 * \param[in] inMarkerTP - bound marker.
		 * \return Value of the isentropic outlet Enthalpy.
		 */
		su2double GetEnthalpyOutIs(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the Entropy generated for turbomachinery performance.
		 * \param[in] inMarkerTP - bound marker.
		 * \return Value of the Entropy generated.
		 */
		su2double GetEntropyGen(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the absolute inlet flow angle to check convergence of BC and turbomachinery performance.
		 * \param[in] inMarkerTP - bound marker.
		 * \return Value of the absolute inlet flow angle.
		 */
		su2double GetAbsFlowAngleIn(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the Total outlet Enthalpy for turbomachinery performance.
		 * \param[in] inMarkerTP - bound marker.
		 * \return Value of the Total outlet Enthalpy.
		 */
		su2double GetTotalEnthalpyOut(unsigned short inMarkerTP);

		/*!
		 * \brief Provide Total outlet isentropic Enthalpy for turbomachinery performance.
		 * \param[in] inMarkerTP - bound marker.
		 * \return Value of the Total outlet isentropic Enthalpy.
		 */
		su2double GetTotalEnthalpyOutIs(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the Total outlet Enthalpy for turbomachinery performance.
		 * \param[in] inMarkerTP - bound marker.
		 * \return Value of the Total outlet Enthalpy.
		 */
		su2double GetTotalRothalpyIn(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the Total outlet Enthalpy for turbomachinery performance.
		 * \param[in] inMarkerTP - bound marker.
		 * \return Value of the Total outlet Enthalpy.
		 */
		su2double GetTotalRothalpyOut(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the absolute outlet flow angle for turbomachinery performance.
		 * \param[in] inMarkerTP - bound marker.
		 * \return Value of the absolute outlet flow angle.
		 */
		su2double GetAbsFlowAngleOut(unsigned short inMarkerTP);

		/*!
		 * \brief Provide the outlet pressure from BC for convergence monitoring.
		 * \param[in] inMarkerTP - bound marker.
		 * \return Value of the Outlet Pressure.
		 */
		su2double GetPressureOut_BC(unsigned short inMarkerTP);

		/*!
		 * \brief Set Total Pressure Losses (turbomachinery performance).
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetTotalPressureLoss(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set Kinetic Energy Losses (turbomachinery performance).
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetKineticEnergyLoss(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set Total-Total Efficiency (turbomachinery performance).
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetTotalTotalEfficiency(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set Total-Static Efficiency (turbomachinery performance).
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetTotalStaticEfficiency(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set the Eulerian Work (turbomachinery performance).
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetEulerianWork(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set the Inlet Total Enthalpy (turbomachinery performance).
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetTotalEnthalpyIn(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set the Inlet Flow Angle (turbomachinery performance).
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetFlowAngleIn(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set the Outlet Flow Angle (turbomachinery performance).
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetFlowAngleOut(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set the Inlet Mass Flow (turbomachinery performance).
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetMassFlowIn(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set the Outlet Mass Flow (turbomachinery performance).
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetMassFlowOut(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set the Inlet Mach number (turbomachinery performance).
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetMachIn(su2double* value, unsigned short inMarkerTP);

		/*!
		 * \brief Set the Outlet Mach number (turbomachinery performance).
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetMachOut(su2double* value, unsigned short inMarkerTP);

		/*!
		 * \brief Set the Outlet Enthalpy (turbomachinery performance).
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetEnthalpyOut(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set the Isentropic Outlet Velocity (turbomachinery performance).
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetVelocityOutIs(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set the total pressure at the inlet for convergence monitoring.
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetTotalPresureIn(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set the total temperature at the inlet for convergence monitoring.
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetTotalTemperatureIn(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set the inlet flow angle from BC for convergence monitoring.
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetFlowAngleIn_BC(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set the inlet entropy for convergence monitoring.
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetEntropyIn(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set the inlet entropy from BC for convergence monitoring.
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetEntropyIn_BC(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set the total inlet enthalpy from BC for convergence monitoring.
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetTotalEnthalpyIn_BC(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set Pressure Ratio (tubomachinery performance).
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetPressureRatio(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set inlet density.
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetDensityIn(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set inlet pressure.
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetPressureIn(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set inlet normal velocity.
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetTurboVelocityIn(su2double* value, unsigned short inMarkerTP);

		/*!
		 * \brief Set outlet density.
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetDensityOut(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set outlet pressure.
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetPressureOut(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set outlet normal velocity.
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetTurboVelocityOut(su2double* value, unsigned short inMarkerTP);

		/*!
		 * \brief Set outlet isentropic Enthalpy.
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetEnthalpyOutIs(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set entropy generated.
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetEntropyGen(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set abslote inlet flow angle.
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetAbsFlowAngleIn(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set total outlet enthalpy.
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetTotalEnthalpyOut(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set total inlet rothalpy.
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetTotalRothalpyIn(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set total outlet rothalpy.
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetTotalRothalpyOut(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set total isentropic outlet enthalpy.
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetTotalEnthalpyOutIs(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set outlet flow angle.
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetAbsFlowAngleOut(su2double value, unsigned short inMarkerTP);

		/*!
		 * \brief Set outlet BC pressure.
		 * \param[in] value      - turboperformance value to set.
		 * \param[in] inMarkerTP - turboperformance marker.
		 */
		void SetPressureOut_BC(su2double value, unsigned short inMarkerTP);

};

#include "turbomachinery.inl"
