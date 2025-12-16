/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    rocketMotor

Description
    Solver for simulating gas-particle flow in solid rocket motors with
    optional propellant surface regression.

    The solver treats the propellant, gas, and particle phases as three
    separate phases and models their coupled interaction through
    momentum, mass, and energy exchange.

    Can also be used without the propellant phase, making it applicable
    for steady flow simulations of two-phase (gas-particle) flows in 
    rocket motors and convergent-divergent nozzles. 

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "GeometricField.H"
#include "processorFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
  #include "postProcess.H"

  #include "addCheckCaseOptions.H"
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"
  #include "createControl.H"
  #include "createTimeControls.H"
  #include "createFields.H"

  #include "CourantNo.H"
  #include "setInitialDeltaT.H"
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  Info<< "\nStarting time loop\n" << endl;

  while (runTime.run())
  {
    #include "readTimeControls.H"

    #include "CourantNo.H"
    // #include "setDeltaTFactor.H"

    runTime++;
    Info<< "Time = " << runTime.timeName() << nl << endl;

    // Regression code goes here --------


    runTime.write();
    runTime.printExecutionTime(Info);
  }
  Info<< "End\n" << endl;

  return 0;
}


// ************************************************************************* //
