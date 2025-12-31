{
    Info << "RegressFluid: Solving for liquid volume fraction" << endl;
    
    const scalarField& V(mesh.V());
    const scalar dt(runTime.deltaT().value());

    // Calculate alpha source
    alphaSource.primitiveFieldRef() = rb*(mag(surf.normal().internalField())/mesh.V());
    volScalarField distributedAlphaSource(0*alphaSource);
    
    forAll(alphaSource, celli)
    {
        if (alphaSource[celli]*dt > alphaOld[celli])
        {
            label adjacent = findAdjacent(celli, alphaOld);
            if (adjacent == -1)
            {
                if (meshObject::debug)
                {
                    Pout << " Adjacent cell not found in present processor to pass on alpha source. \n "
                     << " Celli: " << celli << " \n "
                     << " alpha: " << alphaOld[celli] << "\n "
                     << " alphaSource*dt: " << alphaSource[celli]*dt << endl;
                }
                
                distributedAlphaSource[celli] += alphaSource[celli];
            }
            else
            {
                distributedAlphaSource[adjacent] = (alphaSource[celli]*dt - alphaOld[celli])*V[celli]/(V[adjacent]*dt);
                distributedAlphaSource[celli] = alphaOld[celli]/dt;
                if (meshObject::debug)
                {
                    Pout << "Moving the source from " << celli << " to " << adjacent << endl;
                    Pout << " Celli: " << celli << " \n "
                    << " alpha: " << alphaOld[celli] << "\n "
                    << " alphaSource*dt: " << alphaSource[celli]*dt << endl;
                    Pout << "alphaSource[celli]: " << distributedAlphaSource[celli]*dt << " alphaSource[adjacent]: " << distributedAlphaSource[adjacent]*dt << endl;
                }
                Pout << "~~~~~~~~~~~~ CELL TRANSITION ~~~~~~~~~~~~~~ " << endl;
                
                if (distributedAlphaSource[adjacent]*dt > alphaOld[adjacent])
                {
                    if (meshObject::debug)
                    {
                        Warning << " cell: " << adjacent << " is going to become negative. "
                            << " alphaSource*dt: " << distributedAlphaSource[adjacent]*dt
                            << " alpha: " << alpha[adjacent]
                            << " ----- Please Reduce the time step!" << exit(FatalError);
                    }
                    Pout << "Time step is too large! " << " Please reduce the time step. " << exit(FatalError);
                }
            }
    
         }
        else
        {
            distributedAlphaSource[celli] += alphaSource[celli];
    
            // check for negative alpha in the adjacent cell
            if (distributedAlphaSource[celli]*dt > alphaOld[celli])
            {
                
                if (meshObject::debug)
                {
                    Warning << " cell: " << celli << " is going to become negative. "
                        << " alphaSource: " << distributedAlphaSource[celli]
                        << " alpha: " << alpha[celli] 
                        << " ----- Please Reduce the time step!" << exit(FatalError);
                }
                Pout << "Time step is too large! " << " Please reduce the time step. " << exit(FatalError);
            }
        }
    }
    
    // For processor boundary patches
    distributedAlphaSource.correctBoundaryConditions(); // Send/Recv processor boundary information
    
    {
        const volScalarField::Boundary& alphaOldbF(alphaOld.boundaryField());
        const volScalarField::Boundary& VbF(cellV.boundaryField());
        volScalarField::Boundary& distributedAlphaSourcebF(distributedAlphaSource.boundaryFieldRef());
        forAll(distributedAlphaSourcebF, bFi)
        {
            word BCtype = mesh.boundaryMesh().types()[bFi];
            if ( BCtype == "processor" )    // For Processor Boundary
            {
                const processorPolyPatch& pp = refCast<const processorPolyPatch>( mesh.boundaryMesh()[bFi]);
                scalarField VNF(VbF[bFi].patchNeighbourField());
                const scalarField alphaOldNF(alphaOldbF[bFi].patchNeighbourField());
                const UList<label>& fC(mesh.boundaryMesh()[bFi].faceCells());
    
                scalarField distributedAlphaSourceNF(distributedAlphaSourcebF[bFi].patchNeighbourField());
                forAll(fC, i)
                {
                    // Receiving patch
                    if (distributedAlphaSourceNF[i]*dt > alphaOldNF[i])
                    {
                        if ((alphaOldNF[i] < 1.0) && (alphaOld[fC[i]] == 1.0))
                        {
                            distributedAlphaSource[fC[i]] = (distributedAlphaSourceNF[i]*dt - alphaOldNF[i])*VNF[i]/(V[fC[i]]*dt);
                            if (meshObject::debug)
                            {
                                Pout << "Receiving source from: " << pp.neighbProcNo() << endl;
                                Pout << "AlphaSource-Neighbour[celli]: " << fC[i] << " -> "<< distributedAlphaSourceNF[i] << endl;
                                Pout << "New Alphasource: " << distributedAlphaSource[fC[i]] << " added to " << fC[i] << " alpha here is: " << alphaOld[fC[i]] << endl;
                            }
                            Pout << "\n~~~~~~~~~~~~ INTER PROCESSOR -> CELL TRANSITION ~~~~~~~~~~~~~~ " << endl;
    
                            if (distributedAlphaSource[fC[i]]*dt > alphaOld[fC[i]])
                            {
                                Pout << "Time step is too large! " << " Please reduce the time step. " << exit(FatalError);
                            }
                        }
                    }
    
                    // Sending Patch
                    if (distributedAlphaSource[fC[i]]*dt > alphaOld[fC[i]])
                    {
                        if (alphaOldNF[i] == 1.0)
                        {
                            if (meshObject::debug)
                            {
                                Pout << "Sending source from: " << pp.myProcNo() << endl;
                                Pout << "AlphaSource[i]: " << i << " -> "<< distributedAlphaSource[fC[i]] << endl;
                            }
                            distributedAlphaSource[fC[i]] = alphaOld[fC[i]]/dt;
                        }
    
                    }
                }
            }
        }
    }
    distributedAlphaSource.correctBoundaryConditions();
    
    // For doamin boundary patches
    {
        const volScalarField::Boundary& alphaOldbF(alphaOld.boundaryField());
        const volScalarField::Boundary& VbF(cellV.boundaryField());
        volScalarField::Boundary& distributedAlphaSourcebF(distributedAlphaSource.boundaryFieldRef());
        forAll(distributedAlphaSourcebF, bFi)
        {
            word BCtype = mesh.boundaryMesh().names()[bFi];
            if ( BCtype == bottomPatchName )    // For Non-processor Boundary
            {
                const UList<label>& fC(mesh.boundaryMesh()[bFi].faceCells());
                forAll(fC, i)
                {
                    if (distributedAlphaSource[fC[i]]*dt > alphaOld[fC[i]])
                    {
                        // Info Statements
                        Info << " Boundary patch name: " << BCtype << endl;
                        Pout << "Interface is at the domain boundary! So emptying the cell! " << endl;
                        Pout << "Celli: " << fC[i] << endl << "alphaSource*dt: " << distributedAlphaSource[fC[i]]*dt << endl;
                        distributedAlphaSource[fC[i]] = alphaOld[fC[i]]/dt;
                    }
                }
            }
        }
        
    }
    distributedAlphaSource.correctBoundaryConditions();
    
    forAll(alpha, celli)
    {
        alpha[celli] = alphaOld[celli] - distributedAlphaSource[celli]*dt;
        if (mag(alpha[celli]) < SMALL)
        {
            alpha[celli] = 0.0;
        }
    }
    
    alpha.correctBoundaryConditions();
    distributedAlphaSource.correctBoundaryConditions();
    
    Info<< alpha.name() << " Propellant volume fraction = "
                    << alpha.weightedAverage(mesh.V()).value()
                    << "  Min(alpha) = " << min(alpha).value()
                    << "  Max(alpha) = " << max(alpha).value()
                    << endl;
}
    