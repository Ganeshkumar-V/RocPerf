if (adjustTimeStep)
{
    scalar factor = 1.2;
    scalar translimit = pimple.dict().get<scalar>("translimit");
    {
        // Propellant burning surface transition - Stabilization
        const scalarField iF(phases[propellantIndex].internalField());
        const tmp<scalarField> tinterface((iF - 0.001)*(0.999 - iF));
        const scalarField& interface(tinterface());

        if (gSum(interface*pos0(interface))/surfaces > translimit)
        {
            maxCo = pimple.dict().get<scalar>("maxCo");
        }
        else
        {
            maxCo = pimple.dict().get<scalar>("transCo");
        }
    }
    scalar maxDeltaTFact = maxCo/(CoNum + SMALL);
    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), factor);

    runTime.setDeltaT
    (
        min
        (
            deltaTFact*runTime.deltaTValue(),
            maxDeltaT
        )
    );

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}
