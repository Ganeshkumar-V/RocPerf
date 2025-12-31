#include "createRDeltaT.H"
#include "readGravitationalAcceleration.H"
#include "readhRef.H"

Info<< "Reading propellant field\n" << endl;
volScalarField alpha
(
    IOobject
    (
        "alpha.propellant",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField alphaOld
(
    IOobject
    (
        "alpha.propellant",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha
);

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("", dimVelocity, vector(0, 0, 0))
);

surfaceScalarField phi(fvc::flux(U));

volScalarField alphaSource
(
    IOobject
    (
        "alphaSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("", dimVelocity/dimLength, 0.0)
);

//- Pointer to reconstruction scheme
autoPtr<reconstructionSchemes> surf_
(
    reconstructionSchemes::New
    (
        alpha,
        phi,
        U,
        mesh.solverDict(alpha.name())
    )
);
reconstructionSchemes& surf(surf_());

dimensionedScalar rb("rb", dimVelocity, mesh.solverDict(alpha.name()).get<scalar>("rb"));

volScalarField cellV
(
    IOobject
    (
        "cellVolume",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("", dimVolume, 0.0)
);
cellV.primitiveFieldRef() = mesh.V();
cellV.correctBoundaryConditions();

word bottomPatchName("Wall");
