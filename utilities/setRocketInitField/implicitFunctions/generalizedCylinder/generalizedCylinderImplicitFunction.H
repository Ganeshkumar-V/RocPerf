/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2019-2020 DLR
    Copyright (C) 2025 Ganeshkumar V
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

Class
    Foam::implicitFunctions::generalizedCylinderImplicitFunction

Description
    creates a infintite long generalizedCylinderImplicitFunction

SourceFiles
    generalizedCylinderImplicitFunction.C

\*---------------------------------------------------------------------------*/

#ifndef implicitFunction_generalizedCylinderImplicitFunction_H
#define implicitFunction_generalizedCylinderImplicitFunction_H

#include "implicitFunction.H"
#include "point.H"
#include "tensor.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunctions
{

/*---------------------------------------------------------------------------*\
                   Class generalizedCylinderImplicitFunction Declaration
\*---------------------------------------------------------------------------*/

class generalizedCylinderImplicitFunction
:
    public implicitFunction
{
private:

    // Private Data

        //- Origin point
        const point origin_;

        //- Radius - scalar function
        autoPtr<Function1<scalar>> radius_;

        const vector direction_;

        const tensor project_;


public:

    //- Runtime type information
    TypeName("generalizedCylinder");


    // Constructors

        //- Construct from components
        generalizedCylinderImplicitFunction
        (
            const point& origin,
            const vector& direction
        );

        //- Construct from dictionary
        explicit generalizedCylinderImplicitFunction(const dictionary& dict);


    //- Destructor
    virtual ~generalizedCylinderImplicitFunction() = default;


    // Member Functions
    scalar radius(const vector& p) const
    {
        return radius_->value(p & direction_);
    }

    virtual scalar value(const vector& p) const
    {
        return -mag(project_ & (p - origin_)) + radius(p);
    }

    virtual vector grad(const vector& p) const
    {
        return -(project_ & (p - origin_));
    }

    virtual scalar distanceToSurfaces(const vector& p) const
    {
        return mag(mag(project_ & (p - origin_)) - radius(p));
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace implicitFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
