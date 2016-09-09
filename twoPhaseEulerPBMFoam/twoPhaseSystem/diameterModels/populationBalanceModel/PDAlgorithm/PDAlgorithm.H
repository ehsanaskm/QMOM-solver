/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2004 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    PDAlgorithm 
    rachid Bannari: Class method CM
    brahim SELMA: QMOM and DQMOM

SourceFiles
    PDAlgorithm.C

\*---------------------------------------------------------------------------*/

#ifndef PDAlgorithm_H
#define PDAlgorithm_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class PDAlgorithm
{
        // Private data
        dictionary dict_;

        //- Method name selection
        word name_;

        //- Mean diameter
        dimensionedScalar dmoy_;

        //- Quadrature nodes number (max 3)
        dimensionedScalar noeuds_;

		

public:

    // Constructors

        PDAlgorithm
        (
            const fvMesh& mesh,
            const dictionary& populationBlanceProperties,
            const word& ModelName
        );


    // Selectors

        //- Return a reference to the selected population  model
        static autoPtr<PDAlgorithm> New
        (
            const fvMesh& mesh,
            const dictionary& populationBlanceProperties,
            const word& ModelName
        );


    // Destructor

        ~PDAlgorithm();


    // Member Functions

        const word& name() const
        {
            return name_;
        }
        const dimensionedScalar& dmoy() const
        {
            return dmoy_;
        }   
        const dimensionedScalar& noeuds() const
        {
            return noeuds_;
        }
// QMOM
void WeightsAbscissas
( 
     PtrList<volScalarField>&,  
     PtrList<volScalarField>&, 
     PtrList<volScalarField>&,  
     const volScalarField 
) 
const;       


}; //End of Class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
