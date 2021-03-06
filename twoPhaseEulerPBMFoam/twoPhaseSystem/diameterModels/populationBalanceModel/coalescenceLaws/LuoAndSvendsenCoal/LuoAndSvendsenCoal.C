/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "LuoAndSvendsenCoal.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "incompleteGammaFunction.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LuoAndSvendsenCoal, 0);
    addToRunTimeSelectionTable(coalescenceLaw, LuoAndSvendsenCoal, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::LuoAndSvendsenCoal::LuoAndSvendsenCoal
(
    const word& name,
    const volScalarField& alpha,
    const dictionary& dict
)
:
    coalescenceLaw(name, alpha, dict),
    dictLuoAndSvendsenCoal_
    (
        dict.subDict("LuoAndSvendsenCoeffs")
    ),
    
    // Luo and Svendsen model parameters
    rhoa_(dictLuoAndSvendsenCoal_.lookup("rhoa")),
    rhob_(dictLuoAndSvendsenCoal_.lookup("rhob")),
    sigma_(dictLuoAndSvendsenCoal_.lookup("sigma"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LuoAndSvendsenCoal::~LuoAndSvendsenCoal()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::LuoAndSvendsenCoal::rhoa() const
{
    return rhoa_; 
}

Foam::dimensionedScalar Foam::LuoAndSvendsenCoal::rhob() const
{
    return rhob_;
}

Foam::dimensionedScalar Foam::LuoAndSvendsenCoal::coalescenceRate(scalar& di, scalar& dj, volScalarField& epsf) const
{

	scalar omega_coal, xi, p_coal, weber, num, den, u_i, u_j, u_i_j,coalRate;
	xi=di/(dj+1e-15);
       dimensionedScalar result("coalRate",dimMass/pow3(dimLength)/dimTime,0.0); 
	u_i = 1.43*pow(max(epsf).value()*di, 1./3);

	u_j = 1.43*pow(max(epsf).value()*dj, 1./3);
	u_i_j = sqrt(pow(u_i, 2.)+pow(u_j, 2.));
	weber = rhob_.value()*di*pow(u_i_j,2.0)/sigma_.value();
	num = pow((0.75*(1.+pow(xi,2.))*(1.+pow(xi,3.))), 0.5);
   den = pow((0.5+(rhoa_/rhob_).value()),0.5)*pow((1.+xi),3.);
   p_coal = exp(-(num/den)*pow(weber, 0.5));
 
   omega_coal = (M_PI/4.)*pow((di+dj),2.)*u_i_j;
   coalRate = omega_coal*p_coal;

 
  
   if(coalRate < 0.0) 
   {
	  result = dimensionedScalar
				  (
						"zero",
						dimMass/pow3(dimLength)/dimTime,
						0.0
				  );
   } 
   else
   {
	  result = dimensionedScalar
				  (
						"value",
						dimMass/pow3(dimLength)/dimTime,
						coalRate
				  );
   }

   return result;

}

// ************************************************************************* //
