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

Class
    qmom

SourceFiles
    qmom.C

Description
	 solve population balance equation

References

    [1]
    [2]

Author
    Ehsan Askari
	 Sep. 2016

\*---------------------------------------------------------------------------*/

#ifndef qmom_H
#define qmom_H

#include "populationBalanceReturn.H"
#include "breakupLaw.H"
//Added by Ehsan
#include "diameterModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  class populationBalanceModel;

/*---------------------------------------------------------------------------*\
                         Class qmom Declaration
\*---------------------------------------------------------------------------*/

class qmom
:
    public populationBalanceReturn
{
    // Private data


			// Reference to populationBalance model
			populationBalanceModel& populationBalanceModel_;

			// volume fraction
			volScalarField alpha_;

			//- Menétrey-willam dictionary to define all constant
			dictionary dict_;

			// Sauter diameter d32
			volScalarField d32_;



                       volScalarField dsauter;


			//- mean diameter
			dimensionedScalar dm_;


                       // correction parameters
                       dimensionedScalar ResidualAlphaForDsauter_;
                       dimensionedScalar ResidualAlphaForAdjust_;
                       dimensionedScalar ResidualAlphaForCorrect_;

      dimensionedScalar m0_;
      dimensionedScalar m1_;
      dimensionedScalar m2_;
      dimensionedScalar m3_;
      dimensionedScalar m4_;
      dimensionedScalar m5_;
      dimensionedScalar w1_;
      dimensionedScalar w2_;
      dimensionedScalar w3_;
      dimensionedScalar L1_;
      dimensionedScalar L2_;
      dimensionedScalar L3_;

			//- number of classes  (is a label not scalar!!)
//			label Nc_;

//- number of moments  (is a label not scalar!!)
			label Nm_;

      //- number of nodes  (is a label not scalar!!)
      			label nodes_;

			//- coalescence: initial thickness
			dimensionedScalar h0_;

			//- coalescence: rupture thickness
			dimensionedScalar hf_;

			//-
			dimensionedScalar beta1_;

			//-
			dimensionedScalar beta2_;

			//-
			dimensionedScalar betaPB_;

			//-
			dimensionedScalar k1_;

			//- surface tension
			dimensionedScalar sigma_;

			//- breakup coefficient
			dimensionedScalar breakCoeff_;

			//- coalescence coefficient
			dimensionedScalar coalCoeff_;

			//- Ramkrishna model paramter: equivalent to V(i+1) = S*V(i) for geometric evolution
			dimensionedScalar S_;

			//-
			dimensionedScalar dMin_;
			dimensionedScalar dMax_;

			//- Max of iterations
			label maxIters_;

			//- tolerance
			scalar loopTolerance_;

            //           PtrList<volScalarField> f;
            /*    volScalarField m0;
                volScalarField m1;
                volScalarField m2;
                volScalarField m3;
                volScalarField m4;
                volScalarField m5;
                volScalarField w1;
                volScalarField w2;
                volScalarField w3;
                volScalarField L1;
                volScalarField L2;
                volScalarField L3;*/





			//- class method selector
			//word populationBalanceMethodName_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        qmom(const qmom&);

        //- Disallow default bitwise assignment
        void operator=(const qmom&);

public:



    //- Runtime type information
    TypeName("qmom");

    // Static data members


    // Constructors

        //- Construct from dictionary
        qmom
        (
            const word& name,
            populationBalanceModel& populationBalanceModel
        );


    // Destructor

        virtual ~qmom();





    // Member Functions
/////////////////////////////////////////////////////////////////////
    //    virtual tmp<volScalarField> dSauter() const;
//////////////////////////////////////////////////////////////////////

        //- Correct the rheological model
        virtual void correct();

		  scalar funcTest(scalar&, scalar&, scalar&) const;

		  scalar kronecker(label&, label&) const;

		  //------------------------------ bubbles breakup kernel -----------------------//
/*		  virtual void breakupKernel
		  (
				PtrList<volScalarField>&,
				PtrList<volScalarField>&,
				const volScalarField&
		  ) const;

		  //------------------------------ bubbles coalescence kernel -------------------//
		  virtual void coalescenceKernel
		  (
				PtrList<volScalarField>&,
				PtrList<volScalarField>&,
				const volScalarField&,
				const volScalarField&
		  ) const;
*/

        //------------------------------ PD algorithm -------------------//
        virtual void weightsAbscissas
        (
          PtrList<volScalarField>&,
          PtrList<volScalarField>&,
          PtrList<volScalarField>&,
          const volScalarField
        ) const;

	     //------------------------------ Sauter diameter function ---------------------//
	     tmp<volScalarField> SauterDiameter(PtrList<volScalarField>&, const volScalarField&) const ;



	     //---
	     void adjust
			(
				PtrList<volScalarField>&,
        PtrList<volScalarField>&,
				const volScalarField&
			)
			const;


        //- Return slave friction value fraction
        virtual void updateSauterDiameter();


        //- Write law dictionary
        virtual void writeDict(Ostream& os) const {};
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
