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

sourceFiles
    qmom.H

Author
    Ehsan Askari, M.Sc
    ehsan.askari@usherbrooke.ca

\*---------------------------------------------------------------------------*/

#include "qmom.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "populationBalanceModel.H"
#include "bound.H"
#include "fvm.H"
//Added by Ehsan
#include "twoPhaseSystem.H"
#include "scalar.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(qmom, 0);
    addToRunTimeSelectionTable(populationBalanceReturn, qmom, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
qmom::qmom
(
	const word& name,
	populationBalanceModel& populationBalanceModel
)
:
	populationBalanceReturn(name, populationBalanceModel),
	populationBalanceModel_(populationBalanceModel),
	alpha_(populationBalanceModel.alpha()),

	dict_
    (
        populationBalanceModel_.subDict("qmomCoeffs")
    ),

	dsauter
	(
		IOobject
		(
			"dsauter",
			alpha_.time().timeName(),
			alpha_.db(),
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		alpha_.mesh()
	),

	// Sauter diameter d32
	d32_
	(
		IOobject
		(
			"d32",
			alpha_.time().timeName(),
			alpha_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
         	alpha_.mesh(),
		dimensionedScalar("SMALL", dimLength, SMALL)
	),

  /////////////////////////////////////////////////////////////////
  	m0
  	(
  		IOobject
  		(
  			"m0",
  			alpha_.time().timeName(),
  			alpha_.db(),
  			IOobject::MUST_READ,
  			IOobject::AUTO_WRITE
  		),
  		alpha_.mesh()
  	),

  	m1
  	(
  		IOobject
  		(
  			"m1",
  			alpha_.time().timeName(),
  			alpha_.db(),
  			IOobject::MUST_READ,
  			IOobject::AUTO_WRITE
  		),
  		alpha_.mesh()
  	),
  	m2
  	(
  		IOobject
  		(
  			"m2",
  			alpha_.time().timeName(),
  			alpha_.db(),
  			IOobject::MUST_READ,
  			IOobject::AUTO_WRITE
  		),
  		alpha_.mesh()
  	),
  	m3
  	(
  		IOobject
  		(
  			"m3",
  			alpha_.time().timeName(),
  			alpha_.db(),
  			IOobject::MUST_READ,
  			IOobject::AUTO_WRITE
  		),
  		alpha_.mesh()
  	),
  	m4
  	(
  		IOobject
  		(
  			"m4",
  			alpha_.time().timeName(),
  			alpha_.db(),
  			IOobject::MUST_READ,
  			IOobject::AUTO_WRITE
  		),
  		alpha_.mesh()
  	),
  	m5
  	(
  		IOobject
  		(
  			"m5",
  			alpha_.time().timeName(),
  			alpha_.db(),
  			IOobject::MUST_READ,
  			IOobject::AUTO_WRITE
  		),
  		alpha_.mesh()
  	),
    w0
  	(
  		IOobject
  		(
  			"w0",
  			alpha_.time().timeName(),
  			alpha_.db(),
  			IOobject::MUST_READ,
  			IOobject::AUTO_WRITE
  		),
  		alpha_.mesh()
  	),
  	w1
  	(
  		IOobject
  		(
  			"w1",
  			alpha_.time().timeName(),
  			alpha_.db(),
  			IOobject::MUST_READ,
  			IOobject::AUTO_WRITE
  		),
  		alpha_.mesh()
  	),
  	w2
  	(
  		IOobject
  		(
  			"w2",
  			alpha_.time().timeName(),
  			alpha_.db(),
  			IOobject::MUST_READ,
  			IOobject::AUTO_WRITE
  		),
  		alpha_.mesh()
  	),
    L0
  	(
  		IOobject
  		(
  			"L0",
  			alpha_.time().timeName(),
  			alpha_.db(),
  			IOobject::MUST_READ,
  			IOobject::AUTO_WRITE
  		),
  		alpha_.mesh()
  	),
  	L1
  	(
  		IOobject
  		(
  			"L1",
  			alpha_.time().timeName(),
  			alpha_.db(),
  			IOobject::MUST_READ,
  			IOobject::AUTO_WRITE
  		),
  		alpha_.mesh()
  	),
  	L2
  	(
  		IOobject
  		(
  			"L2",
  			alpha_.time().timeName(),
  			alpha_.db(),
  			IOobject::MUST_READ,
  			IOobject::AUTO_WRITE
  		),
  		alpha_.mesh()
  	),

	ResidualAlphaForDsauter_(dict_.lookup("ResidualAlphaForDsauter")),
  ResidualAlphaForAdjust_(dict_.lookup("ResidualAlphaForAdjust")),
  ResidualAlphaForCorrect_(dict_.lookup("ResidualAlphaForCorrect")),
  dm_(dict_.lookup("dm")),
	Nm_(readLabel(dict_.lookup("Nm"))),
  nodes_(readLabel(dict_.lookup("nodes"))),
	h0_(dict_.lookup("h0")),
	hf_(dict_.lookup("hf")),
	beta1_(dict_.lookup("beta1")),
	beta2_(dict_.lookup("beta2")),
	betaPB_(dict_.lookup("betaPB")),
	k1_(dict_.lookup("k1")),
	sigma_(dict_.lookup("sigma")),
	breakCoeff_(dict_.lookup("breakCoeff")),
	coalCoeff_(dict_.lookup("coalCoeff")),
	S_(dict_.lookup("S")),
	dMin_(dict_.lookup("dMin")),
	dMax_(dict_.lookup("dMax")),
	maxIters_(readLabel(dict_.lookup("maxIters"))),
	loopTolerance_(readScalar(dict_.lookup("loopTolerance")))



{
	Info << "Selecting method for population balance model" << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

qmom::~qmom()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar qmom::funcTest(scalar& x, scalar& y,	scalar& z) const
{
	 if ( y >= x && y <= z)
    {
      return 1.0;
    }
	 else
    {
      return 0.0;
    }
}

scalar  qmom::kronecker(label& i, label& j) const
{
  if (i == j) return 1.0;

  else return 0.0;
}
/*
void qmom::breakupKernel(PtrList<volScalarField>& S, PtrList<volScalarField>& f, const volScalarField& alpha) const
{
	scalar d_i, x_i, d_i1, x_i1;
   scalar  d_j, x_j, d_i_1, x_i_1;
   scalar *diam = new scalar[Nm_+ 1];
   scalar fbv25 = 0.25, fbv50 = 0.50, fbv75 = 0.75;
  diam[0] = dm_.value()/10.0;

   for(int II=1;II<=Nm_;II++)
    {
      diam[II]= Foam::pow(S_.value(), scalar(II-Nm_/2)/3)*dm_.value();
    }

   volScalarField vf = 1.-alpha;

	for (int i = 1; i<= Nm_-1; i++)
	{
			d_i = diam[i];
			d_i_1 = diam[i-1];
			x_i =  M_PI*pow(d_i, 3.0)/6.0;
			x_i_1 =  M_PI*pow(d_i_1, 3.0)/6.0;

			for(int j = i+1 ; j <= Nm_ ; j++)
			{
				d_j = diam[j];
				x_j =  M_PI*pow(d_j, 3.0)/6.0;

				S[i] +=   populationBalanceModel_.rhoa().value()
									* populationBalanceModel_.breakupModel().breakupRate(vf, d_j, d_i)
                           * f[j]*f[j]*(x_i-x_i_1)*sqr(alpha)
                           * x_i/(x_j+SMALL)/(x_j+SMALL)/2.;
			}
	}

	for (int i = 0; i<= Nm_-1; i++)
	{
			d_i = diam[i];
			d_i1 = diam[i+1];
			x_i =  M_PI*pow(d_i, 3.0)/6.0;
			x_i1 =  M_PI*pow(d_i1, 3.0)/6.0;

			for(int j = i+1 ; j <= Nm_; j++)
			{
				d_j = diam[j];
				x_j = M_PI*pow(d_j, 3.0)/6.0;
				S[i] += populationBalanceModel_.rhoa().value()
				                  * populationBalanceModel_.breakupModel().breakupRate(vf, d_j, d_i)
                              * f[j]*f[j]*(x_i1-x_i)*sqr(alpha)
                              * x_i/(x_j+SMALL)/(x_j+SMALL)/2.;
	      }
	}

	for (int i = 1; i<= Nm_; i++)
	{
			d_i = diam[i];
			x_i =  M_PI*pow(d_i, 3.0)/6.0;
			S[i] -= populationBalanceModel_.rhoa()
			       *(
			            populationBalanceModel_.breakupModel().breakupFrequency(vf, d_i, fbv25)
                   + populationBalanceModel_.breakupModel().breakupFrequency(vf, d_i, fbv50)
                   + populationBalanceModel_.breakupModel().breakupFrequency(vf, d_i, fbv75)
                 )
                 * f[i]*f[i]*sqr(alpha)
                 * x_i*x_i/(x_i+SMALL)/(x_i+SMALL)/4.;
	}

	return;
}

void qmom::coalescenceKernel
(
	PtrList<volScalarField>& S,
	PtrList<volScalarField>& f,
	const volScalarField& alpha,
	const volScalarField& epsilon
)
const
{
	scalar *diam=new scalar[Nm_+1];
	scalar d_i, d_i1, d_i_1, x_i, x_i1;
	scalar  d_j, d_k, x_k,  x_j, x_i_1, v;
	scalar d_0, d_1, d_n, d_n_1, x_0, x_1, x_n, x_n_1;

	diam[0] = dm_.value()/10.0;

	for(int II=1;II<=Nm_;II++)
	{
      diam[II]= Foam::pow(S_.value(), scalar(II-Nm_/2)/3)*dm_.value();
	}


      d_0 = diam[0];
      d_1 = diam[1];
      d_n = diam[Nm_];
      d_n_1 = diam[Nm_-1];
      x_0 =  M_PI*::pow(d_0, 3.0)/6.0;
      x_1 =  M_PI*::pow(d_1, 3.0)/6.0;
      x_n =  M_PI*::pow(d_n, 3.0)/6.0;
      x_n_1 =  M_PI*::pow(d_n_1, 3.0)/6.0;
		volScalarField epsf = epsilon;

	for (int ki = 0; ki<=Nm_; ki++)
	{
			d_k = diam[ki];
			x_k =  M_PI*::pow(d_k, 3.0)/6.0;

			for(int  j = ki ; j <= Nm_ ; j++)
			{
				d_j = diam[j];
				x_j = M_PI*Foam::pow(d_j, 3.0)/6.0;
				v = x_j+x_k;

				S[0] += populationBalanceModel_.rhoa().value()
				                  * (funcTest(x_0,v,x_1)*(x_1-v)/(x_1-x_0 + SMALL))
                              * (1.0-kronecker(ki,j)/2.)
                              * populationBalanceModel_.coalescenceModel().coalescenceRate(d_k,d_j,epsf)
                              * sqr(alpha)*x_0/(x_j+SMALL)/(x_k+SMALL)*f[ki]*f[j] ;

				S[Nm_] += populationBalanceModel_.rhoa().value()
				                  * (funcTest(x_n_1,v,x_n)*(v-x_n_1)/(x_n-x_n_1+SMALL))
                              * (1.0-kronecker(ki,j)/2.)
                              * populationBalanceModel_.coalescenceModel().coalescenceRate(d_k,d_j,epsf)
                              * sqr(alpha)*x_n/(x_j+SMALL)/(x_k+SMALL)*f[ki]*f[j] ;
			}

			S[0] -= populationBalanceModel_.rhoa().value()
			             * populationBalanceModel_.coalescenceModel().coalescenceRate(d_0,d_k,epsf)
			             * f[ki]*f[0]
                      * sqr(alpha)*x_0/(x_k+SMALL)/(x_0+SMALL);
	}

	for (int i = 1; i<= Nm_-1; i++)
	{
			d_i = diam[i];
			d_i_1 = diam[i-1];
			d_i1 = diam[i+1];
			x_i =  M_PI*::pow(d_i, 3.0)/6.0;
			x_i1 =  M_PI*::pow(d_i1, 3.0)/6.0;
			x_i_1 =  M_PI*::pow(d_i_1, 3.0)/6.0;

			for(int  k = 0 ; k <= Nm_ ; k++)
			{
				d_k = diam[k];
				x_k =  M_PI*::pow(d_k, 3.0)/6.0;
				for(int j = k ; j <= Nm_ ; j++)
				{
					d_j = diam[j];
					x_j =  M_PI*::pow(d_j, 3.0)/6.0;
					v = x_j+x_k;

					S[i] += populationBalanceModel_.rhoa().value()
					                  * ((funcTest(x_i_1, v, x_i)*(v-x_i_1)/(x_i-x_i_1+SMALL))
                                 + (funcTest(x_i,v,x_i1)*(x_i1-v)/(x_i1-x_i+SMALL)))
                                 * (1.0-kronecker(k,j)/2.)
                                 * populationBalanceModel_.coalescenceModel().coalescenceRate(d_k,d_j,epsf)
                                 * f[j]*f[k]*sqr(alpha)*x_i/(x_j+SMALL)/(x_k+SMALL);
				}

				   S[i] -= populationBalanceModel_.rhoa().value()
				                     * populationBalanceModel_.coalescenceModel().coalescenceRate(d_i,d_k,epsf)
				                     * f[k]*f[i]*sqr(alpha)*x_i/(x_k+SMALL)/(x_i+SMALL);
			}
	}

	return ;
}
*/

void qmom::weightsAbscissas(PtrList<volScalarField>& w, PtrList<volScalarField>& L, PtrList<volScalarField>& m, const volScalarField alpha) const
{
  label N=lrint(nodes_);
  scalar norm_position[N*N];
  scalar matriceP[N*N][N*N],ai[N*N],bi[N*N];
  scalar d[N*N],a[N*N][N*N],v[N*N][N*N],alphaii[N*N];


  forAll(alpha, II)
  {
	//- construct the Jacobi matrix P[][] calculation
	for(label i=1; i<=2*N+1; i++)
    {
      if(i==1) {matriceP[i][1]=1.0;}
      else     {matriceP[i][1]=0.0;}
    }
  //  m[0][II] = 1.0 ;  // forced to be unit
      for(label i=1; i<2*N+1; i++)
    {
     norm_position[i-1] = m[i-1][II]/max(m[0][II], 1e-10);
    }

    // second column
    for(label i=1; i<=2*N+1; i++)
    {
      //   matriceP[i][2] = pow(-1,i-1.0)*mag(m[i-1][II]);
         matriceP[i][2] = pow(-1,i-1.0)*norm_position[i-1];
    }

    //rest of the column
    for(label j=3; j<= 2*N+1; j++)
    {
       for(label i=1; i<= 2*N+2-j; i++)
       {
			matriceP[i][j] = matriceP[1][j-1]*matriceP[i+1][j-2]
                           - matriceP[1][j-2]*matriceP[i+1][j-1];
       }
    }

    // calculation of alpha[i]
    alphaii[1]=0.0;
    for(label i=2; i<= 2*N; i++)
    {
       	alphaii[i] = matriceP[1][i+1]/(matriceP[1][i]*matriceP[1][i-1]);
    }

    //- calculation of coefficients a[i] and b[i]
    for(label i=1; i<=N; i++)
    {
		ai[i] = alphaii[2*i] + alphaii[2*i-1];
    }

    for(label i=1; i<= N-1; i++)
    {
		bi[i] = Foam::sqrt(mag(alphaii[2*i+1]*alphaii[2*i])); // positive root
    }

	//- eigenValues and eigenVector for the Jacobi matrice
	// Diagonal of symmetric tridiagonal matrice
	for (int i=1;i<=N;i++)
	{
      a[i][i]   = ai[i];
	}

	// Subdiagonal of symmetric tridiagonal matrix
	for (int i=1;i<N;i++)
	{
      a[i+1][i] = a[i][i+1] = bi[i];
	}
      a[3][1]=a[1][3]=0.;

	// tensor definition
	tensor ABC;
	tensor evs;

	// filling the tensor ABC
	// WARNING ! max N=3 (nodes of quadrature)
	ABC.xx()=a[1][1]; ABC.xy()=a[1][2]; ABC.xz()=a[1][3];
	ABC.yx()=a[2][1]; ABC.yy()=a[2][2]; ABC.yz()=a[2][3];
	ABC.zx()=a[3][1]; ABC.zy()=a[3][2]; ABC.zz()=a[3][3];

	// calculating eigenvalues
	vector e_values=eigenValues(ABC);

	// calculating eigenvector
	evs = eigenVectors(ABC);

	// allocation of the eigenvector
	v[1][1] = evs.xx();
	v[1][2] = evs.yx();
	v[1][3] = evs.zx();

	// allocation of the eigenvalues
	d[1] = e_values.x();
	d[2] = e_values.y();
	d[3] = e_values.z();

	//- final allocation
	for(int i=0; i<=N-1; i++)
	{
      // weights are the first component of eigenVectors
      w[i][II] = mag(m[0][II]*sqr(v[1][i+1]));
      if(w[i][II] <= 0 || w[i][II] > 1) w[i][II]=alpha[II];

      // abscissas are the eigenValues
      L[i][II] = mag(d[i+1]);

      // limitation to 10 mm of diameter
//      if(L[i][II] > 0.01)
	//  {
		//L[i][II] = 0.01;
	  //}
	} // end loop i
  } // end loop forALL

  return;
}

tmp<volScalarField> qmom::SauterDiameter(PtrList<volScalarField>& m, const volScalarField& alpha) const
{
  const fvMesh& mesh = alpha_.mesh();

//  label N=lrint(Nm_) - 1;
//  scalar *diameter=new scalar[N+1];

  volScalarField result
  (
			IOobject
			(
				"result",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			  alpha
                         // dsauter

  );

/*
  for(label II=0; II<=N; II++)
    {
      diameter[II] = Foam::pow(S_.value(), scalar(II-N/2)/3)*dm_.value();
   //   Info<<"Diameter\n"<<diameter[II]<<"\n"<<endl;
    }
*/
/*

  		for(label i=0;i<=N;i++)
  		{


	   		forAll(f[i], cI)
				{
 		   			if(f[i][cI] < 0.0) f[i][cI] = 0.0;
		   			if(f[i][cI] > 1.0) f[i][cI] = 1.0;
				}
	   	sum_v += f[i];
	   	sum_s += (f[i]/diameter[i]);
   		}
*/

  // partial result of Sauter diameter d32
    result = m[3]/(m[2]+SMALL);

  forAll(result, cellI)
  {
	  if(alpha[cellI] > ResidualAlphaForDsauter_.value())
	  {
			if (result[cellI] > 1.e-6 )
			{
				result[cellI] = max(dMin_.value(), min(result[cellI], dMax_.value()));
			}
	  }
	  else
	  {
			   result[cellI] = dMin_.value();
	  }
  }



  return tmp<volScalarField>
  (
		new volScalarField
		(
			IOobject
			(
				"d32",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			result//*
			//dimensionedScalar
			//(
			  // "one",
			   //dimLength,
			   //1.
			//)
		)
  );



}



void qmom::adjust(PtrList<volScalarField>& S, PtrList<volScalarField>& L,const volScalarField& alpha) const
{

	for(int j=0; j<=Nm_-1; j++)
   {
      forAll(S[j], II)
		{
			S[j][II]=0;
		}
   }


	for(int i=0; i<=nodes_-1; i++)
   {
      forAll(L[i], I)
		{
///////////////////////////////////////////////////////////////
            if(alpha[I] > ResidualAlphaForDsauter_.value())
            {


                L[i][I] = max(dMin_.value(), min(L[i][I], dMax_.value()));

            }
            else
            {
                 L[I] = dMin_.value();
            }
///////////////////////////////////////////////////////////////

		 }
   }

}


// *******************************************************************************************************************************

// calculacte and correct variables
void qmom::correct()
{
  const fvMesh& mesh = alpha_.mesh();

 PtrList<volScalarField> source(Nm_);
 PtrList<volScalarField> m(Nm_);
 PtrList<volScalarField> w(nodes_);
 PtrList<volScalarField> L(nodes_);


  volScalarField source_ini =
    mesh.objectRegistry::lookupObject<volScalarField>("source_ini");

  volScalarField m_ini =
    mesh.objectRegistry::lookupObject<volScalarField>("m_ini");

  volScalarField w_ini =
    mesh.objectRegistry::lookupObject<volScalarField>("w_ini");

  volScalarField L_ini =
    mesh.objectRegistry::lookupObject<volScalarField>("L_ini");



  for(label j=0;j<=Nm_-1;j++)
	 {
		 word sourceName = "Sb_" + Foam::name(j);
		 word fName = "m_" + Foam::name(j);


		 source.set
		 ( j,
	      volScalarField
	      (
	         IOobject
	         (
				    sourceName,
					 mesh.time().timeName(),
					 mesh,
					 IOobject::NO_READ,
				         IOobject::AUTO_WRITE
			   ),
	         source_ini
	      )
	    );

		 m.set
		 ( j,
	      volScalarField
	      (
	         IOobject
	         (
				    fName,
					 mesh.time().timeName(),
					 mesh,
					 IOobject::NO_READ,
				         IOobject::AUTO_WRITE
			   ),
                m_ini
	      )
	    );


	 }

   for(label j=0;j<=nodes_-1;j++)
 	 {
 		 word wName = "w_" + Foam::name(j);
 		 word LName = "L_" + Foam::name(j);


 		 w.set
 		 ( j,
 	      volScalarField
 	      (
 	         IOobject
 	         (
 				    wName,
 					 mesh.time().timeName(),
 					 mesh,
 					 IOobject::NO_READ,
 				         IOobject::AUTO_WRITE
 			   ),
 	         w_ini
 	      )
 	    );

 		 L.set
 		 ( j,
 	      volScalarField
 	      (
 	         IOobject
 	         (
 				    LName,
 					 mesh.time().timeName(),
 					 mesh,
 					 IOobject::NO_READ,
 				         IOobject::AUTO_WRITE
 			   ),
                 L_ini
 	      )
 	    );


 	 }

   m.set
 ( 0,
    volScalarField
    (
       IOobject
       (
        "m0",

                   alpha_.time().timeName(),
                   alpha_.db(),
       IOobject::NO_READ,
             IOobject::NO_WRITE
     ),
            m0
    )
  );
  m.set
  ( 1,
     volScalarField
     (
        IOobject
        (
         "m1",

                    alpha_.time().timeName(),
                    alpha_.db(),
        IOobject::NO_READ,
              IOobject::NO_WRITE
      ),
             m0
     )
   );
   m.set
 ( 2,
    volScalarField
    (
       IOobject
       (
        "m2",

                   alpha_.time().timeName(),
                   alpha_.db(),
       IOobject::NO_READ,
             IOobject::NO_WRITE
     ),
            m2
    )
  );
  m.set
( 3,
   volScalarField
   (
      IOobject
      (
       "m3",

                  alpha_.time().timeName(),
                  alpha_.db(),
      IOobject::NO_READ,
            IOobject::NO_WRITE
    ),
           m3
   )
 );
 m.set
( 4,
  volScalarField
  (
     IOobject
     (
      "m4",

                 alpha_.time().timeName(),
                 alpha_.db(),
     IOobject::NO_READ,
           IOobject::NO_WRITE
   ),
          m4
  )
);
m.set
( 5,
 volScalarField
 (
    IOobject
    (
     "m5",

                alpha_.time().timeName(),
                alpha_.db(),
    IOobject::NO_READ,
          IOobject::NO_WRITE
  ),
         m5
 )
);
L.set
( 0,
 volScalarField
 (
    IOobject
    (
     "L0",

                alpha_.time().timeName(),
                alpha_.db(),
    IOobject::NO_READ,
          IOobject::NO_WRITE
  ),
         L0
 )
);
L.set
( 1,
 volScalarField
 (
    IOobject
    (
     "L1",

                alpha_.time().timeName(),
                alpha_.db(),
    IOobject::NO_READ,
          IOobject::NO_WRITE
  ),
         L1
 )
);
L.set
( 2,
 volScalarField
 (
    IOobject
    (
     "L2",

                alpha_.time().timeName(),
                alpha_.db(),
    IOobject::NO_READ,
          IOobject::NO_WRITE
  ),
         L2
 )
);
w.set
( 0,
 volScalarField
 (
    IOobject
    (
     "w0",

                alpha_.time().timeName(),
                alpha_.db(),
    IOobject::NO_READ,
          IOobject::NO_WRITE
  ),
         w0
 )
);
w.set
( 1,
 volScalarField
 (
    IOobject
    (
     "w1",

                alpha_.time().timeName(),
                alpha_.db(),
    IOobject::NO_READ,
          IOobject::NO_WRITE
  ),
         w1
 )
);
w.set
( 2,
 volScalarField
 (
    IOobject
    (
     "w2",

                alpha_.time().timeName(),
                alpha_.db(),
    IOobject::NO_READ,
          IOobject::NO_WRITE
  ),
         w2
 )
);

m[0]=m0;m[1]=m1;m[2]=m2;m[3]=m3;m[4]=m4;m[5]=m5;
w[0]=w0;w[1]=w1;w[2]=w2;
L[0]=L0;L[1]=L1;L[2]=L2;


//////////////////////////////////////////////////////////////////////////////////////////////////////////


  volScalarField alpha =
    mesh.objectRegistry::lookupObject<volScalarField>("alpha.air");

  volScalarField epsilon =
   mesh.objectRegistry::lookupObject<volScalarField>("epsilon");


  surfaceScalarField phia =
    mesh.objectRegistry::lookupObject<surfaceScalarField>("phi.air");




		qmom::adjust(source,L,alpha);




 		for(label i=0;i<=Nm_-1;i++)
		{
			// update bubbles breakage/coalescence

	//		qmom::coalescenceKernel(source,f,alpha,epsilon);
	//		qmom::breakupKernel(source,f,alpha) ;
	//		volScalarField Sb = source[i]/
			                        //   populationBalanceModel_.breakupModel().rhoa();



	//		if(i != (Nm_-1)/2)

		//	{
     word fScheme("div(phia,mi)");


     fvScalarMatrix miEqn
     (
         fvm::ddt(alpha, m[i])
        + fvm::div(fvc::flux(phia, alpha, fScheme), m[i], fScheme)
        - fvm::Sp(fvc::div(phia), m[i])
     );

	//			solve(miEqn == Sb);
          miEqn.solve();


	//		}
		}



m[0]=mag(m[0]);
m[1]=mag(m[1]);
m[2]=mag(m[2]);
m[3]=mag(m[3]);
m[4]=mag(m[4]);
m[5]=mag(m[5]);

    qmom::weightsAbscissas(w,L,m,alpha);
  Info << "banana" << "\t" << endl; //moshkel injast

		qmom::adjust(source,L,alpha);


    w[0]=mag(w[0])/(mag(w[0])+mag(w[1])+mag(w[2])+SMALL);
    w[1]=mag(w[1])/(mag(w[0])+mag(w[1])+mag(w[2])+SMALL);
    w[2]=mag(w[2])/(mag(w[0])+mag(w[1])+mag(w[2])+SMALL);

    L[0]=0.012*mag(L[0])/(mag(L[0])+mag(L[1])+mag(L[2])+SMALL);
    L[1]=0.012*mag(L[1])/(mag(L[0])+mag(L[1])+mag(L[2])+SMALL);
    L[2]=0.012*mag(L[2])/(mag(L[0])+mag(L[1])+mag(L[2])+SMALL);


/*
   // bounding f
   volScalarField alphai=f[0];

	for(int i=1; i<Nm_; i++)
	{
		if(i!=(Nm_-1)/2) alphai = alphai + f[i];

	}

	f[(Nm_-1)/2] = 1. - alphai;


	f[(Nm_-1)/2] = max(f[(Nm_-1)/2], scalar(0));
	f[(Nm_-1)/2] = min(f[(Nm_-1)/2], scalar(1));

	alphai = alphai + f[(Nm_-1)/2];



	for(int i=0; i<Nm_; i++)
	{
		f[i]=f[i]/(alphai + SMALL);
	}


	forAll(f[(Nm_-1)/2], cellI)
	{
		if(alpha[cellI] < ResidualAlphaForCorrect_.value())
		{
			f[(Nm_-1)/2][cellI]=0;
		}
	}


*/


	qmom::adjust(source,m,alpha);




	// update Sauter diameter d32, relax and correct BC
        dsauter = SauterDiameter(m,alpha);




	d32_ = SauterDiameter(m,alpha)*dimensionedScalar
			(
			   "one",
			   dimLength,
			   1.
			);

        dsauter.relax();
        dsauter.correctBoundaryConditions();


        d32_.relax();
        d32_.correctBoundaryConditions();

m0=m[0];m1=m[1];m2=m[2];m3=m[3];m4=m[4];m5=m[5];
w0=w[0];w1=w[1];w2=w[2];
L0=L[0];L1=L[1];L2=L[2];


}

// update final calculation






void qmom::updateSauterDiameter()
{
   // mesh declaration




   Info << nl << " *********************************************************************** "<< endl;
   Info << nl << "           Updating Sauter diameter d32 by using " << populationBalanceReturn::name() << endl;
   Info << nl << "           min/avg/max d32: " << min(dsauter).value()*1000
              << "/" << average(dsauter).value()*1000<< "/" << max(dsauter).value()*1000<< "mm." << endl;
   Info << nl << " *********************************************************************** "<< endl;
	Info << " " << endl;
}


} // end of namespace
// **************************************************************************************************************************** //
