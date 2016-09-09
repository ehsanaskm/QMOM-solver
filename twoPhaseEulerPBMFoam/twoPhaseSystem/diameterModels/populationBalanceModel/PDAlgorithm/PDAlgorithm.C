/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2004 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------- -----------------------------------------------
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
    
\*----------------------------------------------------------------------------*/

#include "PDAlgorithm.H"
#include "scalar.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components


Foam::PDAlgorithm::PDAlgorithm
(
    const fvMesh& mesh,
    const dictionary& populationBlanceProperties, 
    const word& ModelName  
)
:
    dict_
    (
        populationBlanceProperties.subDict("Ramkrishna")
    ),
    name_(ModelName),
    
    dmoy_
    (
        dict_.lookup("dmoy")
    ),    
    
    noeuds_
    (
        dict_.lookup("noeuds")
    )
    
{}

Foam::autoPtr<Foam::PDAlgorithm> Foam::PDAlgorithm::New
(
    const fvMesh& mesh,
    const dictionary& populationBlanceProperties, 
    const word& ModelName
   
)
{
    return autoPtr<PDAlgorithm>
    (
        new PDAlgorithm(mesh, populationBlanceProperties, ModelName)
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //



Foam::PDAlgorithm::~PDAlgorithm()
{}

void Foam::PDAlgorithm::WeightsAbscissas
(
     PtrList<volScalarField>& W, // Array
     PtrList<volScalarField>& L, // Array
     PtrList<volScalarField>& position,  // Array
     const volScalarField moment0

)
const
{
 
  label N=lrint (noeuds().value());
  scalar norm_position[N*N];
  scalar matriceP[N*N][N*N],ai[N*N],bi[N*N];
  scalar d[N*N],a[N*N][N*N],v[N*N][N*N],alphai[N*N]; 

Info << "POSITION0JUSTBEFOREPD" << average(position[0]) << endl;

  forAll(moment0, II)
  {
	//- construct the Jacobi matrix P[][] calculation
	for(label i=1; i<=2*N+1; i++)
    {
      if(i==1) {matriceP[i][1]=1.0;}
      else     {matriceP[i][1]=0.0;}
    } 

   // position[0][II] = 1.0 ;  // position[0][II] is moment0, which is forced to be unit here

    for(label i=1; i<2*N+1; i++)
    {
     norm_position[i-1] = position[i-1][II]/max(position[0][II], 1e-10);
    } 

    // second column
    for(label i=1; i<=2*N; i++)
    {
         //  matriceP[i][2] = pow(-1,i-1.0)*mag(position[i-1][II]);  
         //  matriceP[i][2] = pow(-1,i-1.0)*(position[i-1][II]);
          matriceP[i][2] = pow(-1,i-1.0)*norm_position[i-1];
    }
     


    //rest of the column   
    for(label j=3; j<= 2*N+1; j++)
    {
       for(label i=1; i<= 2*N+2-j; i++)   
       {
	matriceP[i][j] = matriceP[1][j-1]*matriceP[i+1][j-2]-matriceP[1][j-2]*matriceP[i+1][j-1];
                                 
       }
                     

    } 
       

// calculation of alphai[i]  
    alphai[1]=0.0; // Zeta1 =0.0
    for(label i=2; i<= 2*N; i++)
    {
       	alphai[i] = matriceP[1][i+1]/(matriceP[1][i]*matriceP[1][i-1]+SMALL);
      // alphai[i] = matriceP[1][i+1]/(matriceP[1][i]*matriceP[1][i-1]+SMALL); 
       // Info << "matriceP[1][i]" << matriceP[1][i] << endl;
       // Info << "matriceP[1][i-1]" << matriceP[1][i-1] << endl;
       // Info << "teste" << matriceP[1][i]*matriceP[1][i-1]+SMALL << endl;
    }
    
    //- calculation of coefficients a[i] and b[i]
    for(label i=1; i<=N; i++)
    {
		ai[i] = alphai[2*i] + alphai[2*i-1]; 
    }
          /* Info << "ai1" << "\t" << ai[1] << endl;
           Info << "ai2" << "\t" << ai[2] << endl;
           Info << "ai3" << "\t" << ai[3] << endl;*/


    for(label i=1; i<= N-1; i++)
    {
		bi[i] = Foam::sqrt(mag(alphai[2*i+1]*alphai[2*i])); // positive root
                //bt[i] = -Foam::sqrt(alphai[2*i+1]*alphai[2*i]); // positive root
    } 

          
          // Info << "bi1" << "\t" << bi[1] << endl;
          // Info << "bi2" << "\t" << bi[2] << endl;
  
	//- eigenValues and eigenVector for the Jacobi matrice 
	// Diagonal of symmetric tridiagonal matrice
	for (label i=1;i<=N;i++)
	{
      a[i][i]   = ai[i];
	}  
	
	// Subdiagonal of symmetric tridiagonal matrix
	for (label i=1;i<N;i++)
	{
      a[i+1][i] = a[i][i+1] = bi[i];
	}  
      a[3][1]=a[1][3]=1e-10; // before it was 0.
  
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

/*Info << "v11" << "\t" << v[1][1] << endl;
Info << "v12" << "\t" << v[1][2] << endl;
Info << "v13" << "\t" << v[1][3] << endl;*/

	// allocation of the eigenvalues
	d[1] = e_values.x(); 
	d[2] = e_values.y();
	d[3] = e_values.z();

/*Info << "EigenVa1" << "\t" << d[1] << endl;
Info << "EigenVa2" << "\t" << d[2] << endl;
Info << "EigenVa3" << "\t" << d[3] << endl;*/

   //- final allocation
	for(int i=1; i<=N; i++)
	{
      // Weights (W) are the first component of eigenVectors 
          //   W[i][II] = mag(position[0][II]*sqr(v[1][i]));
         //  W[i][II] = mag(sqr(v[1][i]));
       //   if(W[i][II] <= 0 || W[i][II] > 1) W[i][II]=moment0[II];     

         W[i][II] = position[0][II]*sqr(v[1][i]);

             if(W[i][II] < 0)
             {
                W[i][II]=0.0;
             }           
         
      // Abscissas (L) are the eigenValues
      L[i][II] = fabs(d[i]);
      

  
	} // end loop i
    
    
  } // end loop forALL

Info << "position0PD" << average(position[0]) << endl;

/*Info << "matriceP[1][2]" << "\t" << matriceP[1][3-1] << endl;
        Info << "matriceP[2][1]" << "\t" << matriceP[1+1][3-2] << endl;   
        Info << "matriceP[1][1]" << "\t" << matriceP[1][3-2] << endl;
        Info << "matriceP[2][2]" << "\t" << matriceP[1+1][3-1] << endl;


  Info << "W1PD" << average(W[1]) << endl;
Info << "W2PD" << average(W[2]) << endl;
Info << "W3PD" << average(W[3]) << endl;
  Info << "L1PD" << average(L[1]) << endl;
Info << "L2PD" << average(L[2]) << endl;
Info << "L3PD" << average(L[3]) << endl;

Info << "Pmatrice11" << "\t" << matriceP[1][1] << endl;
Info << "Pmatrice12" << "\t" << matriceP[1][2] << endl;
Info << "Pmatrice13" << "\t" << matriceP[1][3] << endl;
Info << "Pmatrice14" << "\t" << matriceP[1][4] << endl;
Info << "Pmatrice15" << "\t" << matriceP[1][5] << endl;
Info << "Pmatrice16" << "\t" << matriceP[1][6] << endl;
Info << "Pmatrice17" << "\t" << matriceP[1][7] << endl;


Info << "Pmatrice21" << "\t" << matriceP[2][1] << endl;
Info << "Pmatrice22" << "\t" << matriceP[2][2] << endl;
Info << "Pmatrice23" << "\t" << matriceP[2][3] << endl;
Info << "Pmatrice24" << "\t" << matriceP[2][4] << endl;
Info << "Pmatrice25" << "\t" << matriceP[2][5] << endl;
Info << "Pmatrice26" << "\t" << matriceP[2][6] << endl;
Info << "Pmatrice27" << "\t" << matriceP[2][7] << endl;

Info << "Pmatrice31" << "\t" << matriceP[3][1] << endl;
Info << "Pmatrice32" << "\t" << matriceP[3][2] << endl;
Info << "Pmatrice33" << "\t" << matriceP[3][3] << endl;
Info << "Pmatrice34" << "\t" << matriceP[3][4] << endl;
Info << "Pmatrice35" << "\t" << matriceP[3][5] << endl;
Info << "Pmatrice36" << "\t" << matriceP[3][6] << endl;
Info << "Pmatrice37" << "\t" << matriceP[3][7] << endl;

Info << "Pmatrice41" << "\t" << matriceP[4][1] << endl;
Info << "Pmatrice42" << "\t" << matriceP[4][2] << endl;
Info << "Pmatrice43" << "\t" << matriceP[4][3] << endl;
Info << "Pmatrice44" << "\t" << matriceP[4][4] << endl;
Info << "Pmatrice45" << "\t" << matriceP[4][5] << endl;
Info << "Pmatrice46" << "\t" << matriceP[4][6] << endl;
Info << "Pmatrice47" << "\t" << matriceP[4][7] << endl;

Info << "Pmatrice51" << "\t" << matriceP[5][1] << endl;
Info << "Pmatrice52" << "\t" << matriceP[5][2] << endl;
Info << "Pmatrice53" << "\t" << matriceP[5][3] << endl;
Info << "Pmatrice54" << "\t" << matriceP[5][4] << endl;
Info << "Pmatrice55" << "\t" << matriceP[5][5] << endl;
Info << "Pmatrice56" << "\t" << matriceP[5][6] << endl;
Info << "Pmatrice57" << "\t" << matriceP[5][7] << endl;


Info << "Pmatrice61" << "\t" << matriceP[6][1] << endl;
Info << "Pmatrice62" << "\t" << matriceP[6][2] << endl;
Info << "Pmatrice63" << "\t" << matriceP[6][3] << endl;
Info << "Pmatrice64" << "\t" << matriceP[6][4] << endl;
Info << "Pmatrice65" << "\t" << matriceP[6][5] << endl;
Info << "Pmatrice66" << "\t" << matriceP[6][6] << endl;
Info << "Pmatrice67" << "\t" << matriceP[6][7] << endl;


Info << "Pmatrice71" << "\t" << matriceP[7][1] << endl;
Info << "Pmatrice72" << "\t" << matriceP[7][2] << endl;
Info << "Pmatrice73" << "\t" << matriceP[7][3] << endl;
Info << "Pmatrice74" << "\t" << matriceP[7][4] << endl;
Info << "Pmatrice75" << "\t" << matriceP[7][5] << endl;
Info << "Pmatrice76" << "\t" << matriceP[7][6] << endl;
Info << "Pmatrice77" << "\t" << matriceP[7][7] << endl;

*/

        


  return;
 
}










