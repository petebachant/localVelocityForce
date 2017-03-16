/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "localVelocityForce.H"
#include "fvMatrices.H"
#include "DimensionedField.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(localVelocityForce, 1);

    addToRunTimeSelectionTable
    (
        option,
        localVelocityForce,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::fv::localVelocityForce::writePropslocal
(
    scalarField gradPlocal_
) const
{
    // Only write on output time
    if (mesh_.time().outputTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                name_ + "_gradPlocal",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            )
        );
        propsDict.add("gradient", gradPlocal_);
        propsDict.regIOobject::write();
     }        
}


void Foam::fv::localVelocityForce::writeProps
(
    const scalar gradP
) const
{
    // Only write on output time
    if (mesh_.time().outputTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                name_ + "gradP",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            )
        );
        propsDict.add("gradient", gradP);
        propsDict.regIOobject::write();
     }        
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::fv::localVelocityForce::localVelocityForce
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(sourceName, modelType, dict, mesh),
    Ubar_(coeffs_.lookup("Ubar")),
    gradP0_(0.0),
    dGradP_(0.0),
    flowDir_(Ubar_/mag(Ubar_)),
    relaxation_(coeffs_.lookupOrDefault<scalar>("relaxation", 1.0)),
    rAPtr_(NULL)

{
    coeffs_.lookup("fields") >> fieldNames_;
	coeffs_.lookup("freestream") >> freestreamVel_;
	coeffs_.lookup("position") >> position_;
	coeffs_.lookup("gaussFactor") >> gaussFactor_;
	coeffs_.lookup("gaussType") >> gaussType_;
	coeffs_.lookup("diameter") >> diameter_;
	coeffs_.lookup("height") >> height_;

    if (fieldNames_.size() != 1)
    {
        FatalErrorIn
        (
            "Foam::fv::localVelocityForce::"
            "localVelocityForce"
            "("
                "const word&, "
                "const word&, "
                "const dictionary&, "
                "const fvMesh&"
            ")"
        )   << "Source can only be applied to a single field.  Current "
            << "settings are:" << fieldNames_ << exit(FatalError);
    }

    applied_.setSize(fieldNames_.size(), false);

    // Read the initial pressure gradient from file if it exists
    IFstream propsFile
    (
        mesh_.time().timePath()/"uniform"/(name_ + "Properties")
    );

    if (propsFile.good())
    {
        Info<< "    Reading pressure gradient from file" << endl;
        dictionary propsDict(dictionary::null, propsFile);
        propsDict.lookup("gradient") >> gradP0_;
    }

    Info<< "    Initial pressure gradient = " << gradP0_ << nl << endl;
    
    
        // Read nu from object registry
    const dictionary& transportProperties = mesh_.lookupObject<IOdictionary>
    (
        "transportProperties"
    );   
    dimensionedScalar nu;
    transportProperties.lookup("nu") >> nu;
    nu_ = nu.value();
    
    
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fv::localVelocityForce::magUbarAve
(
    const volVectorField& U
) const
{
    scalar magUbarAve = 0.0;

    const scalarField& cv = mesh_.V();
    forAll(cells_, i)
    {
        label cellI = cells_[i];
        scalar volCell = cv[cellI];
        magUbarAve += (flowDir_ & U[cellI])*volCell;
    }

    // reduce(magUbarAve, sumOp<scalar>());

    magUbarAve /= V_;    

    return magUbarAve;
}


Foam::scalarField Foam::fv::localVelocityForce::magUbarlocal
(
    const volVectorField& U
) const
{
    scalarField magUbarlocal_(cells_.size(), 0.0);

    const scalarField& cv = mesh_.V();
    forAll(cells_, i)
    {
        label cellI = cells_[i];
        magUbarlocal_[i] = (flowDir_ & U[cellI]);
    }

    // reduce(magUbarlocal_, sumOp<scalarField>());

    return magUbarlocal_;
}


void Foam::fv::localVelocityForce::correct(volVectorField& U)
{
    const scalarField& rAU = rAPtr_().internalField(); 
    //1.0/UEqn.A() Diagonal part of matrix U [s/m]

    // Integrate flow variables over cell set
    scalar rAUave = 0.0;

    const scalarField& cv = mesh_.V();
    
    forAll(cells_, i)
    {
        label cellI = cells_[i];
        scalar volCell = cv[cellI];
        rAUave += rAU[cellI]*volCell;
    }

    // Collect across all processors
    // reduce(rAUave, sumOp<scalar>());

    // Volume averages
    rAUave /= V_;
    
    // Integrate flow variables over cell set
    scalarField rAUlocal_(cells_.size(), 0.0);
    //_(cells_.size(),0);

    forAll(cells_, i)
    {
        label cellI = cells_[i];
        rAUlocal_[i] = rAU[cellI];
    }

    // Collect across all processors
    // reduce(rAUlocal_, sumOp<scalarField>());

 
    scalar magUbarAve = this->magUbarAve(U);   
    scalarField magUbarlocal_ = this->magUbarlocal(U);  

    // Calculate the pressure gradient increment needed to adjust the average
    // flow-rate to the desired value

    dGradP_ = relaxation_*(mag(Ubar_) - magUbarAve)/rAUave;   
    dGradPlocal_ = relaxation_*(mag(Ubar_) - magUbarlocal_)/(1e-12+rAUlocal_);

    // Apply correction to velocity field          
    forAll(cells_, i)
    {
        label cellI = cells_[i];
        //U[cellI] += flowDir_*rAU[cellI]*dGradP_;
        U[cellI] =  flowDir_ * rAU[cellI] * dGradPlocal_[cellI];
    }

	scalar gradP = gradP0_ + dGradP_; 
	scalarField gradPlocal_ = gradP0_+ dGradPlocal_; 
	vectorField ibSource_(cells_.size(),vector::zero);
	ibSource_ = (flowDir_*gradPlocal_);

    Info<< "Pressure gradient source: uncorrected Ubar = " << magUbarAve
        << ", pressure gradient = " << gradP << endl;
        
    //Info<< "Pressure gradient source: uncorrected U local = " << magUlocal_
        //<< ", pressure gradient local = " << dGradPlocal_ << endl;

    writePropslocal(gradPlocal_);
	
    if (debug == 2)
    {
        //Info << "dGradP_ :" << dGradP_ << endl;
        //Info << "dGradPlocal_ :" << dGradPlocal_ << endl;
        Info << "ibSource_ :" << ibSource_ << endl;
    }   	
}


void Foam::fv::localVelocityForce::addSup
(
    const vectorField& ibSource_,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{   
    DimensionedField<vector, volMesh> Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldI] + "_Su",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", eqn.dimensions()/dimVolume, vector::zero)
    );

    UIndirectList<vector>(Su,cells_) = ibSource_;

    eqn += Su;
}


void Foam::fv::localVelocityForce::addSup
(
    const vectorField& ibSource_,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{

    this->addSup(ibSource_, eqn, fieldI);
}


void Foam::fv::localVelocityForce::constrain
(
    fvMatrix<vector>& eqn,
    const label
)
{
    if (rAPtr_.empty())
    {
        rAPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    name_ + ":rA",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                1.0/eqn.A()
            )
        );
    }
    else
    {
        rAPtr_() = 1.0/eqn.A();
    }

    gradP0_ += dGradP_;
    dGradP_ = 0.0;
}


// ************************************************************************* //
