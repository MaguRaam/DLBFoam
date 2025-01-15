/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "cemaChemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::cemaChemistryModel<ThermoType>::cemaChemistryModel
(
    const fluidMulticomponentThermo& thermo
)
:
    loadBalanced_pyJacChemistryModel<ThermoType>(thermo),
    nElements_(this->template lookup<label>("nElements")),
    lambdaExp_
    (
        IOobject
        (
            "lambdaExp",
            this->mesh().time().name(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, 0)
    )
{
    Info<< "cemaChemistryModel: Number of elements = " << nElements_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::cemaChemistryModel<ThermoType>::~cemaChemistryModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::cemaChemistryModel<ThermoType>::cema
(
    const label celli,
    const scalarField& YTp,
    const scalarSquareMatrix& J
) const
{
    // Compute eigen values
    eigendecomposition eigen(J);
    scalarField EValsRe(eigen.d());
    const scalarField& EValsIm(eigen.e());

    // Get the size of the matrix
    const label m = EValsRe.size();

    // Compute magnitude square of eigen values
    const scalarField EValsMag(EValsRe*EValsRe + EValsIm*EValsIm);

    // Get the indices of magnitude of eigenvalues in increasing order
    labelList order;
    sortedOrder(EValsMag, order);

    // Skip conservation modes for elements and temperature
    const scalar smallestEVal = -vGreat;

    for (label i = 0; i < nElements_ + 1; ++i)
    {
        EValsRe[order[i]] = smallestEVal;
    }

    // Find the index of maximum real part of eigenvalue
    label iMax = findMax(EValsRe);

    // Get the eigenvector corresponding to the maximum real part of eigenvalue
    Field<scalar> cem = eigen.V().col(m, 0, iMax);

    // Compute the maximum real part of eigenvalue
    lambdaExp_[celli] = EValsRe[iMax];
}


template<class ThermoType>
Foam::scalar Foam::cemaChemistryModel<ThermoType>::solve
(
    const scalar deltaT
)
{
    return chemistryModel<ThermoType>::solve(deltaT);
}


template<class ThermoType>
void Foam::cemaChemistryModel<ThermoType>::jacobian
(
    const scalar t,
    const scalarField& YTp,
    const label celli,
    scalarField& dYTpdt,
    scalarSquareMatrix& J
) const
{
    chemistryModel<ThermoType>::jacobian(t, YTp, celli, dYTpdt, J);

    // Calculate cema fields
    cema(celli, YTp, J);
}


// ************************************************************************* //