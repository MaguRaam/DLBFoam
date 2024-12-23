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

Class
    Foam::cemaChemistryModel

Description

SourceFiles
    cemaChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef cemaChemistryModel_H
#define cemaChemistryModel_H


#include "EigenMatrix.H"
#include "loadBalanced_pyJacChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class cemaChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class cemaChemistryModel
:
    public loadBalanced_pyJacChemistryModel<ThermoType>
{
    // Private Data

        //- The number of elements
        label nElements_;

        //- Chemical explosive mode
        mutable volScalarField cem_;


    // Private Member Functions


public:

    //- Runtime type information
    TypeName("cema");


    // Constructors

        //- Construct from thermo
        cemaChemistryModel(const fluidMulticomponentThermo& thermo);

        //- Disallow default bitwise copy construction
        cemaChemistryModel(const cemaChemistryModel<ThermoType>&) = delete;


    //- Destructor
    virtual ~cemaChemistryModel();


    // Member Functions

        //- Calculate maximum Jacobian eigenvalue
        virtual scalar cema(const scalarSquareMatrix& J) const;

        // Overrides to loadBalanced_pyJacChemistryModel

            //- Calculate the ODE jacobian
            virtual void jacobian
            (
                const scalar t,
                const scalarField& YTp,
                const label li,
                scalarField& dYTpdt,
                scalarSquareMatrix& J
            ) const override;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cemaChemistryModel<ThermoType>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cemaChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
