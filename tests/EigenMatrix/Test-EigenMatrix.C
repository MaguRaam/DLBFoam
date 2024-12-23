/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "EigenMatrix.H"

using namespace Foam;

void testMatrix(const SquareMatrix<doubleScalar>& matrix, word testName)
{
    Info << "Test: " << testName << "\n" << endl;
    Info << "Matrix:\n" << matrix << endl;

    const EigenMatrix<doubleScalar> EM(matrix, true);
    Info << "Eigenvalues (Real): " << EM.EValsRe() << endl;
    Info << "Eigenvalues (Imaginary): " << EM.EValsIm() << endl;
    Info << "\n------------------------\n" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    // Initialize and test Identity Matrix
    SquareMatrix<doubleScalar> identityMatrix(3, Zero);
    identityMatrix(0, 0) = 1; identityMatrix(1, 1) = 1; identityMatrix(2, 2) = 1;
    testMatrix(identityMatrix, "Identity Matrix");

    // Initialize and test Diagonal Matrix
    SquareMatrix<doubleScalar> diagonalMatrix(3, Zero);
    diagonalMatrix(0, 0) = 3; diagonalMatrix(1, 1) = 5; diagonalMatrix(2, 2) = 7;
    testMatrix(diagonalMatrix, "Diagonal Matrix");

    // Initialize and test Symmetric Matrix
    SquareMatrix<doubleScalar> symmetricMatrix(3, Zero);
    symmetricMatrix(0, 0) = 6; symmetricMatrix(0, 1) = 2; symmetricMatrix(1, 0) = 2; symmetricMatrix(1, 1) = 3; symmetricMatrix(2, 2) = 5;
    testMatrix(symmetricMatrix, "Symmetric Matrix");

    // Initialize and test Random Symmetric Matrix
    SquareMatrix<doubleScalar> randomSymmetricMatrix(3, Zero);
    randomSymmetricMatrix(0, 0) = 4; randomSymmetricMatrix(0, 1) = -1; randomSymmetricMatrix(1, 0) = -1;
    randomSymmetricMatrix(1, 1) = 2; randomSymmetricMatrix(0, 2) = 0.5; randomSymmetricMatrix(2, 0) = 0.5; randomSymmetricMatrix(2, 2) = 3;
    testMatrix(randomSymmetricMatrix, "Random Symmetric Matrix");

    Info << "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
