/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "minMaxLength.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(minMaxLength, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        minMaxLength,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
// useless:
bool Foam::functionObjects::minMaxLength::calc(){return true;}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::minMaxLength::minMaxLength
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, dict.lookup("fieldName")),
    position(dict.lookup("position")),
    direction(dict.lookup("direction")),
    criterion(readScalar(dict.lookup("criterion"))),
    threshold(readScalar(dict.lookup("threshold"))),
    ResultOutPut("minMaxLength")
{
    Info<<"criterion = "<<criterion<<endl;
    Info<<"threshold = "<<threshold<<endl;
    ResultOutPut << "time" << ",minLength" << ",maxLength" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::minMaxLength::~minMaxLength()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// useless:
bool Foam::functionObjects::minMaxLength::clear(){return true;}


bool Foam::functionObjects::minMaxLength::execute()
{
    if (foundObject<volScalarField>(fieldName_))
    {
        const volScalarField& field = lookupObject<volScalarField>(fieldName_);
        const scalar field_max = max(field).value();
        Info << "field_max = " << field_max << endl;
        scalar maxLength = 0;
        scalar minLength = great;
        if (field_max > threshold)
        {
            forAll (field, cellI)
            {
                if (field[cellI] >= criterion)
                {
                    vector raw = position - mesh_.C()[cellI];
                    scalar length = mag(raw&direction);

                    if ( length > maxLength)
                    {
                        maxLength = length;
                    }
                    if ( length < minLength)
                    {
                        minLength = length;
                    }
                }
            }
        }
        reduce(maxLength, maxOp<scalar>());
        reduce(minLength, minOp<scalar>());

        Info << "maxLength = " << maxLength << endl;
        Info << "minLength = " << minLength << endl;
        ResultOutPut << mesh_.time().value() << "," << minLength << "," << maxLength << endl;
    }
    else
    {
        Info << "Sorry! I cannot found " << fieldName_ << endl;
    }

    return true;
}

bool Foam::functionObjects::minMaxLength::write()
{
    return true;
}


// ************************************************************************* //