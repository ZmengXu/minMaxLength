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
#include "IOstreams.H"
#include "dictionaryEntry.H"
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
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::minMaxLength::writeFileHeader(const label i)
{

    OFstream& file = this->file();

    writeHeader(file, "maximum and minimum length of a field reaching a criterion value to a cutting plane ");
    writeCommented(file, "Time");

    forAll(fieldSet_, fieldi )
    {
        writeTabbed(file, "minLength(" + fieldSet_[fieldi] + ')');
        writeTabbed(file, "maxLength(" + fieldSet_[fieldi] + ')');
    }

    file<< endl;

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::minMaxLength::minMaxLength
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
//    fvMeshFunctionObject(name, runTime, dict),
    fieldExpression(name, runTime, dict, dict.lookup("fieldName")),
    logFiles(obr_, name),
    position(dict.lookup("position")),
    direction(dict.lookup("direction")),
    ResultOutPut("minMaxLength"),
    fieldSet_(dict.lookup("fields")),
    criterion(fieldSet_.size()),
    threshold(fieldSet_.size())
{
    Info<<"criterion = "<<criterion<<endl;
    Info<<"threshold = "<<threshold<<endl;


//OFstream& file = this->file();
//Istream& is();
//const dictionaryEntry entry(dict.subDict("aa"), is);
/*
fieldName_ = entry.keyword();
criterion[i] = readScalar(dict.lookup("criterion"))
threshold[i] = entry.lookupOrDefault<scalar>("threshold", 0.01*criterion[i]);
*/
    forAll( fieldSet_, fieldi)
    {
        const word& fieldName = fieldSet_[fieldi];
        criterion[fieldi] = readScalar(dict.subDict(fieldName).lookup("criterion"));
        threshold[fieldi] = dict.lookupOrDefault<scalar>("threshold", 0.01*criterion[fieldi]);
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::minMaxLength::~minMaxLength()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::minMaxLength::clear()
{
    return true;
}

bool Foam::functionObjects::minMaxLength::read(const dictionary&)
{
    return true;
}

bool Foam::functionObjects::minMaxLength::execute()
{
    return true;
}

bool Foam::functionObjects::minMaxLength::write()
{
    logFiles::write();

    if (Pstream::master()) writeTime(file());
    Log << type() << " " << name() <<  " write:" << nl;

    forAll(fieldSet_, fieldi)
    {
        const word& fieldName = fieldSet_[fieldi];
        const volScalarField& fieldValue = lookupObject<volScalarField>(fieldName);
        const scalar field_max = max(fieldValue).value();
        Info << "field_max = " << field_max << endl;
        scalar maxLength = 0;
        scalar minLength = great;
        if (field_max > threshold[fieldi])
        {
            forAll (fieldValue, cellI)
            {
                if (fieldValue[cellI] >= criterion[fieldi])
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

        file()<< tab << minLength << tab << maxLength;

        Log << "    min/max(" << fieldName << ") = "
            << minLength << ' ' << maxLength;

        Info << fieldName << "maxLength = " << maxLength << endl;
        Info << fieldName << "minLength = " << minLength << endl;
    }
    if (Pstream::master()) file() << endl;
    Log << endl;

    return true;
}


// ************************************************************************* //