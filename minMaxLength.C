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
#include "addToRunTimeSelectionTable.H"
#include "dictionaryEntry.H"
#include "fieldTypes.H"
#include "fvCFD.H"
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

// * * * * * * * * * * * * * * * fieldItem Class  * * * * * * * * * * * * * //

Foam::functionObjects::fieldItem::~fieldItem()
{}

Foam::functionObjects::fieldItem::fieldItem()
:
    fieldName_("unknown"),
    criterion_(0),
    threshold_(0)
{}

Foam::Istream& Foam::functionObjects::operator>>
(
    Istream& is,
    fieldItem& faItem
)
{
    is.check
    (
        "Foam::Istream& Foam::operator>>"
        "(Foam::Istream&, Foam::functionObjects::fieldItem&)"
    );

    const dictionaryEntry entry(dictionary::null, is);

    faItem.fieldName_ = entry.keyword();

    faItem.criterion_ = readScalar(entry.lookup("criterion"));

    faItem.threshold_ = entry.lookupOrDefault<scalar>("threshold", 0.01*faItem.criterion_);

    return is;
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::minMaxLength::writeFileHeader(const label i)
{

    OFstream& file = this->file();

    writeHeader(file, "maximum and minimum length of a field reaching a criterion value to a cutting plane");

    writeCommented(file, "The cutting plane position and direction :");
    file << position_ << ", and " << direction_ << endl;
    
    writeCommented(file, "The fieldName, criterion and threshold are :");
    forAll(fieldItems_, fieldi)
    {
        file<<"fieldName = "<<fieldItems_[fieldi].fieldName()<<",";
        file<<"criterion = "<<fieldItems_[fieldi].criterion()<<",";
        file<<"threshold = "<<fieldItems_[fieldi].threshold()<<";";
    }
    file<< endl;
            
    writeCommented(file, "Time");
    
    forAll(fieldItems_, fieldi )
    {
        const word& fieldName = fieldItems_[fieldi].fieldName();
        writeTabbed(file, "minLength(" + fieldName + ')');
        writeTabbed(file, "maxLength(" + fieldName + ')');
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
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    position_(dict.lookup("position")),
    direction_(dict.lookup("direction")),
    fieldItems_()
{
    // This is important, or "Requested single file, but multiple files are present"
    resetName(typeName);

    dict.lookup("fields") >> fieldItems_;

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

    forAll(fieldItems_, fieldi)
    {
        const word& fieldName = fieldItems_[fieldi].fieldName();
        const scalar& criterion = fieldItems_[fieldi].criterion();
        const scalar& threshold = fieldItems_[fieldi].threshold();
        const volScalarField& fieldValue = lookupObject<volScalarField>(fieldName);
        const scalar field_max = max(fieldValue).value();
        Info << "field_max = " << field_max << endl;
        scalar maxLength = 0;
        scalar minLength = great;

        if (field_max > threshold)
        {
            forAll (fieldValue, cellI)
            {
                if (fieldValue[cellI] >= criterion)
                {
                    vector raw = position_ - mesh_.C()[cellI];
                    scalar length = mag(raw&direction_);

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

        Log << "    minLength(" << fieldName << ") = " << minLength << nl
            << "    maxLength(" << fieldName << ") = " << maxLength << nl;
    }
    
    if (Pstream::master()) file() << endl;
    
    Log << endl;
    
    return true;
}


// ************************************************************************* //