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

Class
    Foam::functionObjects::minMaxLength

Description
    Calculates the min and max distance between a cutting plane, determined by a position and a direction, to a criterion value of the specified scalar field.

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: fieldMinMax  | yes         |
        log          | write min/max data to standard output | no | no
        location     | write location of the min/max value | no | yes
        mode         | calculation mode: magnitude or component | no | magnitude
    \endtable
    Output data is written to the file \<timeDir\>/fieldMinMax.dat

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    minMaxLength.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_minMaxLength_H
#define functionObjects_minMaxLength_H

#include "fieldExpression.H"
#include "fvCFD.H"
#include "logFiles.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class minMaxLength Declaration
\*---------------------------------------------------------------------------*/

class minMaxLength
:
//    public fvMeshFunctionObject,
    public fieldExpression,
    public logFiles
{
    // Private Data
        vector position;
        vector direction;

        OFstream ResultOutPut;

        //- Fields to assess min/max
        wordList fieldSet_;

        scalarList criterion;
        scalarList threshold;

    // Private Member Functions

        //- Calculate the minMaxLength field and return true if successful
        virtual bool calc();
        //- Output file header information
        virtual void writeFileHeader(const label i);

public:

    //- Runtime type information
    TypeName("minMaxLength");


    // Constructors

        //- Construct from Time and dictionary
        minMaxLength
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~minMaxLength();


    // Member Functions

        //- Read the randomise data
        virtual bool read(const dictionary&);

        //- currently does nothing
        virtual bool execute();
        virtual bool write();
        virtual bool clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //