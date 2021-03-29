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
        type         | type name: minMaxLength  | yes         |
        position     | The point on the cutting plane | yes | -
        direction    | The normal vector of the cutting plane | yes | -        
        criterion    | The criterion to determine the max and min length | yes | -
        threshold    | When the max value of the processed fields is less than the threshold, skip the time step. | no | 0.01*criterion
    \endtable
    Output data is written to the file \<timeDir\>/minMaxLength.dat
    
    This is an example:
    
    \verbatim
    {
        type minMaxLength;
        libs ("libSJfieldFunctionObjects.so");
        position        ( 0 0 0 );//The point on the cutting plane
        direction       ( 0 0 1 );//The normal vector of the cutting plane
        fields
        (
            p   //The fileds to be processed
            {
                criterion       10000; // The criterion to determine the max and min length
                threshold       1000;  // When the max value of the processed fields is less than the threshold, skip the time step.
            }
        );

    }
    \endverbatim

SourceFiles
    minMaxLength.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_minMaxLength_H
#define functionObjects_minMaxLength_H

#include "fvMeshFunctionObject.H"
#include "vector.H"
#include "logFiles.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

// Forward declaration of friend functions and operators
class fieldItem;
Istream& operator>>(Istream&, fieldItem&);

/*---------------------------------------------------------------------------*\
                          Class fieldItem Declaration
\*---------------------------------------------------------------------------*/


class fieldItem
{
    
public:
    
    word fieldName_;
    
    scalar criterion_;
    
    scalar threshold_;

    //- Construct null
    fieldItem();

    //- Destructor
    ~fieldItem();

    // Member Functions

    // Access
    //- Return const access to the field name
    const word& fieldName() const
    {
        return fieldName_;
    }
    const scalar& criterion() const
    {
        return criterion_;
    }
    const scalar& threshold() const
    {
        return threshold_;
    }
    
    // IOstream Operators
    friend Istream& operator>>(Istream&, fieldItem&);
};

/*---------------------------------------------------------------------------*\
                          Class minMaxLength Declaration
\*---------------------------------------------------------------------------*/

class minMaxLength
:
    public fvMeshFunctionObject,
    public logFiles
{

private:
    // Private Data
        vector position_;
        
        vector direction_;

        List<fieldItem> fieldItems_;

    // Private Member Functions

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