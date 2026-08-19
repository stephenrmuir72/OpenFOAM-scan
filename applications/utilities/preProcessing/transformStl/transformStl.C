/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Kuralamuthan Veerapandiyan, IIT Madras
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

Application
    transformStl

Description
    Utility to transform STL surface geometries.
    Supports writing all solids in a single STL or separate transformed
       STL file(s).

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fileName.H"
#include "IFstream.H"
#include "List.H"
#include "Time.H"
#include "vector.H"
#include "word.H"

#include "IOdictionary.H"
#include "mathematicalConstants.H"

#include "OFstream.H"
#include "OSspecific.H"

#include <sstream>
#include <string>

using namespace Foam;
using namespace Foam::constant::mathematical;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Struct for ASCII STL triangle
struct stlFacet
{
    vector normal;
    FixedList<vector, 3> vertices;

    stlFacet()
    :
	    normal(vector::zero),
	    vertices(vector::zero)
    {}
};

//- Struct for geomtry associated with single 'solid ... endsolid' STL block
struct stlSolid
{
    word name;
    List<stlFacet> facets;

    stlSolid()
    :
	    name(word::null),
	    facets()
    {}
};

//- Rigid-body transformation assigned to single STL solid
struct solidTransform
{
    vector translation;
    vector rotationAxis;
    vector rotationCentre;
    scalar rotationAngle;

    solidTransform()
    :
        translation(vector::zero),
	rotationAxis(0,0,1),
	rotationCentre(vector::zero),
	rotationAngle(0.0)
    {}
};


//- Read a scalar triplet from the remaining stream contents
bool readVector
(
    std::istringstream& stream,
    vector& value
)
{
    scalar x = 0.0;
    scalar y = 0.0;
    scalar z = 0.0;

    if (!(stream >> x >> y >> z))
    {
        return false;
    }

    value = vector(x, y, z);
    return true;
}

//- Return true when line contains only whitespace
bool emptyLine(const std::string& line)
{
    return line.find_first_not_of(" \t\r\n") == std::string::npos;
}

//- Read an ACSII multi-solid file
List<stlSolid> readAsciiStl(const fileName& stlFileName)
{
    IFstream stlFile(stlFileName);
    if (!stlFile.good())
    {
        FatalErrorInFunction
		<< "Cannot open STL file: " << stlFileName
		<< exit(FatalError);
    }

    DynamicList<stlSolid> solids;

    std::string line;
    label lineNumber = 0;

    bool inSolid  = false;
    bool inFacet  = false;
    label vertexI = 0;

    stlSolid currentSolid;
    stlFacet currentFacet;

    while (std::getline(stlFile.stdStream(), line))
    {
        ++lineNumber;

	if (emptyLine(line))
	{
	    continue;
	}

	std::istringstream lineStream(line);
	std::string keyword;
	lineStream >> keyword;

	if (keyword == "solid")
	{
            if (inSolid)
	    {
		FatalErrorInFunction
		    << "Encountered nested 'solid' declaration at line "
		    << lineNumber << " in file " << stlFileName << nl
		    << "Expected 'endsolid " << currentSolid.name
		    << "' before the next 'solid' declaration."
		    << exit(FatalError);
	    }
	    std::string solidName;
	    lineStream >> solidName;

	    if (solidName.empty())
	    {
                solidName = "solid_" + Foam::name(solids.size());
	    }
	    currentSolid = stlSolid();
	    currentSolid.name = word(solidName);

    	    inSolid = true;
	    inFacet = false;
	    vertexI = 0;
        }

	else if (keyword == "facet")
	{
            if (!inSolid)
	    {
                FatalErrorInFunction
		    << "Found 'facet' outside a solid block at line "
		    << lineNumber << " in file " << stlFileName
		    << exit(FatalError);
	    }

	    if (inFacet)
	    {
                FatalErrorInFunction
		    << "Found a new 'facet' before completing the previous "
		    << "facet at line " << lineNumber << " in file "
		    << stlFileName << exit(FatalError);
	    }

	    std::string normalKeyword;
	    lineStream >> normalKeyword;

	    if (normalKeyword != "normal")
	    {
                FatalErrorInFunction
		    << "Expected 'facet normal' at line " << lineNumber
		    << " in file " << stlFileName
		    << exit(FatalError);
	    }
	    currentFacet = stlFacet();

	    if (!readVector(lineStream, currentFacet.normal))
	    {
                FatalErrorInFunction
		    << "Cannot read facet normal at line " << lineNumber
		    << " in file " << stlFileName
		    << exit(FatalError);
	    }
	    inFacet = true;
	    vertexI = 0;
        }

	else if (keyword == "vertex")
	{
            if (!inFacet)
	    {
                FatalErrorInFunction
		    << "Found 'vertex' outside a facet block at line "
		    << lineNumber << " in file " << stlFileName
		    << exit(FatalError);
	    }

	    if (vertexI >= 3)
	    {
                FatalErrorInFunction
		    << "Facet has more than three vertices at line "
		    << lineNumber << " in file " << stlFileName
		    << exit(FatalError);
	    }


	    if (!readVector(lineStream, currentFacet.vertices[vertexI]))
	    {
                FatalErrorInFunction
		    << "Cannot read vertex coordinates at line "
		    << lineNumber << " in file " << stlFileName
		    << exit(FatalError);
	    }
            ++vertexI;
        }

	else if (keyword == "endfacet")
	{
            if (!inFacet)
	    {
                FatalErrorInFunction
		<< "Found 'endfacet' without a matching 'facet' at line "
		<< lineNumber << " in file " << stlFileName
		<< exit(FatalError);
	    }

	    if (vertexI != 3)
	    {
                FatalErrorInFunction
		    << "Facet at line " << lineNumber
		    << " does not contain exactly three vertices."
		    << exit(FatalError);
	    }
	    currentSolid.facets.append(currentFacet);
	    inFacet = false;
	    vertexI = 0;
	}

	else if (keyword == "endsolid")
	{
            if (!inSolid)
	    {
                FatalErrorInFunction
		    << "Found 'endsolid' without a matching 'solid' at line "
		    << lineNumber << " in file " << stlFileName
		    << exit(FatalError);
	    }

	    if (inFacet)
	    {
                FatalErrorInFunction
		    << "Solid " << currentSolid.name
		    << " ends with an incomplete facet at line "
		    << lineNumber
		    << exit(FatalError);
	    }
	    solids.append(currentSolid);
	    inSolid = false;
	    vertexI = 0;
	}
    }

    if (inFacet)
    {
        FatalErrorInFunction
	    << "Reached end of file with an incomplete facet in "
	    << stlFileName 
	    << exit(FatalError);
    }

    if (inSolid)
    {
        FatalErrorInFunction
	    << "Reached end of file without 'endsolid' "
	    << currentSolid.name << " in " << stlFileName
	    << exit(FatalError);
    }

    if (solids.empty())
    {
        FatalErrorInFunction
	    << "No STL solids found in " << stlFileName << nl
	    << "transformStl currently supports ASCII STL files with "
	    << "'solid' ... 'endsolid' blocks"
	    << exit(FatalError);
    }
    return List<stlSolid>(solids);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /

//- Rotate the given vector
vector rotateVector
(
    const vector& value,
    const vector& axis,
    const scalar angleRad
)
{
    if(mag(angleRad) < SMALL)
    {
        return value;
    }
    const scalar axisMagnitude = mag(axis);

    if (axisMagnitude < SMALL)
    {
        FatalErrorInFunction
            << "Cannot apply a non-zero rotation with rotationAxis "
            << axis << nl
            << "Specify a non-zero rotationAxis."
            << exit(FatalError);
    }

    const vector e(axis/axisMagnitude);

    const scalar c = Foam::cos(angleRad);
    const scalar s = Foam::sin(angleRad);

    return
        c*value
      + s*(e ^ value)
      + (1.0 - c)*(e & value)*e;
}


//- Read the transformation settings for one STL solid
solidTransform readTransform
(
    const dictionary& transformDict,
    const word& solidName
)
{
    solidTransform transform;
    if (!transformDict.found(solidName))
    {
        return transform;
    }

    const dictionary& solidDict(transformDict.subDict(solidName));
    transform.translation    =
	solidDict.lookupOrDefault<vector>("translation", vector::zero);
    transform.rotationAxis   =
	solidDict.lookupOrDefault<vector>("rotationAxis",vector(0, 0, 1));
    transform.rotationCentre =
	solidDict.lookupOrDefault<vector>("rotationCentre",vector::zero);
    const scalar angleDegrees =
	solidDict.lookupOrDefault<scalar>("rotationAngle",0.0);
    transform.rotationAngle = angleDegrees*pi/180.0;

    return transform;
}

//- Apply a rigid-body transform to all facets of single solid
//- Apply a rigid-body transform to all facets belonging to one solid
void transformSolid
(
    stlSolid& solid,
    const solidTransform& transform
)
{
    forAll(solid.facets, facetI)
    {
        stlFacet& facet = solid.facets[facetI];
        forAll(facet.vertices, vertexI)
        {
            const vector relativePosition
            (
                facet.vertices[vertexI] - transform.rotationCentre
            );

            facet.vertices[vertexI] =
                transform.rotationCentre
              + rotateVector
                (
                    relativePosition,
                    transform.rotationAxis,
                    transform.rotationAngle
                )
              + transform.translation;
        }

        facet.normal =
            rotateVector
            (
                facet.normal,
                transform.rotationAxis,
                transform.rotationAngle
            );

        if (mag(facet.normal) > SMALL)
        {
            facet.normal /= mag(facet.normal);
        }
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /

void writeSolid
(
    OFstream& outputFile,
    const stlSolid& solid
)
{
    outputFile << "solid " << solid.name << nl;
    forAll(solid.facets, facetI)
    {
        const stlFacet& facet = solid.facets[facetI];
	outputFile
	    << "  facet normal "
	    << facet.normal.x() << ' '
	    << facet.normal.y() << ' '
	    << facet.normal.z() << nl;
	outputFile << "    outer loop" << nl;

	forAll(facet.vertices, vertexI)
	{
            const vector& vertex = facet.vertices[vertexI];
	    outputFile
		<< "        vertex "
		<< vertex.x() << ' '
		<< vertex.y() << ' '
		<< vertex.z() << nl;
	}
	outputFile << "    endloop" << nl;
	outputFile << "  endfacet" << nl;
    }
    outputFile << "endsolid " << solid.name << nl;
}

//- Write all solids into one ASCII multi-solid STL file
void writeStl
(
    const fileName& outputFileName,
    const List<stlSolid>& solids
)
{
    OFstream outputFile(outputFileName);
    if(!outputFile.good())
    {
        FatalErrorInFunction
	    << "Cannot open output STL file: " << outputFileName
	    << exit(FatalError);
    }

    forAll(solids, solidI)
    {
        writeSolid(outputFile, solids[solidI]);
    }
    Info<< "Wrote multi-solid STL: " << outputFileName << nl;
}

//-Write every solid into an individual ASCII STL file
void writeSeparateStls
(
    const fileName& outputDirectory,
    const List<stlSolid>& solids
)
{
    forAll(solids, solidI)
    {
        const stlSolid& solid = solids[solidI];
	const fileName outputFileName
	(
            outputDirectory/(solid.name + ".stl")
	);
	OFstream outputFile(outputFileName);

	if (!outputFile.good())
	{
            FatalErrorInFunction
		<< "Cannot open output STL file: " << outputFileName
		<< exit(FatalError);
	}

	writeSolid(outputFile, solid);
	Info<< "Wrote STL for solid " << solid.name
	    << " : " << outputFileName << nl; 
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Read an ASCII single- or multi-solid STL file.\n"
	"This initial implementation validates and reports STL facets."
    );

    argList::addArgument
    (
        "stlFile",
	"ASCII STL input file"
    );

    argList::addOption
    (
        "dict",
	"file",
	"Alternative transformation dictionary",
	"(default: constant/transformStlDict)"
    );

    argList::addOption
    (
        "output",
	"file",
	"Output STL filename. Relative paths are place in "
	"constant/triSurface. Default: transformed.stl"
    );

    argList::addBoolOption
    (
        "separateSTLs",
	"Write one STL file for each of the input solid"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    fileName stlFileName(args[1]);

    if(!stlFileName.isAbsolute())
    {
        stlFileName = runTime.constant()/stlFileName;
    }

    Info << "Reading ASCII STL: " << stlFileName << nl << endl;
    List<stlSolid> solids(readAsciiStl(stlFileName));
    IOdictionary transformStlDict
    (
        IOobject
	(
            args.getOrDefault<word>("dict", "transformStlDict"),
	    runTime.constant(),
	    runTime,
	    IOobject::MUST_READ,
	    IOobject::NO_WRITE
	)
    );

    const dictionary& transformsDict
    (
        transformStlDict.subDict("transforms")
    );

    forAll(solids, solidI)
    {
        stlSolid& solid = solids[solidI];
        const solidTransform transform
	(
	    readTransform(transformsDict, solid.name)
	);
        Info<< "Transforming solid : " << solid.name << nl
	    << "    translation    : " << transform.translation    << nl
	    << "    rotationAxis   : " << transform.rotationAxis   << nl
	    << "    rotationCentre : " << transform.rotationCentre << nl
	    << "    rotationAngle  : " << transform.rotationAngle*180.0/pi
	    << " degrees" << nl << nl;
	transformSolid(solid, transform);
    }

    const bool separateSTLs = args.found("separateSTLs");
    if (separateSTLs && args.found("output"))
    {
        FatalErrorInFunction
	    << "Options '-output' and '-separateSTLs' cannot be used "
	    << "together."
	    << exit(FatalError);
    }

    const fileName triSurfaceDir(runTime.constant()/"triSurface");
    if (!isDir(triSurfaceDir))
    {
        Info<< "Creating directory: " << triSurfaceDir << nl;
	mkDir(triSurfaceDir);
    }

    if (separateSTLs)
    {
        writeSeparateStls(triSurfaceDir, solids);
    }
    else
    {
        fileName outputFileName
	(
            args.getOrDefault<fileName>("output", "transformed.stl")
	);
	if (!outputFileName.isAbsolute())
	{
            outputFileName = triSurfaceDir/outputFileName;
	}
	writeStl(outputFileName, solids);
    }


    label totalFacets = 0;
    forAll(solids, solidI)
    {
        const stlSolid& solid = solids[solidI];
        Info<< "Solid: " << solid.name
	    << ", facets: " << solid.facets.size() << nl;
	totalFacets += solid.facets.size();
    }

    Info<< nl
        << "Processed " << solids.size() << " solids(s) and "
        << totalFacets << " facet(s)." << nl
        << endl;

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /


