#include "CaseProcessor.H"
#include <fstream>
#include <sstream>
#include <string>
#include "fvCFD.H"
#include "IOobjectList.H"
#include "ReadFields.H"

// Constructor
CaseProcessor::CaseProcessor(
        const Foam::fileName& rootPath, 
        const Foam::fileName & caseName, 
        Foam::label nrow, 
        Foam::label ncolx,
        Foam::label ncoly) :
    runTime(Foam::Time::controlDictName, rootPath, caseName),
    mesh 
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    ),
    datax(nrow, ncolx),
    datay(nrow, ncoly),
    casePath(rootPath + '/' + caseName)
{
}

void CaseProcessor::run(
        const Foam::pointField & pointCloud, 
        const Foam::wordList & fieldNamesx,
        const Foam::wordList & fieldNamesy)
{
    //- looping over all timeDir
    const Foam::instantList &timeDirs = Foam::Time::findTimes(casePath);

    Foam::IOobjectList objects(mesh, runTime.timeName());
    Info << objects.size() << endl;
    Info << volScalarField::typeName << endl;
    Foam::IOobjectList fields = objects.lookupClass(volScalarField::typeName);
    Info << fields.size() << endl;

    wordList selectedTimes(timeDirs.size()-1);

    forAll(timeDirs, timei)
    {
        if (timei == 0) continue;
        selectedTimes[timei-1] = timeDirs[timei].name();
    }

    word registryName = "myCache";
    ReadFields<volScalarField>("T",mesh,selectedTimes,registryName);
    const objectRegistry & objectR = mesh.lookupObject<objectRegistry>(registryName);
    Info << objectR.names() << nl;
    const volScalarField &T = objectR.lookupObject<volScalarField>("T");

    forAll(timeDirs, timei)
    {
        if (timei == 0) continue;
        runTime.setTime(timeDirs[timei], timei);
        Info<< "Time = " << runTime.timeName() << endl;
        ReadFields<volScalarField>("T",mesh,selectedTimes,registryName);
    }

}
