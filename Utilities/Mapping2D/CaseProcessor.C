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
        const Foam::HashSet<Foam::word> & fieldNames)
{
    //- looping over all timeDir
    const Foam::instantList &timeDirs = Foam::Time::findTimes(casePath);

    forAll(timeDirs, timei)
    {
        if (timei == 0) continue;
        runTime.setTime(timeDirs[timei], timei);

        Foam::IOobjectList objects(mesh, runTime.timeName());

        LIFOStack<regIOobject*> storedObjects;
        readFields<volScalarField>(mesh, objects, fieldNames, storedObjects);
        Info<< "Time = " << runTime.timeName() <<" stored object size: "<< storedObjects.size() << endl;

        while (!storedObjects.empty())
        {
            volScalarField &fd = refCast<volScalarField>(*(storedObjects.pop()));
            Info << fd.name() << tab << average(fd).value() << endl;
        }

        readFields<volVectorField>(mesh, objects, fieldNames, storedObjects);
        volVectorField &U = refCast<volVectorField>(*(storedObjects.pop()));
        volScalarField magU = mag(U);
        Info << average(magU).value() << endl;
    }

}
