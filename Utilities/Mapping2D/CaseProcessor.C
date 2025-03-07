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
        Info<< "Time = " << runTime.timeName() << endl;

        Foam::IOobjectList objects(mesh, runTime.timeName());

        //IOobjectList fields(objects.lookupClass(volScalarField::typeName));

        //forAllConstIter(IOobjectList, fields, fieldIter)
        //{
        //    const IOobject& io = *fieldIter();
        //    Info << io.name() << tab << io.instance() << endl;
        //}

        LIFOStack<regIOobject*> storedObjects;
        readFields<volScalarField>(mesh, objects, fieldNames, storedObjects);

        const volScalarField &T = mesh.lookupObject<volScalarField>("T");
        Info << T.mesh().time().value() << tab << average(T).value() << nl;
        // readFields<volVectorField>(mesh, objects, fieldNames, storedObjects);
        // const volVectorField &U = mesh.lookupObject<volVectorField>("U");
        // volScalarField magU = mag(U);
        // Info << magU << endl;
    }

}
