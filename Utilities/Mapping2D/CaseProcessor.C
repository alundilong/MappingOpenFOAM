#include "CaseProcessor.H"
#include <fstream>
#include <sstream>
#include <string>
#include "fvCFD.H"

// Constructor
CaseProcessor::CaseProcessor(const Foam::fileName& rootPath, const Foam::fileName & caseName, Foam::label nrow, Foam::label ncol) :
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
    data(nrow, ncol)
{
}

const Foam::scalarRectangularMatrix& CaseProcessor::run(const Foam::pointField & pointCloud, const Foam::wordList & fieldNames)
{
    return data;
}
