#ifndef CASEPROCESSOR_H
#define CASEPROCESSOR_H

#include <string>
#include <vector>
#include <iostream>
#include "Time.H"
#include "fvMesh.H"
#include "fileName.H"

namespace Foam {
    class argList;  // Forward declaration
}

class CaseProcessor {
public:

    // Constructor
    explicit CaseProcessor(const Foam::fileName &rootPath, const Foam::fileName &caseName, Foam::label nrow, Foam::label ncol);

    // do the dirty work
    const Foam::scalarRectangularMatrix& run(const Foam::pointField & pointCloud, const Foam::wordList & fieldNames);

private:
    Foam::Time runTime;
    Foam::fvMesh mesh;
    Foam::scalarRectangularMatrix data;
};

#endif // CASEPROCESSOR_H

