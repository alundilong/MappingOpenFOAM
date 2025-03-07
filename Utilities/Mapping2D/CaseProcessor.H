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
    explicit CaseProcessor(
            const Foam::fileName &rootPath, 
            const Foam::fileName &caseName, 
            Foam::label nrow, 
            Foam::label ncolx,
            Foam::label ncoly
            );

    // do the dirty work
    void run(
            const Foam::pointField & pointCloud, 
            const Foam::wordList & fieldNamesx,
            const Foam::wordList & fieldNamesy
            );

private:
    Foam::Time runTime;
    const Foam::fvMesh mesh;
    Foam::scalarRectangularMatrix datax;
    Foam::scalarRectangularMatrix datay;
    const Foam::fileName casePath;
};

#endif // CASEPROCESSOR_H

