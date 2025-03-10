#ifndef CASEPROCESSOR_H
#define CASEPROCESSOR_H

#include <string>
#include <vector>
#include <iostream>
#include "Time.H"
#include "fvMesh.H"
#include "fileName.H"
#include "HashSet.H"
#include "word.H"

namespace Foam {
    class argList;  // Forward declaration
}

// X refers to input
// Y refers to output
class XYPairData {
public:
    explicit XYPairData(Foam::scalarRectangularMatrix &x, Foam::scalarRectangularMatrix &y):
        datax(x),datay(y)
    {
    }
    inline const Foam::scalarRectangularMatrix& x() const { return datax; }
    inline const Foam::scalarRectangularMatrix& y() const { return datay; }
private:
    Foam::scalarRectangularMatrix &datax;
    Foam::scalarRectangularMatrix &datay;
};

class CaseProcessor {
public:

    // Constructor
    explicit CaseProcessor(
                const Foam::fileName &rootPath, 
                const Foam::fileName &caseName
            );

    // Destructor
    ~CaseProcessor()
    {
        if (ptrDatax != nullptr)
        {
            delete ptrDatax;
        }
        if (ptrDatay != nullptr)
        {
            delete ptrDatay;
        }
        if (ptrData != nullptr)
        {
            delete ptrData;
        }
    }

    // do the dirty work
    const XYPairData & run(
                const Foam::pointField & pointCloud, 
                const Foam::HashSet<Foam::word> & fieldNamesX,
                const Foam::HashSet<Foam::word> & fieldNamesY
            );

private:
    Foam::Time runTime;
    Foam::fvMesh mesh;
    const Foam::fileName casePath;
    Foam::scalarRectangularMatrix *ptrDatax;
    Foam::scalarRectangularMatrix *ptrDatay;
    XYPairData *ptrData;
};

#endif // CASEPROCESSOR_H

