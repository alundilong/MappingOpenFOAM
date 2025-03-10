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
        const Foam::fileName & caseName) :
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
    casePath(rootPath + '/' + caseName),
    ptrDatax(nullptr),
    ptrDatay(nullptr),
    ptrData(nullptr)
{
}

const XYPairData & CaseProcessor::run(
        const Foam::pointField & pointCloud, 
        const Foam::HashSet<Foam::word> & fieldNamesX,
        const Foam::HashSet<Foam::word> & fieldNamesY)
{
    //- looping over all timeDir
    const Foam::instantList &timeDirs = Foam::Time::findTimes(casePath);
    label nrow = (timeDirs.size()-1)*pointCloud.size();
    label ncolx = fieldNamesX.size() + 4; // Cx, Cy, Time, Ux and Uy
    label ncoly = fieldNamesY.size() + 1; // Ux and Uy

    ptrDatax = new Foam::scalarRectangularMatrix(nrow,ncolx);
    ptrDatay = new Foam::scalarRectangularMatrix(nrow,ncoly);

    Foam::scalarRectangularMatrix &datax = *ptrDatax;
    Foam::scalarRectangularMatrix &datay = *ptrDatay;

    ptrData = new XYPairData(datax,datay);
    XYPairData &data = *ptrData;

    //-load
    forAll(timeDirs, timei)
    {
        if (timei == 0) continue;
        runTime.setTime(timeDirs[timei], timei);
        Info<< "Time = " << runTime.timeName() << endl;

        //- load X
        if (runTime.timeName() == "0")
        {
            Foam::IOobjectList objects(mesh, runTime.timeName());

            LIFOStack<regIOobject*> storedObjects;
            readFields<volScalarField>(mesh, objects, fieldNamesX, storedObjects);

            labelList cellIds(mesh.nCells());
            forAll(pointCloud, pI)
            {
                cellIds[pI] = mesh.findCell(point(pointCloud[pI]));
            }

            //- add Cx, Cy, and time
            label row = 0;
            label colx = 0;
            forAll(cellIds, celli)
            {
                const vector & v = mesh.C().internalField()[celli];
                datax[(timei-1)*mesh.nCells()+row][colx] = v.x();
                datax[(timei-1)*mesh.nCells()+row][colx+1] = v.y();
                datax[(timei-1)*mesh.nCells()+row][colx+2] = runTime.value();
                row++;
            }

            colx = colx + 3;
            while (!storedObjects.empty())
            {
                volScalarField &fd = refCast<volScalarField>(*(storedObjects.pop()));
                label row = 0;
                forAll(cellIds, celli)
                {
                    datax[(timei-1)*mesh.nCells()+row][colx] = fd.internalField()[celli];
                    row++;
                }
                colx++;
            }

            readFields<volVectorField>(mesh, objects, fieldNamesX, storedObjects);
            while (!storedObjects.empty())
            {
                volVectorField &fd = refCast<volVectorField>(*(storedObjects.pop()));
                if (fd.name() == "U")
                {
                    label row = 0;
                    forAll(cellIds, celli)
                    {
                        const vector & v = fd.internalField()[celli];
                        datax[(timei-1)*mesh.nCells()+row][colx] = v.x();
                        datax[(timei-1)*mesh.nCells()+row][colx+1] = v.y();
                        row++;
                    }
                }
            }
        }
        else
        {
            //- repeat except time channel
            for(int i = 0; i < pointCloud.size(); i++)
                for(int j = 0; j < ncolx; j++)
                {
                    if (j == 2)
                    {
                        datax[(timei-1)*mesh.nCells()+i][j] = runTime.value();
                    }
                    else
                    {
                        datax[(timei-1)*mesh.nCells()+i][j] = datax[i][j];
                    }
                }
        }

        //- load Y
        {
            Foam::IOobjectList objects(mesh, runTime.timeName());

            LIFOStack<regIOobject*> storedObjects;
            readFields<volScalarField>(mesh, objects, fieldNamesY, storedObjects);

            labelList cellIds(mesh.nCells());
            forAll(pointCloud, pI)
            {
                cellIds[pI] = mesh.findCell(point(pointCloud[pI]));
            }

            label coly = 0;
            while (!storedObjects.empty())
            {
                volScalarField &fd = refCast<volScalarField>(*(storedObjects.pop()));
                label row = 0;
                forAll(cellIds, celli)
                {
                    datay[(timei-1)*mesh.nCells()+row][coly] = fd.internalField()[celli];
                    row++;
                }
                coly++;
            }

            readFields<volVectorField>(mesh, objects, fieldNamesY, storedObjects);
            while (!storedObjects.empty())
            {
                volVectorField &fd = refCast<volVectorField>(*(storedObjects.pop()));
                if (fd.name() == "U")
                {
                    label row = 0;
                    forAll(cellIds, celli)
                    {
                        const vector & v = fd.internalField()[celli];
                        datay[(timei-1)*mesh.nCells()+row][coly] = v.x();
                        datay[(timei-1)*mesh.nCells()+row][coly+1] = v.y();
                        row++;
                    }
                }
            }
        }
    }
    return data;
}
