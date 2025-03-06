IOdictionary mapping2DDict
(
    IOobject
    (
         "mapping2DDict",
         runTime.constant(),
         mesh,
         IOobject::MUST_READ_IF_MODIFIED,
         IOobject::NO_WRITE
    )
);

const dictionary & subDictGrids = mapping2DDict.subDict("Grids");

dimensionedScalar scale
(
    subDictGrids.lookup("scale")
);

dimensionedScalar xLow
(
    subDictGrids.lookup("xLow")*scale
);

dimensionedScalar xHig
(
    subDictGrids.lookup("xHig")*scale
);

dimensionedScalar dx
(
    subDictGrids.lookup("dx")*scale
);

dimensionedScalar yLow
(
    subDictGrids.lookup("yLow")*scale
);

dimensionedScalar yHig
(
    subDictGrids.lookup("yHig")*scale
);

dimensionedScalar dy
(
    subDictGrids.lookup("dy")*scale
);

dimensionedScalar zLow
(
    subDictGrids.lookup("zLow")*scale
);

dimensionedScalar zHig
(
    subDictGrids.lookup("zHig")*scale
);

dimensionedScalar xL = xHig - xLow;
dimensionedScalar yL = yHig - yLow;
dimensionedScalar zL = zHig - zLow;

label Nx = (xL/dx).value();
label Ny = (yL/dy).value();

label nTotal = Nx*Ny;
pointField pClouds(nTotal);

scalar xStart = 0.5*dx.value();
scalar yStart = 0.5*dy.value();
scalar zMiddle = 0.5*(zHig+zLow).value();

for (int j = 0; j < Ny; j++)
{
    for (int i = 0; i < Nx; i++)
    {
        scalar x = xStart + i*dx.value();
        scalar y = yStart + j*dy.value();
        label pointI = i + Nx*j;
        pClouds[pointI] = point(x,y,zMiddle);
    }
}
Info << pClouds << endl;
