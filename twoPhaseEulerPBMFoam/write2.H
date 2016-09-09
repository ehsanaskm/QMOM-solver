    if (runTime.outputTime())
    {
        volVectorField f0
        (
            IOobject
            (
                "f0",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            f0
        );

        runTime.write();
    }
