{
    const volScalarField& psi = thermo.psi();

    const volScalarField rAU("rAU", 1.0/UEqn.A());
    const surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));

    tmp<volScalarField> rAtU
    (
        simple.consistent()
      ? volScalarField::New("rAtU", 1.0/(1.0/rAU - UEqn.H1()))
      : tmp<volScalarField>(nullptr)
    );
    tmp<surfaceScalarField> rhorAtUf
    (
        simple.consistent()
      ? surfaceScalarField::New("rhoRAtUf", fvc::interpolate(rho*rAtU()))
      : tmp<surfaceScalarField>(nullptr)
    );

    const volScalarField& rAAtU = simple.consistent() ? rAtU() : rAU;
    const surfaceScalarField& rhorAAtUf =
        simple.consistent() ? rhorAtUf() : rhorAUf;

    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

    tUEqn.clear();

    bool closedVolume = false;

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(rho)*fvc::flux(HbyA)
    );

    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    closedVolume = adjustPhi(phiHbyA, U, p);

    while (simple.correctNonOrthogonal())
    {
        tmp<fvScalarMatrix> tpEqn;
        {
            tpEqn =
            (
                fvm::laplacian(rhorAAtUf, p)
              + fvModels.source(psi, p, rho.name())
             ==
                fvc::div(phiHbyA)
            );
        }

        fvScalarMatrix& pEqn = tpEqn.ref();

        pEqn.setReference
        (
            pressureReference.refCell(),
            pressureReference.refValue()
        );

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    #include "incompressible/continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    U = HbyA - rAAtU*fvc::grad(p);

    U.correctBoundaryConditions();
    fvConstraints.constrain(U);

    fvConstraints.constrain(p);

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    if (closedVolume)
    {
        p += (initialMass - fvc::domainIntegrate(psi*p))
            /fvc::domainIntegrate(psi);
    }

    rho = thermo.rho();
    rho.relax();
}
