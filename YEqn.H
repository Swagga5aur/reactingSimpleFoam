tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    reaction->correct();
    Qdot = reaction->Qdot();
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
     	if (composition.solve(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->muEff(), Yi)
             ==
                reaction->R(Yi)
              + fvModels.source(rho, Yi)
            );

            YiEqn.relax();

            fvConstraints.constrain(YiEqn);

            YiEqn.solve("Yi");

            fvConstraints.constrain(Yi);

            Yi.max(0.0);
	    Yi.min(1.0);
            Yt += Yi;
        }
    }
    composition.normalise();
}
