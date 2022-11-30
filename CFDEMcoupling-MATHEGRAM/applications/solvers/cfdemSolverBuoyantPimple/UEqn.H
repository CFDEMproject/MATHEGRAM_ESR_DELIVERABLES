    // Solve the Momentum equation


    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho*voidfraction, U)
      + fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + particleCloud.divVoidfractionTau(U, voidfraction)
      + fvm::Sp(Ksl, U)
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        if (modelType=="A")
	{
            solve
            (
                UEqn
             ==
                fvc::reconstruct
                (
                    (
                      - ghf*fvc::snGrad(rho*voidfraction)
                      - fvc::snGrad(p_rgh)
                    )*mesh.magSf()
                )*voidfraction
              + Ksl*Us
            );
        }
	else
	{
            solve
            (
                UEqn
             ==
                fvc::reconstruct
                (
                    (
                      - ghf*fvc::snGrad(rho*voidfraction)
                      - fvc::snGrad(p_rgh)
                    )*mesh.magSf()
                )
              + Ksl*Us
            );
	}

        fvOptions.correct(U);
        K = 0.5*magSqr(U);
    }
