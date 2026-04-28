/*---------------------------------------------------------------------------*\

License
    This file is part of GeoChemFoam, an Open source software using OpenFOAM
    for multiphase multicomponent reactive transport simulation in pore-scale
    geological domain.

    GeoChemFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version. See <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "solidGrowthVelocityPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "polyMesh.H"
#include "primitivePatchInterpolation.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidGrowthVelocityPointPatchVectorField::
    solidGrowthVelocityPointPatchVectorField(
        const pointPatch &p, const DimensionedField<vector, pointMesh> &iF)
    : fixedValuePointPatchField<vector>(p, iF), speciesName_("X_S"),
      UName_("U"), growthCoeff_(1e-5), rhoSolid_(1100.0), Mw_(1.0),
      tauCrit_(1.0), kDetach_(1e-7), detachExponent_(1.0), mu_(1e-3),
      maxVelocity_(1e-6) {}

solidGrowthVelocityPointPatchVectorField::
    solidGrowthVelocityPointPatchVectorField(
        const pointPatch &p, const DimensionedField<vector, pointMesh> &iF,
        const dictionary &dict)
    : fixedValuePointPatchField<vector>(p, iF, dict),
      speciesName_(dict.getOrDefault<word>("speciesName", "X_S")),
      UName_(dict.getOrDefault<word>("UName", "U")),
      growthCoeff_(dict.get<scalar>("growthCoeff")),
      rhoSolid_(dict.get<scalar>("rhoSolid")),
      Mw_(dict.getOrDefault<scalar>("Mw", 1.0)),
      tauCrit_(dict.get<scalar>("tauCrit")),
      kDetach_(dict.getOrDefault<scalar>("kDetach", 1e-7)),
      detachExponent_(dict.getOrDefault<scalar>("detachExponent", 1.0)),
      mu_(dict.get<scalar>("mu")),
      maxVelocity_(dict.getOrDefault<scalar>("maxVelocity", 1e-6)) {
  if (!dict.found("value")) {
    updateCoeffs();
  }
}

solidGrowthVelocityPointPatchVectorField::
    solidGrowthVelocityPointPatchVectorField(
        const solidGrowthVelocityPointPatchVectorField &ptf,
        const pointPatch &p, const DimensionedField<vector, pointMesh> &iF,
        const pointPatchFieldMapper &mapper)
    : fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
      speciesName_(ptf.speciesName_), UName_(ptf.UName_),
      growthCoeff_(ptf.growthCoeff_), rhoSolid_(ptf.rhoSolid_), Mw_(ptf.Mw_),
      tauCrit_(ptf.tauCrit_), kDetach_(ptf.kDetach_),
      detachExponent_(ptf.detachExponent_), mu_(ptf.mu_),
      maxVelocity_(ptf.maxVelocity_) {}

solidGrowthVelocityPointPatchVectorField::
    solidGrowthVelocityPointPatchVectorField(
        const solidGrowthVelocityPointPatchVectorField &ptf,
        const DimensionedField<vector, pointMesh> &iF)
    : fixedValuePointPatchField<vector>(ptf, iF),
      speciesName_(ptf.speciesName_), UName_(ptf.UName_),
      growthCoeff_(ptf.growthCoeff_), rhoSolid_(ptf.rhoSolid_), Mw_(ptf.Mw_),
      tauCrit_(ptf.tauCrit_), kDetach_(ptf.kDetach_),
      detachExponent_(ptf.detachExponent_), mu_(ptf.mu_),
      maxVelocity_(ptf.maxVelocity_) {}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidGrowthVelocityPointPatchVectorField::updateCoeffs() {
  if (this->updated()) {
    return;
  }

  const polyMesh &mesh = this->internalField().mesh()();
  const Time &runTime = mesh.time();

  // Get patch normal vectors (points outward from solid into fluid)
  vectorField n = patch().pointNormals();

  // Set up face-to-point interpolator
  primitivePatchInterpolation patchInterpolator(
      mesh.boundaryMesh()[patch().index()]);

  // === GROWTH / DISSOLUTION TERM ===
  // Get solid-phase species field and compute rate of change dC/dt
  const volScalarField &C =
      this->db().objectRegistry::lookupObject<volScalarField>(speciesName_);
  const fvPatchField<scalar> &Cp = C.boundaryField()[patch().index()];

  // Get old time value for rate calculation
  scalarField Cp_old = Cp;
  if (C.nOldTimes() > 0) {
    Cp_old = C.oldTime().boundaryField()[patch().index()];
  }

  // Rate of change: dC/dt (mol/m^3/s)
  scalar dt = runTime.deltaTValue();
  scalarField dCdt = (Cp - Cp_old) / max(dt, SMALL);

  // Growth velocity: v_growth = growthCoeff * dC/dt * Mw / rhoSolid (m/s)
  // Positive dC/dt -> precipitation/growth (mesh moves inward, clogging)
  // Negative dC/dt -> dissolution (mesh moves outward, pore opening)
  scalarField v_growth = growthCoeff_ * dCdt * Mw_ / rhoSolid_;

  // Interpolate to points
  scalarField v_growth_points =
      patchInterpolator.faceToPointInterpolate<scalar>(v_growth);

  // === DETACHMENT / EROSION TERM (optional, active when kDetach > 0) ===
  scalarField v_detach_points(v_growth_points.size(), scalar(0));

  if (kDetach_ > SMALL) {
    // Get velocity field for shear stress calculation
    const volVectorField &U =
        this->db().objectRegistry::lookupObject<volVectorField>(UName_);
    const fvPatchField<vector> &Up = U.boundaryField()[patch().index()];

    // Wall shear stress: tau = mu * |dU/dn|
    scalarField tauWall = mu_ * mag(Up.snGrad());

    // Detachment velocity: v_detach = kDetach * (tau/tauCrit)^m
    scalarField tauRatio = tauWall / max(tauCrit_, SMALL);
    scalarField v_detach =
        kDetach_ * pow(max(tauRatio, scalar(0)), detachExponent_);

    // Interpolate to points
    v_detach_points =
        patchInterpolator.faceToPointInterpolate<scalar>(v_detach);
  }

  // === NET VELOCITY ===
  // v_mesh = -n * (v_growth - v_detach)
  // Positive v_growth -> mesh moves inward (clogging)
  // Positive v_detach -> mesh moves outward (erosion)
  scalarField v_net = v_growth_points - v_detach_points;

  // Apply velocity limiter
  v_net = max(min(v_net, maxVelocity_), -maxVelocity_);

  // Set mesh velocity (negative n = inward growth)
  Field<vector>::operator=(-n *v_net);

  fixedValuePointPatchField<vector>::updateCoeffs();
}

void solidGrowthVelocityPointPatchVectorField::write(Ostream &os) const {
  pointPatchField<vector>::write(os);
  os.writeEntry("speciesName", speciesName_);
  os.writeEntry("UName", UName_);
  os.writeEntry("growthCoeff", growthCoeff_);
  os.writeEntry("rhoSolid", rhoSolid_);
  os.writeEntry("Mw", Mw_);
  os.writeEntry("tauCrit", tauCrit_);
  os.writeEntry("kDetach", kDetach_);
  os.writeEntry("detachExponent", detachExponent_);
  os.writeEntry("mu", mu_);
  os.writeEntry("maxVelocity", maxVelocity_);
  writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField(pointPatchVectorField,
                        solidGrowthVelocityPointPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
