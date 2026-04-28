/*---------------------------------------------------------------------------*\

License
    This file is part of GeoChemFoam, an Open source software using OpenFOAM
    for multiphase multicomponent reactive transport simulation in pore-scale
    geological domain.

    GeoChemFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version. See <http://www.gnu.org/licenses/>.

Class
    Foam::solidGrowthVelocityPointPatchVectorField

Description
    Boundary condition for mesh point velocity based on solid-phase growth
    or dissolution rate, with optional shear-induced detachment.

    Applicable to biofilm growth, mineral precipitation, and dissolution.

    Net mesh velocity = -n * (v_growth - v_detach)

    Where:
    - v_growth = growthCoeff * dC/dt * Mw / rhoSolid
      (positive dC/dt = precipitation/growth -> clogging)
      (negative dC/dt = dissolution -> pore opening)
    - v_detach = kDetach * (tau_wall / tauCrit)^m (shear erosion, optional)

    Negative sign means mesh grows INTO domain (pore clogging).

SeeAlso
    Foam::dissolSolidVelocityPointPatchVectorField

 SourceFiles
    solidGrowthVelocityPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef solidGrowthVelocityPointPatchVectorField_H
#define solidGrowthVelocityPointPatchVectorField_H

#include "fixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

/*---------------------------------------------------------------------------*\
        Class solidGrowthVelocityPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidGrowthVelocityPointPatchVectorField
    : public fixedValuePointPatchField<vector> {
  // Private data

  //- Name of solid-phase concentration field (e.g. biomass, mineral)
  word speciesName_;

  //- Name of velocity field (for shear stress calculation)
  word UName_;

  //- Volumetric growth coefficient (m^3_solid / mol)
  scalar growthCoeff_;

  //- Solid phase density (kg/m^3)
  scalar rhoSolid_;

  //- Molecular weight of solid phase (kg/mol)
  scalar Mw_;

  //- Critical shear stress for detachment (Pa)
  scalar tauCrit_;

  //- Detachment rate constant (m/s), 0 to disable
  scalar kDetach_;

  //- Detachment exponent
  scalar detachExponent_;

  //- Dynamic viscosity for shear calculation (Pa.s)
  scalar mu_;

  //- Maximum mesh velocity limiter (m/s)
  scalar maxVelocity_;

public:
  //- Runtime type information
  TypeName("solidGrowthVelocity");

  // Constructors

  //- Construct from patch and internal field
  solidGrowthVelocityPointPatchVectorField(
      const pointPatch &, const DimensionedField<vector, pointMesh> &);

  //- Construct from patch, internal field and dictionary
  solidGrowthVelocityPointPatchVectorField(
      const pointPatch &, const DimensionedField<vector, pointMesh> &,
      const dictionary &);

  //- Construct by mapping given patchField onto a new patch
  solidGrowthVelocityPointPatchVectorField(
      const solidGrowthVelocityPointPatchVectorField &, const pointPatch &,
      const DimensionedField<vector, pointMesh> &,
      const pointPatchFieldMapper &);

  //- Construct and return a clone
  virtual autoPtr<pointPatchField<vector>> clone() const {
    return autoPtr<pointPatchField<vector>>(
        new solidGrowthVelocityPointPatchVectorField(*this));
  }

  //- Construct as copy setting internal field reference
  solidGrowthVelocityPointPatchVectorField(
      const solidGrowthVelocityPointPatchVectorField &,
      const DimensionedField<vector, pointMesh> &);

  //- Construct and return a clone setting internal field reference
  virtual autoPtr<pointPatchField<vector>>
  clone(const DimensionedField<vector, pointMesh> &iF) const {
    return autoPtr<pointPatchField<vector>>(
        new solidGrowthVelocityPointPatchVectorField(*this, iF));
  }

  // Member functions

  // Evaluation functions

  //- Update the coefficients associated with the patch field
  virtual void updateCoeffs();

  //- Write
  virtual void write(Ostream &) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
