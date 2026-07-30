/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017,2022 OpenFOAM Foundation
    Copyright (C) 2020-2024 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "slicedVolFields.H"
#include "slicedSurfaceFields.H"
#include "SubField.H"
#include "cyclicFvPatchFields.H"
#include "cyclicAMIFvPatchFields.H"
#include "processorFvPatch.H"
#include "fvcSurfaceIntegrate.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Calculate transform tensor with reference vector (unitAxis1)
// and direction vector (axis2).
//
// This is nearly identical to the meshTools axesRotation
// with an E3_E1 transformation with the following exceptions:
//
// - axis1 (e3 == unitAxis1): is already normalized (unit vector)
// - axis2 (e1 == dirn): no difference
// - transformation is row-vectors, not column-vectors
static inline tensor rotation_e3e1
(
    const vector& unitAxis1,
    vector dirn
)
{
    dirn.removeCollinear(unitAxis1);
    dirn.normalise();

    // Set row vectors
    return tensor
    (
        dirn,
        (unitAxis1^dirn),
        unitAxis1
    );
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMesh::makeSf() const
{
    DebugInFunction << "Assembling face areas" << endl;

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (SfPtr_)
    {
        FatalErrorInFunction
            << "face areas already exist"
            << abort(FatalError);
    }

    SfPtr_ = std::make_unique<slicedSurfaceVectorField>
    (
        IOobject
        (
            "S",
            pointsInstance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        *this,
        // Reduced-order (2D) manifold: face "area" is an edge length
        (nTopologicalD() == 2 ? dimLength : dimArea),
        faceAreas()
    );

    SfPtr_->setOriented();
}


void Foam::fvMesh::makeMagSf() const
{
    DebugInFunction << "Assembling mag face areas" << endl;

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (magSfPtr_)
    {
        FatalErrorInFunction
            << "mag face areas already exist"
            << abort(FatalError);
    }

    // Note: Added stabilisation for faces with exactly zero area.
    // These should be caught on mesh checking but at least this stops
    // the code from producing Nans.
    magSfPtr_ = std::make_unique<surfaceScalarField>
    (
        IOobject
        (
            "magSf",
            pointsInstance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mag(Sf()) + dimensionedScalar("vs", Sf().dimensions(), VSMALL)
    );
}


void Foam::fvMesh::makeC() const
{
    DebugInFunction << "Assembling cell centres" << endl;

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (CPtr_ || faceAreaNormalsPtr_)
    {
        FatalErrorInFunction
            << "cell centres already exist"
            << abort(FatalError);
    }

    // Construct as slices. Only preserve processor (not e.g. cyclic)

    CPtr_ = std::make_unique<slicedVolVectorField>
    (
        IOobject
        (
            "C",
            pointsInstance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        *this,
        dimLength,
        cellCentres(),
        faceCentres(),
        true,               //preserveCouples
        true                //preserveProcOnly
    );


    if (nTopologicalD() == 2)
    {
        DebugInFunction << "Assembling 2D cell normals" << endl;

        //faceAreaNormalsPtr_ = new volVectorField
        //(
        //    IOobject
        //    (
        //        "faceAreaNormals",
        //        pointsInstance(),
        //        meshSubDir,
        //        *this,
        //        IOobject::NO_READ,
        //        IOobject::NO_WRITE,
        //        IOobject::NO_REGISTER
        //    ),
        //    *this,
        //    dimless,
        //    primitiveMesh::cellDir0()
        //);
        //
        //auto& n = *faceAreaNormalsPtr_;
        //const tmp<vectorField> faceDir1(primitiveMesh::faceDir1());
        //for (auto& pfld : n.boundaryFieldRef())
        //{
        //    pfld = pfld.patch().patchSlice(faceDir1());
        //}
        //n.boundaryFieldRef().template evaluateCoupled<processorFvPatch>();
        faceAreaNormalsPtr_ = std::make_unique<slicedVolVectorField>
        (
            IOobject
            (
                "faceAreaNormals",
                pointsInstance(),
                meshSubDir,
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            *this,
            dimless,
            cellDir0(),
            faceDir1(),
            true,               //preserveCouples
            true                //preserveProcOnly
        );

        faceAreaNormalsPtr_->setOriented();
    }
}


void Foam::fvMesh::makeCf() const
{
    DebugInFunction << "Assembling face centres" << endl;

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (CfPtr_ || LePtr_ || edgeAreaNormalsPtr_)
    {
        FatalErrorInFunction
            << "face centres already exist"
            << abort(FatalError);
    }

    CfPtr_ = std::make_unique<slicedSurfaceVectorField>
    (
        IOobject
        (
            "Cf",
            pointsInstance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        *this,
        dimLength,
        faceCentres()
    );

    if (nTopologicalD() == 2)
    {
        DebugInFunction << "Assembling 2D face vectors" << endl;

        LePtr_ = std::make_unique<slicedSurfaceVectorField>
        (
            IOobject
            (
                "Le",
                pointsInstance(),
                meshSubDir,
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            *this,
            dimLength,
            faceVec2()
        );

        LePtr_->setOriented();

        edgeAreaNormalsPtr_ = std::make_unique<slicedSurfaceVectorField>
        (
            IOobject
            (
                "edgeAreaNormals",
                pointsInstance(),
                meshSubDir,
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            *this,
            dimless,
            faceDir1()
        );

        //- Related to neighbouring faces, not the edge itself
        //edgeAreaNormalsPtr_->setOriented();
    }
}


void Foam::fvMesh::makeFaceCurvatures() const
{
    DebugInFunction
        << "Calculating face curvatures" << endl;

    if (faceCurvaturesPtr_)
    {
        FatalErrorInFunction
            << "faceCurvaturesPtr_ already allocated"
            << abort(FatalError);
    }

    faceCurvaturesPtr_ = std::make_unique<volScalarField>
    (
        IOobject
        (
            "faceCurvatures",
            pointsInstance(),
            meshSubDir,
            thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        *this,
        dimless/dimLength
    );
    volScalarField& faceCurvatures = *faceCurvaturesPtr_;


//     faceCurvatures =
//         fac::edgeIntegrate(Le()*edgeLengthCorrection())
//         &faceAreaNormals();

    //- Using fvc instead of fac just to avoid dependency. Dimensions are
    //- handled inside fvc::surfaceIntegrate using nTopologicalD.
    auto kN(fvc::surfaceIntegrate(Le()*edgeLengthCorrection()));

    faceCurvatures = sign(kN()&faceAreaNormals())*mag(kN());
}


void Foam::fvMesh::calcEdgeTransformTensors() const
{
    DebugInFunction
        << "Calculating edge transformation tensors" << endl;

    if (edgeTransformTensorsPtr_)
    {
        FatalErrorInFunction
            << "edgeTransformTensorsPtr_ already allocated"
            << abort(FatalError);
    }

    edgeTransformTensorsPtr_ =
        std::make_unique<FieldField<Field, tensor>>(nFaces());
    auto& edgeTransformTensors = *edgeTransformTensorsPtr_;

    // Initialize all transformation tensors
    for (label edgei = 0; edgei < nFaces(); ++edgei)
    {
        edgeTransformTensors.set(edgei, new Field<tensor>(3, tensor::I));
    }

    const volVectorField& Nf = faceAreaNormals();
    const volVectorField& Cf = areaCentres();

    const surfaceVectorField& Ne = edgeAreaNormals();
    const surfaceVectorField& Ce = edgeCentres();

    // Halo-swapped neighbour-cell geometry for the processor-edge branch below.
    //
    // faceAreaNormals()/areaCentres() are cached *sliced* fields: their coupled
    // (processor) patch value is initialised to the boundary-edge slice and is
    // never evaluate()d, so patchNeighbourField() would return that edge value
    // rather than the neighbour area-cell centre/normal. Build short-lived,
    // halo-swapped copies here so the neighbour transform tensor is built from
    // the *neighbour* cell's tangent plane. This is essential for the
    // non-orthogonal Laplacian correction across processor edges on a curved
    // manifold (it is a no-op on a flat plane, where the transform is identity)
    // and mirrors the legacy faMesh, which evaluateCoupled()s these exact
    // geometry fields (finiteArea/faMesh/faMeshDemandDrivenData.C). Local copies
    // are used so the shared cached cell-centre/normal fields are left intact.
    autoPtr<volVectorField> NfSwapPtr;
    autoPtr<volVectorField> CfSwapPtr;
    if (Pstream::parRun())
    {
        NfSwapPtr.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "faceAreaNormalsProcSwap",
                    pointsInstance(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                Nf
            )
        );
        NfSwapPtr->boundaryFieldRef()
            .template evaluateCoupled<processorFvPatch>();

        CfSwapPtr.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "areaCentresProcSwap",
                    pointsInstance(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                Cf
            )
        );
        CfSwapPtr->boundaryFieldRef()
            .template evaluateCoupled<processorFvPatch>();
    }

    // FA-safe skew-correction vectors for the reduced-order (2-manifold) mesh.
    //
    // The generic finiteVolume skewCorrectionVectors forms the correction as
    // Cpf - ((Sf & Cpf)/(Sf & d))*d, which is singular here: on the reduced
    // mesh a "face" is a 2-point edge whose area vector Sf is zero, so (Sf & d)
    // vanishes and the division is a floating-point exception. Instead derive
    // the correction directly from the edge geometry (area/edge centres and the
    // edge vector), mirroring the legacy finiteArea
    // edgeInterpolation::makeSkewCorrectionVectors and guarding degenerate
    // edges (alpha < SMALL -> zero correction). On an orthogonal/flat mesh the
    // correction is ~0, so this is a no-op there and only activates on skewed
    // curved manifolds.
    surfaceVectorField skewCorrVecs
    (
        IOobject
        (
            "edgeSkewCorrectionVectors",
            pointsInstance(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        *this,
        dimensionedVector(dimLength, Zero)
    );
    {
        const pointField& pts = points();
        const faceList& fcs = faces();

        // Internal edges
        vectorField& scvIn = skewCorrVecs.primitiveFieldRef();
        for (label edgei = 0; edgei < nInternalFaces(); ++edgei)
        {
            const vector& P = Cf.internalField()[owner()[edgei]];
            const vector& N = Cf.internalField()[neighbour()[edgei]];
            const face& f = fcs[edgei];
            const vector& S = pts[f[0]];

            const vector d(N - P);
            const vector e(pts[f[1]] - pts[f[0]]);
            const vector de(d ^ e);
            const scalar alpha = magSqr(de);

            if (alpha < SMALL)
            {
                continue;
            }

            const scalar beta = -((d ^ (S - P)) & de)/alpha;
            const vector E(S + beta*e);

            scvIn[edgei] = Ce.internalField()[edgei] - E;
        }

        // Coupled boundary edges: use the halo-swapped neighbour area-centres
        // (CfSwapPtr) built above. Non-coupled boundary edges keep zero skew,
        // matching the legacy finiteArea treatment.
        if (Pstream::parRun())
        {
            auto& scvBf = skewCorrVecs.boundaryFieldRef();

            forAll(boundary(), patchi)
            {
                if (!boundary()[patchi].coupled())
                {
                    continue;
                }

                const labelUList& edgeFaces = boundary()[patchi].faceCells();
                const label startFace = boundary()[patchi].start();

                vectorField nbrCf
                (
                    CfSwapPtr->boundaryField()[patchi].patchNeighbourField()
                );

                auto& pscv = scvBf[patchi];

                forAll(edgeFaces, bndEdgei)
                {
                    const label facei = startFace + bndEdgei;
                    const vector& P = Cf.internalField()[edgeFaces[bndEdgei]];
                    const vector& N = nbrCf[bndEdgei];
                    const face& f = fcs[facei];
                    const vector& S = pts[f[0]];

                    const vector d(N - P);
                    const vector e(pts[f[1]] - pts[f[0]]);
                    const vector de(d ^ e);
                    const scalar alpha = magSqr(de);

                    if (alpha < SMALL)
                    {
                        continue;
                    }

                    const scalar beta = -((d ^ (S - P)) & de)/alpha;
                    const vector E(S + beta*e);

                    pscv[bndEdgei] =
                        Ce.boundaryField()[patchi][bndEdgei] - E;
                }
            }
        }
    }

    // Internal edges transformation tensors
    for (label edgei = 0; edgei < nInternalFaces(); ++edgei)
    {
        const label ownFacei = owner()[edgei];
        const label neiFacei = neighbour()[edgei];
        auto& tensors = edgeTransformTensors[edgei];

        vector edgeCtr = Ce.internalField()[edgei];
        edgeCtr -= skewCorrVecs.internalField()[edgei];

        // Edge transformation tensor
        tensors[0] =
            rotation_e3e1
            (
                Ne.internalField()[edgei],
                (edgeCtr - Cf.internalField()[ownFacei])
            );

        // Owner transformation tensor
        tensors[1] =
            rotation_e3e1
            (
                Nf.internalField()[ownFacei],
                (edgeCtr - Cf.internalField()[ownFacei])
            );

        // Neighbour transformation tensor
        tensors[2] =
            rotation_e3e1
            (
                Nf.internalField()[neiFacei],
                (Cf.internalField()[neiFacei] - edgeCtr)
            );
    }

    // Boundary edges transformation tensors
    forAll(boundary(), patchi)
    {
        const labelUList& edgeFaces = boundary()[patchi].faceCells();

        if (boundary()[patchi].coupled())
        {
            // Neighbour area-cell centre/normal from the halo-swapped copies
            // (see the NfSwap/CfSwap construction above): the cached sliced
            // fields' patchNeighbourField() would otherwise return the local
            // boundary-edge slice, not the neighbour cell value.
            vectorField nbrCf
            (
                CfSwapPtr->boundaryField()[patchi].patchNeighbourField()
            );
            vectorField nbrNf
            (
                NfSwapPtr->boundaryField()[patchi].patchNeighbourField()
            );

            forAll(edgeFaces, bndEdgei)
            {
                const label ownFacei = edgeFaces[bndEdgei];
                const label meshEdgei(boundary()[patchi].start() + bndEdgei);

                auto& tensors = edgeTransformTensors[meshEdgei];

                vector edgeCtr = Ce.boundaryField()[patchi][bndEdgei];
                edgeCtr -= skewCorrVecs.boundaryField()[patchi][bndEdgei];

                // Edge transformation tensor
                tensors[0] =
                    rotation_e3e1
                    (
                        Ne.boundaryField()[patchi][bndEdgei],
                        (edgeCtr - Cf.internalField()[ownFacei])
                    );

                // Owner transformation tensor
                tensors[1] =
                    rotation_e3e1
                    (
                        Nf.internalField()[ownFacei],
                        (edgeCtr - Cf.internalField()[ownFacei])
                    );

                // Neighbour transformation tensor
                tensors[2] =
                    rotation_e3e1
                    (
                        nbrNf[bndEdgei],
                        (nbrCf[bndEdgei] - edgeCtr)
                    );
            }
        }
        else
        {
            forAll(edgeFaces, bndEdgei)
            {
                const label ownFacei = edgeFaces[bndEdgei];
                const label meshEdgei(boundary()[patchi].start() + bndEdgei);

                auto& tensors = edgeTransformTensors[meshEdgei];

                vector edgeCtr = Ce.boundaryField()[patchi][bndEdgei];
                edgeCtr -= skewCorrVecs.boundaryField()[patchi][bndEdgei];

                // Edge transformation tensor
                tensors[0] =
                    rotation_e3e1
                    (
                        Ne.boundaryField()[patchi][bndEdgei],
                        (edgeCtr - Cf.internalField()[ownFacei])
                    );

                // Owner transformation tensor
                tensors[1] =
                    rotation_e3e1
                    (
                        Nf.internalField()[ownFacei],
                        (edgeCtr - Cf.internalField()[ownFacei])
                    );

                // Neighbour transformation tensor
                tensors[2] = tensor::I;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volScalarField::Internal& Foam::fvMesh::V() const
{
    if (!VPtr_)
    {
        DebugInFunction
            << "Constructing from primitiveMesh::cellVolumes()" << endl;

        VPtr_ = std::make_unique<SlicedDimensionedField<scalar, volMesh>>
        (
            IOobject
            (
                "V",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            *this,
            // Reduced-order (2D) manifold: cell "volume" is a face area
            (nTopologicalD() == 2 ? dimArea : dimVolume),
            cellVolumes()
        );
    }

    return *VPtr_;
}


const Foam::volScalarField::Internal& Foam::fvMesh::V0() const
{
    if (!V0Ptr_)
    {
        FatalErrorInFunction
            << "V0 is not available"
            << abort(FatalError);
    }

    return *V0Ptr_;
}


Foam::volScalarField::Internal& Foam::fvMesh::setV0()
{
    if (!V0Ptr_)
    {
        FatalErrorInFunction
            << "V0 is not available"
            << abort(FatalError);
    }

    return *V0Ptr_;
}


const Foam::volScalarField::Internal& Foam::fvMesh::V00() const
{
    if (!V00Ptr_)
    {
        DebugInFunction << "Constructing from V0" << endl;

        V00Ptr_ = std::make_unique<DimensionedField<scalar, volMesh>>
        (
            IOobject
            (
                "V00",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            V0()
        );


        // If V00 is used then V0 should be stored for restart
        V0Ptr_->writeOpt(IOobject::AUTO_WRITE);
    }

    return *V00Ptr_;
}


Foam::tmp<Foam::volScalarField::Internal> Foam::fvMesh::Vsc() const
{
    if (!steady() && moving() && time().subCycling())
    {
        const TimeState& ts = time();
        const TimeState& ts0 = time().prevTimeState();

        scalar tFrac =
        (
            ts.value() - (ts0.value() - ts0.deltaTValue())
        )/ts0.deltaTValue();

        if (tFrac < (1 - SMALL))
        {
            return V0() + tFrac*(V() - V0());
        }
    }

    return V();
}


Foam::tmp<Foam::volScalarField::Internal> Foam::fvMesh::Vsc0() const
{
    if (!steady() && moving() && time().subCycling())
    {
        const TimeState& ts = time();
        const TimeState& ts0 = time().prevTimeState();

        scalar t0Frac =
        (
            (ts.value() - ts.deltaTValue())
          - (ts0.value() - ts0.deltaTValue())
        )/ts0.deltaTValue();

        if (t0Frac > SMALL)
        {
            return V0() + t0Frac*(V() - V0());
        }
    }

    return V0();
}


const Foam::surfaceVectorField& Foam::fvMesh::Sf() const
{
    if (!SfPtr_)
    {
        makeSf();
    }

    return *SfPtr_;
}


const Foam::surfaceScalarField& Foam::fvMesh::magSf() const
{
    if (!magSfPtr_)
    {
        makeMagSf();
    }

    return *magSfPtr_;
}


Foam::tmp<Foam::surfaceVectorField> Foam::fvMesh::unitSf() const
{
    auto tunitVectors = tmp<surfaceVectorField>::New
    (
        IOobject
        (
            "unit(Sf)",
            pointsInstance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        *this,
        dimless,
        (this->Sf() / this->magSf())
    );

    tunitVectors.ref().oriented() = this->Sf().oriented();
    return tunitVectors;
}


const Foam::volVectorField& Foam::fvMesh::C() const
{
    if (!CPtr_)
    {
        makeC();
    }

    return *CPtr_;
}


const Foam::volVectorField& Foam::fvMesh::faceAreaNormals() const
{
    if (!faceAreaNormalsPtr_)
    {
        if (nTopologicalD() == 3)
        {
            FatalErrorInFunction
                << "undefined cell normal direction in 3 dimensions"
                << exit(FatalError);
        }
        makeC();
    }

    return *faceAreaNormalsPtr_;
}


const Foam::surfaceVectorField& Foam::fvMesh::Cf() const
{
    if (!CfPtr_)
    {
        makeCf();
    }

    return *CfPtr_;
}


const Foam::surfaceVectorField& Foam::fvMesh::edgeAreaNormals() const
{
    if (!edgeAreaNormalsPtr_)
    {
        if (nTopologicalD() == 3)
        {
            FatalErrorInFunction
                << "undefined face tangential direction in 3 dimensions"
                << exit(FatalError);
        }
        makeCf();
    }

    return *edgeAreaNormalsPtr_;
}


const Foam::surfaceVectorField& Foam::fvMesh::Le() const
{
    if (!LePtr_)
    {
        if (nTopologicalD() == 3)
        {
            FatalErrorInFunction
                << "undefined face tangential direction in 3 dimensions"
                << exit(FatalError);
        }
        makeCf();
    }

    return *LePtr_;
}


Foam::tmp<Foam::surfaceScalarField> Foam::fvMesh::magLe() const
{
    return mag(Le());
}


Foam::tmp<Foam::surfaceVectorField> Foam::fvMesh::unitLe() const
{
    return Le()/mag(Le());
}


const Foam::vectorField& Foam::fvMesh::pointAreaNormals() const
{
    return pointDir0();
}


const Foam::volScalarField& Foam::fvMesh::faceCurvatures() const
{
    if (!faceCurvaturesPtr_)
    {
        if (nTopologicalD() == 3)
        {
            FatalErrorInFunction
                << "undefined cell tangential direction in 3 dimensions"
                << exit(FatalError);
        }
        makeFaceCurvatures();
    }

    return *faceCurvaturesPtr_;
}


const Foam::FieldField<Foam::Field, Foam::tensor>&
Foam::fvMesh::edgeTransformTensors() const
{
    if (!edgeTransformTensorsPtr_)
    {
        calcEdgeTransformTensors();
    }

    return *edgeTransformTensorsPtr_;
}


Foam::tmp<Foam::surfaceScalarField> Foam::fvMesh::edgeLengthCorrection() const
{
    DebugInFunction
        << "Calculating edge length correction" << endl;

    auto tcorrection = tmp<surfaceScalarField>::New
    (
        IOobject
        (
            "edgeLengthCorrection",
            pointsInstance(),
            meshSubDir,
            thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        *this,
        dimless
    );
    auto& correction = tcorrection.ref();

    const vectorField& pointNormals = pointAreaNormals();

    const auto angleCorrection =
        [](const vector& a, const vector& b) -> scalar
        {
            return Foam::cos(0.5*Foam::asin(Foam::mag(a ^ b)));
        };


    // Internal
    {
        auto& fld = correction.primitiveFieldRef();

        forAll(fld, facei)
        {
            const auto& f = faces()[facei];
            fld[facei] = angleCorrection
            (
                pointNormals[f[0]],
                pointNormals[f[1]]
            );
        }
    }

    // Boundary
    {
        auto& bfld = correction.boundaryFieldRef();

        forAll(boundary(), patchi)
        {
            const fvPatch& fap = boundary()[patchi];
            scalarField& pfld = bfld[patchi];

            forAll(pfld, patchFacei)
            {
                const auto& f = fap.patch()[patchFacei];
                pfld[patchFacei] = angleCorrection
                (
                    pointNormals[f[0]],
                    pointNormals[f[1]]
                );
            }
        }
    }

    return tcorrection;
}


Foam::tmp<Foam::surfaceVectorField> Foam::fvMesh::delta() const
{
    DebugInFunction << "Calculating face deltas" << endl;

    auto tdelta = tmp<surfaceVectorField>::New
    (
        IOobject
        (
            "delta",
            pointsInstance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        *this,
        dimLength
    );
    auto& delta = tdelta.ref();
    delta.setOriented();

    const volVectorField& C = this->C();
    const labelUList& owner = this->owner();
    const labelUList& neighbour = this->neighbour();

    forAll(owner, facei)
    {
        delta[facei] = C[neighbour[facei]] - C[owner[facei]];
    }

    auto& deltabf = delta.boundaryFieldRef();

    forAll(deltabf, patchi)
    {
        deltabf[patchi] = boundary()[patchi].delta();
    }

    return tdelta;
}


const Foam::surfaceScalarField& Foam::fvMesh::phi() const
{
    if (!phiPtr_)
    {
        FatalErrorInFunction
            << "mesh flux field does not exist, is the mesh actually moving?"
            << abort(FatalError);
    }

    // Set zero current time
    // mesh motion fluxes if the time has been incremented
    if (!time().subCycling() && phiPtr_->timeIndex() != time().timeIndex())
    {
        (*phiPtr_) = Zero;
    }

    phiPtr_->setOriented();

    return *phiPtr_;
}


Foam::refPtr<Foam::surfaceScalarField> Foam::fvMesh::setPhi()
{
    refPtr<surfaceScalarField> phiref;

    // Return non-const reference, or nullptr if not available
    phiref.ref(phiPtr_.get());

    return phiref;
}


// ************************************************************************* //
