// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0
//
/// @file    Blend.h
///
/// @author  Andre Pradhana
///
/// @brief   Define methods to blend two level sets together. One such approach
///          is to carve an excess fillet so that the resulting blend appears
///          smoother than a regular union.
///
/// @details The algorithm used by unionFillet is based on the 2007 SIGGRAPH
///          talk "Levelsets in production: Spider-man 3" by Allen et al.
///          https://dl.acm.org/doi/10.1145/1278780.1278815

#ifndef OPENVDB_TOOLS_BLEND_HAS_BEEN_INCLUDED
#define OPENVDB_TOOLS_BLEND_HAS_BEEN_INCLUDED

#include <openvdb/Grid.h>
#include <openvdb/Types.h>
#include <openvdb/math/Math.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/FastSweeping.h> // for maskSdf
#include <openvdb/tools/Merge.h>
#include <openvdb/tools/Morphology.h>
#include <openvdb/tools/SignedFloodFill.h>

#include <stdexcept>
#include <type_traits>

namespace openvdb {
OPENVDB_USE_VERSION_NAMESPACE
namespace OPENVDB_VERSION_NAME {
namespace tools {

/// @brief Threaded VDB union with fillet that produces a new grid from
/// immutable inputs.
///
/// @param lhs Level set grid to be combined with a second input.
/// @param rhs Level set grid to be combined with the first input.
/// @param mask Optional float grid that controls the strength of the blend.
/// @param alpha Controls the blend radius.
/// @param beta Controls the exponent used to shape the blend falloff.
/// @param gamma Controls the strength of the blend.
/// @param supportDilation Optional number of voxels by which to dilate the
/// local blend support region before extending both inputs into that region.
/// A value of zero preserves the default behavior.
///
/// @return The filleted union of the @a lhs and @a rhs level set inputs.
///
/// @throw std::runtime_error if the transforms of @a lhs, @a rhs, and @a mask do not match.
template<typename GridT,
    typename MaskT = typename GridT::template ValueConverter<float>::Type>
typename GridT::Ptr
unionFillet(const GridT& lhs,
    const GridT& rhs,
    typename MaskT::ConstPtr mask,
    typename GridT::ValueType alpha,
    typename GridT::ValueType beta,
    typename GridT::ValueType gamma,
    int supportDilation = 0);


/// @cond OPENVDB_DOCS_INTERNAL

template<typename GridT,
    typename MaskT = typename GridT::template ValueConverter<float>::Type>
struct UnionWithFillet
{
    using TreeT = typename GridT::TreeType;
    using ValueType = typename TreeT::ValueType;
    using TreePtrType = typename TreeT::Ptr;
    using LeafNodeType = typename TreeT::LeafNodeType;
    using MaskTreeType = typename MaskT::TreeType;
    using MaskLeafNodeType = typename MaskTreeType::LeafNodeType;
    using MaskValueType = typename MaskT::ValueType;

    UnionWithFillet(const GridT& lhsGrid,
        const GridT& rhsGrid,
        typename MaskT::ConstPtr mask,
        const ValueType& bandwidth,
        const ValueType& exponent,
        const ValueType& multiplier)
        : mRhsGrid(&rhsGrid)
        , mLhsTree(&lhsGrid.tree())
        , mRhsTree(&rhsGrid.tree())
        , mMaskTree(mask ? mask->treePtr() : nullptr)
        , mBandwidth(bandwidth)
        , mExponent(exponent)
        , mMultiplier(multiplier)
    {
        static_assert(std::is_floating_point<ValueType>::value,
            "UnionWithFillet requires a scalar floating-point grid.");
    }

    typename GridT::Ptr blend();

private:
    struct FilletParms
    {
        ValueType mAlpha;
        ValueType mBeta;
        ValueType mGamma;
    };

    /// Combine overlapping leaves while the generic union merge operator
    /// handles sparse topology, signed tiles, and non-overlapping branches.
    struct FilletLeafOp
    {
        FilletLeafOp(typename MaskTreeType::ConstPtr mask, const FilletParms& parms)
            : mMaskTree(mask)
            , mParms(parms)
        {
        }

        void operator()(LeafNodeType& lhsLeaf, const LeafNodeType& rhsLeaf,
            bool /*pruneCancelledTiles*/, const ValueType& /*background*/) const
        {
            ValueType* lhsData = lhsLeaf.buffer().data();
            const ValueType* rhsData = rhsLeaf.buffer().data();
            typename LeafNodeType::NodeMaskType& lhsMask = lhsLeaf.getValueMask();
            const typename LeafNodeType::NodeMaskType& rhsMask = rhsLeaf.getValueMask();

            // DynamicNodeManager shares this policy between workers. Probe the
            // immutable mask tree directly instead of sharing a mutable accessor.
            const MaskLeafNodeType* maskLeaf = mMaskTree
                ? mMaskTree->template probeConstNode<MaskLeafNodeType>(lhsLeaf.origin())
                : nullptr;
            const MaskValueType* maskData = maskLeaf ? maskLeaf->buffer().data() : nullptr;
            const ValueType maskBackground = mMaskTree
                ? static_cast<ValueType>(mMaskTree->background())
                : ValueType(1);

            const ValueType alpha = mParms.mAlpha;
            const ValueType beta = mParms.mBeta;
            const ValueType gamma = mParms.mGamma;

            for (Index pos = 0; pos < LeafNodeType::SIZE; ++pos) {
                const ValueType A = lhsData[pos];
                const ValueType B = rhsData[pos];
                const bool lhsActive = lhsMask.isOn(pos);
                const bool rhsActive = rhsMask.isOn(pos);

                const ValueType blend =
                    math::Clamp((alpha - A) / alpha, ValueType(0), ValueType(1)) *
                    math::Clamp((alpha - B) / alpha, ValueType(0), ValueType(1));
                const ValueType offset = lhsActive && rhsActive
                    ? math::Pow(blend, beta) * gamma
                    : ValueType(0);
                const ValueType maskValue = maskData
                    ? static_cast<ValueType>(maskData[pos])
                    : maskBackground;

                // Preserve the existing tie behavior: when A equals B, select
                // B's active state rather than the standard union preference for A.
                const bool useA = A < B;
                lhsData[pos] = (useA ? A : B) - offset * maskValue;
                lhsMask.set(pos, useA ? lhsActive : rhsActive);
            }
        }

        typename MaskTreeType::ConstPtr mMaskTree;
        FilletParms mParms;
    };

    TreePtrType mSegment;
    GridT const * const mRhsGrid;
    TreeT const * const mLhsTree;
    TreeT const * const mRhsTree;
    typename MaskTreeType::ConstPtr mMaskTree;
    ValueType mBandwidth;
    ValueType mExponent;
    ValueType mMultiplier;
};


template<typename GridT>
MaskGrid::Ptr
createBlendSupportMask(const GridT& lhs,
    const GridT& rhs,
    typename GridT::ValueType alpha,
    int supportDilation)
{
    MaskGrid::Ptr supportMask = MaskGrid::create();
    supportMask->setTransform(lhs.transform().copy());

    const typename GridT::ValueType supportRadius =
        alpha + typename GridT::ValueType(supportDilation * lhs.voxelSize()[0]);

    typename GridT::ConstAccessor lhsAcc = lhs.getConstAccessor();
    typename GridT::ConstAccessor rhsAcc = rhs.getConstAccessor();
    MaskGrid::Accessor supportAcc = supportMask->getAccessor();

    for (typename GridT::ValueOnCIter iter = lhs.cbeginValueOn(); iter; ++iter) {
        const Coord& ijk = iter.getCoord();
        // Extend only regions where both SDFs can contribute to the fillet.
        if (iter.getValue() < supportRadius && rhsAcc.getValue(ijk) < supportRadius) {
            supportAcc.setValueOn(ijk);
        }
    }

    for (typename GridT::ValueOnCIter iter = rhs.cbeginValueOn(); iter; ++iter) {
        const Coord& ijk = iter.getCoord();
        if (lhsAcc.getValue(ijk) < supportRadius && iter.getValue() < supportRadius) {
            supportAcc.setValueOn(ijk);
        }
    }

    if (supportMask->activeVoxelCount() > 0 && supportDilation > 0) {
        tools::dilateActiveValues(supportMask->tree(), supportDilation,
            tools::NN_FACE_EDGE_VERTEX, tools::IGNORE_TILES);
    }

    return supportMask;
}


template<typename GridT, typename MaskT>
typename GridT::Ptr
UnionWithFillet<GridT, MaskT>::blend()
{
    const FilletParms parms{mBandwidth, mExponent, mMultiplier};
    const FilletLeafOp leafOp(mMaskTree, parms);

    // The public operation keeps both inputs immutable. Copy A into the
    // mutable destination and copy B before allowing the merge to steal its
    // branches. Using Steal avoids concurrent writes to TreeToMerge's tracking
    // mask when the top-down merge processes a const source in parallel.
    mSegment.reset(new TreeT(*mLhsTree));
    TreePtrType rhsCopy(new TreeT(*mRhsTree));
    CsgUnionOp<TreeT, FilletLeafOp> mergeOp(*rhsCopy, Steal(), leafOp);
    tree::DynamicNodeManager<TreeT> nodeManager(*mSegment);
    nodeManager.foreachTopDown(mergeOp);

    // Leaf signs are set by the fillet policy. Reconstruct signed inactive
    // values above the leaf level after the fillet has moved samples.
    tools::signedFloodFill(
        *mSegment, /*threaded=*/true, /*grainSize=*/1, /*minLevel=*/1);

    typename GridT::Ptr result = GridT::create(mSegment);
    result->setTransform(mRhsGrid->transform().copy());
    result->setGridClass(GRID_LEVEL_SET);
    return result;
}

/// @endcond


template<typename GridT, typename MaskT>
typename GridT::Ptr
unionFillet(const GridT& lhs,
    const GridT& rhs,
    typename MaskT::ConstPtr mask,
    typename GridT::ValueType alpha,
    typename GridT::ValueType beta,
    typename GridT::ValueType gamma,
    int supportDilation)
{
    static_assert(std::is_floating_point<typename GridT::ValueType>::value,
        "unionFillet requires a scalar floating-point grid.");

    const math::Transform& lhsXform = lhs.constTransform();
    const math::Transform& rhsXform = rhs.constTransform();
    if (lhsXform != rhsXform) {
        throw std::runtime_error("The two grids need to have the same transforms.");
    }
    if (mask && lhsXform != mask->constTransform()) {
        throw std::runtime_error("The grids and the mask need to have the same transforms.");
    }

    if (supportDilation > 0) {
        MaskGrid::Ptr supportMask =
            createBlendSupportMask(lhs, rhs, alpha, supportDilation);
        if (supportMask->activeVoxelCount() > 0) {
            typename GridT::Ptr lhsExtended = tools::maskSdf(lhs, *supportMask);
            typename GridT::Ptr rhsExtended = tools::maskSdf(rhs, *supportMask);
            UnionWithFillet<GridT, MaskT> op(
                *lhsExtended, *rhsExtended, mask, alpha, beta, gamma);
            return op.blend();
        }
    }

    UnionWithFillet<GridT, MaskT> op(lhs, rhs, mask, alpha, beta, gamma);
    return op.blend();
}

} // namespace tools
} // namespace OPENVDB_VERSION_NAME
} // namespace openvdb

#endif // OPENVDB_TOOLS_BLEND_HAS_BEEN_INCLUDED
