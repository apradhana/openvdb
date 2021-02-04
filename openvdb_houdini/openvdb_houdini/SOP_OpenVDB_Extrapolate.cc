// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0
//
// @file SOP_OpenVDB_Extrapolate.cc
//
// @author FX R&D OpenVDB team
//
// @brief Extrapolate SDF or attributes off a level set surface

#include <houdini_utils/ParmFactory.h>
#include <openvdb_houdini/Utils.h>
#include <openvdb_houdini/SOP_NodeVDB.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/FastSweeping.h>
#include <stdexcept>
#include <string>

namespace hvdb = openvdb_houdini;
namespace hutil = houdini_utils;


namespace {

struct FastSweepingParms {
    FastSweepingParms() :
        mTime(0.0f),
        mFSGroup(nullptr),
        mExtGroup(nullptr),
        mMode(""),
        mNeedExt(false),
        mNSweeps(1),
        mIgnoreTiles(false),
        mPattern(""),
        mDilate(1),
        mFSPrimName(""),
        mExtPrimName(""),
        mExtFieldProcessed(false)
    { }

    fpreal mTime;
    const GA_PrimitiveGroup* mFSGroup;
    const GA_PrimitiveGroup* mExtGroup;
    UT_String mMode;
    bool mNeedExt;
    int mNSweeps;
    bool mIgnoreTiles;
    UT_String mPattern;
    int mDilate;
    std::string mFSPrimName;
    std::string mExtPrimName;
    bool mExtFieldProcessed;
};


/// @brief Helper class for use with GEOvdbApply and calling
///        openvdb::tools::maskSdf. The mask VDB is allowed to be
///        of any grid type.
template <typename GridT>
struct FastSweepingMaskOp
{
    FastSweepingMaskOp(const FastSweepingParms& parms, typename GridT::ConstPtr inGrid)
        : mParms(parms), mInGrid(inGrid), mOutGrid(nullptr) {}

    template<typename MaskGridType>
    void operator()(const MaskGridType& mask)
    {
        mOutGrid = openvdb::tools::maskSdf(*mInGrid, mask, mParms.mIgnoreTiles, mParms.mNSweeps);
    }

    const FastSweepingParms& mParms;
    typename GridT::ConstPtr mInGrid;
    hvdb::Grid::Ptr mOutGrid;
}; // FastSweepingMaskOp


/// @brief Sampler functor for the Dirichlet boundary condition on the
///        surface of the field to be extended.
template<typename ExtGridT>
struct DirichletSamplerOp
{
    using ExtValueT = typename ExtGridT::ValueType;
    using SamplerT = openvdb::tools::GridSampler<ExtGridT, openvdb::tools::BoxSampler>;

    DirichletSamplerOp(typename ExtGridT::ConstPtr functorGrid, SamplerT sampler)
        : mFunctorGrid (functorGrid),
          mSampler(sampler)
    {}

    ExtValueT operator()(const openvdb::Vec3d& xyz) const
    {
        return static_cast<ExtValueT>(mSampler.isSample(xyz));
    }

    typename ExtGridT::ConstPtr mFunctorGrid;
    SamplerT mSampler;
}; // DirichletSamplerOp

} // unnamed namespace


////////////////////////////////////////


class SOP_OpenVDB_Extrapolate: public hvdb::SOP_NodeVDB
{
public:

    SOP_OpenVDB_Extrapolate(OP_Network*, const char* name, OP_Operator*);

    ~SOP_OpenVDB_Extrapolate() override {}

    static OP_Node* factory(OP_Network*, const char* name, OP_Operator*);

    int isRefInput(unsigned input) const override { return (input == 1); }

    class Cache: public SOP_VDBCacheOptions
    {
    protected:
        OP_ERROR cookVDBSop(OP_Context&) override;
        OP_ERROR evalFastSweepingParms(OP_Context&, FastSweepingParms&);
    private:
        /// @brief Resolve the type of the Fast Sweeping grid, i.e. UT_VDB_FLOAT
        ///        or UT_VDB_DOUBLE
        template<typename FSGridT>
        bool processHelper(
            FastSweepingParms& parms,
            hvdb::GU_PrimVDB* lsPrim,
            typename FSGridT::ValueType fsIsoValue = typename FSGridT::ValueType(0),
            const hvdb::GU_PrimVDB* maskPrim = nullptr);

        template<typename FSGridT, typename ExtGridT = FSGridT>
        bool process(
            const FastSweepingParms& parms,
            hvdb::GU_PrimVDB* lsPrim,
            typename FSGridT::ValueType fsIsoValue = typename FSGridT::ValueType(0),
            const hvdb::GU_PrimVDB* maskPrim = nullptr,
            hvdb::GU_PrimVDB* exPrim = nullptr,
            const typename ExtGridT::ValueType& background = typename ExtGridT::ValueType(0),
            hvdb::Grid::Ptr newLSGrid = nullptr);
    }; // class Cache

protected:
    bool updateParmsFlags() override;
};


////////////////////////////////////////


void
newSopOperator(OP_OperatorTable* table)
{
    if (table == nullptr) return;

    hutil::ParmList parms;

    // Level set/Fog grid
    parms.add(hutil::ParmFactory(PRM_STRING, "group", "Source Group")
        .setChoiceList(&hutil::PrimGroupMenuInput1)
        .setTooltip("Specify a subset of the input VDB scalar grid(s).")
        .setDocumentation(
            "A subset of the input VDB scalar grid(s) to be processed"
            " (see [specifying volumes|/model/volumes#group])"));

    // Extension fields
    parms.add(hutil::ParmFactory(PRM_STRING, "extfields", "Extension Group")
        .setChoiceList(&hutil::PrimGroupMenuInput1)
        .setTooltip("Specify a subset of the VDB grid(s) to be extended off"
            " the surface defined by the scalar VDB defined by the Source Group")
        .setDocumentation("Arbitrary VDB fields picked up by this group"
            " will be extended off an iso-surface of a scalar VDB (fog/level set)"
            " as specified by Source Group. The modes that allow the use of this"
            " parameter are Extend Field(s) Off Scalar VDB and Extend Field(s) and"
            " Recompute SDF."));

    // Mask grid
    parms.add(hutil::ParmFactory(PRM_STRING, "mask", "Mask VDB")
        .setChoiceList(&hutil::PrimGroupMenuInput2)
        .setTooltip("Specify a VDB volume whose active voxels are to be used as a mask.")
        .setDocumentation(
            "A VDB volume whose active voxels are to be used as a mask"
            " (see [specifying volumes|/model/volumes#group])"));

    // Sweep
    parms.add(hutil::ParmFactory(PRM_HEADING, "sweep", "General Sweep")
         .setDocumentation(
             "These control the Fast Sweeping parameters."));

    // Modes
    parms.add(hutil::ParmFactory(PRM_STRING, "mode", "Operation")
        .setChoiceListItems(PRM_CHOICELIST_SINGLE, {
            "dilate",       "Dilate",
            "mask",         "Mask",
            "convert",      "Convert Fog VDB To SDF", ///< @todo move to Convert SOP
            "renormalize",  "Renormalize SDF", // by solving the Eikonal equation
            "scalarext",    "Extend Field(s) Off Scalar VDB",
            "scalarsdfext", "Extend Field(s) and Recompute SDF"
        })
        .setDefault("dilate")
        .setDocumentation(
            "The operation to perform\n\n"
            "Dilate:\n"
            "    Dilates an existing signed distance filed by a specified \n"
            "    number of voxels.\n"
            "Mask:\n"
            "    Expand/extrapolate an existing signed distance fild into\n"
            "    a mask.\n"
            "Convert Scalar VDB To SDF:\n"
            "    Converts a scalar fog volume into a signed distance\n"
            "    function. Active input voxels with scalar values above\n"
            "    the given isoValue will have NEGATIVE distance\n"
            "    values on output, i.e. they are assumed to be INSIDE\n"
            "    the iso-surface.\n"
            "Renormalize SDF:\n"
            "    Given an existing approximate SDF it solves the Eikonal\n"
            "    equation for all its active voxels. Active input voxels\n"
            "    with a signed distance value above the given isoValue\n"
            "    will have POSITIVE distance values on output, i.e. they are\n"
            "    assumed to be OUTSIDE the iso-surface.\n"
            "Extend Field(s) Off Scalar VDB:\n"
            "     Computes the extension of a scalar field, defined by the\n"
            "     specified functor, off an iso-surface from an input\n"
            "     Fog or SDF volume.\n"
            "     Note that this mode only uses the first fog/sdf grid\n"
            "     specified by the Source Group parameter.\n"
            "Extend Field(s) and Recompute SDF:\n"
            "    Computes the signed distance field and the extension of a\n"
            "    scalar field, defined by the specified functor, off an\n"
            "    iso-surface from an input Fog or SDF volume.\n"
            "    Note that this mode only uses the first fog/sdf grid\n"
            "    specified by the Source Group parameter."));

    parms.add(hutil::ParmFactory(PRM_INT_J, "sweeps", "Iterations")
         .setDefault(1)
         .setRange(PRM_RANGE_RESTRICTED, 1, PRM_RANGE_UI, 5)
         .setTooltip(
            "The desired number of iterations of the Fast Sweeping algorithm"
            " (one is often enough)"));

    parms.add(hutil::ParmFactory(PRM_FLT_J, "isovalue", "Isovalue")
         .setDefault(0.5)
         .setRange(PRM_RANGE_UI, -3, PRM_RANGE_UI, 3)
         .setTooltip("Isovalue for which the SDF is computed"));

    parms.add(hutil::ParmFactory(PRM_TOGGLE, "ignoretiles", "Ignore Active Tiles")
        .setDefault(PRMzeroDefaults)
        .setTooltip("Ignore active tiles in scalar field and mask VDBs.")
        .setDocumentation(
            "Ignore active tiles in scalar field and mask VDBs.\n\n"
            "This option should normally be disabled, but note that active tiles\n"
            "(sparsely represented regions of constant value) will in that case\n"
            "be densified, which could significantly increase memory usage.\n\n"
            "Proper signed distance fields don't have active tiles.\n"));

    parms.add(hutil::ParmFactory(PRM_INT_J, "dilate", "Dilation")
        .setDefault(3)
        .setRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_UI, 10)
        .setTooltip("The number of voxels by which to dilate the level set narrow band"));

    // Dilation Pattern
    parms.add(hutil::ParmFactory(PRM_STRING, "pattern", "Dilation Pattern")
         .setChoiceListItems(PRM_CHOICELIST_SINGLE, {
            "NN6",  "Faces",
            "NN18", "Faces and Edges",
            "NN26", "Faces, Edges and Vertices"
         })
         .setDefault("NN6")
         .setTooltip("Select the neighborhood pattern for the dilation operation.")
         .setDocumentation(
             "The neighborhood pattern for the dilation operation\n\n"
             "__Faces__ is fastest.  __Faces, Edges and Vertices__ is slowest\n"
             "but can produce the best results for large dilations.\n"
             "__Faces and Edges__ is intermediate in speed and quality.\n"));

    hvdb::OpenVDBOpFactory("VDB Extrapolate", SOP_OpenVDB_Extrapolate::factory, parms, *table)
        .addInput("Source VDB(s)")
        .addOptionalInput("Mask VDB")
        .setVerb(SOP_NodeVerb::COOK_INPLACE, []() { return new SOP_OpenVDB_Extrapolate::Cache; })
        .setDocumentation("\
#icon: COMMON/openvdb\n\
#tags: vdb\n\
\n\
\"\"\"Extrapolate a VDB signed distance field.\"\"\"\n\
\n\
@overview\n\
\n\
This node extrapolates signed distance fields stored as VDB volumes.\n\
Optionally, extrapolation can be masked with another VDB, so that\n\
new distances are computed only where the mask is active.\n\
\n\
@related\n\
- [OpenVDB Convert|Node:sop/DW_OpenVDBConvert]\n\
- [OpenVDB Rebuild Level Set|Node:sop/DW_OpenVDBRebuildLevelSet]\n\
- [Node:sop/isooffset]\n\
\n\
@examples\n\
\n\
See [openvdb.org|http://www.openvdb.org/download/] for source code\n\
and usage examples.\n");
}


bool
SOP_OpenVDB_Extrapolate::updateParmsFlags()
{
    UT_String mode, tmpStr;
    bool changed = false;
    evalString(mode, "mode", 0, 0);
    changed |= enableParm("extfields", mode == "scalarext" || mode == "scalarsdfext");
    changed |= enableParm("mask", mode == "mask");
    changed |= enableParm("dilate", mode == "dilate");
    changed |= enableParm("pattern", mode == "dilate");
    changed |= enableParm("isovalue", !(mode == "mask" || mode == "dilate"));
    changed |= enableParm("ignoretiles", mode == "mask");

    return changed;
}


////////////////////////////////////////


OP_Node*
SOP_OpenVDB_Extrapolate::factory(OP_Network* net,
    const char* name, OP_Operator* op)
{
    return new SOP_OpenVDB_Extrapolate(net, name, op);
}


SOP_OpenVDB_Extrapolate::SOP_OpenVDB_Extrapolate(OP_Network* net,
    const char* name, OP_Operator* op):
    hvdb::SOP_NodeVDB(net, name, op)
{
}


////////////////////////////////////////

template<typename FSGridT>
bool
SOP_OpenVDB_Extrapolate::Cache::processHelper(
    FastSweepingParms& parms,
    hvdb::GU_PrimVDB* lsPrim,
    typename FSGridT::ValueType fsIsoValue,
    const hvdb::GU_PrimVDB* maskPrim)
{
    using namespace openvdb;
    using namespace openvdb::tools;

    if (parms.mNeedExt) {
        // Get the extension primitive
        std::string tmpStr = evalStdString("extfields", parms.mTime);
        parms.mExtGroup = matchGroup(*gdp, tmpStr);
        hvdb::VdbPrimIterator extPrim(gdp, parms.mExtGroup);

        // If we want to extend a field defined by a group but we cannot find a VDB primitive
        // in that group, then throw
        if (!extPrim) {
            std::string msg = "Cannot find the correct VDB primitive named " + tmpStr + ".";
            throw std::runtime_error(msg);
        }

        hvdb::Grid::Ptr newLSGrid = nullptr;

        // Go through the extension fields and extend each one according to the chosen mode
        for (; extPrim; ++extPrim) {
            openvdb::GridBase::Ptr extGridBase = extPrim->getGridPtr();
            UT_VDBType extType = UTvdbGetGridType(*extGridBase);
            parms.mExtPrimName = extPrim.getPrimitiveNameOrIndex().toStdString();
            if (parms.mExtPrimName == parms.mFSPrimName) {
                std::string msg = "Skipping extending VDB primitive " + parms.mExtPrimName + " off the scalar " +
                    "grid " + parms.mFSPrimName + " because they are the same grid";
                addWarning(SOP_MESSAGE, msg.c_str());
                continue;
            }
            switch (extType) {
                case UT_VDB_FLOAT:
                {
                    openvdb::FloatGrid::Ptr extGrid = openvdb::gridPtrCast<openvdb::FloatGrid>(extGridBase);
                    if (extGrid) {
                        float extBg = static_cast<float>(extGrid->background());
                        process<FSGridT, openvdb::FloatGrid>(parms, lsPrim, fsIsoValue, nullptr /*=maskPrim*/, *extPrim, extBg, newLSGrid);
                    } else {
                        std::string msg = "Skipping extending VDB primitive " + parms.mExtPrimName + " because of cast failure.";
                        addWarning(SOP_MESSAGE, msg.c_str());
                    }
                    break;
                } // UT_VDB_FLOAT
                case UT_VDB_DOUBLE:
                {
                    openvdb::DoubleGrid::Ptr extGrid = openvdb::gridPtrCast<openvdb::DoubleGrid>(extGridBase);
                    if (extGrid) {
                        double extBg = static_cast<double>(extGrid->background());
                        process<FSGridT, openvdb::DoubleGrid>(parms, lsPrim, fsIsoValue, nullptr /*=maskPrim*/, *extPrim, extBg, newLSGrid);
                    } else {
                        std::string msg = "Skipping extending VDB primitive " + parms.mExtPrimName + " because of cast failure.";
                        addWarning(SOP_MESSAGE, msg.c_str());
                    }
                    break;
                } // UT_VDB_DOUBLE
                case UT_VDB_INT32:
                {
                    openvdb::Int32Grid::Ptr extGrid = openvdb::gridPtrCast<openvdb::Int32Grid>(extGridBase);
                    if (extGrid) {
                        int extBg = static_cast<int>(extGrid->background());
                        process<FSGridT, openvdb::Int32Grid>(parms, lsPrim, fsIsoValue, nullptr /*=maskPrim*/, *extPrim, extBg, newLSGrid);
                    } else {
                        std::string msg = "Skipping extending VDB primitive " + parms.mExtPrimName + " because of cast failure.";
                        addWarning(SOP_MESSAGE, msg.c_str());
                    }
                    break;
                } // UT_VDB_INT32
                case UT_VDB_VEC3F:
                {
                    openvdb::Vec3SGrid::Ptr extGrid = openvdb::gridPtrCast<openvdb::Vec3SGrid>(extGridBase);
                    if (extGrid) {
                        openvdb::Vec3f extBg = static_cast<openvdb::Vec3f>(extGrid->background());
                        process<FSGridT, openvdb::Vec3SGrid>(parms, lsPrim, fsIsoValue, nullptr /*=maskPrim*/, *extPrim, extBg, newLSGrid);
                    } else {
                        std::string msg = "Skipping extending VDB primitive " + parms.mExtPrimName + " because of cast failure.";
                        addWarning(SOP_MESSAGE, msg.c_str());
                    }
                    break;
                } // UT_VDB_VEC3F
                case UT_VDB_VEC3D:
                {
                    openvdb::Vec3DGrid::Ptr extGrid = openvdb::gridPtrCast<openvdb::Vec3DGrid>(extGridBase);
                    if (extGrid) {
                        openvdb::Vec3d extBg = static_cast<openvdb::Vec3d>(extGrid->background());
                        process<FSGridT, openvdb::Vec3DGrid>(parms, lsPrim, fsIsoValue, nullptr /*=maskPrim*/, *extPrim, extBg, newLSGrid);
                    } else {
                        std::string msg = "Skipping extending VDB primitive " + parms.mExtPrimName + " because of cast failure.";
                        addWarning(SOP_MESSAGE, msg.c_str());
                    }
                    break;
                } // UT_VDB_VEC3D
                case UT_VDB_VEC3I:
                {
                    openvdb::Vec3IGrid::Ptr extGrid = openvdb::gridPtrCast<openvdb::Vec3IGrid>(extGridBase);
                    if (extGrid) {
                        openvdb::Vec3i extBg = static_cast<openvdb::Vec3i>(extGrid->background());
                        process<FSGridT, openvdb::Vec3IGrid>(parms, lsPrim, fsIsoValue, nullptr /*=maskPrim*/, *extPrim, extBg, newLSGrid);
                    } else {
                        std::string msg = "Skipping extending VDB primitive " + parms.mExtPrimName + " because of cast failure.";
                        addWarning(SOP_MESSAGE, msg.c_str());
                    }
                    break;
                } // UT_VDB_VEC3I
                default:
                {
                    addWarning(SOP_MESSAGE, "Unsupported type of VDB grid chosen for extension");
                    break;
                }
            } // switch
        } // vdbprimiterator over extgroup

        // update level set primitive
        if (newLSGrid)
        hvdb::replaceVdbPrimitive(*gdp, newLSGrid, *lsPrim, true);
    } else {
        process<FSGridT>(parms, lsPrim, fsIsoValue, maskPrim);
    }
    return true;
}


////////////////////////////////////////


template<typename FSGridT, typename ExtGridT>
bool
SOP_OpenVDB_Extrapolate::Cache::process(
    const FastSweepingParms& parms,
    hvdb::GU_PrimVDB* lsPrim,
    typename FSGridT::ValueType fsIsoValue,
    const hvdb::GU_PrimVDB* maskPrim,
    hvdb::GU_PrimVDB* exPrim,
    const typename ExtGridT::ValueType& background,
    hvdb::Grid::Ptr newLSGrid)
{
    using namespace openvdb::tools;

    using SamplerT = openvdb::tools::GridSampler<ExtGridT, openvdb::tools::BoxSampler>;
    using ExtValueT = typename ExtGridT::ValueType;

    typename FSGridT::Ptr fsGrid = openvdb::gridPtrCast<FSGridT>(lsPrim->getGridPtr());

    if (parms.mNeedExt) {
        typename ExtGridT::ConstPtr extGrid = openvdb::gridConstPtrCast<ExtGridT>(exPrim->getConstGridPtr());
        if (!extGrid) {
            std::string msg = "Extension grid (" + extGrid->getName() + ") cannot be converted " +
                              "to the explicit type specified.";
            throw std::runtime_error(msg);
        }
        SamplerT sampler(*extGrid);
        DirichletSamplerOp<ExtGridT> op(extGrid, sampler);
        using OpT = DirichletSamplerOp<ExtGridT>;

        if (parms.mMode == "scalarext") {
            hvdb::Grid::Ptr outGrid;
            if (fsGrid->getGridClass() != openvdb::GRID_LEVEL_SET) {
                std::string msg = "Extending " + extGrid->getName() + " grid using " + parms.mFSPrimName + " Fog grid.";
                addWarning(SOP_MESSAGE, msg.c_str());
                outGrid = fogToExt(*fsGrid, op, background, fsIsoValue, parms.mNSweeps);
            }
            else {
                std::string msg = "Extending " + extGrid->getName() + " grid using " + parms.mFSPrimName + " level set grid.";
                addWarning(SOP_MESSAGE, msg.c_str());
                outGrid = sdfToExt(*fsGrid, op, background, fsIsoValue, parms.mNSweeps);
            }

            // Replace the primitive
            outGrid->setTransform(fsGrid->transform().copy());
            hvdb::replaceVdbPrimitive(*gdp, outGrid, *exPrim, true);
        } else if (parms.mMode == "scalarsdfext") {
            std::pair<hvdb::Grid::Ptr, hvdb::Grid::Ptr> outPair;
            if (fsGrid->getGridClass() != openvdb::GRID_LEVEL_SET) {
                std::string msg = "Extending " + extGrid->getName() + " grid using " + parms.mFSPrimName + " Fog grid.";
                addWarning(SOP_MESSAGE, msg.c_str());
                outPair = fogToSdfAndExt(*fsGrid, op, background, fsIsoValue, parms.mNSweeps);
            }
            else {
                std::string msg = "Extending " + extGrid->getName() + " grid using " + parms.mFSPrimName + " level set grid.";
                addWarning(SOP_MESSAGE, msg.c_str());
                outPair = sdfToSdfAndExt(*fsGrid, op, background, fsIsoValue, parms.mNSweeps);
            }

            // Replace the primitives
            outPair.first->setTransform(fsGrid->transform().copy());
            outPair.first->setGridClass(openvdb::GRID_LEVEL_SET);
            outPair.second->setTransform(fsGrid->transform().copy());
            hvdb::replaceVdbPrimitive(*gdp, outPair.second, *exPrim, true);
            newLSGrid = outPair.first;
        }
    } else {
        hvdb::Grid::Ptr outGrid;
        if (parms.mMode == "dilate") {
            if (fsGrid->getGridClass() != openvdb::GRID_LEVEL_SET) {
                std::string msg = "VDB primitive " + parms.mFSPrimName + " was skipped in dilation because it is not a level set.";
                addWarning(SOP_MESSAGE, msg.c_str());
                return false;
            }

            const NearestNeighbors nn =
                (parms.mPattern == "NN18") ? NN_FACE_EDGE : ((parms.mPattern == "NN26") ? NN_FACE_EDGE_VERTEX : NN_FACE);
            outGrid = dilateSdf(*fsGrid, parms.mDilate, nn, parms.mNSweeps);
        } else if (parms.mMode == "convert") {
            if (fsGrid->getGridClass() == openvdb::GRID_LEVEL_SET) {
                std::string msg = "VDB primitive " + parms.mFSPrimName + " was not converted to SDF because it is already a level set.";
                addWarning(SOP_MESSAGE, msg.c_str());
                return false;
            }

            outGrid = fogToSdf(*fsGrid, fsIsoValue, parms.mNSweeps);
            lsPrim->setVisualization(GEO_VOLUMEVIS_ISO, lsPrim->getVisIso(), lsPrim->getVisDensity());
        } else if (parms.mMode == "renormalize") {
            if (fsGrid->getGridClass() != openvdb::GRID_LEVEL_SET) {
                std::string msg = "VDB primitive " + parms.mFSPrimName + " was not renormalized because it is not a level set."
                    " You may want to convert the FOG VDB into a level set before calling this mode.";
                addWarning(SOP_MESSAGE, msg.c_str());
                return false;
            }

            outGrid = sdfToSdf(*fsGrid, fsIsoValue, parms.mNSweeps);
        } else if (parms.mMode == "mask") {
            if (fsGrid->getGridClass() != openvdb::GRID_LEVEL_SET) {
                std::string msg = "VDB primitive " + parms.mFSPrimName + " was skipped in Mask operation because it is not a level set."
                    " You may want to convert the FOG VDB into a level set before calling this mode.";
                addWarning(SOP_MESSAGE, msg.c_str());
                return false;
            }

            FastSweepingMaskOp<FSGridT> op(parms, fsGrid);
            // calling openvdb::tools::maskSdf with mask grid.
            hvdb::GEOvdbApply<hvdb::AllGridTypes>(*maskPrim, op);
            outGrid = op.mOutGrid;
        }
        // Replace the original VDB primitive with a new primitive that contains
        // the output grid and has the same attributes and group membership.
        outGrid->setGridClass(openvdb::GRID_LEVEL_SET);
        hvdb::replaceVdbPrimitive(*gdp, outGrid, *lsPrim, true);
    } // !parms.mNeedExt

    return true;
}


OP_ERROR
SOP_OpenVDB_Extrapolate::Cache::evalFastSweepingParms(OP_Context& context, FastSweepingParms& parms)
{
    const fpreal time = context.getTime();
    parms.mTime = context.getTime();

    // Get the group of level sets to process
    parms.mFSGroup = matchGroup(*gdp, evalStdString("group", time));

    evalString(parms.mMode, "mode", 0, time);

    parms.mNeedExt = (parms.mMode == "scalarext") || (parms.mMode == "scalarsdfext");
    parms.mNSweeps = static_cast<int>(evalInt("sweeps", 0, time));
    parms.mIgnoreTiles = static_cast<bool>(evalInt("ignoretiles", 0, time));

    // For dilate
    evalString(parms.mPattern, "pattern", 0, time);
    parms.mDilate = static_cast<int>(evalInt("dilate", 0, time));
    return error();
}


////////////////////////////////////////


OP_ERROR
SOP_OpenVDB_Extrapolate::Cache::cookVDBSop(OP_Context& context)
{
    try {
        // Evaluate UI parameters
        FastSweepingParms parms;
        if (evalFastSweepingParms(context, parms) >= UT_ERROR_ABORT) return error();

        // Get the mask primitive if the mode is mask
        const fpreal time = context.getTime();
        const GU_Detail* maskGeo = inputGeo(1);
        const GU_PrimVDB* maskPrim = nullptr;
        hvdb::GridCPtr maskGrid = nullptr;
        if (parms.mMode == "mask") {// selected to use a mask
            if (maskGeo) {// second input exists
                const GA_PrimitiveGroup* maskGroup = parsePrimitiveGroups(
                    evalStdString("mask", time).c_str(), GroupCreator(maskGeo));

                hvdb::VdbPrimCIterator maskIt(maskGeo, maskGroup);
                maskPrim = *maskIt;
                if (!maskPrim) {
                     addError(SOP_MESSAGE, "Mask Geometry not found.\n"
                         "Please provide a mask VDB as a second input.");
                     return error();
                }
                if (maskIt) maskGrid = maskIt->getConstGridPtr();// only use the first grid

                if (++maskIt) {
                    addWarning(SOP_MESSAGE, "Multiple Mask grids were found.\n"
                       "Using the first one for reference.");
                }
            } else {
              addError(SOP_MESSAGE, "Mask Geometry not found.\n"
                  "Please provide a mask VDB as a second input");
              return error();
            }
        }

        UT_AutoInterrupt progress("Performing Fast Sweeping");

        // Go through the VDB primitives and process them
        // We only process grids of types UT_VDB_FLOAT and UT_VDB_DOUBLE in Fast Sweeping.
        for (hvdb::VdbPrimIterator it(gdp, parms.mFSGroup); it; ++it) {
            if (progress.wasInterrupted()) {
                throw std::runtime_error("Processing was interrupted");
            }
            hvdb::Grid& inGrid = it->getGrid();
            UT_VDBType inType = UTvdbGetGridType(inGrid);
            parms.mFSPrimName = it.getPrimitiveNameOrIndex().toStdString();

            switch (inType) {
                case UT_VDB_FLOAT:
                {
                    float isoValue = static_cast<float>(evalFloat("isovalue", 0, time));
                    processHelper<openvdb::FloatGrid>(parms, *it /*lsPrim*/, isoValue, maskPrim);
                    parms.mExtFieldProcessed = true;
                    break;
                }
                case UT_VDB_DOUBLE:
                {
                    double isoValue = static_cast<double>(evalFloat("isovalue", 0, time));
                    processHelper<openvdb::DoubleGrid>(parms, *it /*lsPrim*/, isoValue, maskPrim);
                    parms.mExtFieldProcessed = true;
                    break;
                }
                default:
                {
                    std::string msg = "VDB primitive " + parms.mFSPrimName + " was skipped because it is not a floating-point Grid.";
                    addWarning(SOP_MESSAGE, msg.c_str());
                    break;
                }
            }

            // If we need extension, we only process the first grid
            if (parms.mNeedExt && parms.mExtFieldProcessed) {
                break;
            }
        }
    } catch (std::exception& e) {
        addError(SOP_MESSAGE, e.what());
    }

    return error();
}