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

// Add new items to the *end* of this list, and update NUM_DATA_TYPES.
enum DataType {
    TYPE_FLOAT = 0,
    TYPE_DOUBLE,
    TYPE_INT,
    TYPE_BOOL,
    TYPE_VEC3S,
    TYPE_VEC3D,
    TYPE_VEC3I
};

enum { NUM_DATA_TYPES = TYPE_VEC3I + 1 };

std::string
dataTypeToString(DataType dts)
{
    std::string ret;
    switch (dts) {
        case TYPE_FLOAT:  ret = "float"; break;
        case TYPE_DOUBLE: ret = "double"; break;
        case TYPE_INT:    ret = "int"; break;
        case TYPE_BOOL:   ret = "bool"; break;
        case TYPE_VEC3S:  ret = "vec3s"; break;
        case TYPE_VEC3D:  ret = "vec3d"; break;
        case TYPE_VEC3I:  ret = "vec3i"; break;
    }
    return ret;
}

std::string
dataTypeToMenuItems(DataType dts)
{
    std::string ret;
    switch (dts) {
        case TYPE_FLOAT:  ret = "float"; break;
        case TYPE_DOUBLE: ret = "double"; break;
        case TYPE_INT:    ret = "int"; break;
        case TYPE_BOOL:   ret = "bool"; break;
        case TYPE_VEC3S:  ret = "vec3s (float)"; break;
        case TYPE_VEC3D:  ret = "vec3d (double)"; break;
        case TYPE_VEC3I:  ret = "vec3i (int)"; break;
    }
    return ret;
}

DataType
stringToDataType(const std::string& s)
{
    DataType ret = TYPE_FLOAT;
    std::string str = s;
    //hboost::trim(str);
    //hboost::to_lower(str);
    if (str == dataTypeToString(TYPE_FLOAT)) {
        ret = TYPE_FLOAT;
    } else if (str == dataTypeToString(TYPE_DOUBLE)) {
        ret = TYPE_DOUBLE;
    } else if (str == dataTypeToString(TYPE_INT)) {
        ret = TYPE_INT;
    } else if (str == dataTypeToString(TYPE_BOOL)) {
        ret = TYPE_BOOL;
    } else if (str == dataTypeToString(TYPE_VEC3S)) {
        ret = TYPE_VEC3S;
    } else if (str == dataTypeToString(TYPE_VEC3D)) {
        ret = TYPE_VEC3D;
    } else if (str == dataTypeToString(TYPE_VEC3I)) {
        ret = TYPE_VEC3I;
    }
    return ret;
}

struct FastSweepingParms {
    FastSweepingParms() :
        mFSGroup(nullptr),
        mMode(""),
        mNeedExt(false),
        mExtGridNum(0),
        mNSweeps(1),
        mIgnoreTiles(false),
        mPattern(""),
        mDilate(1)
    { }

    const GA_PrimitiveGroup* mFSGroup;
    UT_String mMode;
    bool mNeedExt;
    int mExtGridNum;
    int mNSweeps;
    bool mIgnoreTiles;
    UT_String mPattern;
    int mDilate;
};

template <typename GridT>
struct FastSweepMaskOp
{
    FastSweepMaskOp(typename GridT::ConstPtr inGrid, bool ignoreTiles, int iter)
        : mInGrid(inGrid), mIgnoreActiveTiles(ignoreTiles), mIter(iter) {}

    template<typename MaskGridType>
    void operator()(const MaskGridType& mask)
    {
        mOutGrid = openvdb::tools::maskSdf(*mInGrid, mask, mIgnoreActiveTiles, mIter);
    }

    typename GridT::ConstPtr mInGrid;
    const bool mIgnoreActiveTiles;
    const int mIter;
    typename GridT::Ptr mOutGrid;
};

struct SamplerOp
{
    using SamplerT = openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler>;

    SamplerOp( openvdb::FloatGrid::ConstPtr functorGrid, SamplerT sampler)
        : mFunctorGrid (functorGrid),
          mSampler(sampler)
    {}

    float operator()(const openvdb::Vec3R& xyz)
    {
        return mSampler.wsSample(xyz);
    }

    openvdb::FloatGrid::ConstPtr mFunctorGrid;
    SamplerT mSampler;
};

/// @brief Functor with signature [](const Vec3R &xyz)->ExtValueT that
///        defines the Dirichlet boundary condition, on the iso-surface,
///        of the field to be extended.
template<typename ExtGridT>
struct DirichletSamplerOp
{
    using ExtValueT = typename ExtGridT::ValueType;
    using SamplerT = openvdb::tools::GridSampler<ExtGridT, openvdb::tools::BoxSampler>;

    DirichletSamplerOp(typename ExtGridT::ConstPtr functorGrid, SamplerT sampler)
        : mFunctorGrid (functorGrid),
          mSampler(sampler)
    {}

    ExtValueT operator()(const openvdb::Vec3R& xyz)
    {
        return mSampler.wsSample(xyz);
    }

    typename ExtGridT::ConstPtr mFunctorGrid;
    SamplerT mSampler;
};

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
        template<typename GridT, typename ExtValueT = typename GridT::ValueType>
        bool processOld(hvdb::GridCPtr maskGrid,
            openvdb::FloatGrid::ConstPtr functorGrid,
            hvdb::GU_PrimVDB* lsPrim,
            fpreal time);

        bool process(hvdb::GU_PrimVDB* lsPrim,
            const FastSweepingParms& parms);
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
        .setTooltip("Specify a subset of the input VDB level sets.")
        .setDocumentation(
            "A subset of the input VDB level sets to be processed"
            " (see [specifying volumes|/model/volumes#group])"));

    // Extension grid 
    parms.add(hutil::ParmFactory(PRM_STRING, "extGroup", "Extension Group")
        .setChoiceList(&hutil::PrimGroupMenuInput1)
        .setTooltip("Select a subset of the input OpenVDB grids to be "
                    " extrapolated as Dirichlet boundary condition in the iso-surface.")
        .setDocumentation(
            "A subset of the input VDB to be used as (a) functor(s) defining"
            " the Dirichlet boundary conditions on the iso-surface"
            " of the field to be extended." ));

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
            "dilate",    "Dilate SDF",
            "mask",      "Extrapolate SDF Into Mask",
            "convert",   "Convert Scalar VDB Into SDF", ///< @todo move to Convert SOP
            "correct",   "Correct Approximate SDF", // by solving the Eikonal equation
            "fogext",    "Extend Scalar VDB",
            "sdfext",    "Extend SDF",
            "fogsdfext", "Convert Scalar VDB Into SDF and Compute Extension",
            "sdfsdfext", "Correct Approximate SDF and Compute Extension"
        })
        .setDefault("dilate")
        .setDocumentation(
            "The operation to perform\n\n"
            "Dilate SDF:\n"
            "    Dilates an existing signed distance filed by a specified \n"
            "    number of voxels.\n"
            "Extrapolate SDF Into Mask:\n"
            "    Expand/extrapolate an existing signed distance fild into\n"
            "    a mask.\n"
            "Convert Scalar VDB Into SDF:\n"
            "    Converts a scalar fog volume into a signed distance\n"
            "    function. Active input voxels with scalar values above\n"
            "    the given isoValue will have NEGATIVE distance\n"
            "    values on output, i.e. they are assumed to be INSIDE\n"
            "    the iso-surface.\n"
            "Correct Approximate SDF:\n"
            "    Given an existing approximate SDF it solves the Eikonal\n"
            "    equation for all its active voxels. Active input voxels\n"
            "    with a signed distance value above the given isoValue\n"
            "    will have POSITIVE distance values on output, i.e. they are\n"
            "    assumed to be OUTSIDE the iso-surface.\n"
            "Extend Scalar VDB:\n"
            "     Computes the extension of a scalar field, defined by the\n"
            "     specified functor, off an iso-surface from an input\n"
            "     FOG volume.\n"
            "Extend SDF:\n"
            "    Computes the extension of a scalar field, defined by the\n"
            "    specified functor, off an iso-surface from an input\n"
            "    SDF volume.\n"
            "Convert Scalar VDB Into SDF and Compute Extension:\n"
            "    Computes the signed distance field and the extension of a\n"
            "    scalar field, defined by the specified functor, off an\n"
            "    iso-surface from an input FOG volume.\n"
            "Correct Approximate SDF and Compute Extension:\n"
            "    Computes the signed distance field and the extension of a\n"
            "    scalar field, defined by the specified functor, off an\n"
            "    iso-surface from an input SDF volume."));

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

    parms.add(hutil::ParmFactory(PRM_FLT_J, "isovalue", "Isovalue")
         .setDefault(0.5)
         .setRange(PRM_RANGE_UI, -3, PRM_RANGE_UI, 3)
         .setTooltip("Isovalue for which the SDF is computed"));

    parms.add(hutil::ParmFactory(PRM_INT_J, "sweeps", "Iterations")
         .setDefault(1)
         .setRange(PRM_RANGE_RESTRICTED, 1, PRM_RANGE_UI, 5)
         .setTooltip(
            "The desired number of iterations of the Fast Sweeping algorithm"
            " (one is often enough)"));

    // Extension fields
    parms.add(hutil::ParmFactory(PRM_HEADING, "extensionFields", "Extension Fields")
         .setDocumentation(
             "These supply the fields to be extended."));

    // Dynamic grid menu
    hutil::ParmList gridParms;
    {
// TODO: Get rid of this
//         {   
//             // // TODO: Get rid of this
//             // // Grid class menu
//             // std::vector<std::string> items;
//             // for (int i = 0; i < openvdb::NUM_GRID_CLASSES; ++i) {
//             //     openvdb::GridClass cls = openvdb::GridClass(i);
//             //     items.push_back(openvdb::GridBase::gridClassToString(cls)); // token
//             //     items.push_back(openvdb::GridBase::gridClassToMenuName(cls)); // label
//             // }
// 
// // TODO: Get rid of this
// //             gridParms.add(hutil::ParmFactory(PRM_STRING, "gridClass#", "Class")
// //                 .setChoiceListItems(PRM_CHOICELIST_SINGLE, items)
// //                 .setTooltip("Specify how voxel values should be interpreted.")
// //                 .setDocumentation("\
// // How voxel values should be interpreted\n\
// // \n\
// // Fog Volume:\n\
// //     The volume represents a density field.  Values should be positive,\n\
// //     with zero representing empty regions.\n\
// // Level Set:\n\
// //     The volume is treated as a narrow-band signed distance field level set.\n\
// //     The voxels within a certain distance&mdash;the \"narrow band width\"&mdash;of\n\
// //     an isosurface are expected to define positive (exterior) and negative (interior)\n\
// //     distances to the surface.  Outside the narrow band, the distance value\n\
// //     is constant and equal to the band width.\n\
// // Staggered Vector Field:\n\
// //     If the volume is vector-valued, the _x_, _y_ and _z_ vector components\n\
// //     are to be treated as lying on the respective faces of voxels,\n\
// //     not at their centers.\n\
// // Other:\n\
// //     No special meaning is assigned to the volume's data.\n"));
//         }

        {   // Element type menu
            std::vector<std::string> items;
            for (int i = 0; i < NUM_DATA_TYPES; ++i) {
                items.push_back(dataTypeToString(DataType(i))); // token
                items.push_back(dataTypeToMenuItems(DataType(i))); // label
            }
            gridParms.add(hutil::ParmFactory(PRM_STRING, "elementType#", "Type")
                .setChoiceListItems(PRM_CHOICELIST_SINGLE, items)
                .setTooltip("The type of value stored at each voxel")
                .setDocumentation(
                    "The type of value stored at each voxel\n\n"
                    "VDB volumes are able to store vector values, unlike Houdini volumes,\n"
                    "which require one scalar volume for each vector component."));
        }

        // Optional grid name string
        gridParms.add(hutil::ParmFactory(PRM_STRING, "gridName#", "Name")
            .setTooltip("A name for this VDB")
            .setDocumentation("A value for the `name` attribute of this VDB primitive"));

        // Default background values
        // {
        const char* bgHelpStr = "The \"default\" value for any voxel not explicitly set";
        gridParms.add(hutil::ParmFactory(PRM_FLT_J, "bgFloat#", "Background Value")
            .setTooltip(bgHelpStr)
            .setDocumentation(bgHelpStr));
        gridParms.add(hutil::ParmFactory(PRM_INT_J, "bgInt#", "Background Value")
            .setDefault(PRMoneDefaults)
            .setTooltip(bgHelpStr)
            .setDocumentation(nullptr));
        gridParms.add(hutil::ParmFactory(PRM_INT_J, "bgBool#", "Background Value")
            .setRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_RESTRICTED, 1)
            .setDefault(PRMoneDefaults)
            .setTooltip(bgHelpStr)
            .setDocumentation(nullptr));
        gridParms.add(hutil::ParmFactory(PRM_FLT_J, "bgVec3f#", "Background Value")
            .setVectorSize(3)
            .setTooltip(bgHelpStr)
            .setDocumentation(nullptr));
        gridParms.add(hutil::ParmFactory(PRM_INT_J, "bgVec3i#", "Background Value")
            .setVectorSize(3)
            .setTooltip(bgHelpStr)
            .setDocumentation(nullptr));
        // }
    }

    parms.add(hutil::ParmFactory(PRM_MULTITYPE_LIST, "extGridNum", "VDBs")
        .setMultiparms(gridParms)
        .setDefault(PRMoneDefaults));

    hvdb::OpenVDBOpFactory("OpenVDB Extrapolate", SOP_OpenVDB_Extrapolate::factory, parms, *table)
        .addInput("Level Set VDB")
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
    changed |= enableParm("mask", mode == "mask");
    changed |= enableParm("dilate", mode == "dilate");
    changed |= enableParm("pattern", mode == "dilate");
    changed |= enableParm("isovalue", mode == "convert");
    changed |= enableParm("ignoretiles", mode != "dilate");

    bool needExt = mode == "fogext" || mode == "sdfext" || mode == "fogsdfext" || mode == "sdfsdfext";
    for (int i = 1, N = static_cast<int>(evalInt("extGridNum", 0, 0)); i <= N; ++i) {
        // TODO: Get rid of this
        // evalStringInst("gridClass#", &i, tmpStr, 0, 0);
        // openvdb::GridClass gridClass = openvdb::GridBase::stringToGridClass(tmpStr.toStdString());

        evalStringInst("elementType#", &i, tmpStr, 0, 0);
        DataType eType = stringToDataType(tmpStr.toStdString());
        // TODO: Get rid of this
        // bool isLevelSet = false;

        // TODO: Get rid of this
        // Force a specific data type for some of the grid classes
        // if (gridClass == openvdb::GRID_LEVEL_SET) {
        //     eType = TYPE_FLOAT;
        //     isLevelSet = true;
        // } else if (gridClass == openvdb::GRID_FOG_VOLUME) {
        //     eType = TYPE_FLOAT;
        // } else if (gridClass == openvdb::GRID_STAGGERED) {
        //     eType = TYPE_VEC3S;
        // }

        // Disable unused bg value options
        changed |= enableParmInst("bgFloat#", &i, (eType == TYPE_FLOAT || eType == TYPE_DOUBLE) && needExt);
        changed |= enableParmInst("bgInt#",   &i, (eType == TYPE_INT || eType == TYPE_BOOL) && needExt);
        changed |= enableParmInst("bgVec3f#", &i, (eType == TYPE_VEC3S || eType == TYPE_VEC3D) && needExt);
        changed |= enableParmInst("bgVec3i#", &i, eType == TYPE_VEC3I && needExt);
        changed |= enableParmInst("vecType#", &i, eType >= TYPE_VEC3S && needExt);

        // Hide unused bg value options.
        changed |= setVisibleStateInst("bgFloat#", &i, eType == TYPE_FLOAT || eType == TYPE_DOUBLE);
        changed |= setVisibleStateInst("bgInt#",   &i, eType == TYPE_INT);
        changed |= setVisibleStateInst("bgBool#",  &i, eType == TYPE_BOOL);
        changed |= setVisibleStateInst("bgVec3f#", &i, eType == TYPE_VEC3S || eType == TYPE_VEC3D);
        changed |= setVisibleStateInst("bgVec3i#", &i, eType == TYPE_VEC3I);
        changed |= setVisibleStateInst("vecType#", &i, eType >= TYPE_VEC3S);

        // Disable all these parameters if we don't need extension
        // changed |= enableParmInst("gridClass#", &i, needExt);
        changed |= enableParmInst("elementType#", &i, needExt);
        changed |= enableParmInst("gridName#", &i, needExt);

        //// Enable different data types
        //changed |= enableParmInst("elementType#", &i, gridClass == openvdb::GRID_UNKNOWN);
        //changed |= setVisibleStateInst("elementType#", &i, gridClass == openvdb::GRID_UNKNOWN);
    }
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


bool
SOP_OpenVDB_Extrapolate::Cache::process(
    hvdb::GU_PrimVDB* lsPrim,
    const FastSweepingParms& parms)
{
    using namespace openvdb::tools;

    // typename GridT::ConstPtr inGrid = openvdb::gridConstPtrCast<GridT>(lsPrim->getConstGridPtr());
    // typename GridT::Ptr outGrid;

    if (parms.mNeedExt) {
        using SamplerT = openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler>;
        //if (parms.mMode == "fogext") {
        //    SamplerT sampler(*functorGrid);
        //    SamplerOp op(functorGrid, sampler);
        //    outGrid = fogToExt(*inGrid, op, isoValue, parms.mNSweeps);
        //} else if (parms.mMode == "sdfext") {
        //    SamplerT sampler(*functorGrid);
        //    SamplerOp op(functorGrid, sampler);
        //    outGrid = sdfToExt(*inGrid, op, isoValue, parms.mNSweeps);
        //} else if (parms.mMode == "fogsdfext") {
        //    SamplerT sampler(*functorGrid);
        //    SamplerOp op(functorGrid, sampler);
        //    // std::array<typename GridT::Ptr, 2>
        //    fogToSdfAndExt(*inGrid, op, isoValue, parms.mNSweeps);
        //} else if (parms.mMode == "sdfsdfext") {
        //    SamplerT sampler(*functorGrid);
        //    SamplerOp op(functorGrid, sampler);
        //    // std::array<typename GridT::Ptr, 2>
        //    sdfToSdfAndExt(*inGrid, op, isoValue, parms.mNSweeps);
        //}
    } else {
        if (parms.mMode == "dilate") {
            const NearestNeighbors nn =
                (parms.mPattern == "NN18") ? NN_FACE_EDGE : ((parms.mPattern == "NN26") ? NN_FACE_EDGE_VERTEX : NN_FACE);
            //outGrid = dilateSdf(*inGrid, parms.mDilate, nn, nSweeps);
        } else if (parms.mMode == "convert") {

        } else if (parms.mMode == "correct") {

        } else if (parms.mMode == "mask") {

        }
    //} else if (mode == "dilate") {
    //    UT_String str;
    //    evalString(str, "pattern", 0, time);
    //    const NearestNeighbors nn =
    //        (str == "NN18") ? NN_FACE_EDGE : ((str == "NN26") ? NN_FACE_EDGE_VERTEX : NN_FACE);
    //    outGrid = dilateSdf(*inGrid, static_cast<int>(evalInt("dilate", 0, time)), nn, parms.mNSweeps);
    //} 

    }

    //const float isoValue = evalFloat("isovalue", 0, time);

    //if (parms.mMode == "mask") {
    //    FastSweepMaskOp<GridT> op(inGrid, parms.mIgnoreTiles, parms.mNSweeps);
    //    UTvdbProcessTypedGridTopology(UTvdbGetGridType(*maskGrid), *maskGrid, op);
    //    outGrid = op.mOutGrid;
    //} else if (mode == "dilate") {
    //    UT_String str;
    //    evalString(str, "pattern", 0, time);
    //    const NearestNeighbors nn =
    //        (str == "NN18") ? NN_FACE_EDGE : ((str == "NN26") ? NN_FACE_EDGE_VERTEX : NN_FACE);
    //    outGrid = dilateSdf(*inGrid, static_cast<int>(evalInt("dilate", 0, time)), nn, parms.mNSweeps);
    //} 
    ////else if (mode == "convert") {
    ////    outGrid = fogToSdf(*inGrid, isoValue, parms.mNSweeps);
    ////    lsPrim->setVisualization(GEO_VOLUMEVIS_ISO, lsPrim->getVisIso(), lsPrim->getVisDensity());
    ////} else if (mode == "correct") {
    ////    outGrid = sdfToSdf(*inGrid, isoValue, parms.mNSweeps);
    ////} else if (mode == "fogext") {
    ////    SamplerT sampler(*functorGrid);
    ////    SamplerOp op(functorGrid, sampler);
    ////    outGrid = fogToExt(*inGrid, op, isoValue, parms.mNSweeps);
    ////} else if (mode == "sdfext") {
    ////    SamplerT sampler(*functorGrid);
    ////    SamplerOp op(functorGrid, sampler);
    ////    outGrid = sdfToExt(*inGrid, op, isoValue, parms.mNSweeps);
    ////} else if (mode == "fogsdfext") {
    ////    SamplerT sampler(*functorGrid);
    ////    SamplerOp op(functorGrid, sampler);
    ////    // std::array<typename GridT::Ptr, 2>
    ////    fogToSdfAndExt(*inGrid, op, isoValue, parms.mNSweeps);
    ////} else if (mode == "sdfsdfext") {
    ////    SamplerT sampler(*functorGrid);
    ////    SamplerOp op(functorGrid, sampler);
    ////    // std::array<typename GridT::Ptr, 2>
    ////    sdfToSdfAndExt(*inGrid, op, isoValue, parms.mNSweeps);
    ////}

    //// Replace the original VDB primitive with a new primitive that contains
    //// the output grid and has the same attributes and group membership.
    //// TODO: Should we add a toggle to replace/append the vdb primitive?
    //hvdb::replaceVdbPrimitive(*gdp, outGrid, *lsPrim, true);

    return true;
}

template<typename GridT, typename ExtValueT = typename GridT::ValueType>
bool
SOP_OpenVDB_Extrapolate::Cache::processOld(
    hvdb::GridCPtr maskGrid,
    openvdb::FloatGrid::ConstPtr functorGrid,
    hvdb::GU_PrimVDB* lsPrim,
    fpreal time)
{
    std::cout << "processOld" << std::endl;
    using namespace openvdb::tools;
    using SamplerT = openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler>;

    typename GridT::ConstPtr inGrid = openvdb::gridConstPtrCast<GridT>(lsPrim->getConstGridPtr());
    typename GridT::Ptr outGrid;

    const int nSweeps = static_cast<int>(evalInt("sweeps", 0, time));
    const float isoValue = evalFloat("isovalue", 0, time);
    UT_String mode;
    evalString(mode, "mode", 0, time);

    if (mode == "mask") {
        FastSweepMaskOp<GridT> op(inGrid, evalInt("ignoretiles", 0, time), nSweeps);
        UTvdbProcessTypedGridTopology(UTvdbGetGridType(*maskGrid), *maskGrid, op);
        outGrid = op.mOutGrid;
    } else if (mode == "dilate") {
        UT_String str;
        evalString(str, "pattern", 0, time);
        const NearestNeighbors nn =
            (str == "NN18") ? NN_FACE_EDGE : ((str == "NN26") ? NN_FACE_EDGE_VERTEX : NN_FACE);
        outGrid = dilateSdf(*inGrid, static_cast<int>(evalInt("dilate", 0, time)), nn, nSweeps);
    } 
    //else if (mode == "convert") {
    //    outGrid = fogToSdf(*inGrid, isoValue, nSweeps);
    //    lsPrim->setVisualization(GEO_VOLUMEVIS_ISO, lsPrim->getVisIso(), lsPrim->getVisDensity());
    //} else if (mode == "correct") {
    //    outGrid = sdfToSdf(*inGrid, isoValue, nSweeps);
    //} else if (mode == "fogext") {
    //    SamplerT sampler(*functorGrid);
    //    SamplerOp op(functorGrid, sampler);
    //    outGrid = fogToExt(*inGrid, op, isoValue, nSweeps);
    //} else if (mode == "sdfext") {
    //    SamplerT sampler(*functorGrid);
    //    SamplerOp op(functorGrid, sampler);
    //    outGrid = sdfToExt(*inGrid, op, isoValue, nSweeps);
    //} else if (mode == "fogsdfext") {
    //    SamplerT sampler(*functorGrid);
    //    SamplerOp op(functorGrid, sampler);
    //    // std::array<typename GridT::Ptr, 2>
    //    fogToSdfAndExt(*inGrid, op, isoValue, nSweeps);
    //} else if (mode == "sdfsdfext") {
    //    SamplerT sampler(*functorGrid);
    //    SamplerOp op(functorGrid, sampler);
    //    // std::array<typename GridT::Ptr, 2>
    //    sdfToSdfAndExt(*inGrid, op, isoValue, nSweeps);
    //}

    // Replace the original VDB primitive with a new primitive that contains
    // the output grid and has the same attributes and group membership.
    // TODO: Should we add a toggle to replace/append the vdb primitive?
    hvdb::replaceVdbPrimitive(*gdp, outGrid, *lsPrim, true);

    return true;
}


OP_ERROR
SOP_OpenVDB_Extrapolate::Cache::evalFastSweepingParms(OP_Context& context, FastSweepingParms& parms)
{
    const fpreal time = context.getTime();

    // Get the group of level sets to process
    parms.mFSGroup = matchGroup(*gdp, evalStdString("group", time));

    evalString(parms.mMode, "mode", 0, time);

    parms.mNeedExt = (parms.mMode == "fogext") || (parms.mMode == "sdfext") || (parms.mMode == "fogsdfext") || (parms.mMode == "sdfsdfext");
    parms.mExtGridNum = static_cast<int>(evalInt("extGridNum", 0, 0));
    parms.mNSweeps = static_cast<int>(evalInt("sweeps", 0, time));
    parms.mIgnoreTiles = static_cast<bool>(evalInt("ignoretiles", 0, time));

    // For dilate
    evalString(parms.mPattern, "pattern", 0, time);
    parms.mDilate = static_cast<int>(evalInt("dilate", 0, time));
    return error();
}

OP_ERROR
SOP_OpenVDB_Extrapolate::Cache::cookVDBSop(OP_Context& context)
{
    try {
        // Evaluate UI parameters
        FastSweepingParms parms;
        if (evalFastSweepingParms(context, parms) >= UT_ERROR_ABORT) return error();

        const fpreal time = context.getTime();

        const GU_Detail* maskGeo = inputGeo(1);

        // get a mask grid if the mode is mask
        hvdb::GridCPtr maskGrid;
        if (parms.mMode == "mask") {// selected to use a mask
            if (maskGeo) {// second input exists
                const GA_PrimitiveGroup* maskGroup = parsePrimitiveGroups(
                    evalStdString("mask", time).c_str(), GroupCreator(maskGeo));

                hvdb::VdbPrimCIterator maskIt(maskGeo, maskGroup);

                if (maskIt) maskGrid = maskIt->getConstGridPtr();// only use the first grid

                if (!maskGrid) {
                    addError(SOP_MESSAGE, "mask VDB not found");
                    return error();
                }
            } else {
              addError(SOP_MESSAGE, "Please provide a mask VDB to the second input");
              return error();
            }
        }

        // TODO: Finish this
        // If needExt, pick the first element of the group and go through the list of 
        // extensions
        if (parms.mNeedExt) {
            hvdb::VdbPrimCIterator vdbIter(gdp, parms.mFSGroup);
            UT_String tmpStr;

            for (int i = 1; i <= parms.mExtGridNum; ++i) {
                evalStringInst("elementType#", &i, tmpStr, 0, time);
                DataType eType = stringToDataType(tmpStr.toStdString());

                switch(eType) {
                    case TYPE_FLOAT:
                    {
                        float background = float(evalFloatInst("bgFloat#", &i, 0, time));

                        // createNewGrid<cvdb::FloatGrid>(
                        //     gridNameStr, background, transform, maskGrid, group, gridClass);
                        break;
                    } // TYPE_FLOAT
                    case TYPE_DOUBLE:
                    {
                        double background = double(evalFloatInst("bgFloat#", &i, 0, time));
                        // createNewGrid<cvdb::DoubleGrid>(
                        //     gridNameStr, background, transform, maskGrid, group, gridClass);
                        break;
                    } // TYPE_DOUBLE
                    case TYPE_INT:
                    {
                        int background = static_cast<int>(evalIntInst("bgInt#", &i, 0, time));
                        // createNewGrid<cvdb::Int32Grid>(
                        //     gridNameStr, background, transform, maskGrid, group, gridClass);
                        break;
                    } // TYPE_INT
                    case TYPE_BOOL:
                    {
                        bool background = evalIntInst("bgBool#", &i, 0, time);
                        // createNewGrid<cvdb::BoolGrid>(
                        //     gridNameStr, background, transform, maskGrid, group, gridClass);
                        break;
                    } // TYPE_BOOL
                    case TYPE_VEC3S:
                    {
                        openvdb::Vec3f background(
                            float(evalFloatInst("bgVec3f#", &i, 0, time)),
                            float(evalFloatInst("bgVec3f#", &i, 1, time)),
                            float(evalFloatInst("bgVec3f#", &i, 2, time)));

                        int vecType = static_cast<int>(evalIntInst("vecType#", &i, 0, time));

                        // createNewGrid<cvdb::Vec3SGrid>(
                        //     gridNameStr, background, transform, maskGrid, group, gridClass, vecType);
                        break;
                    } // TYPE_VEC3S
                    case TYPE_VEC3D:
                    {
                        openvdb::Vec3d background(
                            double(evalFloatInst("bgVec3f#", &i, 0, time)),
                            double(evalFloatInst("bgVec3f#", &i, 1, time)),
                            double(evalFloatInst("bgVec3f#", &i, 2, time)));

                        int vecType = static_cast<int>(evalIntInst("vecType#", &i, 0, time));

                        // createNewGrid<cvdb::Vec3DGrid>(
                        //     gridNameStr, background, transform, maskGrid, group, gridClass, vecType);
                        break;
                    } // TYPE_VEC3D
                    case TYPE_VEC3I:
                    {
                        openvdb::Vec3i background(
                            static_cast<openvdb::Int32>(evalIntInst("bgVec3i#", &i, 0, time)),
                            static_cast<openvdb::Int32>(evalIntInst("bgVec3i#", &i, 1, time)),
                            static_cast<openvdb::Int32>(evalIntInst("bgVec3i#", &i, 2, time)));
                        int vecType = static_cast<int>(evalIntInst("vecType#", &i, 0, time));
                        // createNewGrid<cvdb::Vec3IGrid>(
                        //     gridNameStr, background, transform, maskGrid, group, gridClass, vecType);
                        break;
                    } // TYPE_VEC3I
                } // eType switch
            } // end for

            // Add warning if there are multiple fast sweeping grids
            if (++vdbIter) {
                addWarning(SOP_MESSAGE, "Multiple Fast Sweeping grids were found.\n"
                   "Using the first one for reference.");
            }
        } // parms.mNeedExt
        else {
            for (hvdb::VdbPrimIterator it(gdp, parms.mFSGroup); it; ++it) {
                //std::cout << "float" << std::endl;
                //std::cout << "processOld = " << this->processOld<openvdb::FloatGrid >(maskGrid, nullptr, *it, time) << std::endl;
                //std::cout << "double" << std::endl;
                //std::cout << "processOld = " << this->processOld<openvdb::DoubleGrid >(maskGrid, nullptr, *it, time) << std::endl;
                if (!this->processOld<openvdb::FloatGrid >(maskGrid, nullptr, *it, time) &&
                    !this->processOld<openvdb::DoubleGrid>(maskGrid, nullptr, *it, time) ) {
                    std::string s = it.getPrimitiveNameOrIndex().toStdString();
                    s = "VDB primitive " + s + " was skipped because it is not a floating-point Grid.";
                    addWarning(SOP_MESSAGE, s.c_str());
                    continue;
                }
            }
        } // !parms.mNeedExt






        /*    switch(eType) {
                case TYPE_FLOAT:
                {
                    float voxelSize = float(transform->voxelSize()[0]);
                    float background = 0.0;

                    if (gridClass == openvdb::GRID_LEVEL_SET) {
                        background = float(evalFloatInst("width#", &i, 0, time) * voxelSize);
                    } else {
                        background = float(evalFloatInst("bgFloat#", &i, 0, time));
                    }

                    createNewGrid<cvdb::FloatGrid>(
                        gridNameStr, background, transform, maskGrid, group, gridClass);
                    break;
                }
                case TYPE_DOUBLE:
                {
                    double background = double(evalFloatInst("bgFloat#", &i, 0, time));
                    createNewGrid<cvdb::DoubleGrid>(
                        gridNameStr, background, transform, maskGrid, group, gridClass);
                    break;
                }
                case TYPE_INT:
                {
                    int background = static_cast<int>(evalIntInst("bgInt#", &i, 0, time));
                    createNewGrid<cvdb::Int32Grid>(
                        gridNameStr, background, transform, maskGrid, group, gridClass);
                    break;
                }
                case TYPE_BOOL:
                {
                    bool background = evalIntInst("bgBool#", &i, 0, time);
                    createNewGrid<cvdb::BoolGrid>(
                        gridNameStr, background, transform, maskGrid, group, gridClass);
                    break;
                }
                case TYPE_VEC3S:
                {
                    cvdb::Vec3f background(
                        float(evalFloatInst("bgVec3f#", &i, 0, time)),
                        float(evalFloatInst("bgVec3f#", &i, 1, time)),
                        float(evalFloatInst("bgVec3f#", &i, 2, time)));

                    int vecType = static_cast<int>(evalIntInst("vecType#", &i, 0, time));

                    createNewGrid<cvdb::Vec3SGrid>(
                        gridNameStr, background, transform, maskGrid, group, gridClass, vecType);
                    break;
                }
                case TYPE_VEC3D:
                {
                    cvdb::Vec3d background(
                        double(evalFloatInst("bgVec3f#", &i, 0, time)),
                        double(evalFloatInst("bgVec3f#", &i, 1, time)),
                        double(evalFloatInst("bgVec3f#", &i, 2, time)));

                    int vecType = static_cast<int>(evalIntInst("vecType#", &i, 0, time));

                    createNewGrid<cvdb::Vec3DGrid>(
                        gridNameStr, background, transform, maskGrid, group, gridClass, vecType);
                    break;
                }
                case TYPE_VEC3I:
                {
                    cvdb::Vec3i background(
                        static_cast<cvdb::Int32>(evalIntInst("bgVec3i#", &i, 0, time)),
                        static_cast<cvdb::Int32>(evalIntInst("bgVec3i#", &i, 1, time)),
                        static_cast<cvdb::Int32>(evalIntInst("bgVec3i#", &i, 2, time)));
                    int vecType = static_cast<int>(evalIntInst("vecType#", &i, 0, time));
                    createNewGrid<cvdb::Vec3IGrid>(
                        gridNameStr, background, transform, maskGrid, group, gridClass, vecType);
                    break;
                }
            } // eType switch






















        openvdb::FloatGrid::ConstPtr functorGrid = nullptr;
        if (parms.mMode == "fogext" || parms.mMode == "sdfext" || parms.mMode == "fogsdfext" || parms.mMode == "sdfsdfext") {
            const GA_PrimitiveGroup* functorGroup = matchGroup(*gdp, evalStdString("extGroup", time));

            hvdb::VdbPrimCIterator functorIt(gdp, functorGroup);
            if (functorIt) {
                if (functorIt->getStorageType() != UT_VDB_FLOAT) {
                    addError(SOP_MESSAGE, "Unsupported grid type for functor.");
                    return UT_ERROR_ABORT;
                }
                functorGrid = hvdb::Grid::constGrid<openvdb::FloatGrid>(functorIt->getConstGridPtr());
            }

            if (!functorGrid) {
                addError(SOP_MESSAGE, "Missing functor grid");
                return UT_ERROR_ABORT;
            }
        }
            for (hvdb::VdbPrimIterator it(gdp, parms.mFSGroup); it; ++it) {
                if (!this->process<openvdb::FloatGrid >(maskGrid, functorGrid, *it, time) &&
                    !this->process<openvdb::DoubleGrid>(maskGrid, functorGrid, *it, time) ) {
                    std::string s = it.getPrimitiveNameOrIndex().toStdString();
                    s = "VDB primitive " + s + " was skipped because it is not a floating-point Grid.";
                    addWarning(SOP_MESSAGE, s.c_str());
                    continue;
                }
            }*/

    } catch (std::exception& e) {
        addError(SOP_MESSAGE, e.what());
    }

    return error();
}