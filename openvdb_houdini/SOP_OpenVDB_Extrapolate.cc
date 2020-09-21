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
#include <functional> // std::ref

namespace hvdb = openvdb_houdini;
namespace hutil = houdini_utils;


namespace {

// Add new items to the *end* of this list, and update NUM_DATA_TYPES.
enum DataType {
    TYPE_FLOAT = 0,
    TYPE_DOUBLE,
    TYPE_INT,
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
    if (str == dataTypeToString(TYPE_FLOAT)) {
        ret = TYPE_FLOAT;
    } else if (str == dataTypeToString(TYPE_DOUBLE)) {
        ret = TYPE_DOUBLE;
    } else if (str == dataTypeToString(TYPE_INT)) {
        ret = TYPE_INT;
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
        mTime(0.0f),
        mFSGroup(nullptr),
        mExtGroup(nullptr),
        mMode(""),
        mNeedExt(false),
        mExtGridNum(0),
        mNSweeps(1),
        mIgnoreTiles(false),
        mPattern(""),
        mDilate(1),
        mFSPrimName(""),
        mExtPrimName("")
    { }

    fpreal mTime;
    const GA_PrimitiveGroup* mFSGroup;
    const GA_PrimitiveGroup* mExtGroup;
    UT_String mMode;
    bool mNeedExt;
    int mExtGridNum;
    int mNSweeps;
    bool mIgnoreTiles;
    UT_String mPattern;
    int mDilate;
    std::string mFSPrimName;
    std::string mExtPrimName;
};


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
};


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
            const typename ExtGridT::ValueType& background = typename ExtGridT::ValueType(0));
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
            "dilate",      "Dilate",
            "mask",        "Mask",
            "convert",     "Convert Scalar VDB To SDF", ///< @todo move to Convert SOP
            "renormalize", "Renormalize SDF", // by solving the Eikonal equation
            "fogext",      "Extend Off Scalar VDB",
            "sdfext",      "Extend Off SDF",
            "fogsdfext",   "Convert Scalar VDB To SDF and Extend Field",
            "sdfsdfext",   "Renormalize SDF and Extend Field"
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
            "Extend Off Scalar VDB:\n"
            "     Computes the extension of a scalar field, defined by the\n"
            "     specified functor, off an iso-surface from an input\n"
            "     FOG volume.\n"
            "Extend Off SDF:\n"
            "    Computes the extension of a scalar field, defined by the\n"
            "    specified functor, off an iso-surface from an input\n"
            "    SDF volume.\n"
            "Convert Scalar VDB To SDF and Extend Field:\n"
            "    Computes the signed distance field and the extension of a\n"
            "    scalar field, defined by the specified functor, off an\n"
            "    iso-surface from an input FOG volume.\n"
            "Renormalize SDF and Extend Field:\n"
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
    changed |= enableParm("isovalue", !(mode == "mask" || mode == "dilate"));
    changed |= enableParm("ignoretiles", mode == "mask");

    bool needExt = mode == "fogext" || mode == "sdfext" || mode == "fogsdfext" || mode == "sdfsdfext";
    for (int i = 1, N = static_cast<int>(evalInt("extGridNum", 0, 0)); i <= N; ++i) {
        evalStringInst("elementType#", &i, tmpStr, 0, 0);
        DataType eType = stringToDataType(tmpStr.toStdString());

        // Disable unused bg value options
        changed |= enableParmInst("bgFloat#", &i, (eType == TYPE_FLOAT || eType == TYPE_DOUBLE) && needExt);
        changed |= enableParmInst("bgInt#",   &i, (eType == TYPE_INT) && needExt);
        changed |= enableParmInst("bgVec3f#", &i, (eType == TYPE_VEC3S || eType == TYPE_VEC3D) && needExt);
        changed |= enableParmInst("bgVec3i#", &i, eType == TYPE_VEC3I && needExt);
        changed |= enableParmInst("vecType#", &i, eType >= TYPE_VEC3S && needExt);

        // Hide unused bg value options.
        changed |= setVisibleStateInst("bgFloat#", &i, eType == TYPE_FLOAT || eType == TYPE_DOUBLE);
        changed |= setVisibleStateInst("bgInt#",   &i, eType == TYPE_INT);
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
        UT_String tmpStr;

        for (int i = 1; i <= parms.mExtGridNum; ++i) {
            // Get the extension primitive
            evalStringInst("gridName#", &i, tmpStr, 0, parms.mTime);
            parms.mExtGroup = matchGroup(*gdp, tmpStr);
            hvdb::VdbPrimIterator extPrim(gdp, parms.mExtGroup);
            if (!extPrim) {
                std::string msg = "Cannot find the correct VDB primitive (" + tmpStr.toStdString() + ")";
                throw std::runtime_error(msg);
            }
            parms.mExtPrimName = extPrim.getPrimitiveNameOrIndex().toStdString();

            evalStringInst("elementType#", &i, tmpStr, 0, parms.mTime);
            DataType eType = stringToDataType(tmpStr.toStdString());

            switch(eType) {
                case TYPE_FLOAT:
                {
                    float extBg = static_cast<float>(evalFloatInst("bgFloat#", &i, 0, parms.mTime));
                    process<FSGridT, openvdb::FloatGrid>(parms, lsPrim, fsIsoValue, nullptr /*=maskPrim*/, *extPrim, extBg);
                    break;
                } // TYPE_FLOAT
                case TYPE_DOUBLE:
                {
                    double extBg = static_cast<double>(evalFloatInst("bgFloat#", &i, 0, parms.mTime));
                    process<FSGridT, openvdb::FloatGrid>(parms, lsPrim, fsIsoValue, nullptr /*=maskPrim*/, *extPrim, extBg);
                    break;
                } // TYPE_DOUBLE
                case TYPE_INT:
                {
                    int extBg = static_cast<int>(evalIntInst("bgInt#", &i, 0, parms.mTime));
                    process<FSGridT, openvdb::FloatGrid>(parms, lsPrim, fsIsoValue, nullptr /*=maskPrim*/, *extPrim, extBg);
                    break;
                } // TYPE_INT
                case TYPE_VEC3S:
                {
                    openvdb::Vec3f extBg(
                        static_cast<float>(evalFloatInst("bgVec3f#", &i, 0, parms.mTime)),
                        static_cast<float>(evalFloatInst("bgVec3f#", &i, 1, parms.mTime)),
                        static_cast<float>(evalFloatInst("bgVec3f#", &i, 2, parms.mTime)));
                    process<FSGridT, openvdb::Vec3SGrid>(parms, lsPrim, fsIsoValue, nullptr /*=maskPrim*/, *extPrim, extBg);
                    break;
                } // TYPE_VEC3S
                case TYPE_VEC3D:
                {
                    openvdb::Vec3d extBg(
                        static_cast<double>(evalFloatInst("bgVec3f#", &i, 0, parms.mTime)),
                        static_cast<double>(evalFloatInst("bgVec3f#", &i, 1, parms.mTime)),
                        static_cast<double>(evalFloatInst("bgVec3f#", &i, 2, parms.mTime)));
                    process<FSGridT, openvdb::Vec3DGrid>(parms, lsPrim, fsIsoValue, nullptr /*=maskPrim*/, *extPrim, extBg);
                    break;
                } // TYPE_VEC3D
                case TYPE_VEC3I:
                {
                    openvdb::Vec3i extBg(
                        static_cast<openvdb::Int32>(evalIntInst("bgVec3i#", &i, 0, parms.mTime)),
                        static_cast<openvdb::Int32>(evalIntInst("bgVec3i#", &i, 1, parms.mTime)),
                        static_cast<openvdb::Int32>(evalIntInst("bgVec3i#", &i, 2, parms.mTime)));
                    process<FSGridT, openvdb::Vec3IGrid>(parms, lsPrim, fsIsoValue, nullptr /*=maskPrim*/, *extPrim, extBg);
                    break;
                } // TYPE_VEC3I
            } // eType switch
        } // end for
    } else {
        process<FSGridT>(parms, lsPrim, fsIsoValue, maskPrim);
    }
    return true;
}


////////////////////////////////////////


//template<typename FSGridT, typename ExtValueT>
template<typename FSGridT, typename ExtGridT = FSGridT>
bool
SOP_OpenVDB_Extrapolate::Cache::process(
    const FastSweepingParms& parms,
    hvdb::GU_PrimVDB* lsPrim,
    typename FSGridT::ValueType fsIsoValue,
    const hvdb::GU_PrimVDB* maskPrim,
    hvdb::GU_PrimVDB* exPrim,
    const typename ExtGridT::ValueType& background)
{
    using namespace openvdb::tools;

    using SamplerT = openvdb::tools::GridSampler<ExtGridT, openvdb::tools::BoxSampler>;
    using ExtValueT = typename ExtGridT::ValueType;

    typename FSGridT::ConstPtr fsGrid = openvdb::gridConstPtrCast<FSGridT>(lsPrim->getConstGridPtr());

    if (parms.mNeedExt) {
        typename ExtGridT::ConstPtr extGrid = openvdb::gridConstPtrCast<ExtGridT>(exPrim->getConstGridPtr());
        if (!extGrid) {
           std::string msg = "Extension grid (" + extGrid->getName() + ") cannot be converted " +
                             "to the explicit type specified";
           throw std::runtime_error(msg);
        }
        SamplerT sampler(*extGrid);
        DirichletSamplerOp<ExtGridT> op(extGrid, sampler);
        using OpT = DirichletSamplerOp<ExtGridT>;

        if (parms.mMode == "fogext" || parms.mMode == "sdfext") {
            hvdb::Grid::Ptr outGrid;
            if (parms.mMode == "fogext") {
                outGrid = fogToExt(*fsGrid, op, background, fsIsoValue, parms.mNSweeps);
            }
            else {
                outGrid = sdfToExt(*fsGrid, op, background, fsIsoValue, parms.mNSweeps);
            }

            // Create GEO_PrimVDB*
            std::string primName = parms.mExtPrimName + "Ext";
            outGrid->setTransform(fsGrid->transform().copy());
            outGrid->setName(primName);
            hvdb::createVdbPrimitive(*gdp, outGrid, primName.c_str());
        } else if (parms.mMode == "fogsdfext" || parms.mMode == "sdfsdfext") {
            std::pair<hvdb::Grid::Ptr, hvdb::Grid::Ptr> outPair;
            if (parms.mMode == "fogsdfext") {
                outPair = fogToSdfAndExt(*fsGrid, op, background, fsIsoValue, parms.mNSweeps);
            }
            else {
                outPair = sdfToSdfAndExt(*fsGrid, op, background, fsIsoValue, parms.mNSweeps);
            }

            // Create GEO_PrimVDB*
            std::string sdfPrimName = parms.mFSPrimName + "Sdf";
            std::string extPrimName = parms.mExtPrimName + "Ext";
            outPair.first->setName(sdfPrimName);
            outPair.first->setTransform(fsGrid->transform().copy());
            outPair.second->setTransform(fsGrid->transform().copy());
            outPair.second->setName(extPrimName);
            hvdb::createVdbPrimitive(*gdp, outPair.first, sdfPrimName.c_str());
            hvdb::createVdbPrimitive(*gdp, outPair.second, extPrimName.c_str());
        }
    } else {
        hvdb::Grid::Ptr outGrid;
        if (parms.mMode == "dilate") {
            // TODO: Do I need to enforce that the Grid Class to be LEVEL_SET?
            const NearestNeighbors nn =
                (parms.mPattern == "NN18") ? NN_FACE_EDGE : ((parms.mPattern == "NN26") ? NN_FACE_EDGE_VERTEX : NN_FACE);
            outGrid = dilateSdf(*fsGrid, parms.mDilate, nn, parms.mNSweeps);
        } else if (parms.mMode == "convert") {
            outGrid = fogToSdf(*fsGrid, fsIsoValue, parms.mNSweeps);
            lsPrim->setVisualization(GEO_VOLUMEVIS_ISO, lsPrim->getVisIso(), lsPrim->getVisDensity());
        } else if (parms.mMode == "renormalize") {
            if (fsGrid->getGridClass() != openvdb::GRID_LEVEL_SET) {
                throw std::runtime_error("The input grid for sdf to sdf should be a level set.");
            }
            outGrid = sdfToSdf(*fsGrid, fsIsoValue, parms.mNSweeps);
        } else if (parms.mMode == "mask") {
            FastSweepingMaskOp<FSGridT> op(parms, fsGrid);
            hvdb::GEOvdbApply<hvdb::AllGridTypes>(*maskPrim, op);
            outGrid = op.mOutGrid;
        }
        // Replace the original VDB primitive with a new primitive that contains
        // the output grid and has the same attributes and group membership.
        // TODO: Should we add a toggle to replace/append the vdb primitive?
        outGrid->setGridClass(openvdb::GRID_LEVEL_SET);
        hvdb::replaceVdbPrimitive(*gdp, outGrid, *lsPrim, true);
    } // !parms.mNeedExt

    // add various warnings
    // Mode is expecting level-sets, but the grid class is not a level set
    if (fsGrid->getGridClass() != openvdb::GRID_LEVEL_SET &&
        (parms.mMode == "dilate" || parms.mMode == "renormalize" || parms.mMode == "mask" ||
        parms.mMode == "sdfext" || parms.mMode == "sdfsdfext")) {
        std::string msg = "Grid (" + fsGrid->getName() + ")is expected to be a levelset.";
        addWarning(SOP_MESSAGE, msg.c_str());
    }
    // Mode is expecting fog, but the grid class is a level set
    if (fsGrid->getGridClass() == openvdb::GRID_LEVEL_SET &&
        (parms.mMode == "convert" || parms.mMode == "fogext" || parms.mMode == "fogsdfext")) {
        std::string msg = "Grid (" + fsGrid->getName() + ")is not expected to be a levelset.";
        addWarning(SOP_MESSAGE, msg.c_str());
    }
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

    parms.mNeedExt = (parms.mMode == "fogext") || (parms.mMode == "sdfext") || (parms.mMode == "fogsdfext") || (parms.mMode == "sdfsdfext");
    parms.mExtGridNum = static_cast<int>(evalInt("extGridNum", 0, 0));
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
                     addError(SOP_MESSAGE, "Mask Primitive not found.\n"
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

        // Go through the VDB primitives and process them
        for (hvdb::VdbPrimIterator it(gdp, parms.mFSGroup); it;) {
            hvdb::Grid& inGrid = it->getGrid();
            UT_VDBType inType = UTvdbGetGridType(inGrid);
            parms.mFSPrimName = it.getPrimitiveNameOrIndex().toStdString();

            switch (inType) {
                case UT_VDB_FLOAT:
                {
                    float isoValue = static_cast<float>(evalFloat("isovalue", 0, time));
                    processHelper<openvdb::FloatGrid>(parms, *it /*lsPrim*/, isoValue, maskPrim);
                    break;
                }
                case UT_VDB_DOUBLE:
                {
                    double isoValue = static_cast<double>(evalFloat("isovalue", 0, time));
                    processHelper<openvdb::DoubleGrid>(parms, *it /*lsPrim*/, isoValue, maskPrim);
                    break;
                }
                default:
                    std::string s = it.getPrimitiveNameOrIndex().toStdString();
                    s = "VDB primitive " + s + " was skipped because it is not a floating-point Grid.";
                    addWarning(SOP_MESSAGE, s.c_str());
                    break;
            }

            // If we need extension, we only process the first grid
            ++it;
            if (parms.mNeedExt && it) {
                addWarning(SOP_MESSAGE, "Multiple Fast Sweeping grids were found.\n"
                   "Using the first one for reference.");
                break;
            }
        }
    } catch (std::exception& e) {
        addError(SOP_MESSAGE, e.what());
    }

    return error();
}