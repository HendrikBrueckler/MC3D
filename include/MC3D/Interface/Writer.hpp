#ifndef MC3D_WRITER_HPP
#define MC3D_WRITER_HPP

#include "MC3D/Mesh/TetMeshNavigator.hpp"

#include <fstream>
#include <map>
#include <string>

namespace mc3d
{

class Writer : public TetMeshNavigator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        FILE_INACCESSIBLE = 1, // Could not access file
        MISSING_CHART = 4,     // Some or all tets in the mesh did not have a complete parameterization
        MISSING_WALLS = 7,     // No face wall markers found, but writing faces was demanded
    };

    /**
     * @brief Create a writer, that reads a mesh with parametrization \p meshProps
     *  and writes it in .hexex format to \p meshProps.
     *
     * @param meshProps IN: this will contain the read mesh
     * @param fileName IN: file to read
     * @param exactRationalParam IN: whether the parametrization values should be represented as exact rationals
     *                               (may be "numerator/denominator" or "integer"). Rationals are stored
     *                               via the stringification of mpq_class. otherwise double precision floating
     *                               values are used for the parametrization.
     */
    Writer(const TetMeshProps& meshProps, const std::string& fileName, bool exactRationalParam = false);

    /**
     * @brief Write the mesh with parameterization to the given file.
     * Requires properties: CHART
     *
     * @return RetCode SUCCESS or errorcode
     */
    RetCode writeSeamlessParam();

    /**
     * @brief Write the mesh with parameterization to the given file.
     * Requires properties: CHART, IS_WALL
     * Optional properties: WALL_DIST
     *
     * @return RetCode SUCCESS or errorcode
     */
    RetCode writeSeamlessParamAndWalls();

    /**
     * @brief Write integer grid map to specified file in hexex format
     *
     * @return RetCode SUCCESS or errorcode
     */
    RetCode writeIGM();

    /**
     * @brief Write integer grid map to specified file in hexex format
     *
     * @return RetCode SUCCESS or errorcode
     */
    RetCode writeIGMAndWalls();

  private:
    const std::string _fileName;
    std::ofstream _os;
    bool _exact;

    map<int, int> _vtx2idx; // Internal map of vtx mesh idx to vtx output idx

    /**
     * @brief  Check if file is writeable
     *
     * @return RetCode SUCCESS or FILE_INACCESSIBLE
     */
    RetCode checkFile() const;

    /**
     * @brief Write vertex position to internal stream
     *
     * @return RetCode SUCCESS
     */
    RetCode writeVertices();

    /**
     * @brief Write tets and charts to internal stream
     *
     * @param igm IN: whether to write IGM instead of
     * @return RetCode SUCCESS or MISSING_CHART
     */
    RetCode writeTetsAndCharts(bool igm = false);

    /**
     * @brief Write wall faces to internal stream
     *
     * @return RetCode SUCCESS or MISSING_WALLS
     */
    RetCode writeWalls();

    /**
     * @brief Write features faces to internal stream
     *
     * @return RetCode SUCCESS
     */
    RetCode writeFeatures();
};

} // namespace mc3d

#endif
