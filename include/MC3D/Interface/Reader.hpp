#ifndef MC3D_READER_HPP
#define MC3D_READER_HPP

#include "MC3D/Mesh/TetMeshManipulator.hpp"
#include "MC3D/Mesh/TetMeshProps.hpp"

#include <fstream>
#include <string>

namespace mc3d
{

class Reader : public TetMeshManipulator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        FILE_INACCESSIBLE = 1, // Could not access file
        MISSING_VERTICES = 2,  // Less vertices than expected
        MISSING_TETS = 3,      // Less tets than expected
        MISSING_CHART = 4,     // Missing UVW coordinates per tet
        INVALID_CHART = 5,     // Parametrically inverted or degenerate tet
        MISSING_WALLS = 6,     // Less wall faces than expected
        INVALID_WALLS = 7,     // Read wallface is not part of the read mesh
    };

    /**
     * @brief Create a reader, that reads a mesh with parametrization in .hexex format
     *  from \p fileName and writes it into \p meshProps.
     *
     * @param meshProps OUT: this will contain the read mesh
     * @param fileName IN: file to read
     * @param forceSanitization IN: whether sanitization of input should be forced, even for exact rational input
     */
    Reader(TetMeshProps& meshProps, const std::string& fileName, bool forceSanitization = false);

    /**
     * @brief Read only the mesh structure and parametrization from the given file.
     * Allocates properties: CHART
     *
     * @return RetCode SUCCESS or an error code
     */
    RetCode readSeamlessParam();

    /**
     * @brief Read mesh structure, parametrization, walls and walldist from the given file.
     * Allocates properties: CHART, IS_WALL and WALL_DIST
     *
     * Expects the following appendix of the otherwise standard .hexex-File:
     * N_WALLS(int)
     * FACE1_VTX1(int) FACE1_VTX2(int) FACE1_VTX3(int) DIST(double)
     * FACE2_VTX1(int) FACE2_VTX2(int) FACE2_VTX3(int) DIST(double)
     * ...
     *
     * @return RetCode SUCCESS or an error code
     */
    RetCode readSeamlessParamWithWalls();

  private:
    const std::string _fileName;
    std::ifstream _is;
    bool _forceSanitization;
    bool _exactInput = false;

    /**
     * @brief Check if file is readable
     *
     * @return RetCode SUCCESS or FILE_INACCESSIBLE
     */
    RetCode checkFile();

    /**
     * @brief Read vertices from internal stream
     *
     * @return RetCode SUCCESS or MISSING_VERTICES
     */
    RetCode readVertices();

    /**
     * @brief Read tets and their charts from internal stream
     *
     * @return RetCode SUCCESS or MISSING_TETS or MISSING_CHART or INVALID_CHART
     */
    RetCode readTetsAndCharts();

    /**
     * @brief Read feature markers from internal stream
     *
     * @return RetCode SUCCESS
     */
    RetCode readFeatures();

    /**
     * @brief Read walls from internal stream
     *
     * @return RetCode SUCCESS or MISSING_WALLS or INVALID_WALLS
     */
    RetCode readWalls();

    /**
     * @brief Sanitizes the input parametrization (make input truly seamless)
     *
     * @return RetCode SUCCESS or INVALID_MAP
     */
    RetCode sanitizeInput();
};

} // namespace mc3d

#endif
