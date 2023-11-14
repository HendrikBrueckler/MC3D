#ifndef MC3D_MOTORCYCLE_HPP
#define MC3D_MOTORCYCLE_HPP

#include "MC3D/Types.hpp"

#include <queue>

namespace mc3d
{

struct Motorcycle
{
    size_t ID;

    CH tet;  // Tet into which the motorcycle propagates
    EH edge; // edge from which the motorcycle propagates

    Vec3i encodedCoords; // isoCoord and propagationCoord combined (for rotation via transition)

    Q isoValue;   // isovalue in coordinate system of tet
    Q startValue; // startvalue in coordinate system of tet
    Q dist;       // total distance travelled
    Q directDist; // direct distance from origin

    /**
     * @brief Create a motorcycle
     *
     * @param tet_ IN: tet to propagate into
     * @param edge_ IN: edge to propagate from
     * @param encodedCoords_ IN: coordinates as encoded by encodeCoords()
     * @param isoValue_ IN: isovalue in coordinate system of \p tet
     * @param startValue_ IN: startvalue (of propagationCoord) in coordinate system of \p tet
     * @param dist_ IN: starting distance (from origin) of the motorcycle
     */
    Motorcycle(CH tet_, EH edge_, Vec3i encodedCoords_, Q isoValue_, Q startValue_, Q dist_, Q directDist_);

    // Keys used to encode coordinates in a vector
    const static unsigned DEFAULT_KEY = 0;
    const static unsigned ISO_KEY = 1;
    const static unsigned PROPAGATION_KEY = 2;

    /**
     * @brief Get the iso-coordinate and the propagation-coordinate from an encoded vector
     *
     * @param encodedCoords_ IN: encoded vector
     * @return std::pair<int, int> first: iso-coord, second: propagation-coord
     */
    static std::pair<int, int> decodeCoords(const Vec3i& encodedCoords_);

    /**
     * @brief Get the iso coordinate from an encoded vector
     *
     * @param encodedCoords_ IN: encoded vector
     * @return int iso coordinate
     */
    static int isoCoord(const Vec3i& encodedCoords_);

    /**
     * @brief Get the propagation coordinate from an encoded vector
     *
     * @param encodedCoords_ IN: encoded vector
     * @return int propagation coordinate
     */
    static int propagationCoord(const Vec3i& encodedCoords_);

    /**
     * @brief Encode the given iso and propagation coords into \p encodedCoords_
     *
     * @param encodedCoords_ OUT: encoded vector
     * @param isoCoord IN: iso coordinate
     * @param propagationCoord IN: propagation coordinate
     */
    static void encodeCoords(Vec3i& encodedCoords_, int isoCoord, int propagationCoord);

    /**
     * @brief Get the iso-coordinate and the propagation-coordinate of the motorcycle
     *
     * @return std::pair<int, int> first: iso-coord, second: propagation-coord
     */
    std::pair<int, int> decodeCoords() const;

    /**
     * @brief Get the iso coordinate of the motorcycle
     *
     * @param encodedCoords_ IN: encoded vector
     * @return int iso coordinate
     */
    int isoCoord() const;

    /**
     * @brief Get the propagation coordinate from of the motorcycle
     *
     * @param encodedCoords_ IN: encoded vector
     * @return int propagation coordinate
     */
    int propagationCoord() const;

    /**
     * @brief Encode the given iso and propagation coords into this motorcycle
     *
     * @param isoCoord IN: iso coordinate
     * @param propagationCoord IN: propagation coordinate
     */
    void encodeCoords(int isoCoord, int propagationCoord);

    /**
     * @brief Comparator to process first those motorcycles travelled the least far with a FIFO tiebreak
     *
     */
    struct ShorterDistFIFOTieBreak
    {
        bool operator()(const Motorcycle& a, const Motorcycle& b);
    };

    /**
     * @brief Comparator to process first those motorcycles travelled the least far with a LIFO tiebreak
     *
     */
    struct ShorterDistLIFOTieBreak
    {
        bool operator()(const Motorcycle& a, const Motorcycle& b);
    };
};

using MotorcycleQueue = std::priority_queue<Motorcycle, std::deque<Motorcycle>, Motorcycle::ShorterDistFIFOTieBreak>;

} // namespace mc3d

#endif
