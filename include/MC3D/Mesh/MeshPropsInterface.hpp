#ifndef MC3D_MESHWITHPROPSINTERFACE_HPP
#define MC3D_MESHWITHPROPSINTERFACE_HPP

/**
 * @brief Use this to declare a mesh-associated property to be used via \ref mc3d::MeshPropsInterface.
 *        Use this when you want an array-based property, that will allocate space for each entity and its value
 *        and allow for O(1) access.
 *
 * @param NAME the name that will be given to a property. UPPER_CASE_SNAKE_CASE is recommended. Use this as template
 *             parameter to query the property from property managers later
 * @param OVM_ENTITY_T the type of OVM mesh element that is going to be used as key.
 *                     Possible values: Vertex, Edge, HalfEdge, Face, HalfFace, Cell
 * @param VALUE_T the type of the value to be mapped to each entity. Arbitrary types are accepted, but note the
 * necessity to use MC3D_ARG in case there are "," in your type
 */
#define MC3D_PROPERTY(NAME, OVM_ENTITY_T, VALUE_T)                                                                     \
    struct NAME                                                                                                        \
    {                                                                                                                  \
        using entity_t = OVM::Entity::OVM_ENTITY_T;                                                                    \
        using value_t = VALUE_T;                                                                                       \
        using handle_t = OVM::HandleT<entity_t>;                                                                       \
        using prop_t = OVM::PropertyPtr<value_t, entity_t>;                                                             \
        using ref_t = vector<value_t>::reference;                                                                      \
        using const_ref_t = vector<value_t>::const_reference;                                                          \
        const static bool IS_MAPPED = false;                                                                           \
        std::unique_ptr<prop_t> ptr;                                                                                   \
        value_t def = value_t();                                                                                       \
        inline static std::string name()                                                                               \
        {                                                                                                              \
            return "MC3D_" #NAME;                                                                                      \
        }                                                                                                              \
    }

/**
 * @brief Use this to declare a mesh-associated property to be used via \ref mc3d::MeshPropsInterface.
 *        Use this when you want a hashmap-based property, that will allocate space for each entity and its value
 *        and allow for amortized O(1) access. This is best for sparse-(non-default)-valued properties.
 *
 * @param NAME the name that will be given to a property. UPPER_CASE_SNAKE_CASE is recommended. Use this as template
 *             parameter to query the property from property managers later
 * @param OVM_ENTITY_T the type of OVM mesh element that is going to be used as key.
 *                     Possible values: Vertex, Edge, HalfEdge, Face, HalfFace, Cell
 * @param VALUE_T the type of the value to be mapped to each non-default valued entity. Arbitrary types are accepted,
 *                but note the necessity to use MC3D_ARG in case there are "," in your type
 *
 */
#define MC3D_MAP_PROPERTY(NAME, OVM_ENTITY_T, VALUE_T)                                                                 \
    struct NAME                                                                                                        \
    {                                                                                                                  \
        using entity_t = OVM::Entity::OVM_ENTITY_T;                                                                    \
        using value_t = VALUE_T;                                                                                       \
        using handle_t = OVM::HandleT<entity_t>;                                                                       \
        using prop_t = std::unordered_map<handle_t, value_t, handle_hash<entity_t>>;                                   \
        using ref_t = value_t&;                                                                                        \
        using const_ref_t = const value_t&;                                                                            \
        const static bool IS_MAPPED = true;                                                                            \
        std::unique_ptr<prop_t> ptr;                                                                                   \
        value_t def = value_t();                                                                                       \
        inline static std::string name()                                                                               \
        {                                                                                                              \
            return "MC3D_MAPPED_" #NAME;                                                                               \
        }                                                                                                              \
    }

// For passing types containing "," in template parameter lists
#define MC3D_ARG(...) __VA_ARGS__

#include "MC3D/Types.hpp"

#include <memory>
#include <type_traits>
#include <unordered_map>

namespace mc3d
{

/**
 * @brief For use of OVM entities as keys in hashmap
 *
 * @tparam ENTITY_T OVM::Vertex/Edge/HalfEdge/Face/HalfFace/Cell
 */
template <typename ENTITY_T>
struct handle_hash
{
    std::size_t operator()(const mc3d::OVM::HandleT<ENTITY_T>& k) const
    {
        return std::hash<int>()(k.idx());
    }
};

template <typename T, typename... List>
struct is_any_of;

template <typename T, typename Head, typename... Tail>
struct is_any_of<T, Head, Tail...>
{
    enum
    {
        value = std::is_same<T, Head>::value || is_any_of<T, Tail...>::value
    };
};

template <typename T>
struct is_any_of<T>
{
    enum
    {
        value = false
    };
};

template <typename... List>
struct no_dupe;

template <typename Head, typename... Tail>
struct no_dupe<Head, Tail...>
{
    enum
    {
        value = !is_any_of<Head, Tail...>::value && no_dupe<Tail...>::value
    };
};

template <>
struct no_dupe<>
{
    enum
    {
        value = true
    };
};

template <typename T, typename... Ts>
struct Index;

template <typename T, typename... Ts>
struct Index<T, T, Ts...> : std::integral_constant<std::size_t, 0>
{
};

template <typename T, typename U, typename... Ts>
struct Index<T, U, Ts...> : std::integral_constant<std::size_t, 1 + Index<T, Ts...>::value>
{
};

/**
 * @brief Properties are classes that are passed as template parameters, and the mesh property manager is
 *        effectively a std::tuple wrapper, as that enables uniform handling of properties and avoids the
 *        need to code dozens of getters/setters/members. Inherit from instantiations of this templated
 *        interface to create mesh property manager classes.
 *        The reasoning behind using a compile-time fixed set of properties is that it is less error-prone
 *        than properties being loosely defineable at runtime: programming mistakes become obvious during
 *        compilation and not during runtime.
 *        Properties are always specific to the property manager. Multiple managers managing the same property
 *        type for the same mesh always handle unique fully-owned instances of these properties.
 *
 *
 * @tparam Mesh the type of OVM-mesh to associate a set of properties with
 * @tparam Props a list of properties (created with either MC3D_PROPERTY or MC3D_MAP_PROPERTY)
 */
template <typename Mesh, typename... Props>
class MeshPropsInterface
{
    static_assert(no_dupe<Props...>::value, "PROPERTIES NEED TO BE UNIQUE, BUT YOU PASSED DUPLICATES");

  public:
    Mesh& mesh;
    size_t id;

    /**
     * @brief Create a wrapper around \p mesh_ that manages properties of its mesh elements.
     *
     * @param mesh_ IN/OUT: mesh to augment by properties
     */
    MeshPropsInterface(Mesh& mesh_) : mesh(mesh_), id(nextID())
    {
    }

    ~MeshPropsInterface()
    {
        clearRecurse<0>();
    }

    /**
     * @brief Whether this property manager class manages a property of type \p Prop
     *
     * @tparam Prop manageable property
     * @return true if \p Prop is managed by this class
     * @return false else
     */
    template <typename Prop>
    static constexpr bool manages()
    {
        return is_any_of<Prop, Props...>::value;
    }

    /**
     * @brief Query whether a certain property is already allocated
     *
     * @tparam Prop property to query
     * @return true if already allocated in this manager
     * @return false else
     */
    template <typename Prop>
    bool isAllocated() const
    {
        return propIsAllocated<Prop>(std::get<Index<Prop, Props...>::value>(_props).ptr);
    }

    /**
     * @brief Allocate a certain type of property \p prop in this manager.
     *        A property MUST be allocated before it can be used.
     *        A property MUST be released before being allocated again.
     *
     * @tparam Prop property to allocate
     * @return Prop::prop_t& a reference to the allocated property object (map or OVM-prop-handle)
     */
    template <typename Prop>
    typename Prop::prop_t& allocate()
    {
        return allocate<Prop>(std::get<Index<Prop, Props...>::value>(_props).def);
    }

    /**
     * @brief Allocate a certain type of property \p prop in this manager with a default value of \p def.
     *        A property MUST be allocated before it can be used.
     *        A property MUST be released before being allocated again.
     *
     * @tparam Prop property to allocate
     * @param def default property value given to each element unless explicitly overwritten
     * @return Prop::prop_t& a reference to the allocated property object (map or OVM-prop-handle)
     */
    template <typename Prop>
    typename Prop::prop_t& allocate(const typename Prop::value_t& def)
    {
        if (isAllocated<Prop>())
        {
            LOG(WARNING) << "Tried to allocate an already allocated property, resetting to new default";
            release<Prop>();
        }
        setDefault<Prop>(def);
        return allocateProp<Prop>(std::get<Index<Prop, Props...>::value>(_props).ptr, def);
    }

    /**
     * @brief Set the default value of property \p prop in this manager. This MUST be called before
     *        actually allocating the property. You can not change the default value of an already
     *        allocated property!
     *
     * @tparam Prop property to set default value of
     * @param def default property value given to each element unless explicitly overwritten
     */
    template <typename Prop>
    void setDefault(const typename Prop::value_t& def)
    {
        static_assert(is_any_of<Prop, Props...>::value, "NO SUCH PROPERTY MANAGED BY THIS CLASS");
        assert(!isAllocated<Prop>());
        std::get<Index<Prop, Props...>::value>(_props).def = def;
    }

    /**
     * @brief Release a certain type of property \p prop managed by this class.
     *        A property MUST be reallocated before it can be used again.
     *        A property MUST be released before being allocated again.
     *
     * @tparam Prop property to release
     * @return Prop::prop_t& a reference to the allocated property object (map or OVM-prop-handle)
     */
    template <typename Prop>
    void release()
    {
        static_assert(is_any_of<Prop, Props...>::value, "NO SUCH PROPERTY MANAGED BY THIS CLASS");
        if (!isAllocated<Prop>())
            LOG(WARNING) << "Tried to release a non-allocated property";
        else
            clearProp<Prop>(std::get<Index<Prop, Props...>::value>(_props).ptr,
                            std::get<Index<Prop, Props...>::value>(_props).def);
    }

    /**
     * @brief Release all mesh properties and clear all elements of the managed mesh
     */
    void clearAll()
    {
        clearRecurse<0>();
        mesh.clear(false);
    }

    /**
     * @brief Retrieve the property object itself. The property allows for const access
     *        via .at() (in the map case) or operator[] (in the OVM-managed array case).
     *
     * @tparam Prop property to query
     * @return const Prop::prop_t& the property object (a map or an OVM property handle)
     */
    template <typename Prop>
    const typename Prop::prop_t& prop() const
    {
        return getProp<Prop>(std::get<Index<Prop, Props...>::value>(_props).ptr);
    }

    /**
     * @brief Retrieve the property object itself. The property allows for non-const access
     *        via operator[] .
     *
     * @tparam Prop property to query
     * @return Prop::prop_t& the property object (a map or an OVM property handle)
     */
    template <typename Prop>
    typename Prop::prop_t& prop()
    {
        return getProp<Prop>(std::get<Index<Prop, Props...>::value>(_props).ptr);
    }

    /**
     * @brief Get the mapped property value of prop \p Prop for element \p handle
     *
     * @tparam Prop property to query
     * @param handle IN: element for which the property value should be retrieved
     * @return Prop::value_t property value of \p handle
     */
    template <typename Prop, typename std::enable_if<Prop::IS_MAPPED, int>::type = 0>
    typename Prop::value_t get(const typename Prop::handle_t& handle) const
    {
        static_assert(is_any_of<Prop, Props...>::value, "NO SUCH PROPERTY MANAGED BY THIS CLASS");
        assert(handle.is_valid());
        const typename Prop::prop_t& p = prop<Prop>();
        auto it = p.find(handle);
        return it == p.end() ? std::get<Index<Prop, Props...>::value>(_props).def : it->second;
    }

    /**
     * @brief Get the mapped property value of prop \p Prop for element \p handle
     *
     * @tparam Prop property to query
     * @param handle IN: element for which the property value should be retrieved
     * @return Prop::value_t property value of \p handle
     */
    template <typename Prop, typename std::enable_if<!Prop::IS_MAPPED, int>::type = 0>
    typename Prop::value_t get(const typename Prop::handle_t& handle) const
    {
        return prop<Prop>()[handle];
    }

    /**
     * @brief Get the mesh property value of prop \p Prop for the managed mesh
     *
     * @tparam Prop property to query
     * @return Prop::value_t property value of \p handle
     */
    template <typename Prop,
              typename = typename std::enable_if<std::is_same<typename Prop::entity_t, OVM::Entity::Mesh>::value>::type>
    typename Prop::value_t get() const
    {
        return get<Prop>(typename Prop::handle_t(0));
    }

    /**
     * @brief Get the mapped property value of prop \p Prop for element \p handle by reference.
     *        This is NEVER const for map-based properties
     *
     * @tparam Prop property to query
     * @param handle IN: element for which the property value should be retrieved
     * @return Prop::ref_t reference to property value of \p handle
     */
    template <typename Prop>
    typename Prop::ref_t ref(const typename Prop::handle_t& handle)
    {
        assert(handle.is_valid());
        return prop<Prop>()[handle];
    }

    /**
     * @brief Get the mapped property value of prop \p Prop for element \p handle by reference.
     *        This is NEVER const for map-based properties
     *
     * @tparam Prop property to query
     * @param handle IN: element for which the property value should be retrieved
     * @return Prop::const_ref_t reference to property value of \p handle
     */
    template <typename Prop, typename std::enable_if<!Prop::IS_MAPPED, int>::type = 0>
    typename Prop::const_ref_t ref(const typename Prop::handle_t& handle) const
    {
        assert(handle.is_valid());
        return prop<Prop>()[handle];
    }

    /**
     * @brief Get the mesh property value of prop \p Prop for the managed mesh by reference.
     *
     * @tparam Prop property to query
     * @return Prop::ref_t reference to property value of \p handle
     */
    template <typename Prop,
              typename = typename std::enable_if<std::is_same<typename Prop::entity_t, OVM::Entity::Mesh>::value>::type>
    typename Prop::ref_t ref()
    {
        return ref<Prop>(typename Prop::handle_t(0));
    }

    /**
     * @brief Get the mesh property value of prop \p Prop for the managed mesh by reference.
     *
     * @tparam Prop property to query
     * @return Prop::ref_t reference to property value of \p handle
     */
    template <typename Prop,
              typename = typename std::enable_if<std::is_same<typename Prop::entity_t, OVM::Entity::Mesh>::value>::type,
              typename std::enable_if<!Prop::IS_MAPPED, int>::type = 0>
    typename Prop::const_ref_t ref() const
    {
        return ref<Prop>(typename Prop::handle_t(0));
    }

    /**
     * @brief Set the mapped property value of prop \p Prop for element \p handle
     *
     * @tparam Prop property to set
     * @param handle IN: element for which the property value should be set
     * @param val IN: value to set the property of \p handle to
     */
    template <typename Prop, typename std::enable_if<Prop::IS_MAPPED, int>::type = 0>
    void set(const typename Prop::handle_t& handle, const typename Prop::value_t& val)
    {
        static_assert(is_any_of<Prop, Props...>::value, "NO SUCH PROPERTY MANAGED BY THIS CLASS");
        assert(handle.is_valid());
        if (val == std::get<Index<Prop, Props...>::value>(_props).def)
        {
            typename Prop::prop_t& p = prop<Prop>();
            auto it = p.find(handle);
            if (it != p.end())
                p.erase(it);
        }
        else
        {
            prop<Prop>()[handle] = val;
        }
    }

    /**
     * @brief Set the mapped property value of prop \p Prop for element \p handle
     *
     * @tparam Prop property to set
     * @param handle IN: element for which the property value should be set
     * @param val IN: value to set the property of \p handle to
     */
    template <typename Prop, typename std::enable_if<!Prop::IS_MAPPED, int>::type = 0>
    void set(const typename Prop::handle_t& handle, const typename Prop::value_t& val)
    {
        prop<Prop>()[handle] = val;
    }

    /**
     * @brief Set the mesh property value of prop \p Prop for the managed mesh
     *
     * @tparam Prop mesh property to set
     * @param val IN: value to set the mesh property to
     */
    template <typename Prop,
              typename = typename std::enable_if<std::is_same<typename Prop::entity_t, OVM::Entity::Mesh>::value>::type>
    void set(const typename Prop::value_t& val)
    {
        set<Prop>(typename Prop::handle_t(0), val);
    }

    /**
     * @brief Reset the mapped property value of prop \p Prop for element \p handle to the default value
     *
     * @tparam Prop property to set
     * @param handle IN: element for which the property value should be reset to default
     */
    template <typename Prop>
    void reset(const typename Prop::handle_t& handle)
    {
        set<Prop>(handle, std::get<Index<Prop, Props...>::value>(_props).def);
    }

    /**
     * @brief Reset the mesh property value of prop \p Prop for the managed mesh to the default value
     *
     * @tparam Prop property to set
     */
    template <typename Prop,
              typename = typename std::enable_if<std::is_same<typename Prop::entity_t, OVM::Entity::Mesh>::value>::type>
    void reset()
    {
        set<Prop>(typename Prop::handle_t(0), std::get<Index<Prop, Props...>::value>(_props).def);
    }

    /**
     * @brief Copy the mapped property value of prop \p Prop from element \p from to \p to
     *
     * @tparam Prop property to copy
     * @param from IN: element to copy property from
     * @param to IN: element to copy property to
     */
    template <typename Prop, typename std::enable_if<Prop::IS_MAPPED, int>::type = 0>
    void clone(const typename Prop::handle_t& from, const typename Prop::handle_t& to)
    {
        set<Prop>(to, get<Prop>(from));
    }

    /**
     * @brief Copy the mapped property value of prop \p Prop from element \p from to \p to
     *
     * @tparam Prop property to copy
     * @param from IN: element to copy property from
     * @param to IN: element to copy property to
     */
    template <typename Prop, typename std::enable_if<!Prop::IS_MAPPED, int>::type = 0>
    void clone(const typename Prop::handle_t& from, const typename Prop::handle_t& to)
    {
        set<Prop>(to, ref<Prop>(from));
    }

    /**
     * @brief Copy the mapped property values of all managed properties for entity-type \p HANDLE_T
     *        from element \p from to \p to
     *
     * @tparam HANDLE_T handle type
     * @param from IN: element to copy properties from
     * @param to IN: element to copy properties to
     */
    template <typename HANDLE_T>
    void cloneAll(const HANDLE_T& from, const HANDLE_T& to)
    {
        cloneRecurse<HANDLE_T, 0>(from, to);
    }

    /**
     * @brief Reset the mapped property value of all properties for entity-type \p HANDLE_T
     *        of element \p handle to their respective default value
     *
     * @tparam HANDLE_T handle type
     * @param handle IN: element for which the property value should be reset to default
     */
    template <typename HANDLE_T>
    void resetAll(const HANDLE_T& handle)
    {
        resetRecurse<HANDLE_T, 0>(handle);
    }

    /**
     * @brief Log a summary of all managed properties using glog
     */
    void logAllProps()
    {
        logRecurse<0>();
    }

  protected:
    template <size_t I = 0, typename std::enable_if<(I < sizeof...(Props)), int>::type = 0>
    void clearRecurse()
    {
        using Prop = typename std::tuple_element<I, std::tuple<Props...>>::type;

        if (isAllocated<Prop>())
            release<Prop>();

        clearRecurse<I + 1>();
    }

    template <size_t I = 0, typename std::enable_if<(I == sizeof...(Props)), int>::type = 0>
    void clearRecurse()
    {
    }

    template <typename HANDLE_T,
              size_t I = 0,
              typename std::enable_if<(I < sizeof...(Props)), int>::type = 0,
              typename std::enable_if<
                  std::is_same<HANDLE_T, typename std::tuple_element<I, std::tuple<Props...>>::type::handle_t>::value,
                  int>::type
              = 0>
    void cloneRecurse(const HANDLE_T& from, const HANDLE_T& to)
    {
        using Prop = typename std::tuple_element<I, std::tuple<Props...>>::type;

        if (isAllocated<Prop>())
            clone<Prop>(from, to);

        cloneRecurse<HANDLE_T, I + 1>(from, to);
    }

    template <typename HANDLE_T,
              size_t I = 0,
              typename std::enable_if<(I < sizeof...(Props)), int>::type = 0,
              typename std::enable_if<
                  !std::is_same<HANDLE_T, typename std::tuple_element<I, std::tuple<Props...>>::type::handle_t>::value,
                  int>::type
              = 0>
    void cloneRecurse(const HANDLE_T& from, const HANDLE_T& to)
    {
        cloneRecurse<HANDLE_T, I + 1>(from, to);
    }

    template <typename HANDLE_T, size_t I = 0, typename std::enable_if<(I == sizeof...(Props)), int>::type = 0>
    void cloneRecurse(const HANDLE_T& from, const HANDLE_T& to)
    {
        (void)from;
        (void)to;
    }

    template <typename HANDLE_T,
              size_t I = 0,
              typename std::enable_if<(I < sizeof...(Props)), int>::type = 0,
              typename std::enable_if<
                  std::is_same<HANDLE_T, typename std::tuple_element<I, std::tuple<Props...>>::type::handle_t>::value,
                  int>::type
              = 0>
    void resetRecurse(const HANDLE_T& handle)
    {
        using Prop = typename std::tuple_element<I, std::tuple<Props...>>::type;

        if (isAllocated<Prop>())
            reset<Prop>(handle);

        resetRecurse<HANDLE_T, I + 1>(handle);
    }

    template <typename HANDLE_T,
              size_t I = 0,
              typename std::enable_if<(I < sizeof...(Props)), int>::type = 0,
              typename std::enable_if<
                  !std::is_same<HANDLE_T, typename std::tuple_element<I, std::tuple<Props...>>::type::handle_t>::value,
                  int>::type
              = 0>
    void resetRecurse(const HANDLE_T& handle)
    {
        resetRecurse<HANDLE_T, I + 1>(handle);
    }

    template <typename HANDLE_T, size_t I = 0, typename std::enable_if<(I == sizeof...(Props)), int>::type = 0>
    void resetRecurse(const HANDLE_T& handle)
    {
        (void)handle;
    }

    template <size_t I = 0, typename std::enable_if<(I < sizeof...(Props)), int>::type = 0>
    void logRecurse()
    {
        using Prop = typename std::tuple_element<I, std::tuple<Props...>>::type;

        LOG(INFO) << "HAS " << Prop::name() << id << "? " << isAllocated<Prop>();

        logRecurse<I + 1>();
    }

    template <size_t I = 0, typename std::enable_if<(I == sizeof...(Props)), int>::type = 0>
    void logRecurse()
    {
    }

    template <typename Prop>
    bool propIsAllocated(const std::unique_ptr<typename Prop::prop_t>& ptr) const
    {
        static_assert(is_any_of<Prop, Props...>::value, "NO SUCH PROPERTY MANAGED BY THIS CLASS");
        return ptr != nullptr;
    }

    template <typename Prop, typename std::enable_if<Prop::IS_MAPPED, int>::type = 0>
    typename Prop::prop_t& allocateProp(std::unique_ptr<typename Prop::prop_t>& ptr, const typename Prop::value_t& def)
    {
        (void)def;
        static_assert(is_any_of<Prop, Props...>::value, "NO SUCH PROPERTY MANAGED BY THIS CLASS");
        assert(!propIsAllocated<Prop>(ptr));
        ptr = std::unique_ptr<typename Prop::prop_t>(new typename Prop::prop_t);
        return *ptr;
    }

    template <typename Prop, typename std::enable_if<!Prop::IS_MAPPED, int>::type = 0>
    typename Prop::prop_t& allocateProp(std::unique_ptr<typename Prop::prop_t>& ptr, const typename Prop::value_t& def)
    {
        static_assert(is_any_of<Prop, Props...>::value, "NO SUCH PROPERTY MANAGED BY THIS CLASS");
        assert(!propIsAllocated<Prop>(ptr));
        ptr = std::unique_ptr<typename Prop::prop_t>(
            new typename Prop::prop_t(mesh.template request_property<typename Prop::value_t, typename Prop::entity_t>(
                Prop::name() + std::to_string(id), def)));
        mesh.set_persistent(*ptr);
        return *ptr;
    }

    template <typename Prop>
    typename Prop::prop_t& getProp(const std::unique_ptr<typename Prop::prop_t>& ptr) const
    {
        static_assert(is_any_of<Prop, Props...>::value, "NO SUCH PROPERTY MANAGED BY THIS CLASS");
        assert(propIsAllocated<Prop>(ptr));
        return *ptr;
    }

    template <typename Prop, typename std::enable_if<Prop::IS_MAPPED, int>::type = 0>
    void clearProp(std::unique_ptr<typename Prop::prop_t>& ptr, const typename Prop::value_t& def)
    {
        (void)def;
        static_assert(is_any_of<Prop, Props...>::value, "NO SUCH PROPERTY MANAGED BY THIS CLASS");
        assert(propIsAllocated<Prop>(ptr));
        ptr.reset();
    }

    template <typename Prop, typename std::enable_if<!Prop::IS_MAPPED, int>::type = 0>
    void clearProp(std::unique_ptr<typename Prop::prop_t>& ptr, const typename Prop::value_t& def)
    {
        static_assert(is_any_of<Prop, Props...>::value, "NO SUCH PROPERTY MANAGED BY THIS CLASS");
        assert(propIsAllocated<Prop>(ptr));
        mesh.set_persistent(*ptr, false);
        for (auto it = ptr->begin(); it != ptr->end(); it++)
            *it = def;
        ptr.reset();
    }

  private:
    static size_t nextID()
    {
        static size_t MAX_ID = 0;
        return MAX_ID++;
    }

    std::tuple<Props...> _props;
};

} // namespace mc3d

#endif
