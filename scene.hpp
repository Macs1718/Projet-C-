#pragma once
#include "surface.hpp"
#include <memory>
#include <vector>
using namespace geometry;

template <typename Real>
class scene {
  public:
    using value_t        = std::shared_ptr<geometry::surface<Real>>;
    using container      = std::vector<value_t>;
    using iterator       = typename container::iterator;
    using const_iterator = typename container::const_iterator;

    scene() : m_objects(), ambient_color{0.25,0.25,0.25,1.} { m_objects.reserve(1024); }

    unsigned size() const { return m_objects.size(); }

    geometry::surface<Real>& operator [] ( unsigned i )
    {
    	return *m_objects[i];
    }

    const geometry::surface<Real>& operator [] ( unsigned i ) const
    {
    	return *m_objects[i];
    }

    geometry::surface<Real>& get( unsigned i )
    {
    	return *m_objects.get(i);
    }

    const geometry::surface<Real>& get( unsigned i ) const
    {
    	return *m_objects.get(i);
    }

    void push( const surface<Real>& obj )
    {
    	m_objects.emplace_back(obj.clone());
    }

    iterator begin() { return m_objects.begin(); }
    const_iterator cbegin() const { return m_objects.cbegin(); }
    const_iterator begin() const { return cbegin(); }

    iterator end() { return m_objects.end(); }
    const_iterator cend() const { return m_objects.cend(); }
    const_iterator end() const { return cend(); }

    color::rgba ambient_color;

  private:
    container m_objects;
};
