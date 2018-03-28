// [header]
// A very basic raytracer example.
// [/header]
// [compile]
// c++ -o raytracer -O3 -Wall raytracer.cpp
// [/compile]
// [ignore]
// Copyright (C) 2012  www.scratchapixel.com
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// [/ignore]
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <limits>

# include "scene.hpp"
# include "solids.hpp"
using namespace geometry;

//[comment]
// This variable controls the maximum recursion depth
//[/comment]
const double pi = 3.14159265358979323846;
#define MAX_RAY_DEPTH 5

template<typename Real>
Real mix(const Real &a, const Real &b, const Real & mix)
{
    return b * mix + a * (1 - mix);
}


//[comment]
// This is the main trace function. It takes a ray as argument (defined by its origin
// and direction). We test if this ray intersects any of the geometry in the scene.
// If the ray intersects an object, we compute the intersection point, the normal
// at the intersection point, and shade this point using this information.
// Shading depends on the surface property (is it transparent, reflective, diffuse).
// The function returns a color for the ray. If the ray intersects an object that
// is the color of the object at the intersection point, otherwise it returns
// the background color.
//[/comment]
template<typename Real>
color::rgba trace( const rayon<Real>& ray, const scene<Real>& scene, const int& depth )
{
    //if (raydir.length() != 1) std::cerr << "Error " << raydir << std::endl;
    Real nearest = std::numeric_limits<double>::infinity();
    std::shared_ptr<geometry::surface<Real>> first_intersected_obj = nullptr;
    geometry::point<Real> intersect_point;
    geometry::vecteur<Real> intersect_normal;
    // find intersection of this ray with an object in the scene
    for ( const auto& obj : scene ) {
        auto res_intersect = obj->intersect(ray);
        if (std::get<0>(res_intersect)) {
            auto diff = std::get<1>(res_intersect) - ray.origin;
            Real dist2 = diff.square_norm();
            if (dist2 < nearest) {
                nearest = dist2;
                first_intersected_obj = obj;
                intersect_point = std::get<1>(res_intersect);
                intersect_normal = std::get<2>(res_intersect);
            }
        }
    }
    if ( first_intersected_obj == nullptr) return {0.,0.,0.,0.};
    // if there's no intersection return black or background color
    color::rgba surfaceColor; // color of the ray/surfaceof the object intersected by the ray
    // If the normal and the view direction are not opposite to each other
    // reverse the normal direction. That also means we are inside the sphere so set
    // the inside bool to true. Finally reverse the sign of IdotN which we want
    // positive.
    Real bias = std::numeric_limits<double>::epsilon(); // add some bias to the point from which we will be tracing
    bool inside = false;
    if ((ray.direction|intersect_normal) > 0)
    {
        intersect_normal = -intersect_normal;
        inside = true;
    }
    Real transparency = 1. - first_intersected_obj->surface_color[color::rgba::alpha];
    if ((transparency > 0 || first_intersected_obj->reflection > 0) && depth < MAX_RAY_DEPTH) {
        Real facingratio = -(ray.direction|intersect_normal);
        // change the mix value to tweak the effect
        Real fresneleffect = mix(std::pow(1 - facingratio, 3), 1., 0.1);
        // compute reflection direction (not need to normalize because all vectors
        // are already normalized)
        vecteur<Real> refldir = ray.direction - 2. * ( ray.direction | intersect_normal) * intersect_normal;
        refldir.normalize();
        color::rgba reflection = trace({intersect_point + bias * intersect_normal, refldir}, scene, depth+1 );
        color::rgba refraction = {0.,0.,0.,0.};
        // if the sphere is also transparent compute refraction ray (transmission)
        if (transparency > 0.) {
            Real ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface?
            Real cosi = -(intersect_normal|ray.direction);
            Real k = 1 - eta * eta * (1 - cosi * cosi);
            vecteur<Real> refrdir = eta * ray.direction + (eta *  cosi - std::sqrt(k)) * intersect_normal;
            refrdir.normalize();
            refraction = trace({intersect_point - bias * intersect_normal, refrdir}, scene, depth + 1);
        }
        // the result is a mix of reflection and refraction (if the sphere is transparent)
        surfaceColor = ( fresneleffect * reflection + (1 - fresneleffect) * transparency * refraction ) * first_intersected_obj->surface_color;
    }
    else {
        point<Real> bias_intersect = intersect_point + bias * intersect_normal;
        // it's a diffuse object, no need to raytrace any further
        // Search lighting sources :
        for ( unsigned iobj = 0; iobj < scene.size(); ++iobj )
        {
            if ( scene[iobj].emission_color[color::rgba::alpha] > 0) {
                // this is a light
                const surface<Real>& light = scene[iobj];
                Real transmission = 1;
                vecteur<Real> lightDirection = light.origin - intersect_point;
                lightDirection.normalize();

                for (unsigned jobj = 0; jobj < scene.size(); ++jobj) {
                    bool intersect = false;
                    if (iobj != jobj) {
                        std::tie(intersect,std::ignore,std::ignore) = scene[jobj].intersect({bias_intersect, lightDirection});
                        if (intersect) {
                            transmission = 0;
                            break;
                        }
                    }
                }
                surfaceColor += transmission * first_intersected_obj->surface_color  *
                                std::max(Real(0), (intersect_normal|lightDirection)) * light.emission_color;
            }
        }
    }
    return surfaceColor + first_intersected_obj->emission_color + scene.ambient_color;
}

//[comment]
// Main rendering function. We compute a camera ray for each pixel of the image
// trace it and return a color. If the ray hits a sphere, we return the color of the
// sphere at the intersection point, else we return the background color.
//[/comment]
template<typename Real>
void render(const scene<Real> & world )
{
    unsigned width = 1280, height = 1024;
    std::vector<color::rgba> image(width * height);
    color::rgba* pixel = image.data();
    Real invWidth = 1 / Real(width), invHeight = 1 / Real(height);
    Real fov = 30, aspectratio = width / Real(height);
    float angle = std::tan(pi * 0.5 * fov / 180.);
    // Trace rays
    for (unsigned y = 0; y < height; ++y) {
        std::cout << "Rendering line " << y << "\r";
        std::flush(std::cout);
        for (unsigned x = 0; x < width; ++x, ++pixel) {
            Real xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
            Real yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            vecteur<Real> raydir(xx, yy, -1);
            raydir.normalize();
            *pixel = trace({{0.,0.,0.},raydir}, world, 0);
        }
    }
    std::cout << std::endl;
    // Save result to a PPM image (keep these flags if you compile under Windows)
    std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (unsigned i = 0; i < width * height; ++i) {
        ofs << (unsigned char)(std::min(1., image[i][color::rgba::red]  ) * 255) <<
               (unsigned char)(std::min(1., image[i][color::rgba::green]) * 255) <<
               (unsigned char)(std::min(1., image[i][color::rgba::blue] ) * 255);
    }
    ofs.close();
}

//[comment]
// In the main function, we will create the scene which is composed of 5 spheres
// and 1 light (which is also a sphere). Then, once the scene description is complete
// we render that scene, by calling the render() function.
//[/comment]
int main(int argc, char **argv)
{
    using real = double;
    srand48(13);
    scene<real> world;
    world.ambient_color = color::rgba(0.1,0.1,0.1,1.);
    // position, radius, surface color, reflectivity, transparency, emission color
    //world.push(plane<real>{ {0.,-4.,0.}, {0.0,1.0,0.0}, {0.30, 0.20, 0.20, 1.}, 0.0, true});
    world.push(plane<real>{ {0.,0.,-50.}, {0.0,0.0,1.0}, {0.00, 0.70, 0.90, 1.}, 0.0, true});
    world.push( sphere<real>{{0.0, -5004, -20}, 5000, {0.30, 0.20, 0.20, 1.}, 0.0, true});
    
    const int nbSpheres = 100;
    for ( int i = 0; i < nbSpheres; ++i ) {
      real x,y,z,rd,r,b,g,t;
      x = (rand()/(1.*RAND_MAX))*20.-10.;
      y = (rand()/(1.*RAND_MAX))*2.-1.;
      z = (rand()/(1.*RAND_MAX))*10.-25.;
      rd = (rand()/(1.*RAND_MAX))*0.9+0.1;
      r  = (rand()/(4.*RAND_MAX))+0.75;
      g  = (rand()/(4.*RAND_MAX))+0.75;
      b  = (rand()/(4.*RAND_MAX))+0.75;
      t  = (rand()/(4.*RAND_MAX))+0.75;
      world.push( sphere<real>({x, y, z}, rd, {r, g, b, t}, 1.0, true));
    }
    point<real> p0{-0.5,-1,-10}, p1{0.5,-1,-10}, p2{0,-1,-11}, p3{0,   0, -11};
    world.push( triangle<real>(p1, p0, p2, {1., 1., 1., 1.0}, 1.0, true));
    world.push( triangle<real>(p0, p2, p3, {1., 1., 1., 1.0}, 1.0, true));
    world.push( triangle<real>(p1, p3, p2, {1., 1., 1., 1.0}, 1.0, true));
    world.push( triangle<real>(p1, p0, p3, {1., 1., 1., 1.0}, 1.0, true));
    // light
    world.push( sphere<real>({0.0, 20, -3.}, 3, {1., 1., 1., 1.}, 0.0, false));
    render(world);
    
    return 0;
}
