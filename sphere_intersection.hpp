 std::tuple<bool, point<Real>, vecteur<Real>> intersect(const rayon<Real>& ray) const override {
            vecteur<Real> l   = origin - ray.origin;
            Real          tca = (l | ray.direction);
            if (tca < 0) return std::make_tuple(false, point<Real>{0, 0, 0}, vecteur<Real>{0, 0, 0});
            Real d2 = (l | l) - tca * tca;
            if (d2 > radius * radius) return std::make_tuple(false, point<Real>{0, 0, 0}, vecteur<Real>{0, 0, 0});
            Real thc = std::sqrt(radius * radius - d2);
            Real t0  = tca - thc;
            Real t1  = tca + thc;
            if (t0 < 0) t0 = t1;
            vecteur<Real> phit = ray.origin + t0 * ray.direction; // point of intersection
            vecteur<Real> nhit = phit - origin;    // normal at the intersection point
            nhit.normalize();                                     // normalize normal direction
            return std::make_tuple(true, phit, nhit);
        }
