        std::tuple<bool, point<Real>, vecteur<Real>> intersect(const rayon<Real>& ray) const override {
            if ( (ray.direction|normal) == 0 ) return std::make_tuple(false, point<Real>{0, 0, 0}, vecteur<Real>{0, 0, 0});
            vecteur<Real> diff = origin - ray.origin;
            Real alpha =  (diff|normal)/(ray.direction|normal);
            if (alpha <= 0. ) return std::make_tuple(false, point<Real>{0, 0, 0}, vecteur<Real>{0, 0, 0});
            vecteur<Real> phit = ray.origin + alpha * ray.direction;
            return std::make_tuple(true, phit, normal);
        }
