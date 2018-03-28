  std::tuple<bool, point<Real>, vecteur<Real>> intersect(const rayon<Real>& ray) const override {
            if ( (ray.direction|normal) == 0 ) return std::make_tuple(false, point<Real>{0, 0, 0}, vecteur<Real>{0, 0, 0});
            vecteur<Real> diff = origin - ray.origin;
            Real T =  (diff|normal)/(ray.direction|normal);
            if (T <= 0. ) return std::make_tuple(false, point<Real>{0, 0, 0}, vecteur<Real>{0, 0, 0});
            Real e00 = (E0|E0), e01 = (E0|E1), e11 = (E1|E1);
            vecteur<Real> Q = T * ray.direction - diff;
            Real q0 = (E0|Q), q1 = (E1|Q);
            Real delta = e00*e11 - e01*e01;
            Real sigma0 = e11*q0 - e01*q1;
            Real sigma1 = e00*q1 - e01*q0;
            if ( (sigma0 >= 0) and (sigma1 >= 0) and (sigma0+sigma1 <= delta) ) 
                {
                    return std::make_tuple(true, Q+origin, normal);
                }
            else return std::make_tuple(false, point<Real>{0, 0, 0}, vecteur<Real>{0, 0, 0});
        }
