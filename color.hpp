#pragma once
# include <array>
# include <algorithm>
# include <limits>
# include <iostream>

namespace color
{
	class rgba
	{
	public:
		using container = std::array<double,4>;
		enum components {
		    red = 0,
		    green = 1,
		    blue = 2,
		    alpha = 3
		};

		rgba() : m_value({0.,0.,0.,1.})
		{}
	    rgba( double vr, double vg, double vb, double va = 1. ) :
	        m_value({vr,vg,vb,va})
	    {}
	    rgba( const rgba& ) = default;
	    rgba( rgba&& ) = default;
	    ~rgba() = default;

	    rgba& operator = ( const rgba& ) = default;
	    rgba& operator = ( rgba&& ) = default;

	    double operator[] ( components c ) const
	    {
	    	return m_value[c];
	    }
	    double& operator[] ( components c )
	    {
	    	return m_value[c];
	    }

	    const container get_values() const
	    {
	    	return m_value;
	    }

	    container get_values()
	    {
	    	return m_value;
	    }

	    rgba& operator += ( const rgba& col ) {
	    	m_value[red] += col[red]; m_value[blue] += col[blue]; m_value[green] += col[green];
	    	return *this;
	    }

	    rgba operator + ( const rgba& col ) const {
	    	return { m_value[red]+col[red], m_value[green] + col[green], m_value[blue] + col[blue], m_value[alpha] };
	    }

	    rgba operator * ( const double& c ) const {
	    	return { c * m_value[red], c * m_value[green], c * m_value[blue], m_value[alpha]};
	    }

	    rgba operator * ( const rgba& col ) const {
	    	return { m_value[red] * col[red], m_value[green] * col[green], m_value[blue] * col[blue], m_value[alpha] };
	    }

	private:
		container m_value;
	};

	inline rgba operator * ( const double& a, const rgba& col )
	{
		return col * a;
	}

	inline std::ostream& operator << ( std::ostream& out, const rgba& col )
	{
		out << "<" << col[rgba::red] << ", " << col[rgba::green] << ", " << col[rgba::blue] << ", " << col[rgba::alpha] << ">";
		return out;
	}

	class hsva
	{
	public:
		using container = std::array<double,4>;
		enum component {
			hue=0, saturation, value, alpha
		};
		hsva() : m_value({0.,0.,0.,1.})
		{}
		hsva( double vh, double vs, double vv, double va ) : m_value({vh,vs,vv,va})
		{}
		hsva( const rgba& col ) : m_value()
		{
			m_value[alpha] = col[rgba::alpha];
			std::pair<double,double> vminmax;
			vminmax = std::minmax({col[rgba::red], col[rgba::green], col[rgba::blue]});

			m_value[value] = vminmax.second; // v
			double delta = vminmax.second - vminmax.first;
    		if (delta < 0.00001)
    		{
    			m_value[saturation] = 0;
    			m_value[hue]        = 0;  // undefined, maybe nan?
    			return;
    		}
    		if( vminmax.second > 0.0 ) { // NOTE: if Max is == 0, this divide would cause a crash
        		m_value[saturation] = (delta / vminmax.second);                  // s
			} else {
			    // if max is 0, then r = g = b = 0 =>  s = 0, h is undefined
			    m_value[saturation] = 0.0;
			    m_value[hue]        = std::numeric_limits<double>::quiet_NaN();
			    return;
    		}
    		if (col[rgba::red] >= vminmax.second) // > is bogus, just keeps compilor happy
    			m_value[hue] = (col[rgba::green] - col[rgba::blue])/delta; 
		    else if( col[rgba::green] >= vminmax.second )
		        m_value[hue] = 2.0 + (col[rgba::blue] - col[rgba::red]) / delta;  // between cyan & yellow
    		else
        		m_value[hue] = 4.0 + (col[rgba::red] - col[rgba::green]) / delta;  // between magenta & cyan

        	m_value[hue] *= 60.0; // degrees

        	if (m_value[hue] < 0.0)
        		m_value[hue] += 360.0;
		}
		hsva( const hsva& ) = default;
		hsva( hsva&& )      = default;
		~hsva() = default;

		hsva& operator = ( const hsva& ) = default;
		hsva& operator = ( hsva&& )      = default;

		double operator[] ( component c ) const
		{
			return m_value[c];
		}
		double& operator [] ( component c )
		{
			return m_value[c];
		}

		operator rgba() const
		{
			double hh, p, q, t, ff;
    		long        i;
    		rgba out;

    		out[rgba::alpha] = m_value[alpha];
    		if(m_value[saturation]<= 0.0) {       // < is bogus, just shuts up warnings
    			out[rgba::red] = m_value[value];
    			out[rgba::green] = m_value[value];
    			out[rgba::blue] = m_value[value];
		        return out;
    		}
		    hh = m_value[hue];
		    if(hh >= 360.0) hh = 0.0;
    		hh /= 60.0;
    		i = (long)hh;
    		ff = hh - i;
    		p = m_value[value] * (1.0 - m_value[saturation]);
		    q = m_value[value] * (1.0 - (m_value[saturation] * ff));
 			t = m_value[value] * (1.0 - (m_value[saturation] * (1.0 - ff)));

		    switch(i) {
    		case 0:
    			out[rgba::red]   = m_value[value];
    			out[rgba::green] = t;
    			out[rgba::blue]  = p;
        		break;
		    case 1:
		    	out[rgba::red] = q;
		    	out[rgba::green] = m_value[value];
		    	out[rgba::blue] = p;
		        break;
			case 2:
				out[rgba::red]   = p;
				out[rgba::green] = m_value[value];
				out[rgba::blue]  = t;
        		break;

	 	   case 3:
        		out[rgba::red]   = p;
        		out[rgba::green] = q;
        		out[rgba::blue]  = m_value[value];
        		break;
    		case 4:
    			out[rgba::red]   = t;
    			out[rgba::green] = p;
    			out[rgba::blue]  = m_value[value];
        		break;
    		case 5:
    		default:
    			out[rgba::red]   = m_value[value];
    			out[rgba::green] = p;
    			out[rgba::blue]  = q;
        		break;
    		}
    		return out;			
		}
	private:
		container m_value;
	};

}