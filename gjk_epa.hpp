/*
 * gjk_epa.hpp
 *
 *  Created on: 19 февр. 2018 г.
 *      Author: werqa
 */

#ifndef GJK_EPA_HPP_
#define GJK_EPA_HPP_

#include <math.h>
#include <limits>
#include <list>

#define EPSILON 0.0001

class Vector3
{

public:

	//--------------------[ Attributes ]--------------------//

	float x;
	float y;
	float z;


public:

        //----------------[ constructors ]--------------------------

    Vector3()
    : x(0), y(0), z(0)
    {
    }

    Vector3(float nx, float ny, float nz)
     : x(nx), y(ny), z(nz)
    {
    }



    //---------------[ vector arithmetic operator ]--------------

    Vector3 operator+(const Vector3& rhs) const
    {
        return Vector3(x + rhs.x, y + rhs.y, z + rhs.z);
    }


    Vector3 operator-(const Vector3& rhs) const
    {
        return Vector3(x - rhs.x, y - rhs.y, z - rhs.z);
    }

    Vector3 operator*(const Vector3& rhs) const
    {
        return Vector3(x * rhs.x, y * rhs.y, z * rhs.z);
    }

    Vector3 operator/(const Vector3& rhs) const
    {
        return Vector3(x / rhs.x, y / rhs.y, z / rhs.z);
    }

    //--------------[ scalar vector operator ]--------------------

    Vector3 operator*(float rhs) const
    {
        return Vector3(x * rhs, y * rhs, z * rhs);
    }

    Vector3 operator/(float rhs) const
    {
        return Vector3(x / rhs, y / rhs, z / rhs);
    }


    //-------------[ unary operations ]--------------------------

    Vector3 operator-() const
    {
        return Vector3(-x, -y, -z);
    }



    //-------------[ function operations ]---------------------------

    float dot(const Vector3& rhs) const
    {
        return x * rhs.x + y * rhs.y + z * rhs.z;
    }

    Vector3 cross(const Vector3& rhs) const
    {
        return Vector3(y * rhs.z - rhs.y * z,
                       z * rhs.x - rhs.z * x,
                       x * rhs.y - rhs.x * y);
    }

    float length() const
    {
        return sqrt(x * x + y * y + z * z);
    }

    float lengthSquare() const
    {
        return x * x + y * y + z * z;
    }

    void normalize()
    {
        float s = length();
        x /= s;
        y /= s;
        z /= s;
    }

    Vector3 getUnit() const
    {
        float lengthVector = length();
        if (lengthVector < EPSILON)
        {
            return *this;
        }
        // Compute and return the unit vector
        float lengthInv = 1.0 / lengthVector;
        return Vector3( x * lengthInv,
                        y * lengthInv,
                        z * lengthInv);
    }


    void rotateAroundAxis(const Vector3 &axis, const float angle)
    {
    	Vector3 self = *this;
    	Vector3 x = self - axis * (axis * self);
    	Vector3 y = x.cross(axis);
    	Vector3 delta = x * (cos(angle)) + y * (sin(angle)) - x;
    	self = self + delta;
    	*this = self;
    }
};


float dot(const Vector3& a , const Vector3& b)
{
	return a.dot(b);
}

Vector3 cross(const Vector3& a , const Vector3& b)
{
	return a.cross(b);
}



//--------------------------------  gjk-epa ----------------------------------//


// FLOATING POINT VALIDITY
template<typename T> bool is_infinite(const T &value)
{
    T max_value = (std::numeric_limits<T>::max)();
    T min_value = -max_value;
    return !(min_value <= value && value <= max_value);
}

template<typename T> bool is_nan(const T &value)
{
    // True if NAN
    return value != value;
}

template<typename T> bool is_valid(const T &value)
{
    return !is_infinite(value) && !is_nan(value);
}


class ContactBasicData
{
  public:

        Vector3 point; // point is on the surface of A (inside of B if penetration is positive)
        Vector3 normal; // points from B to A
        float penetration;

        // the local support points of the minkowski difference triangle that was used to generate this contact originally
        // these points can be used to update the contact data even after the two objects have been arbitrarily changed
        Vector3 triangleSupports_world[2][3];
};




class GJKEPAGenerator
{
  public:

	GJKEPAGenerator();



	struct SupportPoint
	{
		// WORLD SPACE
		Vector3 v;
		Vector3 Supports[2];

		int uniqueId;

		// id is only used in the expanding polytope algorithm while generating edges
		inline SupportPoint(int uniqueId = 0) : uniqueId(uniqueId) {}
	};



	// Minkowski Differential
	template <typename ConvexTemplate>
	struct	MinkowskiDiff
	{
		const ConvexTemplate* mConvexAPtr;
		const ConvexTemplate* mConvexBPtr;

		bool				  mEnableMargin;


		MinkowskiDiff(const ConvexTemplate* a ,
				      const ConvexTemplate* b )
		:mConvexAPtr(a),
		 mConvexBPtr(b)
		{
		}

		void EnableMargin(bool enable)
		{
			mEnableMargin = enable;
		}

		SupportPoint generateSupport( Vector3 directionWorld , int* uniqueIdCounter) const
		{
			directionWorld.normalize();

			SupportPoint ret(uniqueIdCounter ? (*uniqueIdCounter)++ : 0);
			ret.Supports[0] = mConvexAPtr->support_world( directionWorld );
			ret.Supports[1] = mConvexBPtr->support_world(-directionWorld );

			if(mEnableMargin)
			{
				// .....!?
			}

			ret.v = (ret.Supports[0]) -
					(ret.Supports[1]);

			return ret;
		}

	};



	struct Edge
	{
		SupportPoint points[2];

		inline Edge(const SupportPoint &a,const SupportPoint &b)
		{
			points[0] = a;
			points[1] = b;
		}
	};




	struct Triangle
	{
		SupportPoint points[3];
		Vector3 n;

		inline Triangle(const SupportPoint &a,const SupportPoint &b,const SupportPoint &c)
		{
			points[0] = a;
			points[1] = b;
			points[2] = c;
			n = cross((b.v-a.v),(c.v-a.v)).getUnit();
		}
	};





	struct Simplex
	{

	public:

		SupportPoint _simplex[4];

		int num;
		SupportPoint &a;
		SupportPoint &b;
		SupportPoint &c;
		SupportPoint &d;

		inline Simplex() :
         a(_simplex[0]),
		 b(_simplex[1]),
		 c(_simplex[2]),
		 d(_simplex[3])
		{
			clear();
		}

		inline void clear() { num = 0; }

		inline void set(SupportPoint a,SupportPoint b,SupportPoint c,SupportPoint d)
		{
			num = 4;
			this->a = a;
			this->b = b;
			this->c = c;
			this->d = d;
		}

		inline void set(SupportPoint a,SupportPoint b,SupportPoint c)
		{
			num = 3;
			this->a = a;
			this->b = b;
			this->c = c;
		}

		inline void set(SupportPoint a,SupportPoint b)
		{
			num = 2;
			this->a = a;
			this->b = b;
		}

		inline void set(SupportPoint a)
		{
			num = 1;
			this->a = a;
		}

		inline void push(SupportPoint p)
		{
			num = (std::min)(num+1,4);
			for(int i = num-1; i > 0; i--) _simplex[i] = _simplex[i-1];
			_simplex[0] = p;
		}

	};

	template <typename ConvexTemplate>
	static bool extrapolateContactInformation(ContactBasicData* contactData , const MinkowskiDiff<ConvexTemplate>&  MinkowskiShape)
	{

		SupportPoint supports[3];

		for(uint8_t i = 0; i < 3; i++)
		{
			supports[i].Supports[0] = contactData->triangleSupports_world[0][i];
			supports[i].Supports[1] = contactData->triangleSupports_world[1][i];
			supports[i].v = (supports[i].Supports[0]) -
					        (supports[i].Supports[1]);
		}

		Triangle triangle(supports[0],supports[1],supports[2]);

		// actually extrapolate contact information
		return extrapolateContactInformation(&triangle,contactData , MinkowskiShape);
	}


	template <typename ConvexTemplate>
	static bool extrapolateContactInformation(const Triangle* triangle , ContactBasicData* contactData , const MinkowskiDiff<ConvexTemplate>&  MinkowskiShape)
	{

		// project the origin onto the triangle
		const float distanceFromOrigin = dot(triangle->n,triangle->points[0].v);

		// calculate the barycentric coordinates of the closest triangle with respect to
		// the projection of the origin onto the triangle
		float bary_u,bary_v,bary_w;
		barycentric(triangle->n * distanceFromOrigin,
				triangle->points[0].v,
				triangle->points[1].v,
				triangle->points[2].v,
				&bary_u,
				&bary_v,
				&bary_w);

		// barycentric can fail and generate invalid coordinates, if this happens return false
		if(!is_valid(bary_u) || !is_valid(bary_v) || !is_valid(bary_w)) return false;

		// if any of the barycentric coefficients have a magnitude greater than 1, then the origin is not within the triangular prism described by 'triangle'
		// thus, there is no collision here, return false
		if(fabs(bary_u)>1.0f || fabs(bary_v)>1.0f || fabs(bary_w)>1.0f) return false;

		// collision point on object a in world space
		const Vector3 point(
				((triangle->points[0].Supports[0]) * bary_u)  +
				((triangle->points[1].Supports[0]) * bary_v)  +
				((triangle->points[2].Supports[0]) * bary_w) );

		// collision normal
		const Vector3 normal = -triangle->n;

		// penetration depth
		const float penetrationDepth = distanceFromOrigin;


		contactData->point = point;
		contactData->normal = normal;
		contactData->penetration = penetrationDepth;

		// copy data necessary for contact persistence into the contactData structure
		for(uint8_t i = 0; i < 2; i++)
		{
			for(uint8_t j = 0; j < 3; j++)
			{
				contactData->triangleSupports_world[i][j] = triangle->points[j].Supports[i];
			}
		}

		return true;
	}

	// FUNCTIONS
	static void barycentric(const Vector3 &p,const Vector3 &a,const Vector3 &b,const Vector3 &c,float *u,float *v,float *w)
	{
		// code from Crister Erickson's Real-Time Collision Detection
		Vector3 v0 = b - a,v1 = c - a,v2 = p - a;
		float d00 = dot(v0,v0);
		float d01 = dot(v0,v1);
		float d11 = dot(v1,v1);
		float d20 = dot(v2,v0);
		float d21 = dot(v2,v1);
		float denom = d00 * d11 - d01 * d01;
		*v = (d11 * d20 - d01 * d21) / denom;
		*w = (d00 * d21 - d01 * d20) / denom;
		*u = 1.0f - *v - *w;
	}


	static bool gjk_simtest(const Vector3& v , const Vector3& ao )
	{
		return (v).dot(ao) > 0;
	}


	template <typename ConvexTemplate>
	static bool GJKIntersectionTest(Simplex* inout_simplex , int* uniqueIdCounter , const MinkowskiDiff<ConvexTemplate>&  MinkowskiShape , ContactBasicData* in_contactData )
	{
		Simplex& sim = *inout_simplex;

		/*
		 ********** Note **********
	    There is an issue with with ellipsoids that have very non-uniform scaling.
	    For some reason we fail to pass the origin on the first or second early out test.
		 **************************/
		const unsigned _EXIT_ITERATION_LIMIT = 30;
		unsigned _EXIT_ITERATION_NUM = 0;

		// build the initial simplex, a single support point in an arbitrary position
		sim.clear();
		Vector3 dir =Vector3(1,1,1);


		SupportPoint s = MinkowskiShape.generateSupport(dir,uniqueIdCounter);

		if(fabs(dot(dir,s.v))>=s.v.length()*0.8f)
		{
			// the chosen initial direction is invalid, will produce (0,0,0) for a subsequent direction later, try a different arbitrary direction
			dir = Vector3(1,0,0);
			s = MinkowskiShape.generateSupport(dir,uniqueIdCounter);
		}

		sim.push(s);
		dir = -s.v;


		// iteratively attempt to build a simplex that encloses the origin
		while(true)
		{
			if(_EXIT_ITERATION_NUM++ >= _EXIT_ITERATION_LIMIT) return false;

			// error, for some reason the direction vector is broken
			if(dir.lengthSquare()<=0.0001f) return false;

			// get the next point in the direction of the origin
			SupportPoint a = MinkowskiShape.generateSupport(dir,uniqueIdCounter);


			// early out: if a.v dot d is less than zero that means that the new point did not go past the origin
			// and thus there is no intersection
			const float proj = dot(a.v,dir);

			if(proj < 0) return false;
			sim.push(a);

			const Vector3 ao = -sim.a.v;

			bool _triangleSimplex_needToTestAboveBelow = true; // used to skip a calculation in the tetrahedron case


			// simplex tests
			if(sim.num == 2)
			{
				// simplex is a line, being here means that the early out was passed, and thus
				// the origin must be between point a and b
				// search direction is perpendicular to ab and coplanar with ao

				const Vector3 ab = (sim.b.v-sim.a.v);
				dir = cross(ab,cross(ao,ab));
				continue;
			}
			else if(sim.num == 3)
			{
				jmp_triangleSimplex: // used used in the sim.num == 4 case

				// simplex is a triangle
				const Vector3 ab = (sim.b.v-sim.a.v);
				const Vector3 ac = (sim.c.v-sim.a.v);
				const Vector3 abc = cross(ab,ac);

				if(gjk_simtest(cross(ab,abc),ao))
				{
					// origin is in front of segment ab
					// we already know that the origin is in front of segment bc and behind the projection of segment bc onto point a
					// this is because point a is the support in the direction of bc's normal

					// [a b] ab % ao % ab
					sim.set(sim.a,sim.b);
					dir = cross(ab,cross(ao,ab));
				}
				else
				{
					if(gjk_simtest(cross(abc,ac),ao))
					{
						// origin is in front of segment ac
						// just as above, we know that the origin in in front of segment bc and behind ...

						// [a c] ac % ao % ac
						sim.set(sim.a,sim.c);
						dir = cross(ac,cross(ao,ac));
					}
					else
					{
						// if we have been redirected to the triangle case by the tetrahedron case this test is unnecessary
						// because we know the origin is above the triangle already
						if(!_triangleSimplex_needToTestAboveBelow || gjk_simtest(abc,ao))
						{
							// origin is above the triangle

							// [a b c] abc
							//sim.set(sim.a,sim.b,sim.a); it is already equal to this
							dir = abc;
						}
						else
						{
							// origin is below the triangle

							// [a c b] -abc
							sim.set(sim.a,sim.c,sim.b);
							dir = -abc;
						}
					}
				}
				continue;
			}
			else
			{        // == 4
				// the simplex is a tetrahedron, must check if it is outside any of the side triangles, (abc, acd, adb)
				// if it is then set the simplex equal to that triangle and jump back up to the triangle case, otherwise we know
				// there is an intersection and exit

				// check the triangles (abc,acd,adb), scoped as the temporary variables used here
				// will no longer be valid afterward
				{
					const Vector3 ab = (sim.b.v-sim.a.v);
					const Vector3 ac = (sim.c.v-sim.a.v);

					if(gjk_simtest(cross(ab,ac),ao))
					{
						// origin is in front of triangle abc, simplex is already what it needs to be
						// [a b c] goto jmp_triangleSimplex
						sim.num = 3;
						_triangleSimplex_needToTestAboveBelow = false;
						goto jmp_triangleSimplex;
					}

					const Vector3 ad = (sim.d.v-sim.a.v);

					if(gjk_simtest(cross(ac,ad),ao))
					{
						// origin is in front of triangle acd, simplex is set to this triangle
						// [a c d] goto jmp_triangleSimplex
						sim.set(sim.a,sim.c,sim.d);
						_triangleSimplex_needToTestAboveBelow = false;
						goto jmp_triangleSimplex;
					}

					if(gjk_simtest(cross(ad,ab),ao))
					{
						// origin is in front of triangle adb, simplex is set to this triangle
						// [a d b] goto jmp_triangleSimplex
						sim.set(sim.a,sim.d,sim.b);
						_triangleSimplex_needToTestAboveBelow = false;
						goto jmp_triangleSimplex;
					}

					// intersection confirmed, break from the loop
					break;
				}
			}

		}


		return true;
	}




	template <typename ConvexTemplate>
	static bool EPAManifoldExtraction(const Simplex* in_initialSimplex , int* uniqueIdCounter , MinkowskiDiff<ConvexTemplate>  MinkowskiShape , ContactBasicData* inout_contactData)
	{

	    // the algorithm will iteratively try to grow the polytope to find the face on the minkowski difference that is closest to the origin
	    // it will continue to do this until if fails to push the polytope more than _EXIT_GROWTH_THRESHOLD
	    const float _EXIT_GROWTH_THRESHOLD = 0.002f;
	    const unsigned _EXIT_ITERATION_LIMIT = 30;
	    unsigned _EXIT_ITERATION_CUR = 0;

	    std::list<Triangle> lst_triangles;
	    std::list<Edge>     lst_edges;

	    // process the specified edge, if another edge with the same points in the
	    // opposite order exists then it is removed and the new point is also not added
	    // this ensures only the outermost ring edges of a cluster of triangles remain
	    // in the list
	    auto lam_addEdge = [&](const SupportPoint &a,const SupportPoint &b)->void
	    {
	        for(auto it = lst_edges.begin(); it != lst_edges.end(); it++)
	        {
	            if(it->points[0].uniqueId==b.uniqueId && it->points[1].uniqueId==a.uniqueId)
	            {
	                //opposite edge found, remove it and do not add new one
	                lst_edges.erase(it);
	                return;
	            }
	        }
	        lst_edges.emplace_back(a,b);
	    };

	    // add the GJK simplex triangles to the list
	    lst_triangles.emplace_back(in_initialSimplex->a,in_initialSimplex->b,in_initialSimplex->c);
	    lst_triangles.emplace_back(in_initialSimplex->a,in_initialSimplex->c,in_initialSimplex->d);
	    lst_triangles.emplace_back(in_initialSimplex->a,in_initialSimplex->d,in_initialSimplex->b);
	    lst_triangles.emplace_back(in_initialSimplex->b,in_initialSimplex->d,in_initialSimplex->c);

	    while(true)
	    {
	        if(_EXIT_ITERATION_CUR++ >= _EXIT_ITERATION_LIMIT) return false;

	        // find closest triangle to origin
	        std::list<Triangle>::iterator closestTriangle_it;
	        float closestTriangle_dst = 3.402823466e+38F;
	        for(auto it = lst_triangles.begin(); it != lst_triangles.end(); it++)
	        {
	            const float dst = dot(it->n , it->points[0].v);
	            if(dst < closestTriangle_dst)
	            {
	                closestTriangle_dst = dst;
	                closestTriangle_it = it;

	            }
	        }

	        // get the next support point in front of the current triangle, away from the origin
	        const SupportPoint entry_new_support = MinkowskiShape.generateSupport(closestTriangle_it->n,uniqueIdCounter);




	        // checks how much further this new point will take us from the origin
	        // if it is not far enough then we assume we have found the closest triangle
	        // on the hull from the origin
	        const float newDst = dot(closestTriangle_it->n,entry_new_support.v);
	        const float growth = newDst - closestTriangle_dst;
	        //if(thisItGrowth < 0) break;
	        if((growth < _EXIT_GROWTH_THRESHOLD)) {
	            //printf("EPA: exit growth: %f\n",growth);

	            //////////////////////////////////////////////////////////////////////////
	            //////////////////////////////////////////////////////////////////////////
	            // GENERATE CONTACT INFO AND RETURN
	            return extrapolateContactInformation(&*closestTriangle_it,
	                                                 inout_contactData,MinkowskiShape);


	            break;
	        }

	        for(auto it = lst_triangles.begin(); it != lst_triangles.end();)
	        {
	            // can this face be 'seen' by entry_new_support?
	            if(dot(it->n , (entry_new_support.v - it->points[0].v)) > 0)
	            {
	                lam_addEdge(it->points[0],it->points[1]);
	                lam_addEdge(it->points[1],it->points[2]);
	                lam_addEdge(it->points[2],it->points[0]);
	                it = lst_triangles.erase(it);
	                continue;
	            }
	            it++;
	        }

	        // create new triangles from the edges in the edge list
	        for(auto it = lst_edges.begin(); it != lst_edges.end(); it++)
	        {
	            lst_triangles.emplace_back(entry_new_support,it->points[0],it->points[1]);
	        }

	        lst_edges.clear();
	    }

	    return false;
	}




	template <typename ConvexTemplate>
	static bool CollisionSolveGJK_EPA(const ConvexTemplate* a,
			                          const ConvexTemplate* b,
									  ContactBasicData* input_contactData )
	{
		MinkowskiDiff<ConvexTemplate> MinkowskiShape(a,b);

		Simplex testSimplex;
		int supportPointNextUniqueId = 0;

	    bool isValid = (GJKIntersectionTest(&testSimplex,&supportPointNextUniqueId,MinkowskiShape , input_contactData) &&
	                    EPAManifoldExtraction(&testSimplex,&supportPointNextUniqueId,MinkowskiShape,input_contactData) );

	    return isValid;
	}


};


#endif /* GJK_EPA_HPP_ */
