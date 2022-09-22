#ifndef VECH
#define VECH

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include "rtweekend.h"

class vec3 {
public:
	vec3() {}
	vec3(float e0, float e1, float e2) {e[0]=e0;e[1]=e1;e[2]=e2;}
	
	inline float x() const {return e[0];}
	inline float y() const {return e[1];}
	inline float z() const {return e[2];}
	inline float r() const {return e[0];}
	inline float g() const {return e[1];}
	inline float b() const {return e[2];}
	
	inline const vec3& operator+() const {return *this;}
	inline vec3 operator-() const {return vec3(-e[0],-e[1],-e[2]);}
	inline float operator[](int i) const {return e[i];}
	inline float& operator[](int i) {return e[i];}
	
	inline vec3& operator+=(const vec3 &v2);
	inline vec3& operator-=(const vec3 &v2);
	inline vec3& operator*=(const vec3 &v2);
	inline vec3& operator/=(const vec3 &v2);
	inline vec3& operator*=(const float t);
	//inline vec3& operator/=(const float t);
	
	inline vec3& operator/=(const double t) {return *this *= 1/t;}
	
	inline float length() const {return sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);}
	inline float length_squared() const {return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];}
	inline void make_unit_vector();
	
	float e[3];

	inline static vec3 random() {
        return vec3(random_double(), random_double(), random_double());
    }

    inline static vec3 random(double min, double max) {
        return vec3(random_double(min,max), random_double(min,max), random_double(min,max));
    }

	bool near_zero() const {
        // Return true if the vector is close to zero in all dimensions.
        const auto s = 1e-8;
        return (fabs(e[0]) < s) && (fabs(e[1]) < s) && (fabs(e[2]) < s);
    }
};

//type aliases 
using point3 = vec3;   // 3D point
using color = vec3; 


inline std::istream& operator>>(std::istream &is, vec3 &t) {
 is >> t.e[0] >> t.e[1] >> t.e[2];
 return is;
}

inline std::ostream& operator<<(std::ostream &os, vec3 &t) {
 os << t.e[0] << " " << t.e[1] << " " << t.e[2];
 return os;
}

inline void vec3::make_unit_vector() {
 float k = 1.0 / sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
 e[0] *= k; e[1] *= k; e[2] *= k;
}

inline vec3 operator+(const vec3 &v1, const vec3 &v2) {
 return vec3(v1.e[0]+v2.e[0],v1.e[1]+v2.e[1],v1.e[2]+v2.e[2]);
}
inline vec3 operator-(const vec3 &v1, const vec3 &v2) {
 return vec3(v1.e[0]-v2.e[0],v1.e[1]-v2.e[1],v1.e[2]-v2.e[2]);
}
inline vec3 operator*(const vec3 &v1, const vec3 &v2) {
 return vec3(v1.e[0]*v2.e[0],v1.e[1]*v2.e[1],v1.e[2]*v2.e[2]);
}
inline vec3 operator/(const vec3 &v1, const vec3 &v2) {
 return vec3(v1.e[0]/v2.e[0],v1.e[1]/v2.e[1],v1.e[2]/v2.e[2]);
}
inline vec3 operator*(float t, const vec3 &v) {
 return vec3(t*v.e[0],t*v.e[1],t*v.e[2]);
}
inline vec3 operator*(const vec3 &v,float t) {
 return vec3(t*v.e[0],t*v.e[1],t*v.e[2]);
}
inline vec3 operator/(float t, const vec3 &v) {
 return vec3(t/v.e[0],t/v.e[1],t/v.e[2]);
}
inline vec3 operator/(const vec3 &v, float t) {
 return vec3(v.e[0]/t,v.e[1]/t,v.e[2]/t);
}

inline float dot(const vec3 &v1, const vec3 &v2) {
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

inline vec3 cross(const vec3 &v1, const vec3 &v2) {
	return vec3(v1.e[1]*v2.e[2] - v1.e[2]*v2[1], 
	-(v1.e[0]*v2.e[2] - v1.e[2]*v2[0]), 
	v1.e[0]*v2.e[1] - v1.e[1]*v2[0]);
}

inline vec3& vec3::operator+=(const vec3 &v) {
 e[0] += v.e[0];
 e[1] += v.e[1];
 e[2] += v.e[2];
 return *this;
}
inline vec3& vec3::operator-=(const vec3 &v) {
 e[0] -= v.e[0];
 e[1] -= v.e[1];
 e[2] -= v.e[2];
 return *this;
}
inline vec3& vec3::operator*=(const vec3 &v) {
 e[0] *= v.e[0];
 e[1] *= v.e[1];
 e[2] *= v.e[2];
 return *this;
}
inline vec3& vec3::operator/=(const vec3 &v) {
 e[0] /= v.e[0];
 e[1] /= v.e[1];
 e[2] /= v.e[2];
 return *this;
}
inline vec3& vec3::operator*=(const float t) {
 e[0] *= t;
 e[1] *= t;
 e[2] *= t;
 return *this;
}

inline vec3 unit_vector(vec3 v) {
 return v / v.length();
}

vec3 reflect(const vec3& v, const vec3& n) {
 return v - 2*dot(v,n)*n;
}


vec3 random_in_unit_sphere() {
    while (true) {
        auto p = vec3::random(-1,1);
        if (p.length_squared() >= 1) continue;
        return p;
    }
}

vec3 random_unit_vector() {
    return unit_vector(random_in_unit_sphere());
}

vec3 random_in_hemisphere(const vec3& normal) {
    vec3 in_unit_sphere = random_in_unit_sphere();
    if (dot(in_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
        return in_unit_sphere;
    else
        return -in_unit_sphere;
}

vec3 random_in_unit_disk() {
    while (true) {
        auto p = vec3(random_double(-1,1), random_double(-1,1), 0);
        if (p.length_squared() >= 1) continue;
        return p;
    }
}

vec3 refract(const vec3& uv, const vec3& n, double etai_over_etat) {
    auto cos_theta = fmin(dot(-uv, n), 1.0);
    vec3 r_out_perp =  etai_over_etat * (uv + cos_theta*n);
    vec3 r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.length_squared())) * n;
    return r_out_perp + r_out_parallel;
}

float schlick(float cosine, float ref_idx) {
 float r0 = (1-ref_idx) / (1+ref_idx);
 r0 = r0*r0;
 return r0 + (1-r0)*pow((1-cosine),5);
}

#endif

