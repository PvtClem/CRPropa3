#ifndef CRPROPA_GRID_H
#define CRPROPA_GRID_H

#include "crpropa/Referenced.h"
#include "crpropa/Vector3.h"

#include "kiss/string.h"
#include "kiss/logger.h"

#include <vector>
#include <immintrin.h>
#include <smmintrin.h>

namespace crpropa {

/** If set to TRILINEAR, use trilinear interpolation (standard)
If set to TRICUBIC, use tricubic interpolation instead of trilinear interpolation
If set to NEAREST_NEIGHBOUR , use nearest neighbour interpolation instead of trilinear interpolation */
enum interpolationType {
  TRILINEAR = 0,
  TRICUBIC,
  NEAREST_NEIGHBOUR
};

/** Lower and upper neighbour in a periodically continued unit grid */
inline void periodicClamp(double x, int n, int &lo, int &hi) {
	lo = ((int(floor(x)) % n) + n) % n;
	hi = (lo + 1) % n;
}

inline int periodicBoundary(int index, int n) {
	return ((index % n) + n) % n;
}

/** Lower and upper neighbour in a reflectively repeated unit grid */
inline void reflectiveClamp(double x, int n, int &lo, int &hi) {
	while ((x < 0) or (x > n))
		x = 2 * n * (x > n) - x;
	lo = floor(x);
	hi = lo + (lo < n-1);
}

/** Symmetrical round */
inline double round(double r) {
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}

/**
 * \addtogroup Core
 * @{
 */

/**
 @class GridProperties
 @brief Combines parameters that uniquely define Grid class
 */
class GridProperties: public Referenced {
public:
	size_t Nx, Ny, Nz;
	Vector3d origin;
	Vector3d spacing;
	bool reflective;

	/** Constructor for cubic grid
	 @param	origin	Position of the lower left front corner of the volume
	 @param	N		Number of grid points in one direction
	 @param spacing	Spacing between grid points
	 */
	GridProperties(Vector3d origin, size_t N, double spacing) :
		origin(origin), Nx(N), Ny(N), Nz(N), spacing(Vector3d(spacing)), reflective(false) {
	}

	/** Constructor for non-cubic grid
	 @param	origin	Position of the lower left front corner of the volume
	 @param	Nx		Number of grid points in x-direction
	 @param	Ny		Number of grid points in y-direction
	 @param	Nz		Number of grid points in z-direction
	 @param spacing	Spacing between grid points
	 */
	GridProperties(Vector3d origin, size_t Nx, size_t Ny, size_t Nz, double spacing) :
		origin(origin), Nx(Nx), Ny(Ny), Nz(Nz), spacing(Vector3d(spacing)), reflective(false) {
	}

	/** Constructor for non-cubic grid with spacing vector
	 @param	origin	Position of the lower left front corner of the volume
	 @param	Nx		Number of grid points in x-direction
	 @param	Ny		Number of grid points in y-direction
	 @param	Nz		Number of grid points in z-direction
	 @param spacing	Spacing vector between grid points
	*/
	GridProperties(Vector3d origin, size_t Nx, size_t Ny, size_t Nz, Vector3d spacing) :
		origin(origin), Nx(Nx), Ny(Ny), Nz(Nz), spacing(spacing), reflective(false) {
	}

	virtual ~GridProperties() {
	}

	void setReflective(bool b) {
		reflective = b;
	}
};

/**
 @class Grid
 @brief Template class for fields on a periodic grid with trilinear interpolation

 The grid spacing is constant with diffrent resulution along all three axes.
 Values are calculated by trilinear interpolation of the surrounding 8 grid points.
 The grid is periodically (default) or reflectively extended.
 The grid sample positions are at 1/2 * size/N, 3/2 * size/N ... (2N-1)/2 * size/N.
 */
template<typename T>
class Grid: public Referenced {
	std::vector<T> grid;
	size_t Nx, Ny, Nz; /**< Number of grid points */
	Vector3d origin; /**< Origin of the volume that is represented by the grid. */
	Vector3d gridOrigin; /**< Grid origin */
	Vector3d spacing; /**< Distance between grid points, determines the extension of the grid */
	bool reflective; /**< If set to true, the grid is repeated reflectively instead of periodically */

	//  grid.setInterpolationType(INTERPOLATION_TYPE_LINEAR)

  interpolationType ipolType;

public:
	/** Constructor for cubic grid
	 @param	origin	Position of the lower left front corner of the volume
	 @param	N		Number of grid points in one direction
	 @param spacing	Spacing between grid points
	 */
	Grid(Vector3d origin, size_t N, double spacing) {
		setOrigin(origin);
		setGridSize(N, N, N);
		setSpacing(Vector3d(spacing));
		setReflective(false);
	}

	/** Constructor for non-cubic grid
	 @param	origin	Position of the lower left front corner of the volume
	 @param	Nx		Number of grid points in x-direction
	 @param	Ny		Number of grid points in y-direction
	 @param	Nz		Number of grid points in z-direction
	 @param spacing	Spacing between grid points
	 */
	Grid(Vector3d origin, size_t Nx, size_t Ny, size_t Nz, double spacing) {
		setOrigin(origin);
		setGridSize(Nx, Ny, Nz);
		setSpacing(Vector3d(spacing));
		setReflective(false);
	}

	/** Constructor for non-cubic grid with spacing vector
	 @param	origin	Position of the lower left front corner of the volume
	 @param	Nx		Number of grid points in x-direction
	 @param	Ny		Number of grid points in y-direction
	 @param	Nz		Number of grid points in z-direction
	 @param spacing	Spacing vector between grid points
	*/
	Grid(Vector3d origin, size_t Nx, size_t Ny, size_t Nz, Vector3d spacing) {
	 	setOrigin(origin);
	 	setGridSize(Nx, Ny, Nz);
	 	setSpacing(spacing);
	 	setReflective(false);
	}

	/** Constructor for GridProperties
 	 @param p	GridProperties instance
     */
	Grid(const GridProperties &p) :
		origin(p.origin), spacing(p.spacing), reflective(p.reflective) {
	 	setGridSize(p.Nx, p.Ny, p.Nz);
	}

	void setOrigin(Vector3d origin) {
		this->origin = origin;
		this->gridOrigin = origin + spacing/2;
	}

	/** Resize grid, also enlarges the volume as the spacing stays constant */
	void setGridSize(size_t Nx, size_t Ny, size_t Nz) {
		this->Nx = Nx;
		this->Ny = Ny;
		this->Nz = Nz;
		grid.resize(Nx * Ny * Nz);
		setOrigin(origin);
	}

	void setSpacing(Vector3d spacing) {
		this->spacing = spacing;
		setOrigin(origin);
	}

	void setReflective(bool b) {
		reflective = b;
	}

	void setInterpolationType(interpolationType ipolType) {
	  if (ipolType == TRILINEAR || ipolType == TRICUBIC || ipolType == NEAREST_NEIGHBOUR) {
	    ipolType = ipolType;
	    std::cout << "changed interpolation routine" << std::endl;
	  } else {
	    throw std::runtime_error("InterpolationType: unknown interpolation type");
	  }
	}

	Vector3d getOrigin() const {
		return origin;
	}
	size_t getNx() const {
		return Nx;
	}

	size_t getNy() const {
		return Ny;
	}

	size_t getNz() const {
		return Nz;
	}

	/** Calculates the total size of the grid in bytes */
	size_t getSizeOf() const {
		return sizeof(grid) + (sizeof(grid[0]) * grid.size());
	}

	Vector3d getSpacing() const {
		return spacing;
	}

	bool isReflective() const {
		return reflective;
	}

	T interpolate(const Vector3d &position) {
    if (ipolType == TRICUBIC)
        return tricubicInterpolate(T(), position);
    else if (ipolType == NEAREST_NEIGHBOUR)
      return nearestNeighbourInterpolate(position);
    else
      return trilinearInterpolate(position);
  }

	/** Inspector & Mutator */
	T &get(size_t ix, size_t iy, size_t iz) {
		return grid[ix * Ny * Nz + iy * Nz + iz];
	}

	/** Inspector */
	const T &get(size_t ix, size_t iy, size_t iz) const {
		return grid[ix * Ny * Nz + iy * Nz + iz];
	}
	
	const T &periodicGet(size_t ix, size_t iy, size_t iz) const {
		ix = periodicBoundary(ix, Nx);
		iy = periodicBoundary(iy, Ny);
		iz = periodicBoundary(iz, Nz);
		return grid[ix * Ny * Nz + iy * Nz + iz];
	}

	T getValue(size_t ix, size_t iy, size_t iz) {
		return grid[ix * Ny * Nz + iy * Nz + iz];
	}

	void setValue(size_t ix, size_t iy, size_t iz, T value) {
		grid[ix * Ny * Nz + iy * Nz + iz] = value;
	}

	/** Return a reference to the grid values */
	std::vector<T> &getGrid() {
		return grid;
	}

	/** Position of the grid point of a given index */
	Vector3d positionFromIndex(int index) const {
		int ix = index / (Ny * Nz);
		int iy = (index / Nz) % Ny;
		int iz = index % Nz;
		return Vector3d(ix, iy, iz) * spacing + gridOrigin;
	}

	/** Value of a grid point that is closest to a given position */
	T closestValue(const Vector3d &position) const {
		Vector3d r = (position - gridOrigin) / spacing;
		int ix = round(r.x);
		int iy = round(r.y);
		int iz = round(r.z);
		if (reflective) {
			/** TODO: this is a fix for reflective boundaries, I think
			make sure this gets reviewed appropriately */
			while ((ix < 0) or (ix >= Nx))
				ix = 2 * Nx * (ix >= Nx) - ix - 1;
			while ((iy < 0) or (iy >= Ny))
				iy = 2 * Ny * (iy >= Ny) - iy - 1;
			while ((iz < 0) or (iz >= Nz))
				iz = 2 * Nz * (iz >= Nz) - iz - 1;
		} else {
			ix = ((ix % Nx) + Nx) % Nx;
			iy = ((iy % Ny) + Ny) % Ny;
			iz = ((iz % Nz) + Nz) % Nz;
		}
		return get(ix, iy, iz);
	}



private:	

  __m128 simdPeriodicGet(size_t ix, size_t iy, size_t iz) const {
		ix = periodicBoundary(ix, Nx);
		iy = periodicBoundary(iy, Ny);
		iz = periodicBoundary(iz, Nz);
		return convertVector3fToSimd(grid[ix * Ny * Nz + iy * Nz + iz]);
	}

  __m128 convertVector3fToSimd(const Vector3f v) const {
		__m128 simdVar = _mm_set_ps(0,v.z,v.y,v.x);
		return simdVar;
	}
	
  Vector3d convertSimdToVector3f(__m128 res) const {
		float vec[4];	
		_mm_store_ps(&vec[0], res);
		Vector3f result = Vector3d(vec[0], vec[1], vec[2]);
		return result;
	}

  /** Vectorized cubic Interpolator in 1D */
  __m128 CubicInterpolate(__m128 p0,__m128 p1,__m128 p2,__m128 p3,double position) const {
     __m128 c1 = _mm_set1_ps (1/2.);
     __m128 c2 = _mm_set1_ps (3/2.);
     __m128 c3 = _mm_set1_ps (2.);
     __m128 c4 = _mm_set1_ps (5/2.);

     __m128 pos  = _mm_set1_ps (position);
     __m128 pos2 = _mm_set1_ps (position*position);
     __m128 pos3 = _mm_set1_ps (position*position*position);

     __m128 res = _mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_sub_ps(_mm_add_ps(_mm_mul_ps(c2,p1),_mm_mul_ps(c1,p3)),_mm_add_ps(_mm_mul_ps(c1,p0),_mm_mul_ps(c2,p2))),pos3),
              _mm_mul_ps(_mm_sub_ps(_mm_add_ps(p0,_mm_mul_ps(c3,p2)),_mm_add_ps(_mm_mul_ps(c4,p1),_mm_mul_ps(c1,p3))),pos2)) , _mm_mul_ps(_mm_sub_ps(_mm_mul_ps(c1,p2),_mm_mul_ps(c1,p0)),pos)) ,p1);
     return res;
  }

  /** Interpolate the grid tricubic at a given position */
	Vector3d tricubicInterpolate(Vector3f, const Vector3d &position) const {
		// position on a unit grid
		Vector3d r = (position - gridOrigin) / spacing;

		/** indices of lower and upper neighbours */
		int iX0, iX1, iY0, iY1, iZ0, iZ1;
		if (reflective) {
			reflectiveClamp(r.x, Nx, iX0, iX1);
			reflectiveClamp(r.y, Ny, iY0, iY1);
			reflectiveClamp(r.z, Nz, iZ0, iZ1);
		} else {
			periodicClamp(r.x, Nx, iX0, iX1);
			periodicClamp(r.y, Ny, iY0, iY1);
			periodicClamp(r.z, Nz, iZ0, iZ1);
		}
		
		double fX = r.x - floor(r.x);
		double fY = r.y - floor(r.y);
		double fZ = r.z - floor(r.z);

		int nrCubicInterpolations = 4;
    __m128 interpolateVaryX[nrCubicInterpolations];
    __m128 interpolateVaryY[nrCubicInterpolations];
    __m128 interpolateVaryZ[nrCubicInterpolations];
    for (int iLoopX = -1; iLoopX < nrCubicInterpolations-1; iLoopX++) {
      for (int iLoopY = -1; iLoopY < nrCubicInterpolations-1; iLoopY++) {
        for (int iLoopZ = -1; iLoopZ < nrCubicInterpolations-1; iLoopZ++) {
          interpolateVaryZ[iLoopZ+1] = simdPeriodicGet(iX0+iLoopX, iY0+iLoopY, iZ0+iLoopZ);
        }
        interpolateVaryY[iLoopY+1] = CubicInterpolate(interpolateVaryZ[0], interpolateVaryZ[1], interpolateVaryZ[2], interpolateVaryZ[3], fZ);
      }
      interpolateVaryX[iLoopX+1] = CubicInterpolate(interpolateVaryY[0], interpolateVaryY[1], interpolateVaryY[2], interpolateVaryY[3], fY);
    }
    __m128 result = CubicInterpolate(interpolateVaryX[0], interpolateVaryX[1], interpolateVaryX[2], interpolateVaryX[3], fX);
		return convertSimdToVector3f(result);
	}


  //Vectorized cubic Interpolator in 1D
  double CubicInterpolateScalar(double p0,double p1,double p2,double p3,double pos) const
  {
     return((-0.5*p0+3/2.*p1-3/2.*p2+0.5*p3)*pos*pos*pos+(p0-5/2.*p1+p2*2-0.5*p3)*pos*pos+(-0.5*p0+0.5*p2)*pos+p1);
  }

  /** Interpolate the grid tricubic at a given position */
	double tricubicInterpolate(double, const Vector3d &position) const {
		/** position on a unit grid */
		Vector3d r = (position - gridOrigin) / spacing;

		/** indices of lower and upper neighbors */
		int iX, iX1, iY, iY1, iZ, iZ1;
		if (reflective) {
			reflectiveClamp(r.x, Nx, iX, iX1);
			reflectiveClamp(r.y, Ny, iY, iY1);
			reflectiveClamp(r.z, Nz, iZ, iZ1);
		} else {
			periodicClamp(r.x, Nx, iX, iX1);
			periodicClamp(r.y, Ny, iY, iY1);
			periodicClamp(r.z, Nz, iZ, iZ1);
		}

    double fX = r.x - floor(r.x);
		double fY = r.y - floor(r.y);
		double fZ = r.z - floor(r.z);

    int nrCubicInterpolations = 4;
		double interpolateVaryX[nrCubicInterpolations];
		double interpolateVaryY[nrCubicInterpolations];
		double interpolateVaryZ[nrCubicInterpolations];
		for (int iLoopX = -1; iLoopX < nrCubicInterpolations-1; iLoopX++) {
      for (int iLoopY = -1; iLoopY < nrCubicInterpolations-1; iLoopY++) {
        for (int iLoopZ = -1; iLoopZ < nrCubicInterpolations-1; iLoopZ++) {
          interpolateVaryZ[iLoopZ+1] = periodicGet(iX+iLoopX, iY+iLoopY, iZ+iLoopZ);
        }
        interpolateVaryY[iLoopY+1] = CubicInterpolateScalar(interpolateVaryZ[0], interpolateVaryZ[1], interpolateVaryZ[2], interpolateVaryZ[3], fZ);
      }
      interpolateVaryX[iLoopX+1] = CubicInterpolateScalar(interpolateVaryY[0], interpolateVaryY[1], interpolateVaryY[2], interpolateVaryY[3], fY);
    }
    double result = CubicInterpolateScalar(interpolateVaryX[0], interpolateVaryX[1], interpolateVaryX[2], interpolateVaryX[3], fX);
		return result;
	}


	/** Interpolate the grid trilinear at a given position */
	T trilinearInterpolate(const Vector3d &position) const {
		/** position on a unit grid */
		Vector3d r = (position - gridOrigin) / spacing;

		/** indices of lower (0) and upper (1) neighbours. The neighbours span a grid
		with the origin at [iX0, iY0, iZ0] and the most distant corner [iX1, iY1, iZ1]. */
		int iX0, iX1, iY0, iY1, iZ0, iZ1;
		if (reflective) {
			reflectiveClamp(r.x, Nx, iX0, iX1);
			reflectiveClamp(r.y, Ny, iY0, iY1);
			reflectiveClamp(r.z, Nz, iZ0, iZ1);
		} else {
			periodicClamp(r.x, Nx, iX0, iX1);
			periodicClamp(r.y, Ny, iY0, iY1);
			periodicClamp(r.z, Nz, iZ0, iZ1);
		}

		// linear fraction to lower and upper neighbours
		double fX0 = r.x - floor(r.x);
		double fX1 = 1 - fX0;
		double fY0 = r.y - floor(r.y);
		double fY1 = 1 - fY0;
		double fZ0 = r.z - floor(r.z);
		double fZ1 = 1 - fZ0;

		/** trilinear interpolation (see http://paulbourke.net/miscellaneous/interpolation) */
		T b(0.);
		b += get(iX0, iY0, iZ0) * fX1 * fY1 * fZ1;
		b += get(iX1, iY0, iZ0) * fX0 * fY1 * fZ1;
		b += get(iX0, iY1, iZ0) * fX1 * fY0 * fZ1;
		b += get(iX0, iY0, iZ1) * fX1 * fY1 * fZ0;
		b += get(iX1, iY0, iZ1) * fX0 * fY1 * fZ0;
		b += get(iX0, iY1, iZ1) * fX1 * fY0 * fZ0;
		b += get(iX1, iY1, iZ0) * fX0 * fY0 * fZ1;
		b += get(iX1, iY1, iZ1) * fX0 * fY0 * fZ0;

		return b;
	}

	/** Interpolate the grid at a given position using the nearest neighbour interpolation */
	T nearestNeighbourInterpolate(const Vector3d &position) const {
		return closestValue(position);
  }
};

typedef Grid<Vector3f> Grid3f;
typedef Grid<Vector3d> Grid3d;
typedef Grid<float> Grid1f;
typedef Grid<double> Grid1d;

// DEPRICATED: Will be removed in CRPropa v3.9
class VectorGrid: public Grid3f {
	void printDeprication() const {
		KISS_LOG_WARNING << "VectorGrid is deprecated and will be removed in the future. Replace it with Grid3f (float) or Grid3d (double).";
	}
public:
	VectorGrid(Vector3d origin, size_t N, double spacing) : Grid3f(origin, N, spacing) {
		printDeprication();
	}

	VectorGrid(Vector3d origin, size_t Nx, size_t Ny, size_t Nz, double spacing) : Grid3f(origin, Nx, Ny, Nz, spacing) {
		printDeprication();
	}

	VectorGrid(Vector3d origin, size_t Nx, size_t Ny, size_t Nz, Vector3d spacing) : Grid3f(origin, Nx, Ny, Nz, spacing) {
		printDeprication();
	}
};

// DEPRICATED: Will be removed in CRPropa v3.9
class ScalarGrid: public Grid1f {
	void printDeprication() const {
		KISS_LOG_WARNING << "ScalarGrid is deprecated and will be removed in the future. Replace with Grid1f (float) or Grid1d (double).";
	}
public:
	ScalarGrid(Vector3d origin, size_t N, double spacing) : Grid1f(origin, N, spacing) {
		printDeprication();
	}

	ScalarGrid(Vector3d origin, size_t Nx, size_t Ny, size_t Nz, double spacing) : Grid1f(origin, Nx, Ny, Nz, spacing) {
		printDeprication();
	}

	ScalarGrid(Vector3d origin, size_t Nx, size_t Ny, size_t Nz, Vector3d spacing) : Grid1f(origin, Nx, Ny, Nz, spacing) {
		printDeprication();
	}
};


/** @}*/

} // namespace crpropa

#endif // CRPROPA_GRID_H
