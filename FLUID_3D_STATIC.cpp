//////////////////////////////////////////////////////////////////////
// This file is part of Wavelet Turbulence.
// 
// Wavelet Turbulence is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// Wavelet Turbulence is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Wavelet Turbulence.  If not, see <http://www.gnu.org/licenses/>.
// 
// Copyright 2008 Theodore Kim and Nils Thuerey
// 
// FLUID_3D.cpp: implementation of the static functions of the FLUID_3D class.
//
//////////////////////////////////////////////////////////////////////

#include <zlib.h>
#include "FLUID_3D.h"
#include "IMAGE.h"
#include "WTURBULENCE.h"
#include "INTERPOLATE.h"

//////////////////////////////////////////////////////////////////////
// add a test cube of density to the center
//////////////////////////////////////////////////////////////////////
void FLUID_3D::addSmokeColumn() {
	addSmokeTestCase(_density, _res);
	addSmokeTestCase(_heat, _res);
	if (_wTurbulence) {
		addSmokeTestCase(_wTurbulence->getDensityBig(), _wTurbulence->getResBig());
	}
}

//////////////////////////////////////////////////////////////////////
// generic static version, so that it can be applied to the
// WTURBULENCE grid as well
//////////////////////////////////////////////////////////////////////
void FLUID_3D::addSmokeTestCase(float* field, Vec3Int res)
{
	const int slabSize = res[0]*res[1]; int maxRes = (int)MAX3V(res);
	float dx = 1.0f / (float)maxRes;

	float xTotal = dx * res[0];
	float yTotal = dx * res[1];
	float zTotal = dx * res[2];

  float heighMin = 0.05;
  float heighMax = 0.10;

  for (int z = 0; z < res[2]; z++)
    for (int y = (int)(heighMin*res[1]); y <= (int)(heighMax * res[1]); y++)
      for (int x = 0; x < res[0]; x++)
      {
        float xLength = x * dx - xTotal * 0.4f;
        float zLength = z * dx - zTotal * 0.5f;
        float radius = sqrtf(xLength * xLength + zLength * zLength);

        if (radius < 0.075f * xTotal)
        {
          int index = x + y * res[0] + z * slabSize;
          field[index] = 1.0f;
        }
      }
}

//////////////////////////////////////////////////////////////////////
// add a test cube of density to the center
//////////////////////////////////////////////////////////////////////
void FLUID_3D::addSmokeBlob() {
	addSmokeBlobTestCase(_density, _res);
	addSmokeBlobTestCase(_heat, _res);
	if (_wTurbulence) {
		addSmokeBlobTestCase(_wTurbulence->getDensityBig(), _wTurbulence->getResBig());
	}
}

//////////////////////////////////////////////////////////////////////
// generic static version, so that it can be applied to the
// WTURBULENCE grid as well
//////////////////////////////////////////////////////////////////////
void FLUID_3D::addSmokeBlobTestCase(float* field, Vec3Int res)
{
	const int slabSize = res[0]*res[1]; int maxRes = (int)MAX3V(res);
	float dx = 1.0f / (float)maxRes;

	float xTotal = dx * res[0];
	float yTotal = dx * res[1];
	float zTotal = dx * res[2];

  float heighMin = 0.05;
  float heighMax = 0.10;

  for (int z = 0; z < res[2]; z++)
    for (int y = 0; y < res[1]; y++)
      for (int x = 0; x < res[0]; x++)
      {
        float xLength = x * dx - xTotal * 0.5f;
        float yLength = y * dx - yTotal * 0.5f;
        float zLength = z * dx - zTotal * 0.5f;
        float radius = sqrtf(xLength * xLength + yLength * yLength + zLength * zLength);

        if (radius < 0.03f * xTotal)
        {
          int index = x + y * res[0] + z * slabSize;
          field[index] = 1.0f;
        }
      }
}

//////////////////////////////////////////////////////////////////////
// set x direction to Neumann boundary conditions
//////////////////////////////////////////////////////////////////////
void FLUID_3D::setNeumannX(float* field, Vec3Int res)
{
	const int slabSize = res[0] * res[1];
	int index;
	for (int z = 0; z < res[2]; z++)
		for (int y = 0; y < res[1]; y++)
		{
			// left slab
			index = y * res[0] + z * slabSize;
			field[index] = field[index + 2];

			// right slab
			index += res[0] - 1;
			field[index] = field[index - 2];
		}
 }

//////////////////////////////////////////////////////////////////////
// set y direction to Neumann boundary conditions
//////////////////////////////////////////////////////////////////////
void FLUID_3D::setNeumannY(float* field, Vec3Int res)
{
	const int slabSize = res[0] * res[1];
	int index;
	for (int z = 0; z < res[2]; z++)
		for (int x = 0; x < res[0]; x++)
		{
			// bottom slab
			index = x + z * slabSize;
			field[index] = field[index + 2 * res[0]];

			// top slab
			index += slabSize - res[0];
			field[index] = field[index - 2 * res[0]];
		}

	// fix, force top slab to only allow outwards flux
	for (int z = 0; z < res[2]; z++)
		for (int x = 0; x < res[0]; x++)
		{
			// top slab
			int index = x + z * slabSize;
			index += slabSize - res[0];
			if(field[index]<0.) field[index] = 0.;
			index -= res[0];
			if(field[index]<0.) field[index] = 0.;
		}
}

//////////////////////////////////////////////////////////////////////
// set z direction to Neumann boundary conditions
//////////////////////////////////////////////////////////////////////
void FLUID_3D::setNeumannZ(float* field, Vec3Int res)
{
	const int slabSize = res[0] * res[1];
	const int totalCells = res[0] * res[1] * res[2];
	int index;
	for (int y = 0; y < res[1]; y++)
		for (int x = 0; x < res[0]; x++)
		{
			// front slab
			index = x + y * res[0];
			field[index] = field[index + 2 * slabSize];

			// back slab
			index += totalCells - slabSize;
			field[index] = field[index - 2 * slabSize];
		}
}

//////////////////////////////////////////////////////////////////////
// set x direction to zero
//////////////////////////////////////////////////////////////////////
void FLUID_3D::setZeroX(float* field, Vec3Int res)
{
	const int slabSize = res[0] * res[1];
	int index;
	for (int z = 0; z < res[2]; z++)
		for (int y = 0; y < res[1]; y++)
		{
			// left slab
			index = y * res[0] + z * slabSize;
			field[index] = 0.0f;

			// right slab
			index += res[0] - 1;
			field[index] = 0.0f;
		}
}

//////////////////////////////////////////////////////////////////////
// set y direction to zero
//////////////////////////////////////////////////////////////////////
void FLUID_3D::setZeroY(float* field, Vec3Int res)
{
	const int slabSize = res[0] * res[1];
	int index;
	for (int z = 0; z < res[2]; z++)
		for (int x = 0; x < res[0]; x++)
		{
			// bottom slab
			index = x + z * slabSize;
			field[index] = 0.0f;

			// top slab
			index += slabSize - res[0];
			field[index] = 0.0f;
		}
}

//////////////////////////////////////////////////////////////////////
// set z direction to zero
//////////////////////////////////////////////////////////////////////
void FLUID_3D::setZeroZ(float* field, Vec3Int res)
{
	const int slabSize = res[0] * res[1];
	const int totalCells = res[0] * res[1] * res[2];
	int index;
	for (int y = 0; y < res[1]; y++)
		for (int x = 0; x < res[0]; x++)
		{
			// front slab
			index = x + y * res[0];
			field[index] = 0.0f;

			// back slab
			index += totalCells - slabSize;
			field[index] = 0.0f;
		}
 }

//////////////////////////////////////////////////////////////////////
// copy grid boundary
//////////////////////////////////////////////////////////////////////
void FLUID_3D::copyBorderX(float* field, Vec3Int res)
{
	const int slabSize = res[0] * res[1];
	int index;
	for (int z = 0; z < res[2]; z++)
		for (int y = 0; y < res[1]; y++)
		{
			// left slab
			index = y * res[0] + z * slabSize;
			field[index] = field[index + 1];

			// right slab
			index += res[0] - 1;
			field[index] = field[index - 1];
		}
}
void FLUID_3D::copyBorderY(float* field, Vec3Int res)
{
	const int slabSize = res[0] * res[1];
	const int totalCells = res[0] * res[1] * res[2];
	int index;
	for (int z = 0; z < res[2]; z++)
		for (int x = 0; x < res[0]; x++)
		{
			// bottom slab
			index = x + z * slabSize;
			field[index] = field[index + res[0]]; 
			// top slab
			index += slabSize - res[0];
			field[index] = field[index - res[0]];
		}
}
void FLUID_3D::copyBorderZ(float* field, Vec3Int res)
{
	const int slabSize = res[0] * res[1];
	const int totalCells = res[0] * res[1] * res[2];
	int index;
	for (int y = 0; y < res[1]; y++)
		for (int x = 0; x < res[0]; x++)
		{
			// front slab
			index = x + y * res[0];
			field[index] = field[index + slabSize]; 
			// back slab
			index += totalCells - slabSize;
			field[index] = field[index - slabSize];
		}
}

/////////////////////////////////////////////////////////////////////
// advect field with the semi lagrangian method
//////////////////////////////////////////////////////////////////////
void FLUID_3D::advectFieldSemiLagrange(const float dt, const float* velx, const float* vely,  const float* velz,
		float* oldField, float* newField, Vec3Int res)
{
	const int xres = res[0];
	const int yres = res[1];
	const int zres = res[2];
	static int hits = 0;
	static int total = 0;
	const int slabSize = res[0] * res[1];

	// scale dt up to grid resolution
#if PARALLEL==1 && !_WIN32
#pragma omp parallel
#pragma omp for  schedule(static)
#endif
	for (int z = 0; z < zres; z++)
		for (int y = 0; y < yres; y++)
			for (int x = 0; x < xres; x++)
			{
				const int index = x + y * xres + z * xres*yres;
				
        // backtrace
				float xTrace = x - dt * velx[index];
				float yTrace = y - dt * vely[index];
				float zTrace = z - dt * velz[index];

				// clamp backtrace to grid boundaries
				if (xTrace < 0.5) xTrace = 0.5;
				if (xTrace > xres - 1.5) xTrace = xres - 1.5;
				if (yTrace < 0.5) yTrace = 0.5;
				if (yTrace > yres - 1.5) yTrace = yres - 1.5;
				if (zTrace < 0.5) zTrace = 0.5;
				if (zTrace > zres - 1.5) zTrace = zres - 1.5;

				// locate neighbors to interpolate
				const int x0 = (int)xTrace;
				const int x1 = x0 + 1;
				const int y0 = (int)yTrace;
				const int y1 = y0 + 1;
				const int z0 = (int)zTrace;
				const int z1 = z0 + 1;

				// get interpolation weights
				const float s1 = xTrace - x0;
				const float s0 = 1.0f - s1;
				const float t1 = yTrace - y0;
				const float t0 = 1.0f - t1;
				const float u1 = zTrace - z0;
				const float u0 = 1.0f - u1;

				const int i000 = x0 + y0 * xres + z0 * slabSize;
				const int i010 = x0 + y1 * xres + z0 * slabSize;
				const int i100 = x1 + y0 * xres + z0 * slabSize;
				const int i110 = x1 + y1 * xres + z0 * slabSize;
				const int i001 = x0 + y0 * xres + z1 * slabSize;
				const int i011 = x0 + y1 * xres + z1 * slabSize;
				const int i101 = x1 + y0 * xres + z1 * slabSize;
				const int i111 = x1 + y1 * xres + z1 * slabSize;

				// interpolate
				// (indices could be computed once)
				newField[index] = u0 * (s0 * (t0 * oldField[i000] +
							t1 * oldField[i010]) +
						s1 * (t0 * oldField[i100] +
							t1 * oldField[i110])) +
					u1 * (s0 * (t0 * oldField[i001] +
								t1 * oldField[i011]) +
							s1 * (t0 * oldField[i101] +
								t1 * oldField[i111]));
			}
}

/////////////////////////////////////////////////////////////////////
// advect field with the maccormack method
//
// comments are the pseudocode from selle's paper
//////////////////////////////////////////////////////////////////////
void FLUID_3D::advectFieldMacCormack(const float dt, const float* xVelocity, const float* yVelocity, const float* zVelocity, 
		float* oldField, float* newField, float* temp1, float* temp2, Vec3Int res, const float* obstacles)
{
	float* phiHatN  = temp1;
	float* phiHatN1 = temp2;
	const int sx= res[0];
	const int sy= res[1];
	const int sz= res[2];

	for (int x = 0; x < sx * sy * sz; x++)
		phiHatN[x] = phiHatN1[x] = oldField[x];

	float*& phiN    = oldField;
	float*& phiN1   = newField;

	// phiHatN1 = A(phiN)
	advectFieldSemiLagrange(  dt, xVelocity, yVelocity, zVelocity, phiN, phiHatN1, res);

	// phiHatN = A^R(phiHatN1)
	advectFieldSemiLagrange( -1.0*dt, xVelocity, yVelocity, zVelocity, phiHatN1, phiHatN, res);

	// phiN1 = phiHatN1 + (phiN - phiHatN) / 2
	const int border = 0; 
	for (int z = border; z < sz-border; z++)
		for (int y = border; y < sy-border; y++)
			for (int x = border; x < sx-border; x++) {
				int index = x + y * sx + z * sx*sy;
				phiN1[index] = phiHatN1[index] + (phiN[index] - phiHatN[index]) * 0.50f;
				//phiN1[index] = phiHatN1[index]; // debug, correction off
			}
	copyBorderX(phiN1, res);
	copyBorderY(phiN1, res);
	copyBorderZ(phiN1, res);

	// clamp any newly created extrema
	clampExtrema(dt, xVelocity, yVelocity, zVelocity, oldField, newField, res);

	// if the error estimate was bad, revert to first order
	clampOutsideRays(dt, xVelocity, yVelocity, zVelocity, oldField, newField, res, obstacles, phiHatN1);
} 


//////////////////////////////////////////////////////////////////////
// Clamp the extrema generated by the BFECC error correction
//////////////////////////////////////////////////////////////////////
void FLUID_3D::clampExtrema(const float dt, const float* velx, const float* vely,  const float* velz,
		float* oldField, float* newField, Vec3Int res)
{
	const int xres= res[0];
	const int yres= res[1];
	const int zres= res[2];
	const int slabSize = res[0] * res[1];
	for (int z = 1; z < zres-1; z++)
		for (int y = 1; y < yres-1; y++)
			for (int x = 1; x < xres-1; x++)
			{
				const int index = x + y * xres+ z * xres*yres;
				// backtrace
				float xTrace = x - dt * velx[index];
				float yTrace = y - dt * vely[index];
				float zTrace = z - dt * velz[index];

				// clamp backtrace to grid boundaries
				if (xTrace < 0.5) xTrace = 0.5;
				if (xTrace > xres - 1.5) xTrace = xres - 1.5;
				if (yTrace < 0.5) yTrace = 0.5;
				if (yTrace > yres - 1.5) yTrace = yres - 1.5;
				if (zTrace < 0.5) zTrace = 0.5;
				if (zTrace > zres - 1.5) zTrace = zres - 1.5;

				// locate neighbors to interpolate
				const int x0 = (int)xTrace;
				const int x1 = x0 + 1;
				const int y0 = (int)yTrace;
				const int y1 = y0 + 1;
				const int z0 = (int)zTrace;
				const int z1 = z0 + 1;

				const int i000 = x0 + y0 * xres + z0 * slabSize;
				const int i010 = x0 + y1 * xres + z0 * slabSize;
				const int i100 = x1 + y0 * xres + z0 * slabSize;
				const int i110 = x1 + y1 * xres + z0 * slabSize;
				const int i001 = x0 + y0 * xres + z1 * slabSize;
				const int i011 = x0 + y1 * xres + z1 * slabSize;
				const int i101 = x1 + y0 * xres + z1 * slabSize;
				const int i111 = x1 + y1 * xres + z1 * slabSize;

				float minField = oldField[i000];
				float maxField = oldField[i000];

				minField = (oldField[i010] < minField) ? oldField[i010] : minField;
				maxField = (oldField[i010] > maxField) ? oldField[i010] : maxField;

				minField = (oldField[i100] < minField) ? oldField[i100] : minField;
				maxField = (oldField[i100] > maxField) ? oldField[i100] : maxField;

				minField = (oldField[i110] < minField) ? oldField[i110] : minField;
				maxField = (oldField[i110] > maxField) ? oldField[i110] : maxField;

				minField = (oldField[i001] < minField) ? oldField[i001] : minField;
				maxField = (oldField[i001] > maxField) ? oldField[i001] : maxField;

				minField = (oldField[i011] < minField) ? oldField[i011] : minField;
				maxField = (oldField[i011] > maxField) ? oldField[i011] : maxField;

				minField = (oldField[i101] < minField) ? oldField[i101] : minField;
				maxField = (oldField[i101] > maxField) ? oldField[i101] : maxField;

				minField = (oldField[i111] < minField) ? oldField[i111] : minField;
				maxField = (oldField[i111] > maxField) ? oldField[i111] : maxField;

				newField[index] = (newField[index] > maxField) ? maxField : newField[index];
				newField[index] = (newField[index] < minField) ? minField : newField[index];
			}
}

//////////////////////////////////////////////////////////////////////
// Reverts any backtraces that go into boundaries back to first 
// order -- in this case the error correction term was totally
// incorrect
//////////////////////////////////////////////////////////////////////
void FLUID_3D::clampOutsideRays(const float dt, const float* velx, const float* vely, const float* velz, 
		float* oldField, float* newField, Vec3Int res, const float* obstacles, const float *oldAdvection)
{
	const int sx= res[0];
	const int sy= res[1];
	const int sz= res[2];
	const int slabSize = res[0] * res[1];
	for (int z = 1; z < sz-1; z++)
		for (int y = 1; y < sy-1; y++)
			for (int x = 1; x < sx-1; x++)
			{
				const int index = x + y * sx+ z * slabSize;
				// backtrace
				float xBackward = x + dt * velx[index];
				float yBackward = y + dt * vely[index];
				float zBackward = z + dt * velz[index];
				float xTrace    = x - dt * velx[index];
				float yTrace    = y - dt * vely[index];
				float zTrace    = z - dt * velz[index];

				// see if it goes outside the boundaries
				bool hasObstacle = 
					(zTrace < 1.0f)    || (zTrace > sz - 2.0f) ||
					(yTrace < 1.0f)    || (yTrace > sy - 2.0f) ||
					(xTrace < 1.0f)    || (xTrace > sx - 2.0f) ||
					(zBackward < 1.0f) || (zBackward > sz - 2.0f) ||
					(yBackward < 1.0f) || (yBackward > sy - 2.0f) ||
					(xBackward < 1.0f) || (xBackward > sx - 2.0f);
				// reuse old advection instead of doing another one...
				if(hasObstacle) { newField[index] = oldAdvection[index]; continue; }

				// clamp to prevent an out of bounds access when looking into
				// the _obstacles array
				zTrace = (zTrace < 0.5f) ? 0.5f : zTrace;
				zTrace = (zTrace > sz - 1.5f) ? sz - 1.5f : zTrace;
				yTrace = (yTrace < 0.5f) ? 0.5f : yTrace;
				yTrace = (yTrace > sy - 1.5f) ? sy - 1.5f : yTrace;
				xTrace = (xTrace < 0.5f) ? 0.5f : xTrace;
				xTrace = (xTrace > sx - 1.5f) ? sx - 1.5f : xTrace;

				// locate neighbors to interpolate,
				// do backward first since we will use the forward indices if a
				// reversion is actually necessary
				zBackward = (zBackward < 0.5f) ? 0.5f : zBackward;
				zBackward = (zBackward > sz - 1.5f) ? sz - 1.5f : zBackward;
				yBackward = (yBackward < 0.5f) ? 0.5f : yBackward;
				yBackward = (yBackward > sy - 1.5f) ? sy - 1.5f : yBackward;
				xBackward = (xBackward < 0.5f) ? 0.5f : xBackward;
				xBackward = (xBackward > sx - 1.5f) ? sx - 1.5f : xBackward;

				int x0 = (int)xBackward;
				int x1 = x0 + 1;
				int y0 = (int)yBackward;
				int y1 = y0 + 1;
				int z0 = (int)zBackward;
				int z1 = z0 + 1;
				if(obstacles && !hasObstacle) {
					hasObstacle = hasObstacle || 
						obstacles[x0 + y0 * sx + z0*slabSize] ||
						obstacles[x0 + y1 * sx + z0*slabSize] ||
						obstacles[x1 + y0 * sx + z0*slabSize] ||
						obstacles[x1 + y1 * sx + z0*slabSize] ||
						obstacles[x0 + y0 * sx + z1*slabSize] ||
						obstacles[x0 + y1 * sx + z1*slabSize] ||
						obstacles[x1 + y0 * sx + z1*slabSize] ||
						obstacles[x1 + y1 * sx + z1*slabSize] ;
				}
				// reuse old advection instead of doing another one...
				if(hasObstacle) { newField[index] = oldAdvection[index]; continue; }

				x0 = (int)xTrace;
				x1 = x0 + 1;
				y0 = (int)yTrace;
				y1 = y0 + 1;
				z0 = (int)zTrace;
				z1 = z0 + 1;
				if(obstacles && !hasObstacle) {
					hasObstacle = hasObstacle || 
						obstacles[x0 + y0 * sx + z0*slabSize] ||
						obstacles[x0 + y1 * sx + z0*slabSize] ||
						obstacles[x1 + y0 * sx + z0*slabSize] ||
						obstacles[x1 + y1 * sx + z0*slabSize] ||
						obstacles[x0 + y0 * sx + z1*slabSize] ||
						obstacles[x0 + y1 * sx + z1*slabSize] ||
						obstacles[x1 + y0 * sx + z1*slabSize] ||
						obstacles[x1 + y1 * sx + z1*slabSize] ;
				} // obstacle array
				// reuse old advection instead of doing another one...
				if(hasObstacle) { newField[index] = oldAdvection[index]; continue; }

				// see if either the forward or backward ray went into
				// a boundary
				if (hasObstacle) {
					// get interpolation weights
					float s1 = xTrace - x0;
					float s0 = 1.0f - s1;
					float t1 = yTrace - y0;
					float t0 = 1.0f - t1;
					float u1 = zTrace - z0;
					float u0 = 1.0f - u1;

					const int i000 = x0 + y0 * sx + z0 * slabSize;
					const int i010 = x0 + y1 * sx + z0 * slabSize;
					const int i100 = x1 + y0 * sx + z0 * slabSize;
					const int i110 = x1 + y1 * sx + z0 * slabSize;
					const int i001 = x0 + y0 * sx + z1 * slabSize;
					const int i011 = x0 + y1 * sx + z1 * slabSize;
					const int i101 = x1 + y0 * sx + z1 * slabSize;
					const int i111 = x1 + y1 * sx + z1 * slabSize;

					// interpolate, (indices could be computed once)
					newField[index] = u0 * (s0 * (
								t0 * oldField[i000] +
								t1 * oldField[i010]) +
							s1 * (t0 * oldField[i100] +
								t1 * oldField[i110])) +
						u1 * (s0 * (t0 * oldField[i001] +
									t1 * oldField[i011]) +
								s1 * (t0 * oldField[i101] +
									t1 * oldField[i111])); 
				}
			} // xyz
}

//////////////////////////////////////////////////////////////////////
// image output
//////////////////////////////////////////////////////////////////////
void FLUID_3D::writeImageSliceXY(const float *field, Vec3Int res, int slice, string prefix, int picCnt, float scale) {
  writeProjectedIntern(field, res, 0,1, prefix, picCnt, scale);
}
void FLUID_3D::writeImageSliceYZ(const float *field, Vec3Int res, int slice, string prefix, int picCnt, float scale) {
  writeProjectedIntern(field, res, 1,2, prefix, picCnt, scale);
}
void FLUID_3D::writeImageSliceXZ(const float *field, Vec3Int res, int slice, string prefix, int picCnt, float scale) {
  writeProjectedIntern(field, res, 0,2, prefix, picCnt, scale);
}

//////////////////////////////////////////////////////////////////////
// Helper function for projecting densities along a dimension
//////////////////////////////////////////////////////////////////////
static int getOtherDir(int dir1, int dir2) {
	switch(dir1) {
		case 0:
			switch(dir2) {
				case 1: return 2;
				case 2: return 1; }
			break;
		case 1:
			switch(dir2) {
				case 0: return 2;
				case 2: return 0; }
			break;
		case 2:
			switch(dir2) {
				case 0: return 1;
				case 1: return 0; }
			break;
	}
}

//////////////////////////////////////////////////////////////////////
// average densities along third spatial direction
//////////////////////////////////////////////////////////////////////
void FLUID_3D::writeProjectedIntern(const float *field, Vec3Int res, 
		int dir1, int dir2, string prefix, int picCnt, float scale) {
	const int nitems = res[dir1]*res[dir2];
	const int otherDir = getOtherDir(dir1,dir2);
	float *buf = new float[nitems];
	Vec3Int min = Vec3Int(0);
	Vec3Int max = res;

	min[otherDir] = 0;
	max[otherDir] = res[otherDir];
	float div = 1./(float)MIN3V(res); // normalize for shorter sides, old: res[otherDir];
	div *= 4.; //slightly increase contrast
	for(int i=0; i<nitems; i++) buf[i] = 0.;

	Vec3Int cnt = 0;
	for (cnt[2] = min[2]; cnt[2] < max[2]; cnt[2]++) {
		for (cnt[1] = min[1]; cnt[1] < max[1]; cnt[1]++)
			for (cnt[0] = min[0]; cnt[0] < max[0]; cnt[0]++)
			{
				const int index = cnt[0] + cnt[1] * res[0] + cnt[2] * res[0]*res[1];
				const int bufindex = cnt[dir1] + cnt[dir2] * res[dir1];
				buf[bufindex] += field[index] * scale *div;
			}
	}
	IMAGE::dumpNumberedPNG(picCnt, prefix, buf, res[dir1], res[dir2]);
	delete[] buf;
}
