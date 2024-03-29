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
// WTURBULENCE handling
///////////////////////////////////////////////////////////////////////////////////

#include "WTURBULENCE.h"
#include "INTERPOLATE.h"
#include "IMAGE.h"
#include <MERSENNETWISTER.h>
#include "WAVELET_NOISE.h"
#include "EIGENVALUE_HELPER.h"
#include "LU_HELPER.h"
#include "SPHERE.h"
#include <zlib.h>

// needed to access static advection functions
#include "FLUID_3D.h"

#if PARALLEL==1
#include <omp.h>
#endif // PARALLEL 

// 2^ {-5/6}
static const float persistence = 0.56123f;

//////////////////////////////////////////////////////////////////////
// constructor
//////////////////////////////////////////////////////////////////////
WTURBULENCE::WTURBULENCE(int xResSm, int yResSm, int zResSm, int amplify)
{
  // if noise magnitude is below this threshold, its contribution
  // is negilgible, so stop evaluating new octaves
  _cullingThreshold = 1e-3;

  // factor by which to increase the simulation resolution
  _amplify = amplify;

  // manually adjust the overall amount of turbulence
  _strength = 2.;

  // add the corresponding octaves of noise
  _octaves = log((float)_amplify) / log(2.0f);
  
  // noise resolution
  _xResBig = _amplify * xResSm;
  _yResBig = _amplify * yResSm;
  _zResBig = _amplify * zResSm;
  _resBig = Vec3Int(_xResBig, _yResBig, _zResBig);
  _invResBig = Vec3(1./(float)_resBig[0], 1./(float)_resBig[1], 1./(float)_resBig[2]);
  _slabSizeBig = _xResBig*_yResBig;
  _totalCellsBig = _slabSizeBig * _zResBig;

  // original / small resolution
  _xResSm = xResSm;
  _yResSm = yResSm;
  _zResSm = zResSm;
  _resSm = Vec3Int(xResSm, yResSm, zResSm);
  _invResSm = Vec3(1./(float)_resSm[0], 1./(float)_resSm[1], 1./(float)_resSm[2] );
  _slabSizeSm = _xResSm*_yResSm;
  _totalCellsSm = _slabSizeSm * _zResSm;

  // allocate high resolution density field
  _totalStepsBig = 0;
  _densityBig = new float[_totalCellsBig];
  _densityBigOld = new float[_totalCellsBig]; 

  // allocate high resolution velocity field. Note that this is only
  // necessary because we use MacCormack advection. For semi-Lagrangian
  // advection, these arrays are not necessary.
  _tempBig1 = _tempBig2 = 
  _bigUx = _bigUy = _bigUz = NULL;
  _tempBig1 = new float[_totalCellsBig];
  _tempBig2 = new float[_totalCellsBig];
  _bigUx = new float[_totalCellsBig];
  _bigUy = new float[_totalCellsBig];
  _bigUz = new float[_totalCellsBig]; 
  for(int i = 0; i < _totalCellsBig; i++) {
    _densityBig[i] = 
    _densityBigOld[i] = 
    _bigUx[i] =  
    _bigUy[i] =  
    _bigUz[i] =  
    _tempBig1[i] =  
    _tempBig2[i] =  0.;
  }

  // allocate & init texture coordinates
  _tcU = new float[_totalCellsSm];
  _tcV = new float[_totalCellsSm];
  _tcW = new float[_totalCellsSm];
  _tcTemp = new float[_totalCellsSm];

  // allocate & init energy terms
  _energy = new float[_totalCellsSm];
  _highFreqEnergy = new float[_totalCellsSm];

  // map all 
  const float dx = 1./(float)(_resSm[0]);
  const float dy = 1./(float)(_resSm[1]);
  const float dz = 1./(float)(_resSm[2]);
  int index = 0;
  for (int z = 0; z < _zResSm; z++) 
    for (int y = 0; y < _yResSm; y++) 
      for (int x = 0; x < _xResSm; x++, index++)
      {
        _tcU[index] = x*dx;
        _tcV[index] = y*dy;
        _tcW[index] = z*dz;
        _tcTemp[index] = 0.;
        _energy[index] = 0.;
      }

  // allocate eigenvalue arrays
  _eigMin = new float[_totalCellsSm];
  _eigMax = new float[_totalCellsSm];
  for(int i=0; i < _totalCellsSm; i++) 
    _eigMin[i] = _eigMax[i] =  0.;

  // noise tiles
  int n3 = noiseTileSize * noiseTileSize * noiseTileSize;
  _noiseTile = new float[n3];
  std::string noiseTileFilename = std::string("noise.wavelets");
  //std::string noiseTileFilename = std::string("noise.fft");
  generateTile(_noiseTile, noiseTileFilename);

/*
  unsigned resolution = noiseTileSize;
  float K_ = 1.0;
  float F_0_ = _xResSm / (2 * M_PI);
  float a_ = 4 / F_0_;
  float omega_0_ = M_PI / 4.0;
  float omega_1_ = M_PI / 4.0;
  float number_of_impulses_per_kernel = 64.0;
  unsigned period = noiseTileSize;
  unsigned random_offset = std::time(0);

  _noise = new noise3d(K_, a_, F_0_, omega_0_, omega_1_, number_of_impulses_per_kernel, period, random_offset);
  
  int c = 0;
  for (g_t = 0; g_t < 2*M_PI; g_t += 2*M_PI / 100) {
    //float minVal = 1e10, maxVal = -1e10;
    
    stringstream ss;
    ss << "gabor_noise" << c << ".gabor";
    c++;

    for (int i = 0; i < n3; i++) {
      _noiseTile[i] = (*_noise)(i%noiseTileSize, (i/noiseTileSize)%noiseTileSize, i/noiseTileSize/noiseTileSize);
      //minVal = std::min(minVal, _noiseTile[i]);
      //maxVal = std::max(maxVal, _noiseTile[i]);
    }

    saveTile(_noiseTile, ss.str());
    std::cout << "Generating noise tile " << c << endl;
  }
  //double c = (minVal + maxVal) / 2;
  //for (int i = 0; i < n3; i++) {
    //_noiseTile[i] -= c;
    //_noiseTile[i] /= 0.1 * (maxVal - minVal);
  //}
*/
  //std::cout << "Init wturb " << c << " " << maxVal - minVal << " " << minVal << " " << maxVal << std::endl;
}

//////////////////////////////////////////////////////////////////////
// destructor
//////////////////////////////////////////////////////////////////////
WTURBULENCE::~WTURBULENCE() {
  delete[] _densityBig;
  delete[] _densityBigOld;

  delete[] _bigUx;
  delete[] _bigUy;
  delete[] _bigUz;
  delete[] _tempBig2;
  delete[] _tempBig1;
  delete[] _tempBig2;

  delete[] _tcU;
  delete[] _tcV;
  delete[] _tcW;
  delete[] _tcTemp;

  delete[] _eigMin;
  delete[] _eigMax;
  delete[] _noiseTile;

  delete[] _energy;
  delete[] _highFreqEnergy;
}

//////////////////////////////////////////////////////////////////////
// Get the smallest valid x derivative
//
// Takes the one-sided finite difference in both directions and
// selects the smaller of the two
//////////////////////////////////////////////////////////////////////
static float minDx(int x, int y, int z, float* input, Vec3Int res)
{
  const int index = x + y * res[0] + z * res[0] * res[1];
  const int maxx = res[0]-2;

  // get grid values
  float center = input[index];
  float left  = (x <= 1)    ? FLT_MAX : input[index - 1];
  float right = (x >= maxx) ? FLT_MAX : input[index + 1];

  const float dx = res[0];

  // get all the derivative estimates
  float dLeft   = (x <= 1)     ? FLT_MAX : (center - left) * dx;
  float dRight  = (x >= maxx)  ? FLT_MAX : (right - center) * dx;
  float dCenter = (x <= 1 || x >= maxx) ? FLT_MAX : (right - left) * dx * 0.5f;

  // if it's on a boundary, only one estimate is valid
  if (x <= 1) return dRight;
  if (x >= maxx) return dLeft;

  // if it's not on a boundary, get the smallest one
  float finalD;
  finalD = (fabs(dCenter) < fabs(dRight)) ? dCenter : dRight;
  finalD = (fabs(finalD)  < fabs(dLeft))  ? finalD  : dLeft;

  return finalD;
}

//////////////////////////////////////////////////////////////////////
// get the smallest valid y derivative
//
// Takes the one-sided finite difference in both directions and
// selects the smaller of the two
//////////////////////////////////////////////////////////////////////
static float minDy(int x, int y, int z, float* input, Vec3Int res)
{
  const int index = x + y * res[0] + z * res[0] * res[1];
  const int maxy = res[1]-2;

  // get grid values
  float center = input[index];
  float down  = (y <= 1) ? FLT_MAX : input[index - res[0]];
  float up = (y >= maxy) ? FLT_MAX : input[index + res[0]];

  const float dx = res[1]; // only for square domains

  // get all the derivative estimates
  float dDown   = (y <= 1)  ? FLT_MAX : (center - down) * dx;
  float dUp  = (y >= maxy)  ? FLT_MAX : (up - center) * dx;
  float dCenter = (y <= 1 || y >= maxy) ? FLT_MAX : (up - down) * dx * 0.5f;

  // if it's on a boundary, only one estimate is valid
  if (y <= 1) return dUp;
  if (y >= maxy) return dDown;

  // if it's not on a boundary, get the smallest one
  float finalD = (fabs(dCenter) < fabs(dUp)) ? dCenter : dUp;
  finalD = (fabs(finalD) < fabs(dDown)) ? finalD : dDown;

  return finalD;
}

//////////////////////////////////////////////////////////////////////
// get the smallest valid z derivative
//
// Takes the one-sided finite difference in both directions and
// selects the smaller of the two
//////////////////////////////////////////////////////////////////////
static float minDz(int x, int y, int z, float* input, Vec3Int res)
{
  const int slab = res[0]*res[1];
  const int index = x + y * res[0] + z * slab;
  const int maxz = res[2]-2;

  // get grid values
  float center = input[index];
  float front  = (z <= 1) ? FLT_MAX : input[index - slab];
  float back = (z >= maxz) ? FLT_MAX : input[index + slab];

  const float dx = res[2]; // only for square domains

  // get all the derivative estimates
  float dfront   = (z <= 1)  ? FLT_MAX : (center - front) * dx;
  float dback  = (z >= maxz)  ? FLT_MAX : (back - center) * dx;
  float dCenter = (z <= 1 || z >= maxz) ? FLT_MAX : (back - front) * dx * 0.5f;

  // if it's on a boundary, only one estimate is valid
  if (z <= 1) return dback;
  if (z >= maxz) return dfront;

  // if it's not on a boundary, get the smallest one
  float finalD = (fabs(dCenter) < fabs(dback)) ? dCenter : dback;
  finalD = (fabs(finalD) < fabs(dfront)) ? finalD : dfront;

  return finalD;
}

//////////////////////////////////////////////////////////////////////
// handle texture coordinates (advection, reset, eigenvalues), 
// Beware -- uses big density maccormack as temporary arrays
////////////////////////////////////////////////////////////////////// 
void WTURBULENCE::advectTextureCoordinates (float dtOrg, float* xvel, float* yvel, float* zvel) {
  // advection
  SWAP_POINTERS(_tcTemp, _tcU);
  FLUID_3D::copyBorderX(_tcTemp, _resSm);
  FLUID_3D::copyBorderY(_tcTemp, _resSm);
  FLUID_3D::copyBorderZ(_tcTemp, _resSm);
  FLUID_3D::advectFieldMacCormack(dtOrg, xvel, yvel, zvel, 
      _tcTemp, _tcU, _tempBig1, _tempBig2, _resSm, NULL);

  SWAP_POINTERS(_tcTemp, _tcV);
  FLUID_3D::copyBorderX(_tcTemp, _resSm);
  FLUID_3D::copyBorderY(_tcTemp, _resSm);
  FLUID_3D::copyBorderZ(_tcTemp, _resSm);
  FLUID_3D::advectFieldMacCormack(dtOrg, xvel, yvel, zvel, 
      _tcTemp, _tcV, _tempBig1, _tempBig2, _resSm, NULL);

  SWAP_POINTERS(_tcTemp, _tcW);
  FLUID_3D::copyBorderX(_tcTemp, _resSm);
  FLUID_3D::copyBorderY(_tcTemp, _resSm);
  FLUID_3D::copyBorderZ(_tcTemp, _resSm);
  FLUID_3D::advectFieldMacCormack(dtOrg, xvel, yvel, zvel, 
      _tcTemp, _tcW, _tempBig1, _tempBig2, _resSm, NULL);
}

//////////////////////////////////////////////////////////////////////
// Compute the eigenvalues of the advected texture
////////////////////////////////////////////////////////////////////// 
void WTURBULENCE::computeEigenvalues() {
  // stats
  float maxeig = -1.;
  float mineig = 10.;

  // texture coordinate eigenvalues
  for (int z = 1; z < _zResSm-1; z++) {
    for (int y = 1; y < _yResSm-1; y++) 
      for (int x = 1; x < _xResSm-1; x++)
      {
        const int index = x+ y *_resSm[0] + z*_slabSizeSm;

        // compute jacobian
        float jacobian[3][3] = {
          { minDx(x, y, z, _tcU, _resSm), minDx(x, y, z, _tcV, _resSm), minDx(x, y, z, _tcW, _resSm) } ,
          { minDy(x, y, z, _tcU, _resSm), minDy(x, y, z, _tcV, _resSm), minDy(x, y, z, _tcW, _resSm) } ,
          { minDz(x, y, z, _tcU, _resSm), minDz(x, y, z, _tcV, _resSm), minDz(x, y, z, _tcW, _resSm) }
        };

        // ONLY compute the eigenvalues after checking that the matrix
        // is nonsingular
        JAMA::LU<float> LU = computeLU3x3(jacobian);

        if (LU.isNonsingular())
        {
          // get the analytic eigenvalues, quite slow right now...
          Vec3 eigenvalues = Vec3(1.);
          computeEigenvalues3x3( &eigenvalues[0], jacobian);
          _eigMax[index] = MAX3V(eigenvalues);
          _eigMin[index] = MIN3V(eigenvalues);
          maxeig = MAX(_eigMax[index],maxeig);
          mineig = MIN(_eigMin[index],mineig);
        }
        else
        {
          _eigMax[index] = 10.0f;
          _eigMin[index] = 0.1;
        }
      }
  }
}

//////////////////////////////////////////////////////////////////////
// advect & reset texture coordinates based on eigenvalues
////////////////////////////////////////////////////////////////////// 
void WTURBULENCE::resetTextureCoordinates() 
{
  // allowed deformation of the textures
  const float limit = 2.f;
  const float limitInv = 1./limit;

  // standard reset
  int resets = 0;
  const float dx = 1./(float)(_resSm[0]);
  const float dy = 1./(float)(_resSm[1]);
  const float dz = 1./(float)(_resSm[2]);

  for (int z = 1; z < _zResSm-1; z++)
    for (int y = 1; y < _yResSm-1; y++)
      for (int x = 1; x < _xResSm-1; x++)
      {
        const int index = x+ y *_resSm[0] + z*_slabSizeSm;
        if (_eigMax[index] > limit || _eigMin[index] < limitInv)
        {
          _tcU[index] = (float)x * dx;
          _tcV[index] = (float)y * dy;
          _tcW[index] = (float)z * dz;
          resets++;
        }
      }
}

//////////////////////////////////////////////////////////////////////
// Compute the highest frequency component of the wavelet
// decomposition
//////////////////////////////////////////////////////////////////////
void WTURBULENCE::decomposeEnergy()
{
  // do the decomposition -- the goal here is to have
  // the energy with the high frequency component stomped out
  // stored in _tcTemp when it is done. _highFreqEnergy is only used
  // as an additional temp array
  
  // downsample input
  downsampleXNeumann(_highFreqEnergy, _energy, _xResSm, _yResSm, _zResSm);
  downsampleYNeumann(_tcTemp, _highFreqEnergy, _xResSm, _yResSm, _zResSm);
  downsampleZNeumann(_highFreqEnergy, _tcTemp, _xResSm, _yResSm, _zResSm);

  // upsample input
  upsampleZNeumann(_tcTemp, _highFreqEnergy, _xResSm, _yResSm, _zResSm);
  upsampleYNeumann(_highFreqEnergy, _tcTemp, _xResSm, _yResSm, _zResSm);
  upsampleXNeumann(_tcTemp, _highFreqEnergy, _xResSm, _yResSm, _zResSm);

  // subtract the down and upsampled field from the original field -- 
  // what should be left over is solely the high frequency component
	int index = 0;
  for (int z = 0; z < _zResSm; z++) 
    for (int y = 0; y < _yResSm; y++) {
      for (int x = 0; x < _xResSm; x++, index++) {
        // brute force reset of boundaries
        if(z >= _zResSm - 1 || x >= _xResSm - 1 || y >= _yResSm - 1 || z <= 0 || y <= 0 || x <= 0) 
          _highFreqEnergy[index] = 0.; 
        else 
          _highFreqEnergy[index] = _energy[index] - _tcTemp[index];
    }
  }
}

//////////////////////////////////////////////////////////////////////
// compute velocity from energies and march into obstacles
// for wavelet decomposition
////////////////////////////////////////////////////////////////////// 
void WTURBULENCE::computeEnergy(float* xvel, float* yvel, float* zvel, unsigned char *obstacles) 
{
  // compute everywhere
  for (int x = 0; x < _totalCellsSm; x++) 
    _energy[x] = 0.5f * (xvel[x] * xvel[x] + yvel[x] * yvel[x] + zvel[x] * zvel[x]);

  FLUID_3D::copyBorderX(_energy, _resSm);
  FLUID_3D::copyBorderY(_energy, _resSm);
  FLUID_3D::copyBorderZ(_energy, _resSm);

  // pseudo-march the values into the obstacles
  // the wavelet upsampler only uses a 3x3 support neighborhood, so
  // propagating the values in by 4 should be sufficient
  int index;

  // iterate
  for (int iter = 0; iter < 4; iter++)
  {
    index = _slabSizeSm + _xResSm + 1;
    for (int z = 1; z < _zResSm - 1; z++, index += 2 * _xResSm)
      for (int y = 1; y < _yResSm - 1; y++, index += 2)
        for (int x = 1; x < _xResSm - 1; x++, index++)
          if (obstacles[index] && obstacles[index] != RETIRED)
          {
            float sum = 0.0f;
            int valid = 0;

            if (!obstacles[index + 1] || obstacles[index + 1] == RETIRED)
            {
              sum += _energy[index + 1];
              valid++;
            }
            if (!obstacles[index - 1] || obstacles[index - 1] == RETIRED)
            {
              sum += _energy[index - 1];
              valid++;
            }
            if (!obstacles[index + _xResSm] || obstacles[index + _xResSm] == RETIRED)
            {
              sum += _energy[index + _xResSm];
              valid++;
            }
            if (!obstacles[index - _xResSm] || obstacles[index - _xResSm] == RETIRED)
            {
              sum += _energy[index - _xResSm];
              valid++;
            }
            if (!obstacles[index + _slabSizeSm] || obstacles[index + _slabSizeSm] == RETIRED)
            {
              sum += _energy[index + _slabSizeSm];
              valid++;
            }
            if (!obstacles[index - _slabSizeSm] || obstacles[index - _slabSizeSm] == RETIRED)
            {
              sum += _energy[index - _slabSizeSm];
              valid++;
            }
            if (valid > 0)
            {
              _energy[index] = sum / valid;
              obstacles[index] = MARCHED;
            }
          }
    index = _slabSizeSm + _xResSm + 1;
    for (int z = 1; z < _zResSm - 1; z++, index += 2 * _xResSm)
      for (int y = 1; y < _yResSm - 1; y++, index += 2)
        for (int x = 1; x < _xResSm - 1; x++, index++)
          if (obstacles[index] == MARCHED)
            obstacles[index] = RETIRED;
  }
  index = _slabSizeSm + _xResSm + 1;
  for (int z = 1; z < _zResSm - 1; z++, index += 2 * _xResSm)
    for (int y = 1; y < _yResSm - 1; y++, index += 2)
      for (int x = 1; x < _xResSm - 1; x++, index++)
        if (obstacles[index])
          obstacles[index] = 1;
}

//////////////////////////////////////////////////////////////////////////////////////////
// Evaluate derivatives
//////////////////////////////////////////////////////////////////////////////////////////
Vec3 WTURBULENCE::WVelocity(Vec3 orgPos)
{
  // arbitrarily offset evaluation points
  const Vec3 p1 = orgPos + Vec3(NOISE_TILE_SIZE/2,0,0);
  const Vec3 p2 = orgPos + Vec3(0,NOISE_TILE_SIZE/2,0);
  const Vec3 p3 = orgPos + Vec3(0,0,NOISE_TILE_SIZE/2);

  const float f1y = WNoiseDy(p1, _noiseTile);
  const float f1z = WNoiseDz(p1, _noiseTile);

  const float f2x = WNoiseDx(p2, _noiseTile);
  const float f2z = WNoiseDz(p2, _noiseTile);

  const float f3x = WNoiseDx(p3, _noiseTile);
  const float f3y = WNoiseDy(p3, _noiseTile);

  Vec3 ret = Vec3( 
      f3y - f2z,
      f1z - f3x,
      f2x - f1y ); 
  return ret;
}

//////////////////////////////////////////////////////////////////////////////////////////
// Evaluate derivatives with Jacobian
//////////////////////////////////////////////////////////////////////////////////////////
Vec3 WTURBULENCE::WVelocityWithJacobian(Vec3 orgPos, float* xUnwarped, float* yUnwarped, float* zUnwarped)
{
  // arbitrarily offset evaluation points
  const Vec3 p1 = orgPos + Vec3(NOISE_TILE_SIZE/2,0,0);
  const Vec3 p2 = orgPos + Vec3(0,NOISE_TILE_SIZE/2,0);
  const Vec3 p3 = orgPos + Vec3(0,0,NOISE_TILE_SIZE/2);

  Vec3 final;
  final[0] = WNoiseDx(p1, _noiseTile);
  final[1] = WNoiseDy(p1, _noiseTile);
  final[2] = WNoiseDz(p1, _noiseTile);
  const float f1x = xUnwarped[0] * final[0] + xUnwarped[1] * final[1] + xUnwarped[2] * final[2];
  const float f1y = yUnwarped[0] * final[0] + yUnwarped[1] * final[1] + yUnwarped[2] * final[2];
  const float f1z = zUnwarped[0] * final[0] + zUnwarped[1] * final[1] + zUnwarped[2] * final[2];

  final[0] = WNoiseDx(p2, _noiseTile);
  final[1] = WNoiseDy(p2, _noiseTile);
  final[2] = WNoiseDz(p2, _noiseTile);
  const float f2x = xUnwarped[0] * final[0] + xUnwarped[1] * final[1] + xUnwarped[2] * final[2];
  const float f2y = yUnwarped[0] * final[0] + yUnwarped[1] * final[1] + yUnwarped[2] * final[2];
  const float f2z = zUnwarped[0] * final[0] + zUnwarped[1] * final[1] + zUnwarped[2] * final[2];

  final[0] = WNoiseDx(p3, _noiseTile);
  final[1] = WNoiseDy(p3, _noiseTile);
  final[2] = WNoiseDz(p3, _noiseTile);
  const float f3x = xUnwarped[0] * final[0] + xUnwarped[1] * final[1] + xUnwarped[2] * final[2];
  const float f3y = yUnwarped[0] * final[0] + yUnwarped[1] * final[1] + yUnwarped[2] * final[2];
  const float f3z = zUnwarped[0] * final[0] + zUnwarped[1] * final[1] + zUnwarped[2] * final[2];

  Vec3 ret = Vec3( 
      f3y - f2z,
      f1z - f3x,
      f2x - f1y ); 
  return ret;
}

//////////////////////////////////////////////////////////////////////
// perform an actual noise advection step
//////////////////////////////////////////////////////////////////////
void WTURBULENCE::stepTurbulenceReadable(float dtOrg, float* xvel, float* yvel, float* zvel, unsigned char *obstacles) 
{
  /*
  static int c = 0;
  if (rand()%13 == 0) c = (c + 1) % 100;
  stringstream ss;
  ss << "gabor_noise" << c << ".gabor";
  loadTile(_noiseTile, ss.str());
  std::cout << "loaded tile " << ss.str() << std::endl;
  */
  // enlarge timestep to match grid
  const float dt = dtOrg * _amplify;
  const float invAmp = 1.0f / _amplify;

  // prepare textures
  advectTextureCoordinates(dtOrg, xvel,yvel,zvel);

  // compute eigenvalues of the texture coordinates
  computeEigenvalues();

  // do wavelet decomposition of energy
  computeEnergy(xvel, yvel, zvel, obstacles);
  decomposeEnergy();

  // zero out coefficients inside of the obstacle
  for (int x = 0; x < _totalCellsSm; x++)
    if (obstacles[x]) _energy[x] = 0.f;

  float maxVelocity = 0.;
  for (int z = 1; z < _zResBig - 1; z++) 
    for (int y = 1; y < _yResBig - 1; y++) 
      for (int x = 1; x < _xResBig - 1; x++)
      {
        // get unit position for both fine and coarse grid
        const Vec3 pos = Vec3(x,y,z);
        const Vec3 posSm = pos * invAmp;
        
        // get grid index for both fine and coarse grid
        const int index = x + y *_xResBig + z *_slabSizeBig;
        const int indexSmall = (int)posSm[0] + (int)posSm[1] * _xResSm + (int)posSm[2] * _slabSizeSm;
        
        // get a linearly interpolated velocity and texcoords
        // from the coarse grid
        Vec3 vel = INTERPOLATE::lerp3dVec( xvel,yvel,zvel, 
            posSm[0], posSm[1], posSm[2], _xResSm,_yResSm,_zResSm);
        Vec3 uvw = INTERPOLATE::lerp3dVec( _tcU,_tcV,_tcW, 
            posSm[0], posSm[1], posSm[2], _xResSm,_yResSm,_zResSm);

        // multiply the texture coordinate by _resSm so that turbulence
        // synthesis begins at the first octave that the coarse grid 
        // cannot capture
        Vec3 texCoord = Vec3(uvw[0] * _resSm[0], 
                             uvw[1] * _resSm[1],
                             uvw[2] * _resSm[2]); 

        // retrieve wavelet energy at highest frequency
        float energy = INTERPOLATE::lerp3d(
            _highFreqEnergy, posSm[0],posSm[1],posSm[2], _xResSm, _yResSm, _zResSm);

        // base amplitude for octave 0
        float coefficient = sqrtf(2.0f * fabs(energy));
        const float amplitude = _strength * fabs(0.5 * coefficient) * persistence;

        // add noise to velocity, but only if the turbulence is
        // sufficiently undeformed, and the energy is large enough
        // to make a difference
        const bool addNoise = _eigMax[indexSmall] < 2. && 
                              _eigMin[indexSmall] > 0.5;
        if (addNoise && amplitude > _cullingThreshold) {
          // base amplitude for octave 0
          float amplitudeScaled = amplitude;

          for (int octave = 0; octave < _octaves; octave++)
          {
            // multiply the vector noise times the maximum allowed
            // noise amplitude at this octave, and add it to the total
            vel += WVelocity(texCoord) * amplitudeScaled;

            // scale coefficient for next octave
            amplitudeScaled *= persistence;
            texCoord *= 2.0f;
          }
        }

        // Store velocity + turbulence in big grid for maccormack step
        //
        // If you wanted to save memory, you would instead perform a 
        // semi-Lagrangian backtrace for the current grid cell here. Then
        // you could just throw the velocity away.
        _bigUx[index] = vel[0];
        _bigUy[index] = vel[1];
        _bigUz[index] = vel[2];

        // compute the velocity magnitude for substepping later
        const float velMag = _bigUx[index] * _bigUx[index] + 
                             _bigUy[index] * _bigUy[index] + 
                             _bigUz[index] * _bigUz[index];
        if (velMag > maxVelocity) maxVelocity = velMag;

        // zero out velocity inside obstacles
        float obsCheck = INTERPOLATE::lerp3dToFloat(
            obstacles, posSm[0], posSm[1], posSm[2], _xResSm, _yResSm, _zResSm); 
        if (obsCheck > 0.95) 
          _bigUx[index] = _bigUy[index] = _bigUz[index] = 0.;
      }

  // prepare density for an advection
  SWAP_POINTERS(_densityBig, _densityBigOld);

  // based on the maximum velocity present, see if we need to substep,
  // but cap the maximum number of substeps to 5
  const int maxSubSteps = 5; 
  maxVelocity = sqrt(maxVelocity) * dt;
  int totalSubsteps = (int)(maxVelocity / (float)maxSubSteps);
  totalSubsteps = (totalSubsteps < 1) ? 1 : totalSubsteps;
  totalSubsteps = (totalSubsteps > maxSubSteps) ? maxSubSteps : totalSubsteps;
  const float dtSubdiv = dt / (float)totalSubsteps;

  // set boundaries of big velocity grid
  FLUID_3D::setZeroX(_bigUx, _resBig); 
  FLUID_3D::setZeroY(_bigUy, _resBig); 
  FLUID_3D::setZeroZ(_bigUz, _resBig);

  // do the MacCormack advection, with substepping if necessary
  for(int substep = 0; substep < totalSubsteps; substep++)
  {
    FLUID_3D::advectFieldMacCormack(dtSubdiv, _bigUx, _bigUy, _bigUz, 
        _densityBigOld, _densityBig, _tempBig1, _tempBig2, _resBig, NULL);

    if (substep < totalSubsteps - 1) 
      SWAP_POINTERS(_densityBig, _densityBigOld);
  } // substep
  
  // wipe the density borders
  FLUID_3D::setZeroBorder(_densityBig, _resBig);
    
  // reset texture coordinates now in preparation for next timestep
  // Shouldn't do this before generating the noise because then the 
  // eigenvalues stored do not reflect the underlying texture coordinates
  resetTextureCoordinates();
  
  // output files
  string prefix = string("./amplified.preview/density_bigxy_");
  FLUID_3D::writeImageSliceXY(_densityBig, _resBig, _resBig[2]/2, prefix, _totalStepsBig, 1.0f);
  //string df3Prefix = string("./df3/density_big_");
  //IMAGE::dumpDF3(_totalStepsBig, df3Prefix, _densityBig, _resBig[0],_resBig[1],_resBig[2]);
  string pbrtPrefix = string("./pbrt/density_big_");
  IMAGE::dumpPBRT(_totalStepsBig, pbrtPrefix, _densityBig, _resBig[0],_resBig[1],_resBig[2]);

  _totalStepsBig++;
}

//////////////////////////////////////////////////////////////////////
// perform the full turbulence algorithm, including OpenMP 
// if available
//////////////////////////////////////////////////////////////////////
void WTURBULENCE::stepTurbulenceFull(float dtOrg, float* xvel, float* yvel, float* zvel, unsigned char *obstacles)
{
  // enlarge timestep to match grid
  const float dt = dtOrg * _amplify;
  const float invAmp = 1.0f / _amplify;

  // prepare textures
  advectTextureCoordinates(dtOrg, xvel,yvel,zvel);

  // do wavelet decomposition of energy
  computeEnergy(xvel, yvel, zvel, obstacles);
  decomposeEnergy();

  // zero out coefficients inside of the obstacle
  for (int x = 0; x < _totalCellsSm; x++)
    if (obstacles[x]) _energy[x] = 0.f;

  // parallel region setup
  float maxVelMagThreads[8] = { -1., -1., -1., -1., -1., -1., -1., -1. };
#if PARALLEL==1 && !_WIN32
#pragma omp parallel
#endif
  { float maxVelMag1 = 0.;
#if PARALLEL==1 && !_WIN32
    const int id  = omp_get_thread_num(), num = omp_get_num_threads();
#else
    const int id  = 0; 
#endif

  // vector noise main loop
#if PARALLEL==1 && !_WIN32
#pragma omp for  schedule(static)
#endif
  for (int zSmall = 0; zSmall < _zResSm; zSmall++) 
  for (int ySmall = 0; ySmall < _yResSm; ySmall++) 
  for (int xSmall = 0; xSmall < _xResSm; xSmall++)
  {
    const int indexSmall = xSmall + ySmall * _xResSm + zSmall * _slabSizeSm;

    // compute jacobian
    float jacobian[3][3] = {
      { minDx(xSmall, ySmall, zSmall, _tcU, _resSm), minDx(xSmall, ySmall, zSmall, _tcV, _resSm), minDx(xSmall, ySmall, zSmall, _tcW, _resSm) } ,
      { minDy(xSmall, ySmall, zSmall, _tcU, _resSm), minDy(xSmall, ySmall, zSmall, _tcV, _resSm), minDy(xSmall, ySmall, zSmall, _tcW, _resSm) } ,
      { minDz(xSmall, ySmall, zSmall, _tcU, _resSm), minDz(xSmall, ySmall, zSmall, _tcV, _resSm), minDz(xSmall, ySmall, zSmall, _tcW, _resSm) }
    };

    // get LU factorization of texture jacobian and apply 
    // it to unit vectors
    JAMA::LU<float> LU = computeLU3x3(jacobian);
    float xUnwarped[] = {1.0f, 0.0f, 0.0f};
    float yUnwarped[] = {0.0f, 1.0f, 0.0f};
    float zUnwarped[] = {0.0f, 0.0f, 1.0f};
    float xWarped[] = {1.0f, 0.0f, 0.0f};
    float yWarped[] = {0.0f, 1.0f, 0.0f};
    float zWarped[] = {0.0f, 0.0f, 1.0f};
    bool nonSingular = LU.isNonsingular();
    float eigMax = 10.0f;
    float eigMin = 0.1f;
    if (nonSingular)
    {
      solveLU3x3(LU, xUnwarped, xWarped);
      solveLU3x3(LU, yUnwarped, yWarped);
      solveLU3x3(LU, zUnwarped, zWarped);

      // compute the eigenvalues while we have the Jacobian available
      Vec3 eigenvalues = Vec3(1.);
      computeEigenvalues3x3( &eigenvalues[0], jacobian);
      _eigMax[indexSmall] = MAX3V(eigenvalues);
      _eigMin[indexSmall] = MIN3V(eigenvalues);
    }
    
    // make sure to skip one on the beginning and end
    int xStart = (xSmall == 0) ? 1 : 0;
    int xEnd   = (xSmall == _xResSm - 1) ? _amplify - 1 : _amplify;
    int yStart = (ySmall == 0) ? 1 : 0;
    int yEnd   = (ySmall == _yResSm - 1) ? _amplify - 1 : _amplify;
    int zStart = (zSmall == 0) ? 1 : 0;
    int zEnd   = (zSmall == _zResSm - 1) ? _amplify - 1 : _amplify;
      
    for (int zBig = zStart; zBig < zEnd; zBig++) 
    for (int yBig = yStart; yBig < yEnd; yBig++) 
    for (int xBig = xStart; xBig < xEnd; xBig++)
    {
      const int x = xSmall * _amplify + xBig;
      const int y = ySmall * _amplify + yBig;
      const int z = zSmall * _amplify + zBig;
      
      // get unit position for both fine and coarse grid
      const Vec3 pos = Vec3(x,y,z);
      const Vec3 posSm = pos * invAmp;
      
      // get grid index for both fine and coarse grid
      const int index = x + y *_xResBig + z *_slabSizeBig;
      
      // get a linearly interpolated velocity and texcoords
      // from the coarse grid
      Vec3 vel = INTERPOLATE::lerp3dVec( xvel,yvel,zvel, 
          posSm[0], posSm[1], posSm[2], _xResSm,_yResSm,_zResSm);
      Vec3 uvw = INTERPOLATE::lerp3dVec( _tcU,_tcV,_tcW, 
          posSm[0], posSm[1], posSm[2], _xResSm,_yResSm,_zResSm);

      // multiply the texture coordinate by _resSm so that turbulence
      // synthesis begins at the first octave that the coarse grid 
      // cannot capture
      Vec3 texCoord = Vec3(uvw[0] * _resSm[0], 
                           uvw[1] * _resSm[1],
                           uvw[2] * _resSm[2]); 

      // retrieve wavelet energy at highest frequency
      float energy = INTERPOLATE::lerp3d(
          _highFreqEnergy, posSm[0],posSm[1],posSm[2], _xResSm, _yResSm, _zResSm);

      // base amplitude for octave 0
      float coefficient = sqrtf(2.0f * fabs(energy));
      const float amplitude = _strength * fabs(0.5 * coefficient) * persistence;

      // add noise to velocity, but only if the turbulence is
      // sufficiently undeformed, and the energy is large enough
      // to make a difference
      const bool addNoise = _eigMax[indexSmall] < 2. && 
                            _eigMin[indexSmall] > 0.5;
      if (addNoise && amplitude > _cullingThreshold) {
        // base amplitude for octave 0
        float amplitudeScaled = amplitude;

        for (int octave = 0; octave < _octaves; octave++)
        {
          // multiply the vector noise times the maximum allowed
          // noise amplitude at this octave, and add it to the total
          vel += WVelocityWithJacobian(texCoord, &xUnwarped[0], &yUnwarped[0], &zUnwarped[0]) * amplitudeScaled;

          // scale coefficient for next octave
          amplitudeScaled *= persistence;
          texCoord *= 2.0f;
        }
      }

      // Store velocity + turbulence in big grid for maccormack step
      //
      // If you wanted to save memory, you would instead perform a 
      // semi-Lagrangian backtrace for the current grid cell here. Then
      // you could just throw the velocity away.
      _bigUx[index] = vel[0];
      _bigUy[index] = vel[1];
      _bigUz[index] = vel[2];

      // compute the velocity magnitude for substepping later
      const float velMag = _bigUx[index] * _bigUx[index] + 
                           _bigUy[index] * _bigUy[index] + 
                           _bigUz[index] * _bigUz[index];
      if (velMag > maxVelMag1) maxVelMag1 = velMag;

      // zero out velocity inside obstacles
      float obsCheck = INTERPOLATE::lerp3dToFloat(
          obstacles, posSm[0], posSm[1], posSm[2], _xResSm, _yResSm, _zResSm); 
      if (obsCheck > 0.95) 
        _bigUx[index] = _bigUy[index] = _bigUz[index] = 0.;
    } // xyz

#if PARALLEL==1 && !_WIN32
    maxVelMagThreads[id] = maxVelMag1;
#else
    maxVelMagThreads[0] = maxVelMag1;
#endif
  }
  } // omp
  
  // compute maximum over threads
  float maxVelMag = maxVelMagThreads[0];
#if PARALLEL==1 && !_WIN32
  for (int i = 1; i < 8; i++) 
    if (maxVelMag < maxVelMagThreads[i]) 
      maxVelMag = maxVelMagThreads[i];
#endif

  // prepare density for an advection
  SWAP_POINTERS(_densityBig, _densityBigOld);

  // based on the maximum velocity present, see if we need to substep,
  // but cap the maximum number of substeps to 5
  const int maxSubSteps = 5; 
  maxVelMag = sqrt(maxVelMag) * dt;
  int totalSubsteps = (int)(maxVelMag / (float)maxSubSteps);
  totalSubsteps = (totalSubsteps < 1) ? 1 : totalSubsteps;
  totalSubsteps = (totalSubsteps > maxSubSteps) ? maxSubSteps : totalSubsteps;
  const float dtSubdiv = dt / (float)totalSubsteps;

  // set boundaries of big velocity grid
  FLUID_3D::setZeroX(_bigUx, _resBig); 
  FLUID_3D::setZeroY(_bigUy, _resBig); 
  FLUID_3D::setZeroZ(_bigUz, _resBig);

  // do the MacCormack advection, with substepping if necessary
  for(int substep = 0; substep < totalSubsteps; substep++)
  {
    FLUID_3D::advectFieldMacCormack(dtSubdiv, _bigUx, _bigUy, _bigUz, 
        _densityBigOld, _densityBig, _tempBig1, _tempBig2, _resBig, NULL);

    if (substep < totalSubsteps - 1) 
      SWAP_POINTERS(_densityBig, _densityBigOld);
  } // substep
  
  // wipe the density borders
  FLUID_3D::setZeroBorder(_densityBig, _resBig);
    
  // reset texture coordinates now in preparation for next timestep
  // Shouldn't do this before generating the noise because then the 
  // eigenvalues stored do not reflect the underlying texture coordinates
  resetTextureCoordinates();
  
  // output files
  string prefix = string("./amplified.preview/density_bigxy_");
  FLUID_3D::writeImageSliceXY(_densityBig, _resBig, _resBig[2]/2, prefix, _totalStepsBig, 1.0f);
  //string df3prefix = string("./df3/density_big_");
  //IMAGE::dumpDF3(_totalStepsBig, df3prefix, _densityBig, _resBig[0],_resBig[1],_resBig[2]);
  string pbrtPrefix = string("./pbrt/density_big_");
  IMAGE::dumpPBRT(_totalStepsBig, pbrtPrefix, _densityBig, _resBig[0],_resBig[1],_resBig[2]);
  
  _totalStepsBig++;
}
