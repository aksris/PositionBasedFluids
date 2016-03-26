#include "macgriddata.h"

float CellSize = 1.f;
int gDimension[3] = {5, 5, 5};
float LERP(float val1, float val2, float t){
    return (1 - t) * val1 + t * val2;
}

MACGridData::MACGridData():
containerBounds(0.f, 0.f, 0.f)
{
}

MACGridData::MACGridData(const MACGridData &val){

    containerBounds = val.containerBounds;
}

MACGridData& MACGridData::operator =(const MACGridData& val){
    if (this == &val){
        return *this;
    }

    data = val.data;
    containerBounds = val.containerBounds;
    return *this;
}

void MACGridData::MACGridDataInitialize(){
    containerBounds[0] = CellSize * gDimension[0];
    containerBounds[1] = CellSize * gDimension[1];
    containerBounds[2] = CellSize * gDimension[2];
    data.resize(gDimension[0]*gDimension[1]*gDimension[2]);
    std::fill(data.begin(), data.end(), 0.f);
}

float& MACGridData::operator ()(int i, int j, int k){

    float ret = 0.f;
    if (i < 0 || j < 0 || k < 0 ||
            i > gDimension[0]-1 ||
            j > gDimension[1]-1 ||
            k > gDimension[2]-1) return ret;

    int x = i;
    int y = k * gDimension[0];
    int z = j * gDimension[0] * gDimension[2];

    return data[x+y+z];
}

void MACGridData::getCell(const vec3 &pt, int &i, int &j, int &k){
    vec3 posL = worldToLocal(pt);

    i = (int) (posL[0]/CellSize);
    j = (int) (posL[1]/CellSize);
    k = (int) (posL[2]/CellSize);
}

int MACGridData::getCellMark(int i, int j, int k){
    float ret = 0.f;
    if (i < 0 || j < 0 || k < 0 ||
            i > gDimension[0]-1 ||
            j > gDimension[1]-1 ||
            k > gDimension[2]-1) return ret;

    int x = i;
    int y = k * gDimension[0];
    int z = j * gDimension[0] * gDimension[2];

    return mData[x+y+z];
}

int MACGridData::getCellIndex(int i, int j, int k){

    int x = i;
    int y = k * gDimension[0];
    int z = j * gDimension[0] * gDimension[2];

    return x+y+z;
}

void MACGridData::setCell(int &i, int &j, int &k, const float val){
    int x = i;
    int y = k * gDimension[0];
    int z = j * gDimension[0] * gDimension[2];
    data[x+y+z] = val;
}

void MACGridData::setCellMark(int &i, int &j, int &k, const int val, bool mark){
    int x = i;
    int y = k * gDimension[0];
    int z = j * gDimension[0] * gDimension[2];
    mData[x+y+z] += val;
}

vec3 MACGridData::worldToLocal(const vec3& pt) const
{
   vec3 ret;
   ret[0] = min(max(0.0f, pt[0] - CellSize*0.5f), containerBounds[0]); //staggering
   ret[1] = min(max(0.0f, pt[1] - CellSize*0.5f), containerBounds[1]);
   ret[2] = min(max(0.0f, pt[2] - CellSize*0.5f), containerBounds[2]);
   return ret;
}

float MACGridData::interpolate(const vec3& pt)
{
   vec3 pos = worldToLocal(pt);

   int i = (int) (pos[0]/CellSize);
   int j = (int) (pos[1]/CellSize);
   int k = (int) (pos[2]/CellSize);

   float scale = 1.0 / CellSize;
   float fract_partx = scale * (pos[0] - i*CellSize);
   float fract_party = scale * (pos[1] - j*CellSize);
   float fract_partz = scale * (pos[2] - k*CellSize);

   float v000 = (*this)(i,j,k);
   float v010 = (*this)(i,j+1,k);
   float v100 = (*this)(i+1,j,k);
   float v110 = (*this)(i+1,j+1,k);
   float v001 = (*this)(i,j,k+1);
   float v011 = (*this)(i,j+1,k+1);
   float v101 = (*this)(i+1,j,k+1);
   float v111 = (*this)(i+1,j+1,k+1);

   float lerp1 = LERP(v000, v010, fract_party);
   float lerp2 = LERP(v100, v110, fract_party);

   float lerp3 = LERP(v001, v011, fract_party);
   float lerp4 = LERP(v101, v111, fract_party);

   float fLerp1 = LERP (lerp1, lerp2, fract_partx);
   float fLerp2 = LERP (lerp3, lerp4, fract_partx);

   float ret = LERP(fLerp1, fLerp2, fract_partz);

   return ret;
}

MACGridDataX::MACGridDataX()
{
}

void MACGridDataX::MACGridDataInitialize(){
    MACGridData::MACGridDataInitialize();
    containerBounds[0] = CellSize * (gDimension[0] + 1);
    containerBounds[1] = CellSize * gDimension[1];
    containerBounds[2] = CellSize * gDimension[2];
    data.resize((gDimension[0] + 1) * gDimension[1] * gDimension[2]);
    std::fill(data.begin(), data.end(), 0.f);
}

float& MACGridDataX::operator ()(int i, int j, int k){

    float ret = 0.f;
    if (i < 0 || i > gDimension[0]) return ret;
    if (j < 0) {
        j = 0;
    }
    if (k < 0){
        k = 0;
    }
    if (j > gDimension[1]){
        j = gDimension[1] - 1;
    }
    if (k > gDimension[2]){
        k = gDimension[2] - 1;
    }

    int x = i;
    int y = k * (gDimension[0] + 1);
    int z = j * (gDimension[0] + 1) * gDimension[2];

    return data[x+y+z];
}

vec3 MACGridDataX::worldToLocal(const vec3 &pt) const {
    vec3 ret;
    ret[0] = min(max(0.0f, pt[0]), containerBounds[0]); //staggering // its already at 0 face of x axis
    ret[1] = min(max(0.0f, pt[1] - CellSize*0.5f), containerBounds[1]);
    ret[2] = min(max(0.0f, pt[2] - CellSize*0.5f), containerBounds[2]);
    return ret;
}

void MACGridDataX::setCell(int &i, int &j, int &k, const float val){
    int x = i;
    int y = k * (gDimension[0] + 1);
    int z = j * (gDimension[0] + 1) * gDimension[2];
    data[x+y+z] = val;
}

MACGridDataY::MACGridDataY()
{
}

void MACGridDataY::MACGridDataInitialize(){
    MACGridData::MACGridDataInitialize();
    containerBounds[0] = CellSize * (gDimension[0]);
    containerBounds[1] = CellSize * (gDimension[1] + 1);
    containerBounds[2] = CellSize * gDimension[2];
    data.resize((gDimension[0]) * (gDimension[1] + 1) * gDimension[2]);
    std::fill(data.begin(), data.end(), 0.f);
}

float& MACGridDataY::operator ()(int i, int j, int k){

    float ret = 0.f;
    if (j < 0 || j > gDimension[1]) return ret;
    if (i < 0) {
        i = 0;
    }
    if (k < 0){
        k = 0;
    }
    if (i > gDimension[0]){
        i = gDimension[0] - 1;
    }
    if (k > gDimension[2]){
        k = gDimension[2] - 1;
    }

    int x = i;
    int y = k * (gDimension[0]);
    int z = j * (gDimension[0]) * gDimension[2];

    return data[x+y+z];
}

vec3 MACGridDataY::worldToLocal(const vec3 &pt) const {
    vec3 ret;
    ret[0] = min(max(0.0f, pt[0]  - CellSize*0.5f), containerBounds[0]);
    ret[1] = min(max(0.0f, pt[1]), containerBounds[1]); //staggering // its already at 0 face of y axis
    ret[2] = min(max(0.0f, pt[2] - CellSize*0.5f), containerBounds[2]);
    return ret;
}

void MACGridDataY::setCell(int &i, int &j, int &k, const float val){
    int x = i;
    int y = k * (gDimension[0]);
    int z = j * (gDimension[0]) * gDimension[2];
    data[x+y+z] = val;
}

MACGridDataZ::MACGridDataZ()
{
}

void MACGridDataZ::MACGridDataInitialize(){
    MACGridData::MACGridDataInitialize();
    containerBounds[0] = CellSize * (gDimension[0]);
    containerBounds[1] = CellSize * gDimension[1];
    containerBounds[2] = CellSize * (gDimension[2] + 1);
    data.resize((gDimension[0]) * gDimension[1] * (gDimension[2] + 1));
    std::fill(data.begin(), data.end(), 0.f);
}

float& MACGridDataZ::operator ()(int i, int j, int k){
    float ret = 0.f;
    if (k < 0 || k > gDimension[2]) return ret;
    if (j < 0) {
        j = 0;
    }
    if (i < 0){
        i = 0;
    }
    if (j > gDimension[1]){
        j = gDimension[1] - 1;
    }
    if (i > gDimension[0]){
        i = gDimension[0] - 1;
    }

    int x = i;
    int y = k * (gDimension[0]);
    int z = j * (gDimension[0]) * (gDimension[2] + 1);

    return data[x+y+z];
}

vec3 MACGridDataZ::worldToLocal(const vec3 &pt) const {
    vec3 ret;
    ret[0] = min(max(0.0f, pt[0] - CellSize*0.5f), containerBounds[0]);
    ret[1] = min(max(0.0f, pt[1] - CellSize*0.5f), containerBounds[1]);
    ret[2] = min(max(0.0f, pt[2]), containerBounds[2]);//staggering // its already at 0 face of z axis
    return ret;
}

void MACGridDataZ::setCell(int &i, int &j, int &k, const float val){
    int x = i;
    int y = k * (gDimension[0]);
    int z = j * (gDimension[0]) * (gDimension[2] + 1);
    data[x+y+z] = val;
}
