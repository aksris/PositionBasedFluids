#ifndef MACGRIDDATA_H
#define MACGRIDDATA_H
#include <glm/glm.hpp>
#include <vector>
#include <algorithm>
using namespace glm;

class MACGridData
{
public:
    MACGridData();
    MACGridData(const MACGridData &val);

    virtual MACGridData& operator=(const MACGridData& val);
    virtual void MACGridDataInitialize();

    virtual float& operator()(int i, int j, int k);

    virtual void getCell(const vec3& pt, int& i, int& j, int& k);
    int getCellMark(int i, int j, int k);
    virtual void setCell(int& i, int& j, int& k, const float val);
    void setCellMark(int &i, int &j, int &k, const int val, bool mark);
    int getCellIndex(int i, int j, int k);
    virtual vec3 worldToLocal(const vec3& pt) const;
    virtual float interpolate(const vec3& pt);
    vec3 containerBounds;
    std::vector<float> data;
    std::vector<int> mData;
protected:
};

class MACGridDataX : public MACGridData
{
public:
   MACGridDataX();

   virtual void MACGridDataInitialize();
   virtual float& operator()(int i, int j, int k);
   virtual void setCell(int& i, int& j, int& k, const float val);
   virtual vec3 worldToLocal(const vec3& pt) const;
};

class MACGridDataY : public MACGridData
{
public:
   MACGridDataY();

   virtual void MACGridDataInitialize();
   virtual float& operator()(int i, int j, int k);
   virtual void setCell(int& i, int& j, int& k, const float val);
   virtual vec3 worldToLocal(const vec3& pt) const;
};

class MACGridDataZ : public MACGridData
{
public:
   MACGridDataZ();

   virtual void MACGridDataInitialize();
   virtual float& operator()(int i, int j, int k);
   virtual void setCell(int& i, int& j, int& k, const float val);
   virtual vec3 worldToLocal(const vec3& pt) const;
};


#endif // MACGRIDDATA_H
