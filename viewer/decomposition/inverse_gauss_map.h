#ifndef INVERSE_GAUSS_MAP_H
#define INVERSE_GAUSS_MAP_H


class InverseGaussMap {

  private:
  public:
    int dimension;
    int resolution;

    int getBin2D(int normal);
    int getBin3D(int normal);
    int getBinFromAngle(float angle, int k);


  public:
    InverseGaussMap(int dimension, int resolution);

    void setResolution(int resolution);
    // TODO - VectorND? 2d or 3d...?
    int getBin(int normal);

};


#endif
