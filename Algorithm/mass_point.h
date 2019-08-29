#ifndef MASS_POINT_H
#define MASS_POINT_H

#include "point.h"
class MassPoint : public Point{
private:
    double value;
public:
    MassPoint() = default;
    MassPoint(double data) :value(data){};
    double get_value() const override;
};

#endif