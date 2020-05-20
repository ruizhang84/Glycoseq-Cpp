#ifndef MODEL_SPECTRUM_SPECTRUM_PEAK_H_
#define MODEL_SPECTRUM_SPECTRUM_PEAK_H_

namespace model {
namespace spectrum {

class Peak
{
public:
    Peak(double mz, double intensity):
        mz_(mz), intensity_(intensity){}

    double MZ() { return mz_; }
    void set_mz(double mz) { mz_ = mz; }
    double Intensity() { return intensity_; }
    void set_intensity(double intensity) 
        { intensity_ = intensity; }

protected:
    double mz_;
    double intensity_;
};

} // namespace spectrum
} // namespace model


#endif