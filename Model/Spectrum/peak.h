#ifndef MODEL_SPECTRUM_PEAK_H_
#define MODEL_SPECTRUM_PEAK_H_

class Peak{
public:
    Peak(double mz, double intensity);
    Peak(const Peak&);
    Peak& operator=(const Peak&);
    double get_mz();
    void set_mz(double mz);
    double get_intensity();
    void set_intensity(double intensity);
private:
    double mz;
    double intensity;
};

#endif