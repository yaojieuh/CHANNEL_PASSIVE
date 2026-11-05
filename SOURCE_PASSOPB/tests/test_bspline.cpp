#include <cmath>
#include <cstdio>
#include <vector>
extern "C" void LABS_fill_full(double*,double*,double*,double*,double*,double*,double*,double*,double*,int*,int*,double*,double*);
int main() {
  for (double stretch : {0.0, 0.99}) {
    int k=7,n=32;
    double ly=2;
    std::vector<double> m[5];
    for (auto &a:m) a.resize(n*n);
    std::vector<double> y(n),knots(n),integral(n),center(n);
    LABS_fill_full(m[0].data(),m[1].data(),m[2].data(),m[3].data(),m[4].data(),y.data(),knots.data(),integral.data(),center.data(),&k,&n,&stretch,&ly);
    double sum=0;
    for (double v:integral) { if(!std::isfinite(v)) return 1; sum+=v; }
    if(std::abs(sum-2)>1e-12) return 2;
    for(int i=0;i<n;++i) for(int d=0;d<5;++d) {
      double value=0,scale=1;
      for(int j=0;j<n;++j) {
        double v=m[d][j+i*n];
        if(!std::isfinite(v)) return 3;
        value+=v; scale+=std::abs(v);
      }
      if(std::abs(value-(d==0?1:0))>1e-11*scale) return 4;
    }
  }
  puts("PASS: B-spline partition, derivative sums, integration, uniform and stretched grids");
}
