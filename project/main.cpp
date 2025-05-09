#include <vector>
#include <random>
#include <iostream>
#include <omp.h>

// Calculate the magnitude of a vector
double magnitude(const std::vector<double>& point) {
  double sum = 0.0;
  for (double v : point) {
    sum += v * v;
  }
  return std::sqrt(sum);
}


double montecarlo(double* elapsed, double (*f)(std::vector<double>), int (*h)(std::vector<double>, double(std::vector<double>)), std::vector<std::pair<double, double>> hypercube, int M, double delta = 0) {
  double stime = omp_get_wtime();
  // Get dimension of the function by the size of the hypercube 
  int dimensions = hypercube.size();

  // Calculate the area of the hypercube
  double hypercube_area = 1;
  for (size_t i = 0; i < dimensions; i++) {
    hypercube_area *= hypercube[i].second - hypercube[i].first;
  }

  std::random_device rd;

  double third = 0, lastRes = 0;
  size_t i, stops = 0, iters = 0;
  int thread_count;

  if (delta == 0) {
    iters = M;
  }


  // Generate random points and 
  bool broken = false;
  if (delta == 0) {
    #pragma omp parallel
    {
      std::vector<double> random_point(dimensions);

      // Seed the random number generator
      thread_count = omp_get_num_threads();
      int tid = omp_get_thread_num();
      std::mt19937 gen(rd() + tid); // unique seed per thread

      #pragma omp for reduction(+:third) schedule(static)
      for (i = 0; i < M; i++) {
        // Generate a random point
        for (size_t j = 0; j < dimensions; j++) {
          std::uniform_real_distribution<> dis(hypercube[j].first, hypercube[j].second);
          random_point[j] = dis(gen);
        }

        // Value of f is added this sum if it passes the indicator function
        third += f(random_point) * h(random_point, f);
      }
    }
  } else {
    double last;
    while (iters < M) {
      #pragma omp parallel
      {
        std::vector<double> random_point(dimensions);

        // Seed the random number generator
        thread_count = omp_get_num_threads();
        int tid = omp_get_thread_num();
        std::mt19937 gen(rd() + tid); // unique seed per thread

        #pragma omp for reduction(+:third) schedule(static)
        for (i = 0; i < 10'000; i++) {
          // Generate a random point
          for (size_t j = 0; j < dimensions; j++) {
            std::uniform_real_distribution<> dis(hypercube[j].first, hypercube[j].second);
            random_point[j] = dis(gen);
          }

          // Value of f is added this sum if it passes the indicator function
          third += f(random_point) * h(random_point, f);
        }
      }
      iters += 10'000;
      if (iters != 10'000) {
        if (fabs(last - (hypercube_area * third / iters)) < delta) {
          break;
        }
      }
      last = hypercube_area * third / iters;
    }
  }

  double ftime = omp_get_wtime();
  *elapsed = ftime - stime;
  
  std::cout << "hypercube area: " << hypercube_area << " | iterations: " << iters << std::endl;
  return hypercube_area * third / iters;
}

int main(int argc, char* argv[]) {
  auto f_pi = [](std::vector<double> point) -> double {
    return 1;
  };
  auto h_pi = [](std::vector<double> point, double (*f)(std::vector<double>)) -> int {
    if (magnitude(point) <= 1) {
      return 1;
    }
    return 0;
  };
  std::vector<std::pair<double, double>> hypercube_pi = {{-1.0, 1.0}, {-1.0, 1.0}};
  double elapsed;
  std::cout << "Approximating pi:" << std::endl;
  double pi_res = montecarlo(&elapsed, f_pi, h_pi, hypercube_pi, 10'000'000);
  std::cout << "result: " << pi_res <<  " | error: " << std::fabs( 3.14159265359 - pi_res) << " | time: " << elapsed << "s" << std::endl;
  // double pi_res_del = montecarlo(&elapsed, f_pi, h_pi, hypercube_pi, 10'000'000, .0000001);
  // std::cout << "result: " << pi_res_del <<  " | error: " << std::fabs( 3.14159265359 - pi_res_del) << " | time: " << elapsed << "s" << std::endl;
  
  // f(x) = x^2
  auto f_quad = [](std::vector<double> point) -> double {
    return point[0] * point[0];
  };

  // Given a random point, (x, y), if y <= f(x), return 1, else, return 0
  auto h_quad = [](std::vector<double> point, double (*f)(std::vector<double>)) -> int {
    if (point[1] <= f(point)) {
      return 1;
    }
    return 0;
  };
  auto h_quad2 = [](std::vector<double> point, double (*f)(std::vector<double>)) -> int {
    return 1;
  };
  // Definite integral from 0 to 2
  std::vector<std::pair<double, double>> hypercube_quad = {{0, 2}};
  std::cout << "\nApproximating integral of x^2 between 0 and 2:" << std::endl;
  double quad_res = montecarlo(&elapsed, f_quad, h_quad2, hypercube_quad, 10'000'000);
  std::cout << "result: " << quad_res <<  " | error: " << std::fabs(2.66666666 - quad_res)  << " | time: " << elapsed << "s" << std::endl;
  // double quad_res_del = montecarlo(&elapsed, f_quad, h_quad2, hypercube_quad, 10'000'000, .000001);
  // std::cout << "result: " << quad_res_del <<  " | error: " << std::fabs(2.66666666 - quad_res_del)  << " | time: " << elapsed << "s" << std::endl;

  // indicator function w radius 2
  auto h_sphere = [](std::vector<double> point, double (*f)(std::vector<double>)) -> int {
    if (magnitude(point) <= 2) {
      return f(point);
    }
    return 0;
  };

  // 3d hypersphere w/ radius 2
  std::vector<std::pair<double, double>> hypercube_hs_3d_r2;
  for (size_t i = 0; i < 3; i++) {
    hypercube_hs_3d_r2.push_back({-2, 2});
  }
  std::cout << "\nApproximating volume of a sphere w/ radius 2:" << std::endl;
  double sphere_res = montecarlo(&elapsed, f_pi, h_sphere, hypercube_hs_3d_r2, 10'000'000);
  std::cout << "result: " << sphere_res <<  " | error: " << std::fabs(33.5103216383 - sphere_res)  << " | time: " << elapsed << "s"  << std::endl;
  // double sphere_res_del = montecarlo(&elapsed, f_pi, h_sphere, hypercube_hs_3d_r2, 10'000'000, .0001);
  // std::cout << "result: " << sphere_res_del <<  " | error: " << std::fabs(33.5103216383 - sphere_res_del)  << " | time: " << elapsed << "s"  << std::endl;

  
  // 4d hypersphere w/ radius 2
  std::vector<std::pair<double, double>> hypercube_hs_4d_r2;
  for (size_t i = 0; i < 4; i++) {
    hypercube_hs_4d_r2.push_back({-2, 2});
  }
  std::cout << "\nApproximating volume of a 4d hypersphere w/ radius 2:" << std::endl;
  double sphere_4d_res = montecarlo(&elapsed, f_pi, h_sphere, hypercube_hs_4d_r2, 10'000'000);
  std::cout << "result: " << sphere_4d_res <<  " | error: " << std::fabs(78.9568352087 - sphere_4d_res)  << " | time: " << elapsed << "s"  << std::endl;
  // double sphere_4d_res_del = montecarlo(&elapsed, f_pi, h_sphere, hypercube_hs_4d_r2, 10'000'000, .0001);
  // std::cout << "result: " << sphere_4d_res_del <<  " | error: " << std::fabs(78.9568352087 - sphere_4d_res_del)  << " | time: " << elapsed << "s"  << std::endl;

  // 5d hypersphere w/ radius 2
  std::vector<std::pair<double, double>> hypercube_hs_5d_r2;
  for (size_t i = 0; i < 5; i++) {
    hypercube_hs_5d_r2.push_back({-2, 2});
  }
  std::cout << "\nApproximating volume of a 5d hypersphere w/ radius 2:" << std::endl;
  double sphere_5d_res = montecarlo(&elapsed, f_pi, h_sphere, hypercube_hs_5d_r2, 10'000'000);
  std::cout << "result: " << sphere_5d_res <<  " | error: " << std::fabs(168.441248445 - sphere_5d_res)  << " | time: " << elapsed << "s"  << std::endl;
  // double sphere_5d_res_del = montecarlo(&elapsed, f_pi, h_sphere, hypercube_hs_5d_r2, 10'000'000, .0001);
  // std::cout << "result: " << sphere_5d_res_del <<  " | error: " << std::fabs(168.441248445 - sphere_5d_res_del)  << " | time: " << elapsed << "s"  << std::endl;
  
  // 10d hypersphere w/ radius 2
  std::vector<std::pair<double, double>> hypercube_hs_10d_r2;
  for (size_t i = 0; i < 10; i++) {
    hypercube_hs_10d_r2.push_back({-2, 2});
  }
  std::cout << "\nApproximating volume of a 10d hypersphere w/ radius 2:" << std::endl;
  double sphere_10d_res = montecarlo(&elapsed, f_pi, h_sphere, hypercube_hs_10d_r2, 10'000'000);
  std::cout << "result: " << sphere_10d_res <<  " | error: " << std::fabs(2611.36797683 - sphere_10d_res)  << " | time: " << elapsed << "s"  << std::endl;
  // double sphere_10d_res_del = montecarlo(&elapsed, f_pi, h_sphere, hypercube_hs_10d_r2, 10'000'000, .0001);
  // std::cout << "result: " << sphere_10d_res_del <<  " | error: " << std::fabs(2611.36797683 - sphere_10d_res_del)  << " | time: " << elapsed << "s"  << std::endl;
  
  return 0;
}