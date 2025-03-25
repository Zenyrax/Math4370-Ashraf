#include <vector>
#include <random>
#include <iostream>

// Calculate the magnitude of a vector
double magnitude(const std::vector<double>& point) {
  double sum = 0.0;
  for (double v : point) {
    sum += v * v;
  }
  return std::sqrt(sum);
}

double montecarlo(double (*f)(std::vector<double>), int (*h)(std::vector<double>, double(std::vector<double>)), std::vector<std::pair<double, double>> hypercube, double M, double real_value = 0) {
  // Get dimension of the function by the size of the hypercube 
  int dimensions = hypercube.size();

  // Calculate the area of the hypercube
  double hypercube_area = 1;
  for (size_t i = 0; i < dimensions; i++) {
    hypercube_area *= hypercube[i].second - hypercube[i].first;
  }

  // Seed the random number generator
  std::random_device rd;
  std::mt19937 gen(rd());

  // Generate random points and 
  int N = 0;
  double third = 0;
  std::vector<double> random_point(dimensions);
  if (real_value == 0) {
    for (size_t i = 0; i < M; i++) {
      // Generate a random point
      for (size_t j = 0; j < dimensions; j++) {
        std::uniform_real_distribution<> dis(hypercube[j].first, hypercube[j].second);
        random_point[j] = dis(gen);
      }

      //Is that random point in our function?
      // if(h(random_point, f) == 1) {
      //   N++;
      //   third += f(random_point);
      // }
      third += f(random_point) * h(random_point, f);
    }
  } else {
    double threshold = M;
    M = 1;
    for (size_t j = 0; j < dimensions; j++) {
      std::uniform_real_distribution<> dis(hypercube[j].first, hypercube[j].second);
      random_point[j] = dis(gen);
    }
    third += f(random_point) * h(random_point, f);
    while (std::fabs(real_value - (hypercube_area * third / M)) > threshold) {
      for (size_t j = 0; j < dimensions; j++) {
        std::uniform_real_distribution<> dis(hypercube[j].first, hypercube[j].second);
        random_point[j] = dis(gen);
      }

      third += f(random_point) * h(random_point, f);
      M++;
    }
  }
  
  std::cout << "hypercube area: " << hypercube_area << " | M: " << M << std::endl;
  return hypercube_area * third / M;
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
  std::cout << montecarlo(f_pi, h_pi, hypercube_pi, 1000000) << std::endl;
  std::cout << montecarlo(f_pi, h_pi, hypercube_pi, .0001, 3.14159265359) << std::endl;
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
  std::cout << montecarlo(f_quad, h_quad2, hypercube_quad, 1000000) << std::endl;
  std::cout << montecarlo(f_quad, h_quad2, hypercube_quad, .0001, 2.66666666) << std::endl;
  return 0;
}