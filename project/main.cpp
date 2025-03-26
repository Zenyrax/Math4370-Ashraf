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

double montecarlo(double (*f)(std::vector<double>), int (*h)(std::vector<double>, double(std::vector<double>)), std::vector<std::pair<double, double>> hypercube, int M, double real_value = 0, double threshold = 0) {
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
  size_t i;
  for (i = 0; i < M; i++) {
    // Generate a random point
    for (size_t j = 0; j < dimensions; j++) {
      std::uniform_real_distribution<> dis(hypercube[j].first, hypercube[j].second);
      random_point[j] = dis(gen);
    }

    // Value of f is added this sum if it passes the indicator function
    third += f(random_point) * h(random_point, f);

    // If a error threshold has been provided, break early if it's met
    if (threshold != 0 && std::fabs(real_value - (hypercube_area * third / i)) < threshold) {
      break;
    }
  }
  
  std::cout << "hypercube area: " << hypercube_area << " | iterations: " << i << std::endl;
  return hypercube_area * third / i;
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
  std::cout << montecarlo(f_pi, h_pi, hypercube_pi, 1000000, 3.14159265359, .0001) << std::endl;
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
  std::cout << montecarlo(f_quad, h_quad2, hypercube_quad, 1'000'000) << std::endl;
  std::cout << montecarlo(f_quad, h_quad2, hypercube_quad, 1'000'000, 2.66666666, .0001) << std::endl;

  // Volume of hyperspheres
  auto f_sphere_r1 = [](std::vector<double> point) -> double {
    return 1;
  };
  auto f_sphere_r2 = [](std::vector<double> point) -> double {
    return 4;
  };
  auto f_sphere_r3 = [](std::vector<double> point) -> double {
    return 9;
  };
  auto h_sphere = [](std::vector<double> point, double (*f)(std::vector<double>)) -> int {
    if (magnitude(point) <= 1) {
      return f(point);
    }
    return 0;
  };
  
  // 4d hypersphere w/ radius 2
  std::vector<std::pair<double, double>> hypercube_hs_4d_r2;
  for (size_t i = 0; i < 4; i++) {
    hypercube_hs_4d_r2.push_back({-2, 2});
  }
  std::cout << montecarlo(f_sphere_r2, h_sphere, hypercube_hs_4d_r2, 1'000'000) << std::endl;
  std::cout << montecarlo(f_sphere_r2, h_sphere, hypercube_hs_4d_r2, 1'000'000, 78.9568352087, .0001) << std::endl;

  return 0;
}