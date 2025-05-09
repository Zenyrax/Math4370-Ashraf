#include <vector>
#include <random>
#include <iostream>
#include <omp.h>
#include "mpi.h"
#include <chrono>

// Calculate the magnitude of a vector
double magnitude(const std::vector<double>& point) {
  double sum = 0.0;
  for (double v : point) {
    sum += v * v;
  }
  return std::sqrt(sum);
}

double montecarlo(double* elapsed, int id, int np, double (*f)(std::vector<double>), int (*h)(std::vector<double>, double(std::vector<double>)), std::vector<std::pair<double, double>> hypercube, int M) {
  int dimensions = hypercube.size();
  double buf[1];
  buf[0] = 0;
  MPI_Status status;
  if (id == 0) {
    auto start = std::chrono::high_resolution_clock::now();
    // Calculate the area of the hypercube

    double third = 0;
    double hypercube_area = 1;
    for (size_t i = 0; i < dimensions; i++) {
      hypercube_area *= hypercube[i].second - hypercube[i].first;
    }

    for (size_t i = 0; i < np - 1; i++) {
      MPI_Recv(buf, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      third += buf[0];
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    *elapsed = std::chrono::duration<double>(end - start).count();
    
    std::cout << "hypercube area: " << hypercube_area << " | iterations: " << M << std::endl;
    return hypercube_area * third / M;
  } else {
    std::random_device rd;
    std::mt19937 gen(rd() + id);
    std::vector<double> random_point(dimensions);
    for (size_t i = 0; i < M/(np - 1); i++) {
      // Generate a random point
      for (size_t j = 0; j < dimensions; j++) {
        std::uniform_real_distribution<> dis(hypercube[j].first, hypercube[j].second);
        random_point[j] = dis(gen);
      }

      // Value of f is added this sum if it passes the indicator function
      buf[0] += f(random_point) * h(random_point, f);
    }
    MPI_Send(buf, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
  }
  return 0;
}

int main(int argc, char* argv[]) {
  int ierr, numprocs, myid;
  ierr = MPI_Init(&argc, &argv);
  if (ierr != MPI_SUCCESS) {
     std::cerr << "Error in calling MPI_Init\n";
     return 1;
  }
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid); 

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
  double pi_res = montecarlo(&elapsed, myid, numprocs, f_pi, h_pi, hypercube_pi, 10'000'000);
  if (myid == 0) std::cout << "result: " << pi_res <<  " | error: " << std::fabs( 3.14159265359 - pi_res) << " | time: " << elapsed << "s" << std::endl;

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
  if (myid == 0) std::cout << "\nApproximating integral of x^2 between 0 and 2:" << std::endl;
  double quad_res = montecarlo(&elapsed, myid, numprocs, f_quad, h_quad2, hypercube_quad, 10'000'000);
  if (myid == 0) std::cout << "result: " << quad_res <<  " | error: " << std::fabs(2.66666666 - quad_res)  << " | time: " << elapsed << "s" << std::endl;

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
  if (myid == 0) std::cout << "\nApproximating volume of a sphere w/ radius 2:" << std::endl;
  double sphere_res = montecarlo(&elapsed, myid, numprocs, f_pi, h_sphere, hypercube_hs_3d_r2, 10'000'000);
  if (myid == 0) std::cout << "result: " << sphere_res <<  " | error: " << std::fabs(33.5103216383 - sphere_res)  << " | time: " << elapsed << "s"  << std::endl;
  
  // 4d hypersphere w/ radius 2
  std::vector<std::pair<double, double>> hypercube_hs_4d_r2;
  for (size_t i = 0; i < 4; i++) {
    hypercube_hs_4d_r2.push_back({-2, 2});
  }
  if (myid == 0) std::cout << "\nApproximating volume of a 4d hypersphere w/ radius 2:" << std::endl;
  double sphere_4d_res = montecarlo(&elapsed, myid, numprocs, f_pi, h_sphere, hypercube_hs_4d_r2, 10'000'000);
  if (myid == 0) std::cout << "result: " << sphere_4d_res <<  " | error: " << std::fabs(78.9568352087 - sphere_4d_res)  << " | time: " << elapsed << "s"  << std::endl;

  // 5d hypersphere w/ radius 2
  std::vector<std::pair<double, double>> hypercube_hs_5d_r2;
  for (size_t i = 0; i < 5; i++) {
    hypercube_hs_5d_r2.push_back({-2, 2});
  }
  if (myid == 0) std::cout << "\nApproximating volume of a 5d hypersphere w/ radius 2:" << std::endl;
  double sphere_5d_res = montecarlo(&elapsed, myid, numprocs, f_pi, h_sphere, hypercube_hs_5d_r2, 10'000'000);
  if (myid == 0) std::cout << "result: " << sphere_5d_res <<  " | error: " << std::fabs(168.441248445 - sphere_5d_res)  << " | time: " << elapsed << "s"  << std::endl;

  // 10d hypersphere w/ radius 2
  std::vector<std::pair<double, double>> hypercube_hs_10d_r2;
  for (size_t i = 0; i < 10; i++) {
    hypercube_hs_10d_r2.push_back({-2, 2});
  }
  if (myid == 0) std::cout << "\nApproximating volume of a 10d hypersphere w/ radius 2:" << std::endl;
  double sphere_10d_res = montecarlo(&elapsed, myid, numprocs, f_pi, h_sphere, hypercube_hs_10d_r2, 10'000'000);
  if (myid == 0) std::cout << "result: " << sphere_10d_res <<  " | error: " << std::fabs(2611.36797683 - sphere_10d_res)  << " | time: " << elapsed << "s"  << std::endl;

  ierr = MPI_Finalize();
  return 0;
}