// g++ -o main main.cpp

#include <iostream>

#include "mul.hpp"

#include <vector>
#include <map>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

#include <memory>

#include <fstream>

/*
class Function {
public:
    // Constructor that takes coefficients for a quadratic function ax^2 + bx + c
    Function(double a, double b, double c) : coeffs(Eigen::Vector3d(a, b, c)) {}

    // Evaluate the function at a given x
    double evaluate(double x) const {
        return coeffs(0) * x * x + coeffs(1) * x + coeffs(2);
    }

    // Evaluate the function at multiple points and store in a map
    std::map<double, double> evaluateAtPoints(const std::vector<double>& points) {
        std::map<double, double> results;
        for (double x : points) {
            results[x] = evaluate(x);
        }
        return results;
    }

    // Print the coefficients
    void printCoefficients() const {
        std::cout << "Coefficients: " << coeffs.transpose() << std::endl;
    }

private:
    Eigen::Vector3d coeffs; // Coefficients of the function
};
*/

/*
class HubbardHamiltonian {
public:
    HubbardHamiltonian(int size, double t, double U)
        : size(size), t(t), U(U), hamiltonian(Eigen::MatrixXd::Zero(size, size)) {
        buildHamiltonian();
    }

    const Eigen::MatrixXd& getHamiltonian() const {
        return hamiltonian;
    }

    void printHamiltonian() const {
        std::cout << "Hubbard Hamiltonian Matrix:\n" << hamiltonian << std::endl;
    }

private:
    int size;                 // Number of sites
    double t;                // Hopping parameter
    double U;                // On-site interaction strength
    Eigen::MatrixXd hamiltonian;

    void buildHamiltonian() {
        for (int i = 0; i < size; ++i) {
            // Hopping terms
            if (i < size - 1) {
                hamiltonian(i, i + 1) = -t; // Hop to the right
                hamiltonian(i + 1, i) = -t; // Hop to the left
            }
            // On-site interaction
            hamiltonian(i, i) += U; // Add on-site interaction
        }
    }
};
*/


class Car
{
public:
    Car(std::string name, double speed)
    {
        name_ = name;
        speed_ = speed;
    }

    void print()
    {
        std::cout << "Car: " << name_ << " Speed: " << speed_ << std::endl;
    }

    double getSpeed()
    {
        return speed_;
    }

    std::string getName()
    {
        return name_;
    }

private:
    std::string name_;
    double speed_;
};

class Truck : public Car
{
public:
    Truck(std::string name, double speed, double load) : Car(name, speed)
    {
        load_ = load;
    }

    void print()
    {
        std::cout << "Truck: " << getName() << " Speed: " << getSpeed() << " Load: " << load_ << std::endl;
    }

private:
    double load_;
};

double add(double a, double b)
{
    return a + b;
}

int main(int argc, char** argv)
{
    std::cout << "Hello World!" << std::endl;

    
    
    
    
    //

    double a = 1.5;
    double b = 2.0;
    double c = add(a, b);
    std::cout << "a + b = " << c << std::endl;

    //

    double d = mul(a, b);
    std::cout << "a * b = " << d << std::endl;

    //

    if (a > b)
    {
        std::cout << "a is greater than b" << std::endl;
    }
    else
    {
        std::cout << "a is not greater than b" << std::endl;
    }

    //

    std::vector<double> vec = {1.0, 2.0, 3.0};
    vec.push_back(4.0);

    // range-based for loop
    for (auto v : vec)
    {
        std::cout << v << std::endl;
    }

    for (int i = 0; i < vec.size(); i++)
    {
        std::cout << vec[i] << std::endl;
    }

    // vfr

    std::map<std::string, double> map;

    map["one"] = 1.0;
    map["two"] = 2.2;
    map["three"] = 3.0;

    for (auto m : map)
    {
        std::cout << m.first << " : " << m.second << std::endl;
    }
    std::cout << map["two"] << std::endl;

    //

    Eigen::MatrixXd mat(2, 2);
    mat << 1, 2, 3, 4;

    std::cout << mat << std::endl;

    mat.transposeInPlace();

    std::cout << mat << std::endl;

    Eigen::EigenSolver<Eigen::MatrixXd> es(mat);

    std::cout << "Eigenvalues: " << es.eigenvalues() << std::endl;

    auto transposed = mat.transpose();

    std::cout << transposed << std::endl;

    //
    Car car("Toyota", 100.0);

    car.print();

    Truck truck("Volvo", 80.0, 1000.0);

    truck.print();

    // introduce smart pointers
    std::shared_ptr<Car> car_ptr{std::make_shared<Car>("Toyota", 100.0)};

    car_ptr->print();
    auto also_car_pointer{car_ptr};

    std::cout << also_car_pointer.use_count() << std::endl;

    // save data to file

        // Sample vector of complex numbers
    std::vector<std::complex<double>> complex_vector = {
        {1, 2}, {3, 4}, {5, -6}
    };

    // Write to a CSV file
    std::ofstream file("complex_numbers.csv");
    if (file.is_open()) {
        file << "Real,Imaginary\n";  // Optional header
        for (const auto& c : complex_vector) {
            file << c.real() << "," << c.imag() << "\n";
        }
        file.close();
    } else {
        std::cerr << "Unable to open file";
    }


    return 0;
}