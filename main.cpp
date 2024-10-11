// g++ -o main main.cpp

#include <iostream>

#include "mul.hpp"

#include <vector>
#include <map>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

#include <memory>

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

int main()
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

    //

    std::map<std::string, double> map;

    map["one"] = 1.0;
    map["two"] = 2.0;
    map["three"] = 3.0;

    for (auto m : map)
    {
        std::cout << m.first << " : " << m.second << std::endl;
    }

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

    return 0;
}