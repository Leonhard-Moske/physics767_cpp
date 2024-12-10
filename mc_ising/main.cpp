// g++ -o main main.cpp

#include <iostream>

#include <string>

#include <vector>
#include <map>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

#include <memory>

#include <fstream>

#include <chrono>

class Ising_MC
{
public:
    Eigen::MatrixXi spins;
    int N;

    Eigen::MatrixXi Lattice(std::string k)
    {
        if (k == "r")
        {
            // generate Matrix with entries -1 or 1 at random
            Eigen::MatrixXi spins = Eigen::MatrixXi::Random(N, N);
            spins = (spins.array() > 0).cast<int>() * 2 - 1;
            return spins;
        }
        else
        {
            // return Matrix with all entries 1
            Eigen::MatrixXi spins = Eigen::MatrixXi::Ones(N, N);
            return spins;
        }
    };

    double energy_change(double j, std::vector<std::vector<int>> Nb, std::vector<int> S, Eigen::MatrixXi Lat)
    {
        double sum = 0;
        for (int i = 0; i < Nb.size(); i++)
        {
            sum += Lat(S[0], S[1]) * Lat(Nb[i][0], Nb[i][1]);
        }
        return 2 * j * sum;
    };

    double acceptance(double de, double beta)
    {
        if (de < 0)
        {
            return -1;
        }
        else
        {
            // draw random number between 0 and 1
            double r = (double)rand() / (RAND_MAX);
            if (r <= exp(-beta * de))
            {
                return -1;
            }
            else
            {
                return 1;
            }
        }
    };

    std::vector<std::vector<int>> neighbor(int i, int j)
    {
        // create empty vector of size 4 with 2 entries
        std::vector<std::vector<int>> Nb(4, std::vector<int>(2));

        if (i == 0)
        {
            Nb[2][0] = N - 1;
            Nb[2][1] = j;
        }
        else
        {
            Nb[2][0] = i - 1;
            Nb[2][1] = j;
        }
        if (i == N - 1)
        {
            Nb[3][0] = 0;
            Nb[3][1] = j;
        }
        else
        {
            Nb[3][0] = i + 1;
            Nb[3][1] = j;
        }
        if (j == 0)
        {
            Nb[0][0] = i;
            Nb[0][1] = N - 1;
        }
        else
        {
            Nb[0][0] = i;
            Nb[0][1] = j - 1;
        }
        if (j == N - 1)
        {
            Nb[1][0] = i;
            Nb[1][1] = 0;
        }
        else
        {
            Nb[1][0] = i;
            Nb[1][1] = j + 1;
        }
        return Nb;
    }

    double energy(double j, Eigen::MatrixXi Lat)
    {
        double E = 0;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (i == 0)
                {
                    E += j * int(Lat(i, j) * Lat(N - 1, j));
                }
                else
                {
                    E += j * int(Lat(i, j) * Lat(i - 1, j));
                }
                if (i == N - 1)
                {
                    E += j * int(Lat(i, j) * Lat(0, j));
                }
                else
                {
                    E += j * int(Lat(i, j) * Lat(i + 1, j));
                }
                if (j == 0)
                {
                    E += j * int(Lat(i, j) * Lat(i, N - 1));
                }
                else
                {
                    E += j * int(Lat(i, j) * Lat(i, j - 1));
                }
                if (j == N - 1)
                {
                    E += j * int(Lat(i, j) * Lat(i, 0));
                }
                else
                {
                    E += j * int(Lat(i, j) * Lat(i, j + 1));
                }
            }
        }
        return E;
    }
};

int main()
{
    int N = 8;
    double j = 1;
    int steps = 10000;
    Eigen::VectorXd E_ar = Eigen::VectorXd::Zero(steps);
    Eigen::VectorXd M_ar = Eigen::VectorXd::Zero(steps);

    Ising_MC Ising;
    Ising.N = N;
    Eigen::MatrixXi lat = Ising.Lattice("r");
    double e = Ising.energy(j, lat);
    double beta = 0.1;
    double acc_rate = 0;

    // start timing here

    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < steps; i++)
    {
        E_ar(i) = e / (N * N);
        M_ar(i) = double(lat.sum()) / double((N * N));

        std::cout << "Lattice: " << lat << std::endl;
        std::cout << "Magn: " << lat.sum() << std::endl;    

        for (int k = 0; k < N * N; k++)
        {
            // draw random site
            int x = rand() % N;
            int y = rand() % N;

            std::vector<std::vector<int>> neighbors = Ising.neighbor(x, y);

            double delta_e = Ising.energy_change(j, neighbors, {x, y}, lat);

            double acc = Ising.acceptance(delta_e, beta);

            if (acc == -1)
            {

                lat(x, y) = -lat(x, y);
                e += delta_e;
                acc_rate += 1;
            }
        }
    }

    // end timing here

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";

    N = 16;
    steps = 10000;
    j = 1;
    std::vector<double> beta_ar = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
    Eigen::MatrixXd E_ar_mat = Eigen::MatrixXd::Zero(beta_ar.size(), steps);
    Eigen::MatrixXd M_ar_mat = Eigen::MatrixXd::Zero(beta_ar.size(), steps);
    for (int t = 0; t < beta_ar.size(); t++)
    {
        Ising_MC Ising_beta;

        Ising_beta.N = N;
        Eigen::MatrixXi lat_beta = Ising_beta.Lattice("f");
        e = Ising_beta.energy(j, lat_beta);
        beta = beta_ar[t];
        acc_rate = 0;

        // start timing here

        start = std::chrono::high_resolution_clock::now();

            for (int i = 0; i < steps; i++)
            {
                E_ar_mat(t, i) = e / (N * N);
                M_ar_mat(t, i) = double(lat_beta.sum()) / double((N * N));
                for (int k = 0; k < N * N; k++)
                {
                    // draw random site
                    int x = rand() % N;
                    int y = rand() % N;

                    std::vector<std::vector<int>> neighbors = Ising_beta.neighbor(x, y);

                    double delta_e = Ising_beta.energy_change(j, neighbors, {x, y}, lat_beta);

                    double acc = Ising_beta.acceptance(delta_e, beta_ar[t]);

                    if (acc == -1)
                    {
                        lat_beta(x, y) = -lat_beta(x, y);
                        e += delta_e;
                        acc_rate += 1;
                    }
                }
            }

            // end timing here

            end = std::chrono::high_resolution_clock::now();

            elapsed_seconds = end - start;

            std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";
    }

    // print magnetization and energy
    std::ofstream file;
    file.open("data.txt");
    file << "E_ar = " << E_ar << std::endl;
    file << "M_ar = " << M_ar << std::endl;
    file << "E_ar_mat = " << E_ar_mat << std::endl;
    file << "M_ar_mat = " << M_ar_mat << std::endl;
    file.close();

    return 0;
}