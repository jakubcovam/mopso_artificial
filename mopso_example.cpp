/**
 * @file mopso_example.cpp
 * @brief Multi-objective Particle Swarm Optimization (MOPSO) for an artificial problem.
 *
 * This program implements the MOPSO algorithm to find the Pareto front for the following two objectives:
 * - **Objective 1 (f1):** \( f1(x1, x2) = x1^2 + x2^2 \)
 * - **Objective 2 (f2):** \( f2(x1, x2) = (x1 - 2)^2 + (x2 - 2)^2 \)
 *
 * **Decision Variables:**
 * - \( x1 \in [0, 2] \)
 * - \( x2 \in [0, 2] \)
 *
 * @details
 * The MOPSO algorithm searches for a set of non-dominated solutions (Pareto front) by optimizing two conflicting objectives.
 * Particles in the swarm represent potential solutions, moving through the decision space guided by personal and global best positions.
 * The problem is inspired by Example 20.2 in a book of Dan Simon (2013): Evolutionary Optimization Algorithms.
 *
 * **Program Output:**
 * - Generates an output file **`all_example_solutions.csv`** containing all evaluated solutions during the optimization process.
 * - The CSV file includes the following columns:
 *   - `f1`: Objective function 1 value.
 *   - `f2`: Objective function 2 value.
 *   - `x1`: Decision variable \( x1 \).
 *   - `x2`: Decision variable \( x2 \).
 *   - `Type`: Indicates whether the solution is part of the final Pareto front (`'Pareto'`) or is a dominated solution (`'Dominated'`).
 *
 * **How to Use the Output:**
 * - The CSV file can be used to visualize the Pareto front and Pareto set.
 * - Use the provided Python script `mopso_example_plot_pareto.py` to generate plots.
 * - Ensure you have the required Python libraries (`matplotlib`, `pandas`) installed.
 *
 * **Compilation and Execution Instructions:**
 * - Compile the program using a C++11 compatible compiler as:
 *   ```
 *   g++ -std=c++11 -o mopso_example mopso_example.cpp
 *   ```
 * - Run the executable:
 *   ```
 *   ./mopso_example
 *   ```
 *
 * @author
 * Michala Jakubcová
 * https://github.com/jakubcovam
 *
 * @date
 * November 12, 2024
 *
 * @version
 * 1.0
 *
 * @note
 * - The program uses a random seed based on the current time. For reproducibility, set a fixed seed using `srand()` with a specific value.
 * - Adjust the PSO parameters (swarm size, number of iterations, coefficients) as needed for your specific use case.
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <limits>
#include <random>
#include <fstream>    // For file operations
#include <tuple>      // For std::tuple
#include <string>     // For std::string

/**
 * @brief Objective function f1: \( f1 = x1^2 + x2^2 \).
 * @param position A vector containing decision variables \( x1 \) and \( x2 \).
 * @return The computed value of objective function f1.
 */
double objective1(const std::vector<double>& position) {
    double x1 = position[0];
    double x2 = position[1];
    return x1 * x1 + x2 * x2;
}

/**
 * @brief Objective function f2: \( f2 = (x1 - 2)^2 + (x2 - 2)^2 \).
 * @param position A vector containing decision variables \( x1 \) and \( x2 \).
 * @return The computed value of objective function f2.
 */
double objective2(const std::vector<double>& position) {
    double x1 = position[0];
    double x2 = position[1];
    return (x1 - 2.0) * (x1 - 2.0) + (x2 - 2.0) * (x2 - 2.0);
}

/**
 * @struct Particle
 * @brief Represents a particle in the swarm with position, velocity, objectives, and personal bests.
 */
struct Particle {
    std::vector<double> position;     ///< Decision variables x1 and x2
    std::vector<double> velocity;     ///< Velocities for x1 and x2
    std::vector<double> objectives;   ///< Objective function values [f1, f2]
    std::vector<double> bestPosition; ///< Personal best position
    std::vector<double> bestObjectives; ///< Objectives at personal best position
};

/**
 * @brief Checks if one solution dominates another in multi-objective optimization.
 * @param obj1 Objectives of the first solution.
 * @param obj2 Objectives of the second solution.
 * @return True if obj1 dominates obj2, false otherwise.
 *
 * @details
 * A solution obj1 dominates obj2 if it is no worse in all objectives and better in at least one.
 */
bool dominates(const std::vector<double>& obj1, const std::vector<double>& obj2) {
    bool betterInAny = false;
    for (size_t i = 0; i < obj1.size(); ++i) {
        if (obj1[i] > obj2[i]) {
            return false; // obj1 is worse in at least one objective
        } else if (obj1[i] < obj2[i]) {
            betterInAny = true;
        }
    }
    return betterInAny;
}

/**
 * @brief Updates the archive of non-dominated solutions (Pareto front).
 * @param particle The particle to consider for inclusion in the archive.
 * @param archive The current archive of non-dominated solutions.
 *
 * @details
 * - Removes any solutions in the archive that are dominated by the new particle.
 * - Adds the new particle to the archive if it is not dominated by any existing solutions.
 */
void updateArchive(const Particle& particle, std::vector<Particle>& archive) {
    // Remove dominated solutions from the archive
    archive.erase(std::remove_if(archive.begin(), archive.end(),
        [&](const Particle& p) { return dominates(particle.objectives, p.objectives); }),
        archive.end());

    // Check if the particle is dominated by any in the archive
    for (const auto& p : archive) {
        if (dominates(p.objectives, particle.objectives)) {
            return; // Particle is dominated, do not add to archive
        }
    }

    // Particle is non-dominated, add to archive
    archive.push_back(particle);
}

/**
 * @brief Selects a global best particle from the archive (random selection).
 * @param archive The archive of non-dominated solutions.
 * @return A randomly selected particle from the archive.
 */
Particle selectGlobalBest(const std::vector<Particle>& archive) {
    size_t index = rand() % archive.size();
    return archive[index];
}

/**
 * @brief Updates a particle's velocity and position based on PSO equations.
 * @param particle The particle to update.
 * @param globalBest The global best particle selected from the archive.
 * @param w Inertia weight.
 * @param c1 Cognitive coefficient.
 * @param c2 Social coefficient.
 *
 * @details
 * - Velocity is updated using the inertia component, cognitive component (personal best), and social component (global best).
 * - Position is updated by adding the new velocity.
 * - Variable bounds are enforced to keep the position within [0, 2].
 */
void updateParticle(Particle& particle, const Particle& globalBest, double w, double c1, double c2) {
    for (size_t i = 0; i < particle.position.size(); ++i) {
        double r1 = static_cast<double>(rand()) / RAND_MAX;
        double r2 = static_cast<double>(rand()) / RAND_MAX;
        particle.velocity[i] = w * particle.velocity[i]
            + c1 * r1 * (particle.bestPosition[i] - particle.position[i])
            + c2 * r2 * (globalBest.position[i] - particle.position[i]);

        // Update position
        particle.position[i] += particle.velocity[i];

        // Apply variable bounds
        if (particle.position[i] < 0.0) particle.position[i] = 0.0;
        if (particle.position[i] > 2.0) particle.position[i] = 2.0;
    }
}

/**
 * @brief Evaluates objectives, updates personal bests, updates the archive, and collects data.
 * @param particle The particle to evaluate.
 * @param archive The archive of non-dominated solutions.
 * @param allSolutions A vector to store all evaluated solutions for visualization.
 *
 * @details
 * - Computes the objective functions for the particle's current position.
 * - Updates the particle's personal best if the current objectives dominate the previous best.
 * - Updates the archive with the new particle if it is non-dominated.
 * - Collects the solution data for later analysis.
 */
void evaluateAndUpdate(Particle& particle, std::vector<Particle>& archive, std::vector<std::tuple<double, double, double, double>>& allSolutions) {
    // Evaluate objectives
    double f1 = objective1(particle.position);
    double f2 = objective2(particle.position);

    particle.objectives = { f1, f2 };

    // Update personal best
    if (dominates(particle.objectives, particle.bestObjectives)) {
        particle.bestPosition = particle.position;
        particle.bestObjectives = particle.objectives;
    }

    // Update archive
    updateArchive(particle, archive);

    // Collect solutions (without 'Type' for now)
    allSolutions.push_back(std::make_tuple(
        particle.objectives[0],  // f1
        particle.objectives[1],  // f2
        particle.position[0],    // x1
        particle.position[1]     // x2
    ));
}

/**
 * @brief Compares two positions with a tolerance to account for floating-point precision.
 * @param pos1 First position vector.
 * @param pos2 Second position vector.
 * @param tol Tolerance for comparison (default is 1e-6).
 * @return True if positions are equal within the tolerance, false otherwise.
 */
bool positionsEqual(const std::vector<double>& pos1, const std::vector<double>& pos2, double tol = 1e-6) {
    if (pos1.size() != pos2.size()) return false;
    for (size_t i = 0; i < pos1.size(); ++i) {
        if (std::fabs(pos1[i] - pos2[i]) > tol) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Main function executing the MOPSO algorithm.
 * @return Exit status.
 *
 * @details
 * - Initializes the swarm and archive.
 * - Runs the main optimization loop, updating particles and the archive.
 * - After optimization, labels solutions based on the final Pareto front.
 * - Writes all solutions to a CSV file for visualization.
 */
int main() {
    srand(static_cast<unsigned>(time(0)));

    // Define PSO parameters
    const int swarmSize = 50;
    const int maxIterations = 100;
    const double w = 0.5;    // Inertia weight
    const double c1 = 1.5;   // Cognitive coefficient
    const double c2 = 1.5;   // Social coefficient

    // Initialize swarm and archive
    std::vector<Particle> swarm(swarmSize);
    std::vector<Particle> archive;

    // Vector to store all solutions for visualization
    std::vector<std::tuple<double, double, double, double>> allSolutions;

    // Initialize particles
    for (auto& particle : swarm) {
        // Initialize positions (x1 and x2)
        particle.position.resize(2);
        particle.velocity.resize(2);
        particle.bestPosition.resize(2);
        particle.bestObjectives = { std::numeric_limits<double>::max(), std::numeric_limits<double>::max() };

        particle.position[0] = static_cast<double>(rand()) / RAND_MAX * 2.0; // x1: 0 to 2
        particle.position[1] = static_cast<double>(rand()) / RAND_MAX * 2.0; // x2: 0 to 2

        // Initialize velocities
        particle.velocity[0] = 0.0;
        particle.velocity[1] = 0.0;

        // Initialize personal bests
        particle.bestPosition = particle.position;

        // Evaluate objectives, update archive, and collect data
        evaluateAndUpdate(particle, archive, allSolutions);
    }

    // Main optimization loop
    for (int iter = 0; iter < maxIterations; ++iter) {
        for (auto& particle : swarm) {
            // Select global best
            Particle globalBest = selectGlobalBest(archive);

            // Update particle
            updateParticle(particle, globalBest, w, c1, c2);

            // Evaluate objectives, update archive, and collect data
            evaluateAndUpdate(particle, archive, allSolutions);
        }

        // Optional: Output iteration information
        std::cout << "Iteration " << iter + 1 << " completed. Archive size: " << archive.size() << std::endl;
    }

    // After optimization, label all solutions based on final Pareto front
    // Create a set of final Pareto front positions for easy lookup
    std::vector<std::vector<double>> finalParetoPositions;
    for (const auto& p : archive) {
        finalParetoPositions.push_back(p.position);
    }

    // Write all solutions to a CSV file
    std::ofstream outFile("all_example_solutions.csv");
    if (!outFile) {
        std::cerr << "Error opening file for writing." << std::endl;
        return 1;
    }

    // Write header
    outFile << "f1,f2,x1,x2,Type\n";
    for (const auto& entry : allSolutions) {
        double f1 = std::get<0>(entry);
        double f2 = std::get<1>(entry);
        double x1 = std::get<2>(entry);
        double x2 = std::get<3>(entry);

        // Check if the position is in the final Pareto front
        std::string solutionType = "Dominated";
        for (const auto& pos : finalParetoPositions) {
            if (positionsEqual(pos, std::vector<double>{x1, x2})) {
                solutionType = "Pareto";
                break;
            }
        }

        outFile << f1 << "," << f2 << ","
                << x1 << "," << x2 << ","
                << solutionType << "\n";
    }
    outFile.close();

    std::cout << "\nAll solutions data saved to 'all_example_solutions.csv'" << std::endl;

    return 0;
}
