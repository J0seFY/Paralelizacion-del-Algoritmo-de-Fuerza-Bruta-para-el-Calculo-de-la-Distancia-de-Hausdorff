#include <vector>
#include <limits>
#include <cmath>
#include <iostream>
#include <chrono>
#include <random>
#include <string>
#include <iomanip> // Para std::setw (formatear la salida)

/**
 * @struct Point
 * @brief Estructura simple para representar un punto 2D.
 */
struct Point{
    double x;
    double y;
};

/**
 * @brief Calcula la distancia euclidiana entre dos puntos 2D.
 *
 * @param p1 Primer punto.
 * @param p2 Segundo punto.
 * @return La distancia euclidiana entre p1 y p2.
 */
inline double euclidean_distance(const Point& p1, const Point& p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return std::sqrt(dx * dx + dy * dy);
}

/**
 * @brief Calcula la distancia de Hausdorff dirigida (h(A, B)).
 *
 * @param A Primer conjunto de puntos.
 * @param B Segundo conjunto de puntos.
 * @return La distancia de Hausdorff dirigida desde A a B.
 */
double directed_hausdorff_distance(const std::vector<Point>& A, const std::vector<Point>& B) {
    double max_min_dist = 0.0;

    // Este es el bucle que vamos a paralelizar
    for (const auto& a : A) {
        double min_dist = std::numeric_limits<double>::max();
        
        // Este es el bucle interno
        for (const auto& b : B) {
            double dist = euclidean_distance(a, b);
            if (dist < min_dist) {
                min_dist = dist;
            }
        }
        
        // Aquí ocurre la posible condición de carrera
        if (min_dist > max_min_dist) {
            max_min_dist = min_dist;
        }
    }

    return max_min_dist;
}

/**
 * @brief Calcula la distancia de Hausdorff completa (simétrica).
 *
 * d_H(A, B) = max(h(A, B), h(B, A))
 *
 * @param A Primer conjunto de puntos.
 * @param B Segundo conjunto de puntos.
 * @return La distancia de Hausdorff completa.
 */
double hausdorff_distance(const std::vector<Point>& A, const std::vector<Point>& B) {
    double h_A_B = directed_hausdorff_distance(A, B);
    double h_B_A = directed_hausdorff_distance(B, A);
    return std::max(h_A_B, h_B_A);
}

/**
 * @brief Genera una nube de puntos aleatorios.
 */
std::vector<Point> generateRandomPointCloud(int numPoints, double minVal = 0.0, double maxVal = 1.0) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(minVal, maxVal);
    std::vector<Point> cloud;
    cloud.reserve(numPoints);
    for (int i = 0; i < numPoints; ++i) {
        cloud.push_back({dis(gen), dis(gen)});
    }
    return cloud;
}

/**
 * @brief Ejecuta un experimento y mide el tiempo promedio.
 *
 * @param testName Nombre del experimento.
 * @param A Conjunto de puntos A.
 * @param B Conjunto de puntos B.
 * @param numRuns Número de veces que se repetirá la prueba.
 */
void run_experiment(const std::string& testName, const std::vector<Point>& A, const std::vector<Point>& B, int numRuns = 5) {
    std::cout << std::left << std::setw(20) << testName 
              << std::setw(12) << A.size() 
              << std::setw(12) << B.size();

    double total_duration_ms = 0.0;
    double dist = 0.0;

    for (int i = 0; i < numRuns; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        
        dist = hausdorff_distance(A, B); // <-- Usamos la función completa

        auto end = std::chrono::high_resolution_clock::now();
        total_duration_ms += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    }

    double avg_duration_ms = total_duration_ms / numRuns;

    std::cout << std::fixed << std::setprecision(2) << std::setw(15) << avg_duration_ms 
              << std::setprecision(6) << "  (Dist: " << dist << ")" << std::endl;
}

int main() {
    
    // --- 1. Definición de los Datos de Prueba ---
    std::cout << "Generando conjuntos de datos de prueba..." << std::endl;
    auto setSmall_1k = generateRandomPointCloud(1000);
    auto setSmall_5k = generateRandomPointCloud(5000); // Set "pequeño" alternativo
    
    auto setMedium_10k = generateRandomPointCloud(10000);
    auto setMedium_20k = generateRandomPointCloud(20000);
    
    auto setLarge_50k = generateRandomPointCloud(50000);
    auto setLarge_70k = generateRandomPointCloud(70000);
    
    std::cout << "¡Datos generados! Iniciando experimentos (promedio de 5 ejecuciones)..." << std::endl;
    std::cout << "--------------------------------------------------------------------------" << std::endl;
    std::cout << std::left << std::setw(20) << "Nombre Test" 
              << std::setw(12) << "Tamaño A" 
              << std::setw(12) << "Tamaño B" 
              << std::setw(15) << "Tiempo Prom (ms)" 
              << "  Resultado (para verificar)" << std::endl;
    std::cout << "--------------------------------------------------------------------------" << std::endl;

    // --- 2. Ejecución Sistemática de Pruebas ---
    
    run_experiment("Pequeño (1k vs 1k)", setSmall_1k, setSmall_1k);
    run_experiment("Pequeño (1k vs 5k)", setSmall_1k, setSmall_5k);
    
    run_experiment("Mediano (10k vs 10k)", setMedium_10k, setMedium_10k);
    run_experiment("Mediano (10k vs 20k)", setMedium_10k, setMedium_20k);
    
    run_experiment("Grande (50k vs 50k)", setLarge_50k, setLarge_50k);
    run_experiment("Grande (50k vs 70k)", setLarge_50k, setLarge_70k);
    
    run_experiment("Asimétrico (1k vs 70k)", setSmall_1k, setLarge_70k);
    
    std::cout << "--------------------------------------------------------------------------" << std::endl;

    return 0;
}