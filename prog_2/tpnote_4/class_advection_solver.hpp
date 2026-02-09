#pragma once
#include "class_mesh.hpp"
#include <cmath>
#include <iostream>
#include <vector>

class advection_solver
{
protected:
    std::vector<double> u_;
    triangulation& mesh_;

public:
    triangulation& mesh() { return mesh_; };

    advection_solver(triangulation& mesh) : mesh_(mesh) {};

    void set_initial_data(point2d Xc, double amp, double sig)
    {
        for (triangulation* scanner = &mesh_; scanner; scanner = (triangulation*)scanner->p_next())
        {
            double x = scanner->item().centroid().x();
            double y = scanner->item().centroid().y();

            // Calcul de la distance au carré par rapport au centre
            double r_squared = std::pow(x - Xc.x(), 2) + std::pow(y - Xc.y(), 2);

            // Valeur de la gaussienne au centre du triangle
            u_.push_back(amp * std::exp(-r_squared / (2.0 * std::pow(sig, 2))));
        }
    }

    void advance(double T_final, point2d velocity)
    {
        double t{0.};
        int step{0};
        std::size_t n_elements = mesh_.n_elements();

        double dt              = 0.5 * mesh().mesh_size() / sqrt(velocity * velocity);
        std::cout << "dt= " << dt << std::endl;

        std::vector<double> next_u(n_elements, 0.0); // Solution at t + dt

        set_initial_data(point2d(0., 0.), 1., 2.);
        export_to_VTK("initial.vtk");

        while (t < T_final)
        {
            if (t + dt > T_final)
                dt = T_final - t;

            int tri_index{0};

            for (triangulation* p_scan = &mesh(); p_scan; p_scan = (triangulation*)p_scan->p_next())
            {
                double flux_total = 0;

                for (auto kk = 0; kk < 3; ++kk)
                {
                    int neighbor_index = p_scan->item().neighbor(kk);

                    if (neighbor_index != -1)
                    {
                        double vn = velocity * p_scan->item().normal(kk);

                        // upwind flux
                        flux_total += vn > 0 ? vn * u_[tri_index] * p_scan->item().faces_length(kk) : vn * u_[neighbor_index] * p_scan->item().faces_length(kk);
                    }
                }

                next_u[tri_index] = u_[tri_index] - (dt / p_scan->item().area()) * flux_total;
                tri_index++;
            }

            // prepare next time iteration
            u_ = next_u;
            t += dt;
            if (++step % 1 == 0)
            {
                std::cout << "Temps---: " << t << std::endl;
            }
        }
        export_to_VTK("final.vtk");
    }

    void export_to_VTK(std::string const& filename)
    {
        std::ofstream file(filename);
        if (!file.is_open())
            return;

        file << "# vtk DataFile Version 3.0" << std::endl;
        file << "Advection Initial Condition" << std::endl;
        file << "ASCII" << std::endl;
        file << "DATASET UNSTRUCTURED_GRID" << std::endl;

        // 1. Écriture des points (sommets)
        file << "POINTS " << mesh_.n_vertices() << " double" << std::endl;
        for (const auto& v : mesh_.vertices())
        {
            file << v->x() << " " << v->y() << " 0.0" << std::endl; // Z=0 pour la 2D
        }

        // 2. Écriture des cellules (triangles)
        // Format : CELLS [nb_cellules] [nb_total_entiers]
        // Pour chaque triangle : 3 (nb de sommets) + index1 + index2 + index3 = 4 entiers
        file << "CELLS " << mesh_.n_elements() << " " << mesh_.n_elements() * 4 << std::endl;

        for (triangulation* scanner = &mesh_; scanner; scanner = (triangulation*)scanner->p_next())
        {
            file << "3 " << scanner->item().vertex(0)->index() << " " << scanner->item().vertex(1)->index() << " " << scanner->item().vertex(2)->index() << std::endl;
        }
        // 3. Type de cellules (5 correspond au triangle dans VTK)
        file << "CELL_TYPES " << mesh_.n_elements() << std::endl;
        for (size_t i = 0; i < mesh_.n_elements(); ++i)
        {
            file << "5" << std::endl;
        }

        // 4. Données aux cellules
        file << "CELL_DATA " << u_.size() << std::endl;
        file << "SCALARS U_Initial double 1" << std::endl;
        file << "LOOKUP_TABLE default" << std::endl;
        for (double val : u_)
        {
            file << val << std::endl;
        }

        file.close();
        std::cout << "Fichier VTK genere : " << filename << std::endl;
    }

    // TP4
    void advance_edge(double T_final, point2d velocity);
};




// tp4
void advection_solver::advance_edge(double T_final, point2d velocity){
    double t{0.};
    int step{0};
    double velocity_n{0.};
    int T_L, T_R;

    

    std::size_t n_elements = mesh_.n_elements();
    std::vector<double> phi(n_elements, 0.0);
    std::vector<double> areas(n_elements, 0.0);

    for (triangulation* scanner = &mesh_; scanner; scanner = (triangulation*)scanner->p_next())
    {
        areas[scanner->item().index()] = scanner->item().area();
    }

    double dt              = 0.5 * mesh().mesh_size() / sqrt(velocity * velocity);
    std::cout << "dt= " << dt << std::endl;

    std::vector<double> next_u(n_elements, 0.0); // Solution at t + dt

    set_initial_data(point2d(0., 0.), 1., 2.);
    export_to_VTK("initial.vtk");
    while (t < T_final)
    {
        if (t + dt > T_final)
            dt = T_final - t;

        for (auto ii = 0; ii < n_elements; ++ii)
            phi[ii] = 0.;

        for (auto& e : mesh_.edges()){
            velocity_n = velocity * e.normal();
            T_L = e.neighbor(0);
            T_R = e.neighbor(1);
            double flux_total = 0;

            if (velocity_n > 0){
                flux_total = velocity_n * u_[T_L] * e.length();
            }
            else{
                flux_total = velocity_n * (u_[T_R]!=-1 ? u_[T_R] : 0) * e.length();
            }
            phi[T_L] += flux_total;
            if (T_R != -1){
                phi[T_R] -= flux_total;
            }
        }

        for (auto ii = 0; ii < n_elements; ++ii){
            u_[ii] -= (dt / areas[ii]) * phi[ii];
        }

        t += dt;
        step++;
        if (step % 1 == 0)
        {
            std::cout << "Temps---: " << t << std::endl;
        }
    }
    export_to_VTK("final.vtk");
}