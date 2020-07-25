#include "wake/export/vtk_exporter.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>

namespace wake {

bool VtkExporter::export_vti_2d(
    const std::string& filepath,
    const LatticeGrid2D& grid,
    const real_t* vorticity
) {
    std::ofstream out(filepath);
    if (!out.is_open()) {
        return false;
    }

    int nx = grid.nx();
    int ny = grid.ny();

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    out << "  <ImageData WholeExtent=\"0 " << (nx - 1) << " 0 " << (ny - 1) << " 0 0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
    out << "    <Piece Extent=\"0 " << (nx - 1) << " 0 " << (ny - 1) << " 0 0\">\n";
    out << "      <PointData Vectors=\"Velocity\" Scalars=\"Speed\">\n";

    // Velocity field
    out << "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            out << grid.ux(x, y) << " " << grid.uy(x, y) << " 0 ";
        }
        out << "\n";
    }
    out << "        </DataArray>\n";

    // Speed scalar
    out << "        <DataArray type=\"Float32\" Name=\"Speed\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            out << grid.speed(x, y) << " ";
        }
        out << "\n";
    }
    out << "        </DataArray>\n";

    // Density scalar
    out << "        <DataArray type=\"Float32\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            out << grid.rho(x, y) << " ";
        }
        out << "\n";
    }
    out << "        </DataArray>\n";

    // Vorticity (if present)
    if (vorticity) {
        out << "        <DataArray type=\"Float32\" Name=\"Vorticity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                out << vorticity[y * nx + x] << " ";
            }
            out << "\n";
        }
        out << "        </DataArray>\n";
    }

    // CellType mask
    out << "        <DataArray type=\"UInt8\" Name=\"CellType\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            out << static_cast<int>(grid.get_cell_type(x, y)) << " ";
        }
        out << "\n";
    }
    out << "        </DataArray>\n";

    out << "      </PointData>\n";
    out << "    </Piece>\n";
    out << "  </ImageData>\n";
    out << "</VTKFile>\n";

    return true;
}

bool VtkExporter::export_vti_3d(
    const std::string& filepath,
    const LatticeGrid3D& grid,
    const real_t* q_criterion,
    const Vec3* vorticity
) {
    std::ofstream out(filepath);
    if (!out.is_open()) {
        return false;
    }

    int nx = grid.nx();
    int ny = grid.ny();
    int nz = grid.nz();

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    out << "  <ImageData WholeExtent=\"0 " << (nx - 1) << " 0 " << (ny - 1) << " 0 " << (nz - 1) << "\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
    out << "    <Piece Extent=\"0 " << (nx - 1) << " 0 " << (ny - 1) << " 0 " << (nz - 1) << "\">\n";
    out << "      <PointData Vectors=\"Velocity\" Scalars=\"Speed\">\n";

    // Velocity field
    out << "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                out << grid.ux(x, y, z) << " " << grid.uy(x, y, z) << " " << grid.uz(x, y, z) << " ";
            }
            out << "\n";
        }
    }
    out << "        </DataArray>\n";

    // Speed scalar
    out << "        <DataArray type=\"Float32\" Name=\"Speed\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                out << grid.speed(x, y, z) << " ";
            }
            out << "\n";
        }
    }
    out << "        </DataArray>\n";

    // Density scalar
    out << "        <DataArray type=\"Float32\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                out << grid.rho(x, y, z) << " ";
            }
            out << "\n";
        }
    }
    out << "        </DataArray>\n";

    // Q-Criterion
    if (q_criterion) {
        out << "        <DataArray type=\"Float32\" Name=\"Q_Criterion\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {
                    out << q_criterion[(z * ny + y) * nx + x] << " ";
                }
                out << "\n";
            }
        }
        out << "        </DataArray>\n";
    }

    // Vorticity vector
    if (vorticity) {
        out << "        <DataArray type=\"Float32\" Name=\"Vorticity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {
                    const Vec3& v = vorticity[(z * ny + y) * nx + x];
                    out << v.x << " " << v.y << " " << v.z << " ";
                }
                out << "\n";
            }
        }
        out << "        </DataArray>\n";
    }

    // CellType mask
    out << "        <DataArray type=\"UInt8\" Name=\"CellType\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                out << static_cast<int>(grid.get_cell_type(x, y, z)) << " ";
            }
            out << "\n";
        }
    }
    out << "        </DataArray>\n";

    out << "      </PointData>\n";
    out << "    </Piece>\n";
    out << "  </ImageData>\n";
    out << "</VTKFile>\n";

    return true;
}

bool VtkExporter::write_pvd_manifest(
    const std::string& manifest_path,
    const std::string& file_prefix,
    int start_step,
    int end_step,
    int step_interval,
    real_t dt
) {
    std::ofstream out(manifest_path);
    if (!out.is_open()) {
        return false;
    }

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <Collection>\n";

    for (int s = start_step; s <= end_step; s += step_interval) {
        std::ostringstream fn;
        fn << file_prefix << "_" << std::setw(6) << std::setfill('0') << s << ".vti";
        real_t t = static_cast<real_t>(s) * dt;
        out << "    <DataSet timestep=\"" << t << "\" group=\"\" part=\"0\" file=\"" << fn.str() << "\"/>\n";
    }

    out << "  </Collection>\n";
    out << "</VTKFile>\n";

    return true;
}

} // namespace wake
