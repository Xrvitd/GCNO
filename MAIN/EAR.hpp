#include <CGAL/Simple_cartesian.h>
#include <CGAL/edge_aware_upsample_point_set.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/IO/write_points.h>
#include <vector>
#include <fstream>
// types
typedef CGAL::Simple_cartesian<double> Kernel;

// Point with normal vector stored in a std::pair.
typedef std::pair<Kernel::Point_3, Kernel::Vector_3> PointVectorPair;
// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;
void EAR()
{
    const std::string input_filename = CGAL::data_file_path("points_3/before_upsample.xyz");
    const char* output_filename = "data/before_upsample_UPSAMPLED.xyz";
    // Reads a .xyz point set file in points[], *with normals*.
    std::vector<PointVectorPair> points;
    if (!CGAL::IO::read_points(input_filename,
        std::back_inserter(points),
        CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())
        .normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>())))
    {
        std::cerr << "Error: cannot read file " << input_filename << std::endl;
        return ;
    }
    //Algorithm parameters
    const double sharpness_angle = 25;   // control sharpness of the result.
    const double edge_sensitivity = 0;    // higher values will sample more points near the edges
    const double neighbor_radius = 0.25;  // initial size of neighborhood.
    const std::size_t number_of_output_points = points.size() * 4;
    //Run algorithm
    CGAL::edge_aware_upsample_point_set<Concurrency_tag>(
        points,
        std::back_inserter(points),
        CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
        normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()).
        sharpness_angle(sharpness_angle).
        edge_sensitivity(edge_sensitivity).
        neighbor_radius(neighbor_radius).
        number_of_output_points(number_of_output_points));
    // Saves point set.
    if (!CGAL::IO::write_points(output_filename, points,
        CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())
        .normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>())
        .stream_precision(17)))
        return ;
    return ;
}