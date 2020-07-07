// redist.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <istream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>

//#include <chrono>
//#include <ctime>
//#include <utility>

#include <boost/mpi.hpp>

// Timing
#include <boost/timer/timer.hpp>
#include <boost/chrono.hpp>

//boost polygon includes
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

// Includes from CGAL
//#include <boost/iterator/zip_iterator.hpp>

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Triangulation_vertex_base_with_info_2.h>
//#include <CGAL/Point_set_2.h>
//#include <CGAL/Polygon_2.h>

/*
struct velpnt {
    double z;
    double vx;
    double vy;
    double vz;
    double diam;
    double ratio;
    int proc;
};
*/

struct PntVel {
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double diam;
    double ratio;
    int proc;
};

//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Triangulation_vertex_base_with_info_2<velpnt, K> vb2D;
//typedef CGAL::Triangulation_data_structure_2<vb2D> Tds2D;
//typedef CGAL::Point_set_2<K, Tds2D>::Vertex_handle Vertex_handle2D;
//typedef CGAL::Point_set_2<K, Tds2D> PointSet2;
//typedef K::Point_2 cgal_point_2;
//typedef CGAL::Polygon_2<K> Polygon_2;


typedef boost::geometry::model::d2::point_xy<double> bpoint;
typedef boost::geometry::model::polygon<bpoint> bpoly;

/*
struct domain {
    std::vector<cgal_point_2> xy;
    //std::vector<double> x;
    //std::vector<double> y;
};
*/

//typedef std::map<int, domain> DomainList;

typedef std::map<int, bpoly> DomainListPoly;

struct inputs {
    int NinitDom;
    int Nzeros;
    std::string prefix;
    std::string suffix;
    std::string ExpadnedDomainFile;
    std::string ActualDomainFile;
    std::string NewPrefix;
};

bool readInputFile(std::string filename, inputs& input) {
    bool outcome = false;
    std::ifstream datafile(filename.c_str());
    if (!datafile.good()) {
        std::cout << "Can't open the file" << filename << std::endl;
    }
    else {
        std::string line;
        {// Number of domains that the initial velocity field is split
            getline(datafile, line);
            std::istringstream inp(line.c_str());
            inp >> input.NinitDom;
            // File name prefix of the existing velocity field
            inp >> input.prefix;
            // Number of zeros
            inp >> input.Nzeros;
            // Ending of file name
            inp >> input.suffix;
        }
        { 
            getline(datafile, line);
            std::istringstream inp(line.c_str());
            inp >> input.ActualDomainFile;
            inp >> input.ExpadnedDomainFile;
        }
        {
            getline(datafile, line);
            std::istringstream inp(line.c_str());
            inp >> input.NewPrefix;
        }
        outcome = true;
        datafile.close();
    }
    return outcome;
}

bool readVelocityFiles(std::string filename, std::vector< PntVel>& points, int myRank, 
    DomainListPoly& expanded, DomainListPoly& actual) {
    bool outcome = false;
    std::ifstream datafile(filename.c_str());
    if (!datafile.good()) {
        std::cout << "Can't open the file" << filename << std::endl;
    }
    else {
        std::cout << "Processor " << myRank << " reads " << filename << std::endl;
        // Start timing
        //auto start = std::chrono::high_resolution_clock::now();
        boost::timer::cpu_timer timer;

        std::string line;
        PntVel pv;
        double x, y;
        int count_lines = 0;
        while (getline(datafile, line)) {
            if (line.size() > 1) {
                count_lines++;
                std::istringstream inp(line.c_str());
                inp >> x;
                inp >> y;
                if (boost::geometry::within(bpoint(x, y), expanded[myRank])) {
                    int proc = -9;
                    if (boost::geometry::within(bpoint(x, y), actual[myRank])) {
                        proc = myRank;
                    }
                    else {
                        DomainListPoly::iterator it = actual.begin();
                        for (; it != actual.end(); ++it) {
                            if (it->first == myRank)
                                continue;
                            if (boost::geometry::within(bpoint(x, y), it->second)) {
                                proc = it->first;
                                break;
                            }
                        }
                    }
                    if (proc == -9) {
                        std::cout << "Cant find actual domain for point (" << x << "," << y << ") in my expanded domain with id " << myRank <<  std::endl;
                    }
                    pv.x = x;
                    pv.y = y;
                    inp >> pv.z;
                    inp >> pv.vx;
                    inp >> pv.vy;
                    inp >> pv.vz;
                    pv.proc = proc;
                    inp >> proc;
                    inp >> pv.diam;
                    inp >> pv.ratio;
                    points.push_back(pv);
                }
            }
        }
        datafile.close();
        
        //auto finish = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> elapsed = finish - start;
        boost::chrono::duration<double> seconds = boost::chrono::nanoseconds(timer.elapsed().user);
        std::cout << ". . . .Proc " << myRank << " spend " << seconds.count() / 60 << " min to read " << count_lines << " lines" << std::endl;
        outcome = true;
    }
    return outcome;

}

/*
bool readVelocityFields(std::string filename, PointSet2& Pset) {
    bool outcome = false;
    std::ifstream datafile(filename.c_str());
    if (!datafile.good()) {
        std::cout << "Can't open the file" << filename << std::endl;
    }
    else {
        std::string line;
        double x, y;
        velpnt vpnt;
        int proc;
        std::vector< std::pair<cgal_point_2, velpnt> > xy_data;
        while (getline(datafile, line)) {
            if (line.size() > 1) {
                std::istringstream inp(line.c_str());
                inp >> x;
                inp >> y;
                inp >> vpnt.z;
                inp >> vpnt.vx;
                inp >> vpnt.vy;
                inp >> vpnt.vz;
                inp >> proc;
                inp >> vpnt.diam;
                inp >> vpnt.ratio;
                vpnt.proc = -9;
                cgal_point_2 p(x, y);
                xy_data.push_back(std::make_pair(p, vpnt));
            }
        }
        Pset.insert(xy_data.begin(), xy_data.end());
        outcome = true;
        datafile.close();
    }
    return outcome;
}
*/

bool readDomain(std::string filename, DomainListPoly& dpoly) {
    bool outcome = false;
    std::ifstream datafile(filename.c_str());
    if (!datafile.good()) {
        std::cout << "Can't open the file" << filename << std::endl;
    }
    else {
        std::string line;
        getline(datafile, line);
        int Npoly, id, nVerts;
        double x, y;
        {
            std::istringstream inp(line.c_str());
            inp >> Npoly;
        }
        for (int i = 0; i < Npoly; ++i) {
            std::vector<bpoint> pnts;
            getline(datafile, line); {
                std::istringstream inp(line.c_str());
                inp >> id;
                inp >> nVerts;
            }
            for (int j = 0; j < nVerts; ++j) {
                getline(datafile, line);
                std::istringstream inp(line.c_str());
                inp >> x;
                inp >> y;
                pnts.push_back(bpoint(x, y));
            }
            bpoly poly;
            boost::geometry::assign_points(poly, pnts);
            boost::geometry::correct(poly);
            dpoly.insert(std::make_pair(id, poly));
        }
        datafile.close();
        outcome = true;
    }
    return outcome;
}

/*
bool readDomain(std::string filename, DomainList& dl) {
    bool outcome = false;
    std::ifstream datafile(filename.c_str());
    if (!datafile.good()) {
        std::cout << "Can't open the file" << filename << std::endl;
    }
    else {
        std::string line;
        getline(datafile, line);
        int Npoly, id, nVerts;
        double x, y;
        {
            std::istringstream inp(line.c_str());
            inp >> Npoly;
        }
        for (int i = 0; i < Npoly; ++i) {
            getline(datafile, line);
            {
                std::istringstream inp(line.c_str());
                inp >> id;
                inp >> nVerts;
            }
            domain d;
            //Polygon_2 poly;
            for (int j = 0; j < nVerts; ++j) {
                getline(datafile, line);
                std::istringstream inp(line.c_str());
                inp >> x;
                inp >> y;
                //poly.push_back(cgal_point_2(x, y));
                //d.x.push_back(x);
                //d.y.push_back(y);
                d.xy.push_back(cgal_point_2(x, y));
            }

            //if (poly.is_clockwise_oriented()) {
            //    poly.reverse_orientation();
            //}
            //for (int j = 0; j < nVerts; ++j) {
            //    d.x.push_back(poly.vertex(i).x());
            //    d.y.push_back(poly.vertex(i).y());
            //}
            dl.insert(std::make_pair(id, d));
        }
        datafile.close();
        outcome = true;
    }
    return outcome;
}
*/
std::string num2Padstr(int i, int n) {
    std::stringstream ss;
    ss << std::setw(n) << std::setfill('0') << i;
    return ss.str();
}

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;

    //std::cout << argc << std::endl;
    //std::cout << argv[0] << std::endl;
    inputs inp;
    bool tf = readInputFile(argv[1], inp);
    if (!tf)
        return 0;

    DomainListPoly actualDom, expandedDom;
    tf = readDomain(inp.ActualDomainFile, actualDom);
    if (!tf)
        return 0;
    tf = readDomain(inp.ExpadnedDomainFile, expandedDom);
    if (!tf)
        return 0;

    std::vector< PntVel> my_data;
    for (int iproc = 0; iproc < inp.NinitDom; ++iproc) {
        std::string filename = inp.prefix + num2Padstr(iproc, inp.Nzeros) + inp.suffix;
        tf = readVelocityFiles(filename, my_data, world.rank(), expandedDom, actualDom);
    }

    std::cout << "____ Processor " << world.rank() << " is printing..." << std::endl;
    std::string outfilename = inp.NewPrefix + num2Padstr(world.rank(), inp.Nzeros) + inp.suffix;
    std::ofstream outstream;
    outstream.open(outfilename.c_str());
    for (std::vector<PntVel>::iterator it = my_data.begin(); it != my_data.end(); ++it) {
        outstream << std::setprecision(2) << std::fixed
            << it->x << " " << it->y << " " << it->z << " "
            << std::setprecision(6) << std::fixed
            << it->vx << " " << it->vy << " " << it->vz << " " << it->proc << " "
            << std::setprecision(1) << std::fixed
            << it->diam << " " << it->ratio << std::endl;
    }
    outstream.close();
    world.barrier();
    return 0;

/*
    DomainList actualDom, expandedDom;
    tf = readDomain(inp.ActualDomainFile, actualDom);
    if (!tf)
        return 0;
    tf = readDomain(inp.ExpadnedDomainFile, expandedDom);
    if (!tf)
        return 0;


    domain myExpandeddomain = expandedDom[world.rank()];
    std::vector< std::pair<cgal_point_2, velpnt> > my_data;
    PointSet2 PsetMyRank; // This contains all the points that exist in the expanded domain of my rank
    // Loop through the velocity files and create Tree structures for each subdomain
    
    for (int iproc = 0; iproc < inp.NinitDom; ++iproc) {
        std::cout << "--Processor " << world.rank() << " is reading data from " << iproc << std::endl;
        // This contains all the points from the processor iproc
        PointSet2 PsetProc;
        
        std::string filename = inp.prefix + num2Padstr(iproc, inp.Nzeros) + inp.suffix;
        tf = readVelocityFields(filename, PsetProc);
        if (!tf) {
            std::cout << "Error while reading " << filename << std::endl;
            return 0;
        }

        // Find the points in my extended domain
        // split the rectangular areas into 2 triangles
        {// First Triangle
            std::list<Vertex_handle2D> LV;
            PsetProc.range_search(myExpandeddomain.xy[0], myExpandeddomain.xy[1], myExpandeddomain.xy[2], std::back_inserter(LV));
            std::list<Vertex_handle2D>::const_iterator it = LV.begin();
            for (; it != LV.end(); ++it) {
                my_data.push_back(std::make_pair((*it)->point(), (*it)->info()));
            }
        }
        {// First Triangle
            std::list<Vertex_handle2D> LV;
            PsetProc.range_search(myExpandeddomain.xy[0], myExpandeddomain.xy[2], myExpandeddomain.xy[3], std::back_inserter(LV));
            std::list<Vertex_handle2D>::const_iterator it = LV.begin();
            for (; it != LV.end(); ++it) {
                my_data.push_back(std::make_pair((*it)->point(), (*it)->info()));
            }
        }
    }
    PsetMyRank.insert(my_data.begin(), my_data.end());

    // Now loop through the new domain boundaries to assign processor ids
    std::cout << "=== Processor " << world.rank() << " is identifying  owners" << std::endl;
    DomainList::iterator itdom = actualDom.begin();
    for (; itdom != actualDom.end(); ++itdom) {
        
        { // Search the first triangle
            std::list<Vertex_handle2D> LV;
            PsetMyRank.range_search(itdom->second.xy[0], itdom->second.xy[1], itdom->second.xy[2], std::back_inserter(LV));
            std::list<Vertex_handle2D>::const_iterator it = LV.begin();
            for (; it != LV.end(); ++it) {
                (*it)->info().proc = itdom->first;
            }
        }
        { // Search the second triangle
            std::list<Vertex_handle2D> LV;
            PsetMyRank.range_search(itdom->second.xy[0], itdom->second.xy[2], itdom->second.xy[3], std::back_inserter(LV));
            std::list<Vertex_handle2D>::const_iterator it = LV.begin();
            for (; it != LV.end(); ++it) {
                (*it)->info().proc = itdom->first;
            }
        }
    }

    // Last print the new domain data
    std::cout << "____ Processor " << world.rank() << " is printing..." << std::endl;
    int count_pnts = 0;
    std::string outfilename = inp.NewPrefix + num2Padstr(world.rank(), inp.Nzeros) + inp.suffix;
    std::ofstream outstream;
    outstream.open(outfilename.c_str());
    {
        CGAL::Triangulation_2<K,Tds2D>::Finite_vertices_iterator itt = PsetMyRank.finite_vertices_begin();
        for (; itt != PsetMyRank.finite_vertices_end(); ++itt) {
            outstream << std::setprecision(2) << std::fixed
                << itt->point().x() << " " << itt->point().y() << " " << itt->info().z << " "
                << std::setprecision(6) << std::fixed
                << itt->info().vx << " " << itt->info().vy << " " << itt->info().vz << " "
                << itt->info().proc << " " << itt->info().diam << " " << itt->info().ratio << std::endl;
            count_pnts++;
        }
    }
    if (count_pnts != my_data.size()) {
        std::cout << "Count Points = " << count_pnts << " but  my_data.size is " << my_data.size() << std::endl;
    }
    outstream.close();

    return 0;
*/

}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
