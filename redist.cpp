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

#include <boost/filesystem.hpp>
#include <boost/mpi.hpp>

// Timing
#include <boost/timer/timer.hpp>
#include <boost/chrono.hpp>

//boost polygon includes
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

/**
 * PntVel is a structure that holds information for a single point.
 */
struct PntVel {
    // x coordinate
    double x;
    // y coordinate
    double y;
    // z coordinate
    double z;
    // x velocity
    //double vx;
    // y velocity
    //double vy;
    // z velocity
    //double vz;
    // the diameter of the mesh element that this point belonged in the original mesh
    //double diam;
    // the ratio of the xy length / z length of the mesh element
    //double ratio;
    // the processor that owns this point in the new distribution
    int proc;
    std::vector<double> doubleInfo;
    std::vector<int> intInfo;
};



// The boost definition for 2D points
typedef boost::geometry::model::d2::point_xy<double> bpoint;
// The boost definition for polygons that use 2D boost points
typedef boost::geometry::model::polygon<bpoint> bpoly;

// The definition for a list of boost 2D polygons
typedef std::map<int, bpoly> DomainListPoly;

/**
 * A structure that contains all the input parameters.
 * 
 * The input parameters are listed in a single file. Each line holds one parameter
 * The format of the input file is not flexible at all. No comments allowed or any literals with spaces or any weird looking characters
 * 
 * They must be listed in the following order
 */
struct inputs {
    // The number of subdomains that the problem was originaly solved. This is the number of processors used in the simulation step
    int NinitDom;
    /**
     * The velocity information is usually printed into files with the following format:
     * 
     * f:\UCDAVIS\ichnos\NPSAT_example\example1__0000.ich.
     * 
     * The above name essentially consists of three parts. 
     * 
     * 1. prefix :  f:\UCDAVIS\ichnos\NPSAT_example\example1__
     * 2. processor number : which must have identical zero padding for all processors. 
     * 3 suffix : .ich
     * 
     * Nzeros is the width of the processor number 
     */
    int Nzeros;
    // This is the file name along with the path if it lives in a different folder than the program
    std::string prefix;
    /** 
     * The characters after the processor id.
     * If the file is example1__0000_test01.ich then the suffix is _test01.ich
     */
    std::string suffix;
    /**
     * The file name that containts the new expanded domain boundaries.
     */
    std::string ExpadnedDomainFile;
    /**
     * The file name that containts the new actual domain boundaries.
     */
    std::string ActualDomainFile;
    // This is the prefix where the redistributed velocity fields will be printed
    std::string NewPrefix;
    std::vector<int> InfoType;
    int Ninfo;
    std::vector<int> prec;

};

/**
 * readInputFile parse the input file.
 * 
 * See the documentation of the #input for details about the format
 * 
 * \param filename
 * \param input
 * \return true if the reading was succefull.
 */
bool readInputFile(std::string filename, inputs& input) {
    bool outcome = false;
    std::ifstream datafile(filename.c_str());
    if (!datafile.good()) {
        std::cout << "Can't open the file " << filename << std::endl;
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

        {// get the information of the points
            getline(datafile, line);
            std::istringstream inp(line.c_str());
            inp >> input.Ninfo;
            for (int i = 0; i < input.Ninfo; i++){
                int t;
                inp >> t;
                input.InfoType.push_back(t);
            }
        }

        {// Get the precision for the floating numbers
            getline(datafile, line);
            std::istringstream inp(line.c_str());
            for (int i = 0; i < input.Ninfo; i++) {
                int t;
                inp >> t;
                input.prec.push_back(t);
            }
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

/**
 * Reads the original velocity files and selects the points that belong to myRank new subdomain.
 * 
 * 
 * \param filename is the original velocity file name.
 * \param points is a vector to holds the points that live on the new subdomain with id myRank
 * \param myRank is the current processor id
 * \param expanded is a list of the expanded domains
 * \param actual is a list of the actual domains
 * \return true if everything was succesfull 
 */
bool readVelocityFiles(std::string filename, std::vector< PntVel>& points, int myRank, 
    DomainListPoly& expanded, DomainListPoly& actual, inputs& input) {
    bool outcome = false;
    std::ifstream datafile(filename.c_str());
    if (!datafile.good()) {
        std::cout << "Can't open the file " << filename << std::endl;
    }
    else {
        std::cout << "Processor " << myRank << " reads " << filename << std::endl;
        // Start timing
        //auto start = std::chrono::high_resolution_clock::now();
        boost::timer::cpu_timer timer;

        std::string line;
        double x, y;
        int count_lines = 0;
        while (getline(datafile, line)) {
            PntVel pv;
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
                    pv.proc = proc;
                    double dd;
                    int ii;
                    for (int i = 0; i < input.Ninfo; i++)
                    {
                        if (input.InfoType[i] == 0) {
                            inp >> dd;
                            pv.doubleInfo.push_back(dd);
                        }
                        else {
                            inp >> ii;
                            pv.intInfo.push_back(ii);
                        }
                    }
                    //inp >> pv.vx;
                    //inp >> pv.vy;
                    //inp >> pv.vz;
                    //pv.proc = proc;
                    //inp >> proc;
                    //inp >> pv.diam;
                    //inp >> pv.ratio;
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

/**
 * Reads the 2D polygons that define the new subdomains.
 * 
 * \param filename is the name of the file
 * \param dpoly is a list of 2D polygons where each polygon is one domain.
 * \return true if reading was wihtout any issues
 * 
 * The fromat of the file is the following:
 * 
 * Line 1: Npoly Number of polygons listed in this file
 * 
 * Line 2: ID Nverts, ID is the processor ID and Nverts is the number of vertices of the polygon
 * 
 * Repeat Nverts times
 * X Y 
 * 
 * Repeat Npoly times the above starting from Line 2
 * 
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

/**
 * This converts an integer to string with specified number of zero padding.
 * 
 * \param i is the integer number to convert to string
 * \param n is the number of zeros to fill in
 * \return the integer in a string format
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
    std::cout << "Current path: " << boost::filesystem::current_path() << std::endl;
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
        tf = readVelocityFiles(filename, my_data, world.rank(), expandedDom, actualDom, inp);
    }

    std::cout << "____ Processor " << world.rank() << " is printing..." << std::endl;
    std::string outfilename = inp.NewPrefix + num2Padstr(world.rank(), inp.Nzeros) + inp.suffix;
    std::ofstream outstream;
    outstream.open(outfilename.c_str());
    for (std::vector<PntVel>::iterator it = my_data.begin(); it != my_data.end(); ++it) {
        outstream << std::setprecision(3) << std::fixed
            << it->x << " " << it->y << " " << it->z << " " << it->proc << " ";
            
        int d_cnt = 0;
        int i_cnt = 0;
        for (int i = 0; i < inp.Ninfo; i++){
            if (inp.InfoType[i] == 0) {
                outstream << std::setprecision(inp.prec[i]) << std::fixed << it->doubleInfo[d_cnt] << " ";
                d_cnt++;
            }
            else {
                outstream << it->intInfo[i_cnt] << " ";
                i_cnt++;
            }
        }
        outstream << std::endl;
    }
    outstream.close();
    world.barrier();
    return 0;
}

