// redist.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <istream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include<string>

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

#include <boost/spirit/include/qi.hpp>

int dgbRank = 0;
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
    std::vector<double> data;
};

struct cellIDint{
    cellIDint();
    std::string to_string();
    cellIDint(int aa, int bb, int cc){
        a = aa;
        b = bb;
        c = cc;
    }
    int a = -1;
    int b = -1;
    int c = -1;

};
std::string cellIDint::to_string(){
    return std::to_string(a) + ":" + std::to_string(b) + "_" + std::to_string(c);
}

cellIDint cellID2int(std::string str){
    std::vector<std::string> vec;
    bool res = boost::spirit::qi::parse(str.begin(), str.end(),
                                        boost::spirit::qi::as_string[*(boost::spirit::qi::char_ - ':' - "_")] %
                                        (boost::spirit::qi::lit(':') | boost::spirit::qi::lit('_')), vec);
    int a, b, c;
    a = stoi(vec[0]);
    b = stoi(vec[1]);
    if (b == 0){
        c = -9;
    }
    else{
        c = stoi(vec[2]);
    }

    return cellIDint(a, b, c);
}

struct cellData{
    std::vector<cellIDint> neighborCells;
    std::vector<PntVel> CellVelocities;
    std::vector<int> VelIds;
    double x;
    double y;
    double z;
    int id;
    //std::vector<int> neighborIds;
};

typedef std::map<int, cellData> L3map;
typedef std::map<int, L3map> L2map;
typedef std::map<int, L2map> L1map;

class cellIdList{
public:
    cellIdList(){};
    bool insert(std::string cell_string, cellData& CData){
        return insert(cellID2int(cell_string), CData);
    };
    bool insert(cellIDint clid, cellData& CData);

    bool addVelocityPoint(cellIDint clid, PntVel& pv);

    bool find(cellIDint clid);
    void clear();

    L1map CellIds;
    L1map::iterator itLev1;
    L2map::iterator itLev2;
    L3map::iterator itLev3;
    int idx = 0;
};

void cellIdList::clear() {
    CellIds.clear();
}

bool cellIdList::find(cellIDint clid) {
    itLev1 = CellIds.find(clid.a);
    if (itLev1 != CellIds.end()){
        itLev2 = itLev1->second.find(clid.b);
        if (itLev2 != itLev1->second.end()){
            itLev3 = itLev2->second.find(clid.c);
            if (itLev3 != itLev2->second.end()){
                return true;
            }
        }
    }
    return false;
}

bool cellIdList::insert(cellIDint clid, cellData& CData) {
    bool out = true;

    //int a = clid.a;
    //int b = clid.b;
    //int c = clid.c;
    itLev1 = CellIds.find(clid.a);
    if (itLev1 != CellIds.end()){
        itLev2 = itLev1->second.find(clid.b);
        if (itLev2 != itLev1->second.end()){
            itLev3 = itLev2->second.find(clid.c);
            if (itLev3 != itLev2->second.end()){
                std::cout << " The cell " << clid.to_string() << " is already in the list" << std::endl;
                out = false;
            }
            else{
                itLev2->second.insert(std::pair<int ,cellData >(clid.c, CData));
            }
        }
        else{
            L3map tmp;
            tmp.insert(std::pair<int, cellData>(clid.c, CData));
            itLev1->second.insert(std::pair<int, L3map>(clid.b, tmp));
        }
    }
    else{// If the a key doesn't exist in the map
        L3map tmp;
        tmp.insert(std::pair<int,cellData>(clid.c, CData));
        L2map tmp1;
        tmp1.insert(std::pair<int, L3map >(clid.b , tmp));
        CellIds.insert(std::pair<int, L2map>(clid.a, tmp1));
    }
    return out;
}

bool cellIdList::addVelocityPoint(cellIDint clid, PntVel &pv) {
    itLev1 = CellIds.find(clid.a);
    if (itLev1 != CellIds.end()){
        itLev2 = itLev1->second.find(clid.b);
        if (itLev2 != itLev1->second.end()){
            itLev3 = itLev2->second.find(clid.c);
            if (itLev3 != itLev2->second.end()){
                itLev3->second.CellVelocities.push_back(pv);
                return true;
            }
        }
    }
    return false;
}


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
    //std::vector<int> InfoType;
    int Ninfo;
    std::vector<int> prec;
    int Nprint;
    std::vector<int> printOrder;

    bool UseGraph;
    int overlapIter = 0;

};

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

int findProcessorOwner(double x, double y, DomainListPoly& actual, int rank){
    if (boost::geometry::within(bpoint(x, y), actual[rank])){
        return rank;
    }
    else{
        for(unsigned int i = 0; i < actual.size(); ++i){
            if ( i == rank)
                continue;
            if (boost::geometry::within(bpoint(x, y), actual[i])){
                return i;
            }
        }
    }
    return -9;
}

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
            // Is NPSAT graph available
            int tmp;
            inp >> tmp;
            if (tmp != 0 ){
                input.UseGraph = true;
                inp >> tmp;
                input.overlapIter = tmp;
            }

        }

        {// get the information of the points
            getline(datafile, line);
            std::istringstream inp(line.c_str());
            inp >> input.Ninfo;
            inp >> input.Nprint;
            for (int i = 0; i < input.Nprint; i++){
                int t;
                inp >> t;
                input.printOrder.push_back(t);
            }
        }

        {// Get the precision for the floating numbers
            getline(datafile, line);
            std::istringstream inp(line.c_str());
            for (int i = 0; i < input.Nprint; i++) {
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

bool readCellGraphFiles(std::string filename, cellIdList& CLact, cellIdList& CLexp,
                        int myrank, DomainListPoly& expanded,DomainListPoly& actual){
    bool outcome = false;
    std::ifstream datafile(filename.c_str());
    if (!datafile.good()) {
        std::cout << "Can't open the file " << filename << std::endl;
    }
    else{
        //if (myrank == dgbRank){
        //    std::cout << "Start HERE " << filename <<  std::endl;
        //}
        std::cout << "Processor " << myrank << " reads graph file " << filename << std::endl;
        boost::timer::cpu_timer timer;
        std::string line;
        double x, y, z;
        int n;
        int nAct = 0;
        int nExp = 0;
        int count_lines = 0;
        std::string celidstr;
        while (getline(datafile, line)){
            if (line.size() > 1){
                std::istringstream inp(line.c_str());
                inp >> x;
                inp >> y;
                inp >> z;
                bool inAct = false;
                bool inExp = false;
                if (boost::geometry::within(bpoint(x, y), actual[myrank])){
                    inAct = true;
                }
                else if (boost::geometry::within(bpoint(x, y), expanded[myrank])){
                    inExp = true;
                }
                if (inExp || inAct){
                    inp >> n;
                    std::vector<cellIDint> neighcells;
                    inp >> celidstr;
                    cellData cdata;
                    cdata.x = x;
                    cdata.y = y;
                    cdata.z = z;

                    cellIDint thisCell = cellID2int(celidstr);
                    for (int i = 0; i < n; ++i){
                        inp >> celidstr;
                        cdata.neighborCells.push_back(cellID2int(celidstr));
                    }
                    if (inAct){
                        CLact.insert(thisCell, cdata);
                        nAct++;
                    }
                    else if (inExp){
                        CLexp.insert(thisCell,cdata);
                        nExp++;
                    }
                }
            }
        }

        //if (myrank == dgbRank){
        //    std::cout << "End HERE " << filename <<  std::endl;
        //}

        boost::chrono::duration<double> seconds = boost::chrono::nanoseconds(timer.elapsed().user);
        std::cout << ". . . .Proc " << myrank << " spend " << seconds.count() / 60 << " min to read " << nAct << " " << nExp  << std::endl;
        outcome = true;
    }
    return outcome;

}

bool readVelocityWithGraph(std::string filename, int myRank, inputs& input, cellIdList& CLact/*, cellIdList& CLexp*/ ){
    bool outcome = false;
    std::ifstream datafile(filename.c_str());
    if (!datafile.good()) {
        std::cout << "Can't open the file " << filename << std::endl;
        return outcome;
    }
    else{
        //if (myRank == dgbRank){
        //    std::cout << "Start HERE " << filename <<  std::endl;
        //}
        std::cout << "Processor " << myRank << " reads velocity file " << filename << std::endl;
        boost::timer::cpu_timer timer;
        std::string line;
        double x, y, z;
        int cnt_data = 0;
        std::string celidstr;
        int cnt_lines = 0;
        while (getline(datafile, line)){
            if (line.size() > 1){
                //if (myRank == dgbRank){
                //std::cout << cnt_lines++ << std::endl;
                //std:: cout << line << std::endl;
                //}

                std::istringstream inp(line.c_str());
                PntVel pv;
                pv.proc = 0;
                inp >> pv.x;
                inp >> pv.y;
                inp >> pv.z;
                double dd;
                for (int i = 0; i < input.Ninfo; i++){
                    inp >> dd;
                    pv.data.push_back(dd);
                }
                inp >> celidstr;
                cellIDint thisCell = cellID2int(celidstr);
                bool tf = CLact.addVelocityPoint(thisCell, pv);
                if (tf) {cnt_data++;}
                //if (!tf){
                //    tf = CLexp.addVelocityPoint(thisCell, pv);
                //    if (tf)
                //        cnt_data++;
                //}
                //else{
                //    cnt_data++;
                //}
            }
        }

        //if (myRank == dgbRank){
        //    std::cout << "Start HERE" << filename <<  std::endl;
        //}

        boost::chrono::duration<double> seconds = boost::chrono::nanoseconds(timer.elapsed().user);
        std::cout << ". . . .Proc " << myRank << " spend " << seconds.count() / 60 << " min to insert " << cnt_data << " data" << std::endl;
        outcome = true;
    }
    return outcome;
}

void printFilesGraphV2(cellIdList& CLact, inputs& input, int rank, DomainListPoly& dpoly){
    boost::timer::cpu_timer timer;
    L1map::iterator itLev1;
    L2map::iterator itLev2;
    L3map::iterator itLev3;
    int cellIdx = 0;
    int velIdx = 0;

    // Add ids to the cells and the velocity points according to the order they are going to be printed
    for (itLev1 = CLact.CellIds.begin(); itLev1 != CLact.CellIds.end(); ++itLev1){
        for (itLev2 = itLev1->second.begin(); itLev2 != itLev1->second.end(); ++itLev2){
            for (itLev3 = itLev2->second.begin(); itLev3 != itLev2->second.end(); ++ itLev3){
                itLev3->second.id = cellIdx;
                cellIdx++;
                for (int i = 0; i < itLev3->second.CellVelocities.size(); ++i){
                    itLev3->second.VelIds.push_back(velIdx);
                    velIdx++;
                }
            }
        }
    }

    // Print the velocity and the Cell graph files
    std::cout << "____ Processor " << rank << " is printing Velocity and Graph files ..." << std::endl;
    std::string graphfilename = input.NewPrefix + num2Padstr(rank, input.Nzeros) + ".grph";
    std::ofstream graphstream;
    graphstream.open(graphfilename.c_str());

    std::string velfilename = input.NewPrefix + num2Padstr(rank, input.Nzeros) + input.suffix;
    std::ofstream velstream;
    velstream.open(velfilename.c_str());

    for (itLev1 = CLact.CellIds.begin(); itLev1 != CLact.CellIds.end(); ++itLev1){
        for (itLev2 = itLev1->second.begin(); itLev2 != itLev1->second.end(); ++itLev2){
            for (itLev3 = itLev2->second.begin(); itLev3 != itLev2->second.end(); ++ itLev3){
                // Print the velocities
                for (unsigned int i = 0; i < itLev3->second.VelIds.size(); ++i){
                    int proc = findProcessorOwner(itLev3->second.CellVelocities[i].x,
                                                  itLev3->second.CellVelocities[i].y,
                                                  dpoly, rank);
                    if (proc < 0){
                        std::cout << "Point "
                                  << itLev3->second.CellVelocities[i].x << " "
                                  << itLev3->second.CellVelocities[i].y
                                  << " is not in any processor polygons" << std::endl;
                    }
                    else{
                        velstream << std::setprecision(3) << std::fixed
                                  << itLev3->second.CellVelocities[i].x << " "
                                  << itLev3->second.CellVelocities[i].y << " "
                                  << itLev3->second.CellVelocities[i].z << " "
                                  << proc << " ";
                        for (unsigned int j = 0; j < input.Nprint; ++j){
                            velstream << std::setprecision(input.prec[j]) << std::fixed
                                      << itLev3->second.CellVelocities[i].data[input.printOrder[j] - 1] << " ";
                        }
                        velstream << std::endl;
                    }
                }

                std::vector<int> tmp;
                for(unsigned int i = 0; i < itLev3->second.neighborCells.size(); ++i){
                    bool tf = CLact.find(itLev3->second.neighborCells[i]);
                    if (tf){
                        tmp.push_back(CLact.itLev3->second.id);
                    }
                }
                // Print the cell graph
                graphstream << std::setprecision(3) << std::fixed
                            << itLev3->second.x << " "
                            << itLev3->second.y << " "
                            << itLev3->second.z << " "
                            << tmp.size() << " "
                            << itLev3->second.CellVelocities.size() << " ";
                for(unsigned int i = 0; i < tmp.size(); ++i){
                    graphstream << tmp[i] << " ";
                }
                for (unsigned int i = 0; i < itLev3->second.VelIds.size(); ++i){
                    graphstream << itLev3->second.VelIds[i] << " ";
                }
                graphstream << std::endl;
            }
        }
    }
    graphstream.close();
    velstream.close();
    boost::chrono::duration<double> seconds = boost::chrono::nanoseconds(timer.elapsed().user);
    std::cout << ". . . .Proc " << rank << " spend " << seconds.count() / 60 << " min to print "
              << cellIdx << " cells and " << velIdx << " velocity points" << std::endl;
}

void printFilesGraph(cellIdList& CLact, inputs& input, int rank, DomainListPoly& dpoly){
    boost::timer::cpu_timer timer;

    L1map::iterator itLev1;
    L2map::iterator itLev2;
    L3map::iterator itLev3;
    int idx = 0;
    //Add ids according to the order they are going to be printed
    for (itLev1 = CLact.CellIds.begin(); itLev1 != CLact.CellIds.end(); ++itLev1){
        for (itLev2 = itLev1->second.begin(); itLev2 != itLev1->second.end(); ++itLev2){
            for (itLev3 = itLev2->second.begin(); itLev3 != itLev2->second.end(); ++ itLev3){
                for (int i = 0; i < itLev3->second.CellVelocities.size(); ++i){
                    itLev3->second.VelIds.push_back(idx);
                    idx++;
                }
            }
        }
    }
    // Add and print the indices of the neighbor cells
    std::cout << "____ Processor " << rank << " is printing Velocity and Graph files ..." << std::endl;
    std::string graphfilename = input.NewPrefix + num2Padstr(rank, input.Nzeros) + ".grph";
    std::ofstream graphstream;
    graphstream.open(graphfilename.c_str());

    std::string velfilename = input.NewPrefix + num2Padstr(rank, input.Nzeros) + input.suffix;
    std::ofstream velstream;
    velstream.open(velfilename.c_str());
    idx = 0;
    for (itLev1 = CLact.CellIds.begin(); itLev1 != CLact.CellIds.end(); ++itLev1){
        for (itLev2 = itLev1->second.begin(); itLev2 != itLev1->second.end(); ++itLev2){
            for (itLev3 = itLev2->second.begin(); itLev3 != itLev2->second.end(); ++ itLev3){
                for (unsigned int i = 0; i < itLev3->second.VelIds.size(); ++i){
                    int proc = findProcessorOwner(itLev3->second.CellVelocities[i].x,
                                       itLev3->second.CellVelocities[i].y,
                                       dpoly, rank);
                    if (proc < 0){
                        std::cout << "Point "
                                << itLev3->second.CellVelocities[i].x << " "
                                << itLev3->second.CellVelocities[i].y
                                << " is not in any processor polygons" << std::endl;
                    }
                    else{
                        velstream << std::setprecision(3) << std::fixed
                                  << itLev3->second.CellVelocities[i].x << " "
                                  << itLev3->second.CellVelocities[i].y << " "
                                  << itLev3->second.CellVelocities[i].z << " "
                                  << proc << " ";
                        for (unsigned int j = 0; j < input.Nprint; ++j){
                            velstream << std::setprecision(input.prec[j]) << std::fixed
                            << itLev3->second.CellVelocities[i].data[input.printOrder[j] - 1] << " ";
                        }
                        velstream << std::endl;
                    }

                    // Make a list of the ids that this velocity node depends on
                    // Add first all the velocity nodes of the
                    std::vector<int> tmp;
                    for (unsigned int j = 0; j < itLev3->second.VelIds.size(); ++j){
                        if (itLev3->second.VelIds[j] == idx)
                            continue;
                        tmp.push_back(itLev3->second.VelIds[j]);
                    }

                    // Add the node ids of the neighboring cells
                    for (unsigned int j = 0; j < itLev3->second.neighborCells.size(); ++j){
                        bool tf = CLact.find(itLev3->second.neighborCells[j]);
                        if (tf){
                            for (int k = 0; k < CLact.itLev3->second.VelIds.size(); ++k){
                                tmp.push_back(CLact.itLev3->second.VelIds[k]);
                            }
                        }
                    }
                    graphstream << tmp.size() << " ";
                    for (unsigned int j = 0; j < tmp.size(); ++j){
                        graphstream << tmp[j] << " ";
                    }
                    graphstream << std::endl;
                    idx++;
                }
            }
        }
    }
    graphstream.close();
    velstream.close();
    boost::chrono::duration<double> seconds = boost::chrono::nanoseconds(timer.elapsed().user);
    std::cout << ". . . .Proc " << rank << " spend " << seconds.count() / 60 << " min to print " << idx << " data" << std::endl;

}

void appendNeighborCells(cellIdList& CLact, cellIdList& CLexp, inputs& input){
    L1map::iterator itLev1;
    L2map::iterator itLev2;
    L3map::iterator itLev3;
    bool foundCell;
    cellIdList CLextra;
    for (int i = 0; i < input.overlapIter; ++i){
        CLextra.clear();
        for (itLev1 = CLact.CellIds.begin(); itLev1 != CLact.CellIds.end(); ++itLev1){
            for (itLev2 = itLev1->second.begin(); itLev2 != itLev1->second.end(); ++itLev2){
                for (itLev3 = itLev2->second.begin(); itLev3 != itLev2->second.end(); ++ itLev3){
                    for (int j = 0; j < itLev3->second.neighborCells.size(); ++j){
                        foundCell = CLact.find(itLev3->second.neighborCells[j]);
                        if (!foundCell){
                            foundCell = CLextra.find(itLev3->second.neighborCells[j]);
                            if (!foundCell){
                                foundCell = CLexp.find(itLev3->second.neighborCells[j]);
                                if (foundCell){
                                    CLextra.insert(itLev3->second.neighborCells[j], CLexp.itLev3->second);
                                }
                                //else{
                                //    std::cout << itLev3->second.neighborCells[j].to_string() << " Not found anywhere" << std::endl;
                                //}
                            }
                        }
                    }
                }
            }
        }

        for (itLev1 = CLextra.CellIds.begin(); itLev1 != CLextra.CellIds.end(); ++itLev1){
            for (itLev2 = itLev1->second.begin(); itLev2 != itLev1->second.end(); ++itLev2){
                for (itLev3 = itLev2->second.begin(); itLev3 != itLev2->second.end(); ++ itLev3){
                    CLact.insert(cellIDint(itLev1->first,itLev2->first, itLev3->first), itLev3->second);
                }
            }
        }
    }

    /*
    std::string outfile = "tmp.dat";
    std::ofstream outstream;
    outstream.open(outfile.c_str());
    for (itLev1 = CLact.CellIds.begin(); itLev1 != CLact.CellIds.end(); ++itLev1){
        for (itLev2 = itLev1->second.begin(); itLev2 != itLev1->second.end(); ++itLev2){
            for (itLev3 = itLev2->second.begin(); itLev3 != itLev2->second.end(); ++ itLev3){
                for (int i = 0; i < itLev3->second.CellVelocities.size(); ++i){
                    outstream << itLev3->second.CellVelocities[i].x << " "
                              << itLev3->second.CellVelocities[i].y << " "
                              << itLev3->second.CellVelocities[i].z << std::endl;
                }
            }
        }
    }
    outstream.close();
    */
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
                        inp >> dd;
                        pv.data.push_back(dd);
                    }
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



int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;

    //if (world.rank() == 0){
        std::cout << "Redist version 1.2" << std::endl;
    //}

    //int DomainId = std::stoi(argv[2]);
    int DomainId = world.rank();
    //std::cout << DomainId << std::endl;

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


    if (inp.UseGraph){
        cellIdList CLactual;
        cellIdList CLbuffer;
        for (int iproc = 0; iproc < inp.NinitDom; ++iproc) {
            std::string filename = inp.prefix + num2Padstr(iproc, inp.Nzeros) + ".grph";
            tf = readCellGraphFiles(filename, CLactual, CLbuffer,
                                    DomainId, expandedDom, actualDom);
        }
        world.barrier();
        appendNeighborCells(CLactual, CLbuffer, inp);
        world.barrier();
        for (int iproc = 0; iproc < inp.NinitDom; ++iproc){
            std::string filename1 = inp.prefix + num2Padstr(iproc, inp.Nzeros) + inp.suffix;
            tf = readVelocityWithGraph(filename1, DomainId, inp, CLactual);
        }
        world.barrier();
        printFilesGraphV2(CLactual, inp, DomainId, actualDom);
    }
    else{
        std::vector< PntVel> my_data;
        for (int iproc = 0; iproc < inp.NinitDom; ++iproc) {
            std::string filename = inp.prefix + num2Padstr(iproc, inp.Nzeros) + inp.suffix;
            tf = readVelocityFiles(filename, my_data, DomainId, expandedDom, actualDom, inp);
        }

        std::cout << "____ Processor " << DomainId << " is printing..." << std::endl;
        std::string outfilename = inp.NewPrefix + num2Padstr(DomainId, inp.Nzeros) + inp.suffix;
        std::ofstream outstream;
        outstream.open(outfilename.c_str());
        for (std::vector<PntVel>::iterator it = my_data.begin(); it != my_data.end(); ++it) {
            outstream << std::setprecision(3) << std::fixed
                      << it->x << " " << it->y << " " << it->z << " " << it->proc << " ";

            for (int i = 0; i < inp.Nprint; i++){
                outstream << std::setprecision(inp.prec[i]) << std::fixed << it->data[ inp.printOrder[i] - 1 ] << " ";
            }
            outstream << std::endl;
        }
        outstream.close();
    }







    world.barrier();
    return 0;
}

