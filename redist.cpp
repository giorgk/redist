// redist.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <istream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <cctype>
#include <stdexcept>

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

int dgbRank = 0;

enum class InputFormat {
    AUTO,
    ASCII,
    BINARY
};

enum class OutputFormat {
    ASCII,
    BINARY
};
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

bool parseCellIdInteger(const std::string& text, int& value) {
    std::istringstream input(text.c_str());
    char extra = 0;
    if (!(input >> value))
        return false;
    return !(input >> extra);
}

cellIDint cellID2int(const std::string& str){
    std::vector<std::string> parts;
    std::string::size_type begin = 0;
    for (std::string::size_type i = 0; i <= str.size(); ++i) {
        if (i == str.size() || str[i] == ':' || str[i] == '_') {
            parts.push_back(str.substr(begin, i - begin));
            begin = i + 1;
        }
    }

    if (parts.size() < 2)
        throw std::runtime_error("Invalid cell ID: " + str);

    int a = -1;
    int b = -1;
    int c = -9;
    const bool haveA = parseCellIdInteger(parts[0], a);
    const bool haveB = parseCellIdInteger(parts[1], b);
    bool haveC = true;
    if (b != 0)
        haveC = parts.size() >= 3 && parseCellIdInteger(parts[2], c);

    if (!haveA || !haveB || !haveC)
        throw std::runtime_error("Invalid cell ID: " + str);

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
    // Optional input format appended to the legacy first configuration line.
    // AUTO preserves legacy files and selects binary input for *.bin suffixes.
    InputFormat inputFormat = InputFormat::AUTO;
    // Redistributed output encoding. ASCII remains the default so that all
    // legacy input files preserve their behavior.
    OutputFormat outputFormat = OutputFormat::ASCII;
    // Output velocity suffix. Its default depends on outputFormat.
    std::string outputSuffix;
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

    bool UseGraph = false;
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

std::string lowerString(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

bool endsWith(const std::string& value, const std::string& ending) {
    return value.size() >= ending.size() &&
           value.compare(value.size() - ending.size(), ending.size(), ending) == 0;
}

InputFormat resolvedInputFormat(const inputs& input) {
    if (input.inputFormat != InputFormat::AUTO)
        return input.inputFormat;
    return endsWith(lowerString(input.suffix), ".bin")
               ? InputFormat::BINARY
               : InputFormat::ASCII;
}

int findProcessorOwner(double x, double y, const DomainListPoly& actual, int rank){
    DomainListPoly::const_iterator rankDomain = actual.find(rank);
    if (rankDomain != actual.end() &&
        boost::geometry::within(bpoint(x, y), rankDomain->second)){
        return rank;
    }
    else{
        for (DomainListPoly::const_iterator it = actual.begin();
             it != actual.end(); ++it){
            if (it->first == rank)
                continue;
            if (boost::geometry::within(bpoint(x, y), it->second)){
                return it->first;
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
            int tmp = 0;
            if (inp >> tmp) {
                if (tmp != 0 ){
                    input.UseGraph = true;
                    if (!(inp >> tmp)) {
                        std::cout << "The graph overlap iteration count is missing."
                                  << std::endl;
                        return false;
                    }
                    input.overlapIter = tmp;
                }
            }
            else
                inp.clear(); // Legacy four-field line: graph input is disabled.

            // Backward compatible optional fields:
            //   ... useGraph [overlapIter]
            //       [inputFormat [outputSuffix] [outputFormat]]
            // If the token following inputFormat is ascii or binary, it is
            // treated as outputFormat and the output suffix is derived.
            std::string format;
            if (inp >> format) {
                format = lowerString(format);
                if (format == "ascii")
                    input.inputFormat = InputFormat::ASCII;
                else if (format == "binary")
                    input.inputFormat = InputFormat::BINARY;
                else if (format == "auto")
                    input.inputFormat = InputFormat::AUTO;
                else {
                    std::cout << "Invalid input format '" << format
                              << "'. Expected ascii, binary, or auto." << std::endl;
                    return false;
                }

                std::string outputToken;
                if (inp >> outputToken) {
                    const std::string lowerOutputToken = lowerString(outputToken);
                    if (lowerOutputToken == "ascii")
                        input.outputFormat = OutputFormat::ASCII;
                    else if (lowerOutputToken == "binary")
                        input.outputFormat = OutputFormat::BINARY;
                    else {
                        input.outputSuffix = outputToken;
                        std::string outputFormat;
                        if (inp >> outputFormat) {
                            outputFormat = lowerString(outputFormat);
                            if (outputFormat == "ascii")
                                input.outputFormat = OutputFormat::ASCII;
                            else if (outputFormat == "binary")
                                input.outputFormat = OutputFormat::BINARY;
                            else {
                                std::cout << "Invalid output format '" << outputFormat
                                          << "'. Expected ascii or binary." << std::endl;
                                return false;
                            }
                        }
                    }
                }
            }

            if (input.outputSuffix.empty()) {
                input.outputSuffix = input.suffix;
                if (input.outputFormat == OutputFormat::BINARY) {
                    if (!endsWith(lowerString(input.outputSuffix), ".bin"))
                        input.outputSuffix += ".bin";
                }
                else if (resolvedInputFormat(input) == InputFormat::BINARY &&
                         endsWith(lowerString(input.outputSuffix), ".bin")) {
                    input.outputSuffix.erase(input.outputSuffix.size() - 4);
                }
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

template <typename T>
bool readPod(std::istream& input, T& value) {
    input.read(reinterpret_cast<char*>(&value), sizeof(T));
    return input.good();
}

template <typename T>
bool writePod(std::ostream& output, const T& value) {
    output.write(reinterpret_cast<const char*>(&value), sizeof(T));
    return output.good();
}

bool writeMagic(std::ostream& output, const char magic[8]) {
    output.write(magic, 8);
    return output.good();
}

bool readBinaryString(std::istream& input, std::string& value) {
    std::uint32_t length = 0;
    if (!readPod(input, length))
        return false;
    value.resize(length);
    if (length > 0)
        input.read(&value[0], static_cast<std::streamsize>(length));
    return input.good();
}

bool readAndCheckMagic(std::istream& input,
                       const char expected[8],
                       const std::string& filename) {
    char actual[8] = {0,0,0,0,0,0,0,0};
    input.read(actual, 8);
    if (!input.good() || std::memcmp(actual, expected, 8) != 0) {
        std::cout << "Invalid binary file magic in " << filename << std::endl;
        return false;
    }
    return true;
}

struct BinaryVelocityHeader {
    std::uint32_t dimension = 0;
    std::uint32_t sourceRank = 0;
    std::uint32_t sourceNproc = 0;
    std::uint32_t pointsPerCell = 0;
    std::uint32_t faceMode = 0;
    std::uint64_t nCells = 0;
    std::uint64_t nPoints = 0;
};

bool readBinaryVelocityHeader(std::ifstream& input,
                              const std::string& filename,
                              const int expectedRank,
                              const int expectedNproc,
                              BinaryVelocityHeader& header) {
    const char expectedMagic[8] = {'N','P','S','I','V','L','1','\0'};
    if (!readAndCheckMagic(input, expectedMagic, filename))
        return false;

    std::uint32_t version = 0;
    if (!readPod(input, version) ||
        !readPod(input, header.dimension) ||
        !readPod(input, header.sourceRank) ||
        !readPod(input, header.sourceNproc) ||
        !readPod(input, header.pointsPerCell) ||
        !readPod(input, header.faceMode) ||
        !readPod(input, header.nCells) ||
        !readPod(input, header.nPoints)) {
        std::cout << "Incomplete binary velocity header in " << filename << std::endl;
        return false;
    }

    if (version != 1 || header.dimension != 3 ||
        header.sourceRank != static_cast<std::uint32_t>(expectedRank) ||
        header.sourceNproc != static_cast<std::uint32_t>(expectedNproc) ||
        (header.pointsPerCell != 1 && header.pointsPerCell != 8) ||
        header.faceMode > 1) {
        std::cout << "Unsupported or inconsistent binary velocity header in "
                  << filename << std::endl;
        return false;
    }
    return true;
}

bool readBinaryVelocityRecord(std::ifstream& input,
                              PntVel& point,
                              std::string& cellId) {
    double vx = 0.0, vy = 0.0, vz = 0.0;
    double diameter = 0.0, ratio = 0.0;
    std::uint32_t sourceRank = 0;
    if (!readPod(input, point.x) ||
        !readPod(input, point.y) ||
        !readPod(input, point.z) ||
        !readPod(input, vx) ||
        !readPod(input, vy) ||
        !readPod(input, vz) ||
        !readPod(input, sourceRank) ||
        !readPod(input, diameter) ||
        !readPod(input, ratio) ||
        !readBinaryString(input, cellId))
        return false;

    point.proc = static_cast<int>(sourceRank);
    point.data.clear();
    // Preserve the legacy NPSAT ASCII logical column order:
    // vx vy vz source_rank diameter ratio.
    point.data.push_back(vx);
    point.data.push_back(vy);
    point.data.push_back(vz);
    point.data.push_back(static_cast<double>(sourceRank));
    point.data.push_back(diameter);
    point.data.push_back(ratio);
    return true;
}

bool readCellGraphFilesBinary(std::string filename,
                              cellIdList& CLact,
                              cellIdList& CLexp,
                              int myrank,
                              DomainListPoly& expanded,
                              DomainListPoly& actual,
                              int expectedSourceRank,
                              int expectedNproc) {
    std::ifstream datafile(filename.c_str(), std::ios::binary);
    if (!datafile.good()) {
        std::cout << "Can't open the file " << filename << std::endl;
        return false;
    }

    const char expectedMagic[8] = {'N','P','S','I','G','R','1','\0'};
    if (!readAndCheckMagic(datafile, expectedMagic, filename))
        return false;

    std::uint32_t version = 0, dimension = 0, sourceRank = 0, sourceNproc = 0;
    std::uint64_t nCells = 0;
    if (!readPod(datafile, version) ||
        !readPod(datafile, dimension) ||
        !readPod(datafile, sourceRank) ||
        !readPod(datafile, sourceNproc) ||
        !readPod(datafile, nCells)) {
        std::cout << "Incomplete binary graph header in " << filename << std::endl;
        return false;
    }
    if (version != 1 || dimension != 3 ||
        sourceRank != static_cast<std::uint32_t>(expectedSourceRank) ||
        sourceNproc != static_cast<std::uint32_t>(expectedNproc)) {
        std::cout << "Unsupported or inconsistent binary graph header in "
                  << filename << std::endl;
        return false;
    }

    std::cout << "Processor " << myrank << " reads binary graph file "
              << filename << std::endl;
    boost::timer::cpu_timer timer;
    int nAct = 0;
    int nExp = 0;

    for (std::uint64_t record = 0; record < nCells; ++record) {
        double x = 0.0, y = 0.0, z = 0.0;
        std::uint32_t nNeighbors = 0;
        std::string cellId;
        if (!readPod(datafile, x) || !readPod(datafile, y) ||
            !readPod(datafile, z) || !readPod(datafile, nNeighbors) ||
            !readBinaryString(datafile, cellId)) {
            std::cout << "Incomplete binary graph record in " << filename << std::endl;
            return false;
        }

        cellData cdata;
        cdata.x = x;
        cdata.y = y;
        cdata.z = z;
        for (std::uint32_t i = 0; i < nNeighbors; ++i) {
            std::string neighborId;
            if (!readBinaryString(datafile, neighborId)) {
                std::cout << "Incomplete binary graph neighbor list in "
                          << filename << std::endl;
                return false;
            }
            cdata.neighborCells.push_back(cellID2int(neighborId));
        }

        if (boost::geometry::within(bpoint(x, y), actual[myrank])) {
            CLact.insert(cellID2int(cellId), cdata);
            ++nAct;
        }
        else if (boost::geometry::within(bpoint(x, y), expanded[myrank])) {
            CLexp.insert(cellID2int(cellId), cdata);
            ++nExp;
        }
    }

    boost::chrono::duration<double> seconds =
        boost::chrono::nanoseconds(timer.elapsed().user);
    std::cout << ". . . .Proc " << myrank << " spend "
              << seconds.count() / 60 << " min to read "
              << nAct << " " << nExp << std::endl;
    return true;
}

bool readVelocityWithGraphBinary(std::string filename,
                                 int myRank,
                                 inputs& input,
                                 cellIdList& CLact,
                                 int expectedSourceRank) {
    if (input.Ninfo != 6) {
        std::cout << "Binary NPSAT velocity files provide exactly 6 logical data columns "
                     "(vx vy vz source_rank diameter ratio); Ndata must be 6."
                  << std::endl;
        return false;
    }

    std::ifstream datafile(filename.c_str(), std::ios::binary);
    if (!datafile.good()) {
        std::cout << "Can't open the file " << filename << std::endl;
        return false;
    }

    BinaryVelocityHeader header;
    if (!readBinaryVelocityHeader(datafile, filename, expectedSourceRank,
                                  input.NinitDom, header))
        return false;

    std::cout << "Processor " << myRank << " reads binary velocity file "
              << filename << std::endl;
    boost::timer::cpu_timer timer;
    int inserted = 0;
    for (std::uint64_t record = 0; record < header.nPoints; ++record) {
        PntVel point;
        std::string cellId;
        if (!readBinaryVelocityRecord(datafile, point, cellId)) {
            std::cout << "Incomplete binary velocity record in " << filename << std::endl;
            return false;
        }
        if (CLact.addVelocityPoint(cellID2int(cellId), point))
            ++inserted;
    }

    boost::chrono::duration<double> seconds =
        boost::chrono::nanoseconds(timer.elapsed().user);
    std::cout << ". . . .Proc " << myRank << " spend "
              << seconds.count() / 60 << " min to insert "
              << inserted << " data" << std::endl;
    return true;
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
        int cnt_data = 0;
        std::string celidstr;
        while (getline(datafile, line)){
            if (line.size() > 1){
                //if (myRank == dgbRank){
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

bool writeRedistributedVelocityHeader(std::ofstream& output,
                                      const inputs& input,
                                      int rank,
                                      int nproc,
                                      std::uint64_t nPoints) {
    const char magic[8] = {'R','D','I','V','E','L','1','\0'};
    const std::uint32_t version = 1;
    const std::uint32_t dimension = 3;
    const std::uint32_t outputRank = static_cast<std::uint32_t>(rank);
    const std::uint32_t outputNproc = static_cast<std::uint32_t>(nproc);
    const std::uint32_t nData = static_cast<std::uint32_t>(input.Nprint);

    if (!writeMagic(output, magic) ||
        !writePod(output, version) ||
        !writePod(output, dimension) ||
        !writePod(output, outputRank) ||
        !writePod(output, outputNproc) ||
        !writePod(output, nData) ||
        !writePod(output, nPoints))
        return false;

    // Preserve the one-based source-column identifiers from line 2 of the
    // configuration file so that a binary reader can identify the payload.
    for (int i = 0; i < input.Nprint; ++i) {
        const std::uint32_t sourceColumn =
            static_cast<std::uint32_t>(input.printOrder[i]);
        if (!writePod(output, sourceColumn))
            return false;
    }
    return true;
}

bool writeRedistributedVelocityRecord(std::ofstream& output,
                                      const PntVel& point,
                                      int owner,
                                      const inputs& input) {
    const std::int32_t outputOwner = static_cast<std::int32_t>(owner);
    if (!writePod(output, point.x) ||
        !writePod(output, point.y) ||
        !writePod(output, point.z) ||
        !writePod(output, outputOwner))
        return false;

    for (int i = 0; i < input.Nprint; ++i) {
        const int sourceIndex = input.printOrder[i] - 1;
        if (sourceIndex < 0 ||
            sourceIndex >= static_cast<int>(point.data.size())) {
            std::cout << "Output column " << input.printOrder[i]
                      << " is outside the available velocity data." << std::endl;
            return false;
        }
        if (!writePod(output, point.data[sourceIndex]))
            return false;
    }
    return true;
}

bool writeRedistributedVelocityBinary(const std::string& filename,
                                      const std::vector<PntVel>& points,
                                      const inputs& input,
                                      int rank,
                                      int nproc) {
    std::ofstream output(filename.c_str(),
                         std::ios::out | std::ios::binary | std::ios::trunc);
    if (!output.good()) {
        std::cout << "Can't open the binary output file " << filename << std::endl;
        return false;
    }

    if (!writeRedistributedVelocityHeader(
            output, input, rank, nproc,
            static_cast<std::uint64_t>(points.size())))
        return false;

    for (std::vector<PntVel>::const_iterator it = points.begin();
         it != points.end(); ++it) {
        if (!writeRedistributedVelocityRecord(output, *it, it->proc, input))
            return false;
    }

    output.close();
    if (!output.good()) {
        std::cout << "Failed while writing binary velocity file "
                  << filename << std::endl;
        return false;
    }
    return true;
}

struct BinaryGraphVelocity {
    const PntVel* point;
    int owner;
};

bool printFilesGraphBinary(cellIdList& CLact,
                           inputs& input,
                           int rank,
                           int nproc,
                           DomainListPoly& dpoly) {
    boost::timer::cpu_timer timer;
    std::vector<BinaryGraphVelocity> velocityRecords;
    std::uint64_t nCells = 0;

    // Assign compact local cell and velocity IDs. Points outside all actual
    // polygons are omitted from both binary files so graph references remain
    // consistent with the velocity record count.
    for (L1map::iterator itLev1 = CLact.CellIds.begin();
         itLev1 != CLact.CellIds.end(); ++itLev1) {
        for (L2map::iterator itLev2 = itLev1->second.begin();
             itLev2 != itLev1->second.end(); ++itLev2) {
            for (L3map::iterator itLev3 = itLev2->second.begin();
                 itLev3 != itLev2->second.end(); ++itLev3) {
                itLev3->second.id = static_cast<int>(nCells++);
                itLev3->second.VelIds.clear();
                for (std::size_t i = 0;
                     i < itLev3->second.CellVelocities.size(); ++i) {
                    PntVel& point = itLev3->second.CellVelocities[i];
                    const int owner =
                        findProcessorOwner(point.x, point.y, dpoly, rank);
                    if (owner < 0) {
                        std::cout << "Point " << point.x << " " << point.y
                                  << " is not in any processor polygons"
                                  << std::endl;
                        continue;
                    }
                    itLev3->second.VelIds.push_back(
                        static_cast<int>(velocityRecords.size()));
                    BinaryGraphVelocity record;
                    record.point = &point;
                    record.owner = owner;
                    velocityRecords.push_back(record);
                }
            }
        }
    }

    const std::string velocityFilename =
        input.NewPrefix + num2Padstr(rank, input.Nzeros) + input.outputSuffix;
    const std::string graphFilename =
        input.NewPrefix + num2Padstr(rank, input.Nzeros) + ".grph.bin";

    std::ofstream velocityStream(
        velocityFilename.c_str(),
        std::ios::out | std::ios::binary | std::ios::trunc);
    if (!velocityStream.good()) {
        std::cout << "Can't open the binary output file "
                  << velocityFilename << std::endl;
        return false;
    }
    if (!writeRedistributedVelocityHeader(
            velocityStream, input, rank, nproc,
            static_cast<std::uint64_t>(velocityRecords.size())))
        return false;
    for (std::size_t i = 0; i < velocityRecords.size(); ++i) {
        if (!writeRedistributedVelocityRecord(
                velocityStream, *velocityRecords[i].point,
                velocityRecords[i].owner, input))
            return false;
    }
    velocityStream.close();
    if (!velocityStream.good()) {
        std::cout << "Failed while writing binary velocity file "
                  << velocityFilename << std::endl;
        return false;
    }

    std::ofstream graphStream(
        graphFilename.c_str(),
        std::ios::out | std::ios::binary | std::ios::trunc);
    if (!graphStream.good()) {
        std::cout << "Can't open the binary output file "
                  << graphFilename << std::endl;
        return false;
    }

    const char graphMagic[8] = {'R','D','I','G','R','P','1','\0'};
    const std::uint32_t version = 1;
    const std::uint32_t dimension = 3;
    const std::uint32_t outputRank = static_cast<std::uint32_t>(rank);
    const std::uint32_t outputNproc = static_cast<std::uint32_t>(nproc);
    const std::uint64_t nVelocityPoints =
        static_cast<std::uint64_t>(velocityRecords.size());
    if (!writeMagic(graphStream, graphMagic) ||
        !writePod(graphStream, version) ||
        !writePod(graphStream, dimension) ||
        !writePod(graphStream, outputRank) ||
        !writePod(graphStream, outputNproc) ||
        !writePod(graphStream, nCells) ||
        !writePod(graphStream, nVelocityPoints))
        return false;

    for (L1map::iterator itLev1 = CLact.CellIds.begin();
         itLev1 != CLact.CellIds.end(); ++itLev1) {
        for (L2map::iterator itLev2 = itLev1->second.begin();
             itLev2 != itLev1->second.end(); ++itLev2) {
            for (L3map::iterator itLev3 = itLev2->second.begin();
                 itLev3 != itLev2->second.end(); ++itLev3) {
                std::vector<std::uint64_t> neighborIds;
                for (std::size_t i = 0;
                     i < itLev3->second.neighborCells.size(); ++i) {
                    if (CLact.find(itLev3->second.neighborCells[i]))
                        neighborIds.push_back(
                            static_cast<std::uint64_t>(CLact.itLev3->second.id));
                }

                const std::uint32_t nNeighbors =
                    static_cast<std::uint32_t>(neighborIds.size());
                const std::uint32_t nCellVelocities =
                    static_cast<std::uint32_t>(itLev3->second.VelIds.size());
                if (!writePod(graphStream, itLev3->second.x) ||
                    !writePod(graphStream, itLev3->second.y) ||
                    !writePod(graphStream, itLev3->second.z) ||
                    !writePod(graphStream, nNeighbors) ||
                    !writePod(graphStream, nCellVelocities))
                    return false;

                for (std::size_t i = 0; i < neighborIds.size(); ++i) {
                    if (!writePod(graphStream, neighborIds[i]))
                        return false;
                }
                for (std::size_t i = 0;
                     i < itLev3->second.VelIds.size(); ++i) {
                    const std::uint64_t velocityId =
                        static_cast<std::uint64_t>(itLev3->second.VelIds[i]);
                    if (!writePod(graphStream, velocityId))
                        return false;
                }
            }
        }
    }

    graphStream.close();
    if (!graphStream.good()) {
        std::cout << "Failed while writing binary graph file "
                  << graphFilename << std::endl;
        return false;
    }

    boost::chrono::duration<double> seconds =
        boost::chrono::nanoseconds(timer.elapsed().user);
    std::cout << ". . . .Proc " << rank << " spend "
              << seconds.count() / 60 << " min to print "
              << nCells << " binary cells and " << velocityRecords.size()
              << " binary velocity points" << std::endl;
    return true;
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
                for (std::size_t i = 0;
                     i < itLev3->second.CellVelocities.size(); ++i){
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

    std::string velfilename = input.NewPrefix + num2Padstr(rank, input.Nzeros) + input.outputSuffix;
    std::ofstream velstream;
    velstream.open(velfilename.c_str());

    for (itLev1 = CLact.CellIds.begin(); itLev1 != CLact.CellIds.end(); ++itLev1){
        for (itLev2 = itLev1->second.begin(); itLev2 != itLev1->second.end(); ++itLev2){
            for (itLev3 = itLev2->second.begin(); itLev3 != itLev2->second.end(); ++ itLev3){
                // Print the velocities
                for (std::size_t i = 0; i < itLev3->second.VelIds.size(); ++i){
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
                        for (int j = 0; j < input.Nprint; ++j){
                            velstream << std::setprecision(input.prec[j]) << std::fixed
                                      << itLev3->second.CellVelocities[i].data[input.printOrder[j] - 1] << " ";
                        }
                        velstream << std::endl;
                    }
                }

                std::vector<int> tmp;
                for(std::size_t i = 0; i < itLev3->second.neighborCells.size(); ++i){
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
                for(std::size_t i = 0; i < tmp.size(); ++i){
                    graphstream << tmp[i] << " ";
                }
                for (std::size_t i = 0; i < itLev3->second.VelIds.size(); ++i){
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
                for (std::size_t i = 0;
                     i < itLev3->second.CellVelocities.size(); ++i){
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

    std::string velfilename = input.NewPrefix + num2Padstr(rank, input.Nzeros) + input.outputSuffix;
    std::ofstream velstream;
    velstream.open(velfilename.c_str());
    idx = 0;
    for (itLev1 = CLact.CellIds.begin(); itLev1 != CLact.CellIds.end(); ++itLev1){
        for (itLev2 = itLev1->second.begin(); itLev2 != itLev1->second.end(); ++itLev2){
            for (itLev3 = itLev2->second.begin(); itLev3 != itLev2->second.end(); ++ itLev3){
                for (std::size_t i = 0; i < itLev3->second.VelIds.size(); ++i){
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
                        for (int j = 0; j < input.Nprint; ++j){
                            velstream << std::setprecision(input.prec[j]) << std::fixed
                            << itLev3->second.CellVelocities[i].data[input.printOrder[j] - 1] << " ";
                        }
                        velstream << std::endl;
                    }

                    // Make a list of the ids that this velocity node depends on
                    // Add first all the velocity nodes of the
                    std::vector<int> tmp;
                    for (std::size_t j = 0; j < itLev3->second.VelIds.size(); ++j){
                        if (itLev3->second.VelIds[j] == idx)
                            continue;
                        tmp.push_back(itLev3->second.VelIds[j]);
                    }

                    // Add the node ids of the neighboring cells
                    for (std::size_t j = 0; j < itLev3->second.neighborCells.size(); ++j){
                        bool tf = CLact.find(itLev3->second.neighborCells[j]);
                        if (tf){
                            for (std::size_t k = 0;
                                 k < CLact.itLev3->second.VelIds.size(); ++k){
                                tmp.push_back(CLact.itLev3->second.VelIds[k]);
                            }
                        }
                    }
                    graphstream << tmp.size() << " ";
                    for (std::size_t j = 0; j < tmp.size(); ++j){
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
                    for (std::size_t j = 0;
                         j < itLev3->second.neighborCells.size(); ++j){
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

bool readVelocityFilesBinary(std::string filename,
                             std::vector<PntVel>& points,
                             int myRank,
                             DomainListPoly& expanded,
                             DomainListPoly& actual,
                             inputs& input,
                             int expectedSourceRank) {
    if (input.Ninfo != 6) {
        std::cout << "Binary NPSAT velocity files provide exactly 6 logical data columns "
                     "(vx vy vz source_rank diameter ratio); Ndata must be 6."
                  << std::endl;
        return false;
    }

    std::ifstream datafile(filename.c_str(), std::ios::binary);
    if (!datafile.good()) {
        std::cout << "Can't open the file " << filename << std::endl;
        return false;
    }

    BinaryVelocityHeader header;
    if (!readBinaryVelocityHeader(datafile, filename, expectedSourceRank,
                                  input.NinitDom, header))
        return false;

    std::cout << "Processor " << myRank << " reads binary velocity file "
              << filename << std::endl;
    boost::timer::cpu_timer timer;
    std::uint64_t selected = 0;

    for (std::uint64_t record = 0; record < header.nPoints; ++record) {
        PntVel point;
        std::string cellId;
        if (!readBinaryVelocityRecord(datafile, point, cellId)) {
            std::cout << "Incomplete binary velocity record in " << filename << std::endl;
            return false;
        }

        if (!boost::geometry::within(bpoint(point.x, point.y), expanded[myRank]))
            continue;

        int proc = -9;
        if (boost::geometry::within(bpoint(point.x, point.y), actual[myRank]))
            proc = myRank;
        else {
            DomainListPoly::iterator it = actual.begin();
            for (; it != actual.end(); ++it) {
                if (it->first == myRank)
                    continue;
                if (boost::geometry::within(bpoint(point.x, point.y), it->second)) {
                    proc = it->first;
                    break;
                }
            }
        }

        if (proc == -9) {
            std::cout << "Cant find actual domain for point ("
                      << point.x << "," << point.y
                      << ") in my expanded domain with id " << myRank << std::endl;
        }
        point.proc = proc;
        points.push_back(point);
        ++selected;
    }

    boost::chrono::duration<double> seconds =
        boost::chrono::nanoseconds(timer.elapsed().user);
    std::cout << ". . . .Proc " << myRank << " spend "
              << seconds.count() / 60 << " min to select "
              << selected << " data" << std::endl;
    return true;
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
        std::cout << "Redist version 1.4" << std::endl;
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

    const InputFormat inputFormat = resolvedInputFormat(inp);
    std::cout << "Input format: "
              << (inputFormat == InputFormat::BINARY ? "binary" : "ascii")
              << "; output format: "
              << (inp.outputFormat == OutputFormat::BINARY ? "binary" : "ascii")
              << "; redistributed Ichnos suffix: " << inp.outputSuffix
              << std::endl;


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
            const std::string graphSuffix =
                inputFormat == InputFormat::BINARY ? ".grph.bin" : ".grph";
            std::string filename = inp.prefix + num2Padstr(iproc, inp.Nzeros) + graphSuffix;
            if (inputFormat == InputFormat::BINARY)
                tf = readCellGraphFilesBinary(filename, CLactual, CLbuffer,
                                              DomainId, expandedDom, actualDom,
                                              iproc, inp.NinitDom);
            else
                tf = readCellGraphFiles(filename, CLactual, CLbuffer,
                                        DomainId, expandedDom, actualDom);
            if (!tf)
                return 0;
        }
        world.barrier();
        appendNeighborCells(CLactual, CLbuffer, inp);
        world.barrier();
        for (int iproc = 0; iproc < inp.NinitDom; ++iproc){
            std::string filename1 = inp.prefix + num2Padstr(iproc, inp.Nzeros) + inp.suffix;
            if (inputFormat == InputFormat::BINARY)
                tf = readVelocityWithGraphBinary(filename1, DomainId, inp,
                                                 CLactual, iproc);
            else
                tf = readVelocityWithGraph(filename1, DomainId, inp, CLactual);
            if (!tf)
                return 0;
        }
        world.barrier();
        if (inp.outputFormat == OutputFormat::BINARY)
            tf = printFilesGraphBinary(CLactual, inp, DomainId,
                                       world.size(), actualDom);
        else {
            printFilesGraphV2(CLactual, inp, DomainId, actualDom);
            tf = true;
        }
        if (!tf)
            return 0;
    }
    else{
        std::vector< PntVel> my_data;
        for (int iproc = 0; iproc < inp.NinitDom; ++iproc) {
            std::string filename = inp.prefix + num2Padstr(iproc, inp.Nzeros) + inp.suffix;
            if (inputFormat == InputFormat::BINARY)
                tf = readVelocityFilesBinary(filename, my_data, DomainId,
                                             expandedDom, actualDom, inp, iproc);
            else
                tf = readVelocityFiles(filename, my_data, DomainId,
                                       expandedDom, actualDom, inp);
            if (!tf)
                return 0;
        }

        std::cout << "____ Processor " << DomainId << " is printing..." << std::endl;
        std::string outfilename = inp.NewPrefix + num2Padstr(DomainId, inp.Nzeros) + inp.outputSuffix;
        if (inp.outputFormat == OutputFormat::BINARY) {
            tf = writeRedistributedVelocityBinary(
                outfilename, my_data, inp, DomainId, world.size());
            if (!tf)
                return 0;
        }
        else {
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
    }







    world.barrier();
    return 0;
}

