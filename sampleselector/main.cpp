//
//  main.cpp
//  sampleselector
//
//  Given a set of geo locations find the maximum set of locations
//  where each location is unique within a defined buffer zone
//
//  Created by Robert Fleischmann on 12/07/2015.
//  Copyright (c) 2015 Robert Fleischmann. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include "math.h"
#include <iomanip>
#include <random>
#include <sys/types.h>
#include <sys/stat.h>
#include <map>
#include <set>
#include <vector>

using namespace std;

// source and destination paths
static const string SOURCE_FILE = "/Users/rendro/Downloads/data.csv";
static const string DEST_FOLDER = "/Users/rendro/Downloads/result";

// min distance for buffer zone
static const double MIN_DISTANCE = 0.15;
// number of sets to create
static const int ITERATIONS = 1e4;
// min set size for exporting
static const int MIN_SET_SIZE = 67;
// write result csv files for all found sets that have at least MIN_SET_SIZE datapoints
static const bool WRITE_FILES = false;

// datapoint with id (line number or generic id), latitude, longitude and a buffer zone
struct datapoint {
    int id;
    string name;
    double lat;
    double lng;
    set<int> buffer_ids;
};

// define dataset as a set of datapoints
typedef map<int, datapoint> dataset;

inline bool operator < (const datapoint& lhs, const datapoint& rhs) {
    return lhs.id < rhs.id;
}

inline bool operator == (const datapoint& lhs, const datapoint& rhs) {
    return lhs.id == rhs.id;
}

// haversine formula to caluclate the distance between two datapoints in km
double haversine(const datapoint& start, const datapoint& end) {
    static const double EARTH_RADIUS = 6371.0;
    static const double DEG2RAD = M_PI / 180;
    double dlng = DEG2RAD * (end.lng - start.lng);
    double dlat = DEG2RAD * (end.lat - start.lat);
    double a = pow(sin(dlat/2), 2) + cos(DEG2RAD * start.lat) * cos(DEG2RAD * end.lat) * pow(sin(dlng/2), 2);
    double b = 2 * atan2(sqrt(a), sqrt(1 - a));
    return EARTH_RADIUS * b;
}

set<int> getBufferIDsForDatapoint(dataset::iterator &item, dataset data) {
    set<int> buffer_ids;
    for (auto it : data) {
        if (item->first != it.first) {
            if (haversine(item->second, it.second) <= MIN_DISTANCE) {
                buffer_ids.insert(it.first);
            }
        }
    }
    return buffer_ids;
}

// read in a csv (columns [0, 1, 2, 3] match [id, name, lat, lng] of a datapoint
dataset csvToDataset(string file) {
    dataset data;
    
    ifstream infile( file );
    
    while ( infile ) {
        string s;
        if ( !getline( infile, s ) ) break;
        
        istringstream ss( s );
        
        vector <string> record;
        
        while (ss) {
            string s;
            if ( !getline( ss, s, ',' ) ) break;
            record.push_back( s );
        }
        
        set<int> buffer_ids;
        int id = stoi( record[0] );
        datapoint datapoint = {
            id,
            record[1],
            stod( record[2] ),
            stod( record[3] ),
            buffer_ids
        };
        
        data.insert({id, datapoint});
    }
    if ( !infile.eof() ) {
        cerr << "NO EOF!" << endl;
    }
    
    for (auto it = data.begin(); it != data.end(); ++it) {
        it->second.buffer_ids = getBufferIDsForDatapoint(it, data);
    }
    
    return data;
}

// generate a result set from two sets of datapoints of which the first set contains all
// datapoints with other datapoints in the buffer zone and the of which the second set
// contains all datapoints without other datapoints in the buffer zone
set<int> generateSet(set<int> &standalonePoints, set<int> &withNearbyPoints, dataset &datapoints) {
    random_device rd;
    mt19937 rng(rd());
    // create copies
    set<int> remainingPoints(withNearbyPoints.begin(), withNearbyPoints.end());
    set<int> result(standalonePoints.begin(), standalonePoints.end());

    while (remainingPoints.size() != 0) {
        auto it = remainingPoints.begin();
        // generate random index
        uniform_int_distribution<int> uni(0, (int)remainingPoints.size());
        int r = uni(rng);
        // pick random datapoint by advancing the iterator to the random position
        advance(it, r % remainingPoints.size());
        // add picked datapoint to result list
        result.insert(*it);
        // remove all datapoints within buffer zone if still in remaining dataset
        set<int> buffer_ids = datapoints.at(*it).buffer_ids;
        for (auto j : buffer_ids) {
            auto tmp = remainingPoints.find(j);
            if (tmp != remainingPoints.end()) {
                remainingPoints.erase(tmp);
            }
        }
        // remove picked datapoint from remaining list
        remainingPoints.erase(remainingPoints.find(*it));
    }
    return result;
}

// main sampleselector
int main(int argc, const char * argv[]) {
    // read in data
    dataset datapoints = csvToDataset(SOURCE_FILE);
    set<int> standalonePoints;
    set<int> withNearbyPoints;
    
    for (auto it : datapoints) {
        it.second.buffer_ids.size() == 0
            ? standalonePoints.insert(it.first)
            : withNearbyPoints.insert(it.first);
    }

    cout << "# Datapoints with datapoints in buffer zone   : " << withNearbyPoints.size() << endl;
    cout << "# Datapoints with no datapoint in buffer zone : " << standalonePoints.size() << endl;
    cout << endl;

    // create result sets
    set<set<int>> results;
    int i = 0;
    double biggestSet = 0;
    while (i < ITERATIONS) {
        i++;
        set<int> resultSet = generateSet(standalonePoints, withNearbyPoints, datapoints);

        if (i % 1000 == 0) {
            cout << "Iteration #" << to_string(i) << endl;
        }

        if (resultSet.size() >= MIN_SET_SIZE) {
            results.insert(resultSet);
//            cout << "Iteration #" << to_string(i) << " => " << resultSet.size() << " datapoints" << endl;
        }
        biggestSet = biggestSet > resultSet.size() ? biggestSet : resultSet.size();
    }

    cout << endl;
    cout << "Unique result sets found: " << results.size() << endl;
    cout << "Biggest result set found: " << biggestSet << endl;
    
    if (WRITE_FILES) {
        cout << "Writing files..." << endl;
        int filenr = 0;
        map<unsigned long, int> setsPerLength;
        for (auto i : results) {
            char path [256];
            sprintf(path, "%s/%lu", DEST_FOLDER.c_str(), i.size());
            mkdir(path, 0775);

            setsPerLength[i.size()]
                ? setsPerLength[i.size()] += 1
                : setsPerLength[i.size()] = 1;

            ofstream outfile(string(path) + "/" + to_string(setsPerLength[i.size()]) + "-" + to_string(filenr++) + ".csv");
            // set precision to max to output coordinates as exact as possible
            outfile << setprecision(numeric_limits<double>::digits10);
            for (auto j : i) {
                datapoint point = datapoints.at(j);
                // write csv row
                outfile << point.id << "," << point.name << "," << point.lat << "," << point.lng << endl;
            }
            outfile.close();
        }

        // write info file
        ofstream outfile(DEST_FOLDER + "/info.txt");
        outfile << "INFO" << endl;
        outfile << "---" << endl;
        outfile << "Num of datapoints: " << datapoints.size() << endl;
        outfile << "Min Distance (m): " << MIN_DISTANCE * 1000 << endl;
        outfile << "Iterations: " << ITERATIONS << endl;
        outfile << "Unique sets found: " << results.size() << endl;
        outfile << "Biggest set found: " << biggestSet << endl;
        outfile << "Found " << results.size() << " sets with " << MIN_SET_SIZE << " or more datapoints:" << endl;

        for (auto kv : setsPerLength) {
            outfile << "  " << kv.first << " -> " << kv.second << " times" << endl;
        }

        outfile.close();
    }
    
    cout << "DONE!" << endl;
    return 0;
}
