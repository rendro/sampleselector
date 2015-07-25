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
static const string DEST_FOLDER = "/Users/rendro/Desktop/result";

// min distance for buffer zone
static const double MIN_DISTANCE = 5.0;
// number of sets to create
static const int ITERATIONS = 1000;
// min set size for exporting
static const int MIN_SET_SIZE = 72;
// write result csv files for all found sets that have at least MIN_SET_SIZE datapoints
static const bool WRITE_FILES = true;

// datapoint with id (line number or generic id), latitude, longitude and a buffer zone
struct datapoint {
    int id;
    double lat;
    double lng;
    set<datapoint> buffer;
};

// define dataset as a set of datapoints
typedef set<datapoint> dataset;

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
    double dlon = DEG2RAD * (end.lng - start.lng);
    double dlat = DEG2RAD * (end.lat - start.lat);
    double a = pow(sin(dlat/2), 2) + cos(DEG2RAD * start.lat) * cos(DEG2RAD * start.lat) * pow(sin(dlon/2), 2);
    double b = 2 * asin(sqrt(a));
    return EARTH_RADIUS * b;
}

// read in a csv (columns [0, 1, 2] match [id, lat, lng] of a datapoint
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
        
        dataset buffer;
        
        datapoint datapoint = {
            stoi( record[0] ),
            stod( record[1] ),
            stod( record[2] ),
            buffer
        };
        
        data.insert( datapoint );
    }
    if ( !infile.eof() ) {
        cerr << "NO EOF!" << endl;
    }
    
    return data;
}

// generate a result set from two sets of datapoints of which the first set contains all
// datapoints with other datapoints in the buffer zone and the of which the second set
// contains all datapoints without other datapoints in the buffer zone
dataset generateSet(dataset& withNearbyDataset, dataset& standaloneDataset) {
    random_device rd;
    mt19937 rng(rd());
    
    dataset remainingDataset(withNearbyDataset.begin(), withNearbyDataset.end());
    dataset resultSet(standaloneDataset.begin(), standaloneDataset.end());
    
    while (remainingDataset.size() != 0) {
        // create iterator
        dataset::iterator it = remainingDataset.begin();
        
        // generate random index
        uniform_int_distribution<int> uni(0, (int)remainingDataset.size());
        int r = uni(rng);
        // pick random datapoint by advancing the iterator to the random position
        advance(it, r % remainingDataset.size());
        
        
        // add picked datapoint to result list
        resultSet.insert(*it);
        
        // remove all datapoints within buffer zone if still in remaining dataset
        for (dataset::iterator j = it->buffer.begin(); j != it->buffer.end(); ++j) {
            dataset::iterator tmp = remainingDataset.find(*j);
            if (tmp != remainingDataset.end()) {
                remainingDataset.erase(tmp);
            }
        }
        
        // remove picked datapoint from remaining list
        remainingDataset.erase(remainingDataset.find(*it));
    }
    
    return resultSet;
}

// main sampleselector
int main(int argc, const char * argv[]) {
    
    dataset withNearbyDataset;
    dataset standaloneDataset;
    
    {
        // read in data
        dataset allDatapoints = csvToDataset(SOURCE_FILE);
        
        // iterate over all datapoints and add nearby points to the buffer zone
        dataset::iterator i;
        for (i = allDatapoints.begin(); i != allDatapoints.end(); ++i) {
            dataset::iterator j;
            dataset buffer;
            for (j = allDatapoints.begin(); j != allDatapoints.end(); ++j) {
                if (i != j) {
                    double dist = haversine(*i, *j);
                    if (dist < MIN_DISTANCE) {
                        buffer.insert(*j);
                    }
                }
            }
            
            if (buffer.size() > 0) {
                datapoint k = {
                    i->id,
                    i->lat,
                    i->lng,
                    buffer
                };
                withNearbyDataset.insert(k);
            } else {
                standaloneDataset.insert(*i);
            }
        }
    }
    
    cout << "# Datapoints with datapoints in buffer zone   : " << withNearbyDataset.size() << endl;
    cout << "# Datapoints with no datapoint in buffer zone : " << standaloneDataset.size() << endl;
    cout << endl;
    
    // create result sets
    set<dataset> results;
    int i = 0;
    double biggestSet = 0;
    while (i < ITERATIONS) {
        i++;
        dataset resultSet = generateSet(withNearbyDataset, standaloneDataset);
        
        results.insert(resultSet);
        biggestSet = biggestSet > resultSet.size() ? biggestSet : resultSet.size();
        
        cout << "Iteration #" << to_string(i) << " => " << resultSet.size() << " datapoints" << endl;
    }
    
    cout << endl;
    cout << "Unique result sets found: " << results.size() << endl;
    cout << "Biggest result set found: " << biggestSet << endl;
    
    // subset with all result sets that have a minimum length of MIN_SET_SIZE
    set<dataset> resultsSubset;
    copy_if(results.begin(), results.end(), inserter(resultsSubset, resultsSubset.begin()), [](const dataset & resultSet){
        return resultSet.size() >= MIN_SET_SIZE;
    });
    
    cout << "Found " << resultsSubset.size() << " sets with " << MIN_SET_SIZE << " or more datapoints" << endl;
    
    if (WRITE_FILES) {
        int filenr = 0;
        
        map<unsigned long, int> setsPerLength;
        
        for (set<dataset>::iterator i = resultsSubset.begin(); i != resultsSubset.end(); ++i) {
            
            char path [256];
            sprintf(path, "%s/%lu", DEST_FOLDER.c_str(), i->size());
            mkdir(path, 0775);
            
            if (setsPerLength[i->size()]) {
                setsPerLength[i->size()] += 1;
            } else {
                setsPerLength[i->size()] = 1;
            }
            
            ofstream outfile(string(path) + "/" + to_string(setsPerLength[i->size()]) + "-" + to_string(filenr++) + ".csv");
            
            for (dataset::iterator j = i->begin(); j != i->end(); ++j) {
                // set precision to max to output coordinates as exact as possible
                outfile << setprecision(numeric_limits<double>::digits10);
                // write csv row
                outfile << j->id << "," << j->lat << "," << j->lng << endl;
            }
            
            outfile.close();
        }
        
        // write info file
        ofstream outfile(DEST_FOLDER + "/info.txt");
        outfile << "INFO" << endl;
        outfile << "---" << endl << endl;
        outfile << "Iterations: " << ITERATIONS << endl;
        outfile << "Unique sets found: " << results.size() << endl;
        outfile << "Biggest set found: " << biggestSet << endl;
        outfile << "Found " << resultsSubset.size() << " sets with " << MIN_SET_SIZE << " or more datapoints:" << endl;
        
        for (auto kv : setsPerLength) {
            outfile << "    " << kv.first << " -> " << kv.second << " times" << endl;
        }
        
        outfile.close();
    }
    
    cout << "DONE!";
    return 0;
}
