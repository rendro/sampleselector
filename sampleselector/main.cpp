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
#include <thread>
#include <chrono>
#include <mutex>

using namespace std;

// source and destination paths
static const string PRESELECTED_FILE = "/Users/rendro/Downloads/preselected.csv";
static const string SOURCE_FILE = "/Users/rendro/Downloads/data.csv";
static const string DEST_FOLDER = "/Users/rendro/Downloads/result";

#define NUM_THREADS 8 

// use preselected set
static const bool WITH_PRESELECTION = true;
// min distance for buffer zone
static const double MIN_DISTANCE = 0.15;
// number of sets to create
static const int ITERATIONS = 1e3;
// min set size for exporting
static const int MIN_SET_SIZE = 63;
// write result csv files for all found sets that have at least MIN_SET_SIZE datapoints
static const bool WRITE_FILES = true;

// set of unique identifiers for points
typedef set<string> keyset;
// datapoint with id (line number or generic id), latitude, longitude and a buffer zone
struct datapoint {
    int id;
    string name;
    double lat;
    double lng;
    keyset buffer_ids;
};
// define dataset as a set of datapoints
typedef map<string, datapoint> dataset;

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

keyset getBufferIDsForDatapoint(dataset::iterator &item, dataset data) {
    keyset buffer_ids;
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
dataset csvToDataset(string file, bool create_buffer = true) {
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
        
        keyset buffer_ids;
        string key = record[1];
        datapoint datapoint = {
            stoi( record[0] ),
            key,
            stod( record[2] ),
            stod( record[3] ),
            buffer_ids
        };
        
        data.insert({key, datapoint});
    }
    if ( !infile.eof() ) {
        cerr << "NO EOF!" << endl;
    }
    
    if (create_buffer) {
        for (auto it = data.begin(); it != data.end(); ++it) {
            it->second.buffer_ids = getBufferIDsForDatapoint(it, data);
        }
    }
    
    return data;
}

// generate a result set from two sets of datapoints of which the first set contains all
// datapoints with other datapoints in the buffer zone and the of which the second set
// contains all datapoints without other datapoints in the buffer zone
keyset generateSet(keyset &standalonePoints, keyset &withNearbyPoints, dataset &datapoints) {
    random_device rd;
    mt19937 rng(rd());
    // create copies
    keyset remainingPoints(withNearbyPoints.begin(), withNearbyPoints.end());
    keyset result(standalonePoints.begin(), standalonePoints.end());

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
        keyset buffer_ids = datapoints.at(*it).buffer_ids;
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

struct thread_data {
    set<keyset> &resultSets;
    mutex &resultsMutex;
    keyset &standalonePoints;
    keyset &withNearbyPoints;
    dataset &datapoints;
    keyset &preselectedKeySet;
    int iterations;
};

void *genAndAddSet(void *threadarg) {
    stringstream ss;
    ss << this_thread::get_id();
    string thread_id = ss.str();
    struct thread_data *data;
    data = (struct thread_data *) threadarg;
    int i = 0;

    char buffer[100];
    sprintf(buffer, "Start thread(%s) with %d iterations\n", thread_id.c_str(), data->iterations);
    cout << buffer;
    memset(buffer, 0, sizeof(buffer));

    while (i < data->iterations) {
        i++;
        keyset set = generateSet(data->standalonePoints, data->withNearbyPoints, data->datapoints);
        if (WITH_PRESELECTION) {
            set.insert(data->preselectedKeySet.begin(), data->preselectedKeySet.end());
        }
        if (set.size() >= MIN_SET_SIZE) {
            lock_guard<mutex> guard(data->resultsMutex);
            data->resultSets.insert(set);
        }
        if (i % 1000 == 0) {
            sprintf(
                buffer,
                "Thread(%s) – Iteration %d – %.2f%%\n",
                thread_id.c_str(),
                i,
                static_cast<float>(i) / static_cast<float>(data->iterations) * 100.0
            );
            cout << buffer;
            memset(buffer, 0, sizeof(buffer));
        }
    }

    sprintf(buffer, "Terminate thread(%s) after %d iterations\n", thread_id.c_str(), i);
    cout << buffer;
    memset(buffer, 0, sizeof(buffer));
    return 0;
}

// main sampleselector
int main(int argc, const char * argv[]) {
    auto time_start = chrono::steady_clock::now();
    // read in data
    dataset preselectedPoints;
    dataset datapoints;
    keyset standalonePoints;
    keyset withNearbyPoints;

    if (WITH_PRESELECTION) {
        preselectedPoints = csvToDataset(PRESELECTED_FILE, false);
        dataset points = csvToDataset(SOURCE_FILE);
        dataset selectedPoints;

        cout << "With preselected sets" << endl;
        cout << "Number of preselected points: " << preselectedPoints.size() << endl;
        cout << "Number of points in source: " << points.size() << endl;

        // only use datapoints that are not in buffer zone of preselected set points
        for (auto point : points) {
            bool use_point = true;
            for (auto preselectedPoint : preselectedPoints) {
                if (haversine(point.second, preselectedPoint.second) < MIN_DISTANCE) {
                    use_point = false;
                }
            }
            if (use_point) {
                cout << "Use point: " << point.first << endl;
                point.second.buffer_ids.size() == 0
                    ? standalonePoints.insert(point.first)
                    : withNearbyPoints.insert(point.first);
            } else {
                cout << "Drop point: " << point.first << endl;
            }
        }
        datapoints.insert(points.begin(), points.end());
        datapoints.insert(preselectedPoints.begin(), preselectedPoints.end());
        
        cout << "Number of points after filter: " << datapoints.size() << endl;
        cout << endl;
    } else {
        datapoints = csvToDataset(SOURCE_FILE);
        for (auto it : datapoints) {
            it.second.buffer_ids.size() == 0
                ? standalonePoints.insert(it.first)
                : withNearbyPoints.insert(it.first);
        }
    }

    cout << "# Datapoints with datapoints in buffer zone   : " << withNearbyPoints.size() << endl;
    cout << "# Datapoints with no datapoint in buffer zone : " << standalonePoints.size() << endl;
    cout << endl;
    
    keyset preselectedKeySet;
    for (auto point : preselectedPoints) {
        preselectedKeySet.insert(point.first);
    }

    thread *threads[NUM_THREADS];
    mutex resultsMutex;
    set<keyset> results;
    thread_data td = {
        results,
        resultsMutex,
        standalonePoints,
        withNearbyPoints,
        datapoints,
        preselectedKeySet,
        ITERATIONS / NUM_THREADS,
    };

    for (int i = 0; i < NUM_THREADS; ++i) {
        threads[i] = new thread(&genAndAddSet, &td);
    }

    // join threads
    for (auto t : threads) {
        t->join();
    }

    double biggestSet = 0;
    for (auto &set : results) {
        biggestSet = biggestSet > set.size() ? biggestSet : set.size();
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
        if (WITH_PRESELECTION) {
            outfile << "Num of preselected points: " << preselectedKeySet.size() << endl;
            outfile << "Num of datapoints total: " << datapoints.size() << endl;
        } else {
            outfile << "Num of datapoints: " << datapoints.size() << endl;
        }
        outfile << "  Datapoints with points in buffer: " << withNearbyPoints.size() << endl;
        outfile << "  Datapoints with empty buffer: " << standalonePoints.size() << endl;
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
    auto time_end = chrono::steady_clock::now();
    cout << "Elapsed time in sec: " << chrono::duration_cast<chrono::seconds>(time_end - time_start).count() << endl;
    return 0;
}
