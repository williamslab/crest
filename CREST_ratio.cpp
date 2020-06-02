//
//  main.cpp
//  Debugcode
//
//  Created by qiaoqiao on 5/21/20.
//  Copyright © 2020 qiaoqiao. All rights reserved.
//
#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <limits.h>
#include <errno.h>
#include <iostream>
#include<vector>
#include<string>
#include<unordered_map>
#include <fstream>
#include <cstring>
#include <array>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cassert>
#include <map>
#include <set>
#include <iterator>
// #include <boost/algorithm/string.hpp>

// #define SIZE 2500
using namespace std;
///////////////////////////////////////////////////////

#define VERSION_NUMBER  "1.0.0"
#define RELEASE_DATE    "30 Apr 2020"

class CmdLineOpts {
  public:

    static bool parseCmdLineOptions(int argc, char **argv);
    static void printUsage(FILE *out, char *programName);

    static char *ibdFile;
    static char *relFile;
    static char *outPrefix;
    static double ibd2;
    static double clusterThres;
};

// define/initialize static members
char  *CmdLineOpts::ibdFile = NULL;
char  *CmdLineOpts::relFile = NULL;
char  *CmdLineOpts::outPrefix = NULL;
double CmdLineOpts::ibd2=0.02;
double CmdLineOpts::clusterThres=10;

// Parses the command line options for the program.
bool CmdLineOpts::parseCmdLineOptions(int argc, char **argv) {
  enum {
    IBD = CHAR_MAX + 1,
    CLS_THRES,
  };

  static struct option const longopts[] =
  {
  
    {"ibd2", required_argument, NULL, IBD},
    {"cluster_thres", required_argument, NULL, CLS_THRES},
    {0, 0, 0, 0}
  };

  // option index for getopt_long()
  int optionIndex = 0;
  int c;
  bool haveGoodArgs = true;

  char optstring[80] = "i:r:o:";
  while ((c = getopt_long(argc, argv, optstring, longopts, &optionIndex))
                                    != -1) {
    errno = 0;    
    char *endptr; 
    switch (c) {
      case 0:
    // flag set by getopt_long()
    break;

      case 'i':
    if (ibdFile != NULL) {
      if (haveGoodArgs)
        fprintf(stderr, "\n");
      fprintf(stderr, "ERROR: multiple definitions of ibd filename\n");
      haveGoodArgs = false;
    }
    ibdFile = optarg;
    break;
      case 'r':
    if (relFile != NULL) {
      if (haveGoodArgs)
        fprintf(stderr, "\n");
      fprintf(stderr, "ERROR: multiple definitions of rel filename\n");
      haveGoodArgs = false;
    }
    relFile = optarg;
    break;
      case 'o':
    outPrefix = optarg;
    break;

      case IBD:
    ibd2 = strtod(optarg, &endptr);
    if (errno != 0 || *endptr != '\0') {
      fprintf(stderr, "ERROR: unable to parse --ibd2 argument as floating point value\n");
      if (errno != 0)
        perror("strtod");
      exit(2);
    }

    if (ibd2 < 0 || ibd2 > 1) {
      fprintf(stderr, "ERROR: --ibd2 value must be between 0 and 1\n");
      exit(5);
    }
    break;

      case CLS_THRES:
    clusterThres = strtod(optarg, &endptr);
    if (errno != 0 || *endptr != '\0') {
      fprintf(stderr, "ERROR: unable to parse --cluster_thres argument as floating point value\n");
      if (errno != 0)
        perror("strtod");
      exit(2);
    }
    
    break;

      case '?':
    // bad option; getopt_long already printed error message
        printUsage(stderr, argv[0]);
    exit(1);
    break;

      default:
    exit(1);
    }
  }
  // Check for errors in command line options

  if (ibdFile == NULL || relFile == NULL || outPrefix == NULL) {
    if (haveGoodArgs)
      fprintf(stderr, "\n");
    fprintf(stderr, "ERROR: ibd, relatives, and output prefix names required\n");
    haveGoodArgs = false;
  }
  
  if (!haveGoodArgs) {
    printUsage(stderr, argv[0]);
  }

  return haveGoodArgs;
}

// Prints usage message to <out>.  <programName> should be argv[0]
void CmdLineOpts::printUsage(FILE *out, char *programName) {
  fprintf(out, "\n");
  fprintf(out, "CREST  v%s    (Released %s)\n\n",
      VERSION_NUMBER, RELEASE_DATE);
  fprintf(out, "Usage:\n");
  fprintf(out, "%s [ARGUMENTS]\n", programName);
  fprintf(out, "\n");
  fprintf(out, "REQUIRED ARGUMENTS:\n");
  fprintf(out, "  -i <filename>\t\tibd file including all ibdsegments among samples\n");
  fprintf(out, "  -r <filename>\t\trel file containing a list of relatives within samples\n");
  fprintf(out, "  -o <prefix>\t\toutput prefix (creates <prefix>.csv.)\n");
  
  fprintf(out, "OPTIONS:\n");
  fprintf(out, "  --ibd2 <#>\tibd2 proportion to exclude some close relative types (default 0.02; 0 disables)\n");
  fprintf(out, "  --cluster_thres <#>\tthreshold in cM to cluster relatives \n");
  fprintf(out, "\n");
}


/////////////////////////////////////////////////////////////////////////////

//used to store mutual relatives for clustering
struct MutRelative {
    public:
    int id;
    double coverage;
    vector<int> rel_list;
    double inst;
    double d1;
    double d2;
};
bool sortByCoverage(const MutRelative &lhs, const MutRelative &rhs) { return lhs.coverage < rhs.coverage; };

// define class block to store data
class block {
    public:
    string id1;
    string id2;
    int chrom;
    string type;
    double length;
    double value1, value2;
    
};

// read in data   ------------- should be consistent with IBIS

istream& operator>>(istream& is, block& IBDseg)
{
    
    string start;
    string end;
    double numMarker;
    double numErr;
    double errDensity;
    is >> IBDseg.id1 >> IBDseg.id2 >> IBDseg.chrom >> start >>  end >>  IBDseg.type >> IBDseg.value1 >> IBDseg.value2 >>IBDseg.length >> numMarker >> numErr >> errDensity;
    return is;
};


// store segment data stucture
struct segment{
    
    int chrom;
    //string type;
    double value1=0;
    double value2=0;
    
};

// structure to store the segments of one sample
struct dataBlock {
    public:
    double total=0;
    double max=0;
    vector<segment> segmentVec;
    bool index=0;
};

struct pairRelated {
    public:
    string id1;
    string id2;
    double kinship;
    double ibd2;
    int segCount;
    int degree;

};

// second degree pair
struct secondPair {
    public:
    string id1;
    string id2;
    double kinship;
    int num=0;
    double total=0;
    double max=0;
    vector<double> lenOfSeg;
    vector<double> numRel;
    vector<double> numCluster;
    vector<double> coverage;
    vector<double> overlapLen1;
    vector<double> overlapLen2;
    //vector<double> overlapt1;
    //vector<double> overlapt2;

    bool operator==(const secondPair& other) { return id1 == other.id1 && id2 == other.id2; }
    
};

// read in relatives
istream& operator>>(istream& is, pairRelated& pair)
{

    is >> pair.id1 >> pair.id2 >> pair.kinship >> pair.ibd2 >> pair.segCount >> pair.degree;
    return is;
};
/////////////////////////////////////////////////////////////////////////////////////////////
// function intersection 
dataBlock intersection( dataBlock x1y, dataBlock x2y ){
    dataBlock result;
    array<vector<int>,22> segByChrom1;
    array<vector<int>,22> segByChrom2;
    
    segment s1,s2; 
    
    // store segments index for different chromosomes
    for(int i=0; i < x1y.segmentVec.size(); i++){
        s1 = x1y.segmentVec[i];
        segByChrom1[s1.chrom-1].push_back(i);
    }
    
    for(int j = 0; j < x2y.segmentVec.size(); j++){
        s2 = x2y.segmentVec[j];
        segByChrom2[s2.chrom-1].push_back(j);
    }
    
    segment seg1;
    segment seg2;
    // double tmp=0;
    
    vector<int> segIndex1, segIndex2;
    // find overlap for each chromosome
    for(int chr=0; chr<22;chr++){
        
        if(segByChrom1[chr].size() != 0 && segByChrom2[chr].size() != 0){
            
            segIndex1 = segByChrom1[chr];
            segIndex2 = segByChrom2[chr];
            
            int m = 0;
            int n = 0;
            int k;
            segment seg;
            while(m < segIndex1.size() && n < segIndex2.size()){
                // if n is the last one, go to the next chrom
                seg1 = x1y.segmentVec[segIndex1[m]];
                seg2 = x2y.segmentVec[segIndex2[n]];
                if(seg1.value1 <= seg2.value1){
                    k = n;
                    while(k < segIndex2.size() && seg1.value2 > x2y.segmentVec[segIndex2[k]].value1){
                        seg = x2y.segmentVec[segIndex2[k]];
                        seg.value2 = min(seg.value2,seg1.value2);
                        result.segmentVec.push_back(seg);
                        result.max = max(result.max,abs(seg.value2-seg.value1));
                        result.total += abs(seg.value2-seg.value1);
                        ++k;
                        
                    }
                    ++m;
                } else{
                    k = m;
                    while(k < segIndex1.size() && seg2.value2 > x1y.segmentVec[segIndex1[k]].value1){
                        seg = x1y.segmentVec[segIndex1[k]];
                        seg.value2 = min(seg.value2,seg2.value2);
                        result.segmentVec.push_back(seg);
                        result.max = max(result.max,abs(seg.value2-seg.value1));
                        result.total += abs(seg.value2-seg.value1);
                        ++k;
                        
                    }
                    ++n;
                    
                }
            }
            
        }
    }
    
    return result;
};
// function combine
dataBlock combine(dataBlock xy, dataBlock tmp){
    dataBlock result;
    
    if(tmp.total == 0){
        result=xy;
        return result;
    }else if(xy.total == 0){
        result=tmp;
        return result;
    }else{
        array<vector<int>,22> segByChrom1;
        array<vector<int>,22> segByChrom2;
        
        segment s1,s2;
        
        // store segments index for different chromosomes
        for(int i=0; i < xy.segmentVec.size(); i++){
            s1=xy.segmentVec[i];
            segByChrom1[s1.chrom-1].push_back(i);
        }
        
        for(int j=0; j < tmp.segmentVec.size(); j++){
            s2=tmp.segmentVec[j];
            segByChrom2[s2.chrom-1].push_back(j);
        }
        
        segment seg1;
        segment seg2;
        // segment seg;
        
        vector<int> segIndex1, segIndex2;
        // find overlap for each chromosome
        for(int chr=0; chr<22;chr++){
            
            // cout<<"Chromosome"<<chr<<endl;
            
            int m = 0;
            int n = 0;
            int num;
            vector<segment> newresult;
            segIndex1 = segByChrom1[chr];
            segIndex2 = segByChrom2[chr];
            while(m < segByChrom1[chr].size() && n < segByChrom2[chr].size()){
                // cout<<m<<" " << n <<endl;
                // cout<<"first while"<<endl;
                seg1 = xy.segmentVec[segIndex1[m]];
                seg2 = tmp.segmentVec[segIndex2[n]];
                num = newresult.size();
                // cout << num<<endl;
                // cout <<"seg1" <<seg1.value1 << " " <<seg1.value2<<endl;
                // cout <<"seg2" <<seg2.value1 << " " <<seg2.value2<<endl;
                if(num != 0){
                    segment& seg = newresult[num-1];
                    // cout <<"seg" <<seg.value1 << " " <<seg.value2<<endl;
                    if(seg.value2 >= seg1.value1){
                        seg.value2 = max(seg.value2,seg1.value2);
                        m++;
                    }else if (seg.value2 >= seg2.value1){
                        seg.value2 = max(seg.value2,seg2.value2);
                        n++;
                    }else if (seg.value2 < seg1.value1 && seg.value2 < seg2.value1){
                        if(seg1.value1 < seg2.value1){
                            newresult.push_back(seg1);
                            m++;
                        }else{
                            newresult.push_back(seg2);
                            n++;
                        }
                    }
                }
                else{
                    if(seg1.value1 < seg2.value1){
                        newresult.push_back(seg1);
                        m++;
                    }else{
                        newresult.push_back(seg2);
                        n++;
                    }
                    
                }
                
            }
            while(m < segByChrom1[chr].size()){
               // cout<<"second  while"<<endl;
                seg1 = xy.segmentVec[segIndex1[m]];
                num = newresult.size();
                if(num != 0){
                    segment& seg = newresult[num-1];
                    if(seg.value2 >= seg1.value1){
                        seg.value2 = max(seg.value2,seg1.value2);
                    }else{
                        newresult.push_back(seg1);
                    }
                }else{
                    newresult.push_back(seg1);
                }
                m++;
                
            }
            
            while(n < segByChrom2[chr].size()){
               // cout <<"third while"<<endl;
                seg2 = tmp.segmentVec[segIndex2[n]];
                num = newresult.size();
                if(num != 0){
                    segment& seg = newresult[num-1];
                    if(seg.value2 >= seg2.value1){
                        seg.value2 = max(seg.value2,seg2.value2);
                    }else{
                        newresult.push_back(seg2);
                    }
                }else{
                    newresult.push_back(seg2);
                }
                n++;
                
            }
            // cout << newresult.size() << endl;
            for (int i = 0; i < newresult.size();i++){
                result.segmentVec.push_back(newresult[i]);
                result.total += abs(newresult[i].value2 - newresult[i].value1);
                result.max = max(result.max,abs(newresult[i].value2 - newresult[i].value1));
            }
            // cout << "total" << result.total << endl;
        }
        return result;
    }
};// end combine
////////////////////////////////////////////////////////////////////////////////

ofstream & operator <<(ofstream& ofs, vector<double>& vec)
{
    for (int i = 0; i < vec.size(); i++)
    {
        ofs << vec[i] << ' ';
    }
    return ofs;
}

//writing the output to a file
ofstream & operator << (ofstream& ofs, secondPair &entry)
{
    ofs << entry.id1 << "," << entry.id2 <<  ",";
    ofs << entry.coverage << ",";
    ofs << entry.overlapLen1 << ",";
    ofs << entry.overlapLen2 << "\n";
    //ofs << entry.overlapt1 << ",";
    //ofs << entry.overlapt2 << "\n";
    return ofs;
}
/////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {
    bool success = CmdLineOpts::parseCmdLineOptions(argc, argv);

    if (!success)
    return -1;

  int outPrefixLen = strlen(CmdLineOpts::outPrefix);
  char *outFile = new char[ outPrefixLen + 4 + 1 ]; //+4 for .csv, + 1 for \0
  if (outFile == NULL) {
    printf("ERROR: out of memory");
    exit(5);
  }

  // open the log file
  sprintf(outFile, "%s.csv", CmdLineOpts::outPrefix);

 // read in second degree pairs list with kinship
    ifstream idfile(CmdLineOpts::relFile);
    string id;
    // map ID into int, find o(logn)
    map<string,int> ID;
    // for tracking back individuals 
    map<int,string> ID_track;
    // used to choose distant relatives 
    pairRelated pair_kinship;
    // a vector of second degree relatives 
    vector<secondPair> setOf2ndPair;
    vector<pairRelated> setOfRelatives;
    secondPair pair;
    // total number of samples in second degree pairs 
    int NUM = 0; 
    double ibd2_thres = CmdLineOpts::ibd2; 
    while(getline(idfile,id)){
        stringstream istring(id);
        
        istring >> pair_kinship;
        // second degree into ID
        if (pair_kinship.degree == 2 && pair_kinship.ibd2<=ibd2_thres){
		  if(ID.find(pair_kinship.id1) == ID.end()){
		          ID[pair_kinship.id1] = ID.size();
		          NUM++;
            }// end if id1
            if(ID.find(pair_kinship.id2) == ID.end()){
                ID[pair_kinship.id2] = ID.size();
                NUM++;
            }// end if id2
            pair.id1 = pair_kinship.id1;
            pair.id2 = pair_kinship.id2;
            pair.kinship = pair_kinship.kinship;
	    setOf2ndPair.push_back(pair);
        }else if(pair_kinship.degree>2 && pair_kinship.degree<7 && pair_kinship.ibd2<=ibd2_thres){
            setOfRelatives.push_back(pair_kinship);
        }
    }// end while read in second pairs list

   //************************************************************************

    const int NUM1=ID.size();
    unordered_map<int,set<int>> relatives;
    map<string,int>::iterator it1;
    map<string,int>::iterator it2; 

    vector<pairRelated>::iterator it3;
    for(it3 = setOfRelatives.begin(); it3 != setOfRelatives.end(); ++it3){
        
        it1 = ID.find((*it3).id1);
        it2 = ID.find((*it3).id2);
        // id1 in the ID and is second degree samples
        
        if( it1 != ID.end() && it1->second <= NUM1){
            // if id2 not in the ID, map it to a number first
            // cout << ID.size() << endl;
            if (it2 == ID.end()){
                ID[(*it3).id2] = ID.size();
            }
        // put id2 in the set of relatives of id1
            relatives[it1->second].insert(ID[(*it3).id2]);
        }// end if id1
        //
        // id1 not in the ID but id2 in the ID
        if(it2 != ID.end() && it2->second <=NUM1){
            // id1 not in the ID
            // cout << ID.size() << endl;
            if (it1 == ID.end()){
            ID[(*it3).id1] = ID.size();
            }
            relatives[it2->second].insert(ID[(*it3).id1]);
        }// end if id2
        
    }

    for(auto id = ID.begin(); id != ID.end(); ++id){
        ID_track[id->second] = id->first;
    }
   
    ifstream IBDdata(CmdLineOpts::ibdFile);
    ofstream out(outFile);
    
    if( ! out || ! IBDdata)    {
        cout << "Error opening file for output" << endl ;
        return -1 ;
    }
    
    string line;
    unordered_map<int,unordered_map<int, dataBlock>> IBDsegment;
    block sample;
    
    // store segments to the array
    while(getline(IBDdata, line)){
        stringstream iss(line);
        
        iss >> sample; 
	assert(sample.value1 <= sample.value2);
    // only read in segments for samples in ID
    if (ID.find(sample.id1) == ID.end() || ID.find(sample.id2) == ID.end()){
        continue;
    }
        int row = ID[sample.id1]-1;
        int col = ID[sample.id2]-1;
        int t;
        // swap the row and col
        if(row > col){
            t = row;
            row = col;
            col = t;
        }
        double length;
        segment oneSegment;
        
        oneSegment.chrom = sample.chrom;
        oneSegment.value1 = sample.value1;
        oneSegment.value2 = sample.value2;

        
        IBDsegment[row][col].segmentVec.push_back(oneSegment);
        length = abs(oneSegment.value2 - oneSegment.value1);
  
        IBDsegment[row][col].total += length;
        
        if(IBDsegment[row][col].max < length){
            IBDsegment[row][col].max = length;
        }
        
    }//end read in lines

    dataBlock x1y;
    dataBlock x2y;
    
    for(auto id = ID.begin(); id != ID.end(); ++id){
        ID_track[id->second] = id->first;
    }
    
    vector<secondPair>::iterator it;
    for(it = setOf2ndPair.begin(); it != setOf2ndPair.end(); ++it){
        int i = ID[(*it).id1]-1;
        int j = ID[(*it).id2]-1;
        int ij;
        if(i > j){
            ij = i;
            i = j;
            j = ij;
        }
        dataBlock tmp1;
        dataBlock tmp;
        dataBlock tmp2;
        dataBlock tmp_union1;
        dataBlock tmp_union2;
                
        (*it).total = IBDsegment[i][j].total;
        (*it).max = IBDsegment[i][j].max;
        (*it).num = IBDsegment[i][j].segmentVec.size();
    
        vector<segment>::iterator len;
        for(len = IBDsegment[i][j].segmentVec.begin(); len != IBDsegment[i][j].segmentVec.end(); ++len){
            (*it).lenOfSeg.push_back(abs((*len).value2-(*len).value1));
        }

        vector<int> commonRelatives;
        
        set_intersection(relatives[ID[(*it).id1]].begin(), relatives[ID[(*it).id1]].end(), relatives[ID[(*it).id2]].begin(), relatives[ID[(*it).id2]].end(), back_inserter(commonRelatives));
        
        //******************************************************************************************************
        // cluster the relatives  
        vector<MutRelative> clusters;
        vector<int> index_rel;
        dataBlock ibd_rel;
        vector<MutRelative> relativeVec;

        if (commonRelatives.size() !=0){
            for (vector<int>::iterator it1 = commonRelatives.begin(); it1!=commonRelatives.end(); ++it1){
                int k = (*it1)-1;
             
                MutRelative rel;
                dataBlock x1y;
                dataBlock x2y;

                rel.id=k+1;
                if(k<i){
                x1y = IBDsegment[k][i];
                x2y = IBDsegment[k][j];
                }else if (k>i && k<j){
                x1y = IBDsegment[i][k];
                x2y = IBDsegment[k][j];
                }else{
                x1y = IBDsegment[i][k];
                x2y = IBDsegment[j][k];
                }
                rel.coverage = max(x1y.total, x2y.total); 
                relativeVec.push_back(rel); 
            }
            sort (relativeVec.begin(), relativeVec.end(), sortByCoverage);

            for (vector<MutRelative>::iterator it = relativeVec.begin(); it!=relativeVec.end(); ++it){
            int h = (*it).id-1;
            if (index_rel.empty()){
                index_rel.push_back(h+1);
                MutRelative rel;
                rel.id = h+1;
                rel.rel_list.push_back(h+1);
                clusters.push_back(rel);
            }else{  
                int indicator = 0;
                for (vector<int>::iterator it1 = index_rel.begin(); it1 !=index_rel.end(); ++it1){
                    int l = (*it1) - 1;
		    dataBlock sharing;
                    if (l < h){
                            sharing = IBDsegment[l][h];
                        } else{
                            sharing = IBDsegment[h][l];
                        }
                    if(sharing.total > CmdLineOpts::clusterThres){
                        indicator = 1;
                            
                    }// end if
                }// end for loop 
                if (indicator==0){
                    index_rel.push_back(h+1);
                    MutRelative rel;
                    rel.id = h+1;
                    rel.rel_list.push_back(h+1);
                    clusters.push_back(rel);
                    }// end if
            } // end else 
        }// end for loop 
        for (vector<MutRelative>::iterator it = relativeVec.begin(); it!=relativeVec.end(); ++it){
            int h = (*it).id-1;
            if (find(index_rel.begin(), index_rel.end(), h+1) == index_rel.end()){
                vector<int> overlappingVec;
                for (vector<MutRelative>::iterator it1 = clusters.begin(); it1!=clusters.end(); ++it1){
                    vector<int> relativesList = (*it1).rel_list; 
                    double overlapping=0;
                    for (vector<int>::iterator it2=relativesList.begin(); it2!=relativesList.end(); ++it2){
                        int id_rel = (*it2)-1;
                        if (id_rel < h){
                            if (overlapping<IBDsegment[id_rel][h].total){
                                overlapping = IBDsegment[id_rel][h].total;
                            }// end if 
                        } else{
                            if (overlapping<IBDsegment[h][id_rel].total){
                                overlapping = IBDsegment[h][id_rel].total;
                            }// end if
                        }// end else 
                    } // end for 
                    overlappingVec.push_back(overlapping); 
                }//end for 
                int maxElementIndex = max_element(overlappingVec.begin(),overlappingVec.end()) - overlappingVec.begin();

                relativeVec[maxElementIndex].rel_list.push_back(h+1);
            }// end if 
        }//end for 
    }// end if 


        (*it).numRel.push_back( commonRelatives.size()/1.0);
        (*it).numCluster.push_back( clusters.size()/1.0);
        double n = 0;
        double d1 = 0;
        double d2 = 0;
        double cov = 0;     

//**********************************************************************************************
        if (clusters.size() !=0){
            dataBlock x1y;
            dataBlock x2y;
            dataBlock tmp1;
            dataBlock tmp;
            dataBlock tmp2;
            dataBlock tmp_union1;
            dataBlock tmp_union2;
            // go through all clusters

            for (vector<MutRelative>::iterator it1 = clusters.begin(); it1!=clusters.end(); ++it1){
                MutRelative c = (*it1);

                //still need to discuss the k,i,j
                for (vector<int>::iterator it2 = c.rel_list.begin(); it2!=c.rel_list.end(); ++it2){
                    int k = (*it2)-1;
                    if(k<i){
                        x1y = IBDsegment[k][i];
                        x2y = IBDsegment[k][j];
                    }else if (k>i && k<j){
                        x1y = IBDsegment[i][k];
                        x2y = IBDsegment[k][j];
                    }else{
                        x1y = IBDsegment[i][k];
                        x2y = IBDsegment[j][k];
                    }
                    if(x1y.segmentVec.size() != 0 && x2y.segmentVec.size() != 0 ){

                        tmp = intersection(x1y,x2y);
                        tmp1 = combine(tmp,tmp1);
                        tmp_union1 = combine(x1y,tmp_union1);
                        tmp_union2 = combine(x2y,tmp_union2);
                    }// end if check segments overlapping
                } // end for loop through the relative list
                
                n += tmp1.total;
                d1 += tmp_union1.total;
                d2 += tmp_union2.total;
                cov += max(tmp_union1.total,tmp_union2.total);
        }
        // to do: store the results to 

        
        (*it).coverage.push_back(cov);
        (*it).overlapLen1.push_back(n / d1); 
        (*it).overlapLen2.push_back(n / d2); 
        //(*it).overlapt1.push_back(d1);
        //(*it).overlapt2.push_back(d2);
            
        }else{
            (*it).coverage.push_back(0);
            (*it).overlapLen1.push_back(0); 
            (*it).overlapLen2.push_back(0); 
            //(*it).overlapt1.push_back(0);
            //(*it).overlapt2.push_back(0);
        }// end if relatives not empty
        out << (*it);
    } // end for loop over 2nd degree pairs
   
    IBDdata.close();
    out.close();
    
    return 0;
}
