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
#include <cmath>
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
    static int maxDegree;
    static double kinship_lw;
    static double kinship_up;
    static double kinship_rel_lw;
    static double kinship_rel_up;
    static int pcflag;
    static double ibd0;
};

// define/initialize static members
char  *CmdLineOpts::ibdFile = NULL;
char  *CmdLineOpts::relFile = NULL;
char  *CmdLineOpts::outPrefix = NULL;
double CmdLineOpts::ibd2 = 0.02;
int CmdLineOpts::maxDegree = 6;
double CmdLineOpts::ibd0 = 0.1;
double CmdLineOpts::kinship_lw = 0.0883883;
double CmdLineOpts::kinship_up = 0.1767767;
double CmdLineOpts::kinship_rel_lw = 0.0055243;
double CmdLineOpts::kinship_rel_up = 0.0883883;
int CmdLineOpts::pcflag=0;

// Parses the command line options for the program.
bool CmdLineOpts::parseCmdLineOptions(int argc, char **argv) {
  enum {
    IBD = CHAR_MAX + 1,
    DEGREE,
    PC,
    IBD0,
    KIN_LW,
    KIN_UP,
    KIN_REL_LW,
    KIN_REL_UP,
    
  };

  static struct option const longopts[] =
  {
  
    {"ibd2", required_argument, NULL, IBD},
    {"max_degree", required_argument, NULL, DEGREE},
    {"pc", no_argument, NULL, PC},
    {"ibd0", required_argument, NULL, IBD0},
    {"kinship_lw", required_argument, NULL, KIN_LW},
    {"kinship_up", required_argument, NULL, KIN_UP},
    {"kinship_rel_lw", required_argument, NULL, KIN_REL_LW},
    {"kinship_rel_up", required_argument, NULL, KIN_REL_UP},
    {0, 0, 0, 0},
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

        case DEGREE:
    maxDegree = strtol(optarg, &endptr, 10);
    if (errno != 0 || *endptr != '\0') {
      fprintf(stderr, "ERROR: unable to parse --max_degree argument as integer\n");
      if (errno != 0)
        perror("strtol");
      exit(2);
    }
    
    break;

    case PC:
        pcflag = 1;
        break;
                
        case IBD0:
    ibd0 = strtod(optarg, &endptr);
    if (errno != 0 || *endptr != '\0') {
      fprintf(stderr, "ERROR: unable to parse --ibd0 argument as floating point value\n");
      if (errno != 0)
        perror("strtod");
      exit(2);
    }

    if (ibd0 < 0 || ibd0 > 1) {
      fprintf(stderr, "ERROR: --ibd0 value must be between 0 and 1\n");
      exit(5);
    }
    break;
            
    case KIN_UP:
    kinship_up = strtod(optarg, &endptr);
    if (errno != 0 || *endptr != '\0') {
      fprintf(stderr, "ERROR: unable to parse --kinship_up argument as floating point value\n");
      if (errno != 0)
        perror("strtod");
      exit(2);
    }

    if (kinship_up < 0 || kinship_up > 1) {
      fprintf(stderr, "ERROR: --kinship_up value must be between 0 and 1\n");
      exit(5);
    }
    break;
            
    case KIN_LW:
       kinship_lw = strtod(optarg, &endptr);
       if (errno != 0 || *endptr != '\0') {
         fprintf(stderr, "ERROR: unable to parse --kinship_lw argument as floating point value\n");
         if (errno != 0)
           perror("strtod");
         exit(2);
       }

       if (kinship_lw < 0 || kinship_lw > 1) {
         fprintf(stderr, "ERROR: --kinship_lw value must be between 0 and 1\n");
         exit(5);
       }
       break;
            
    case KIN_REL_UP:
    kinship_rel_up = strtod(optarg, &endptr);
    if (errno != 0 || *endptr != '\0') {
      fprintf(stderr, "ERROR: unable to parse --kinship_rel_up argument as floating point value\n");
      if (errno != 0)
        perror("strtod");
      exit(2);
    }

    if (kinship_rel_up < 0 || kinship_rel_up > 1) {
      fprintf(stderr, "ERROR: --kinship_rel_up value must be between 0 and 1\n");
      exit(5);
    }
            break;
            
    case KIN_REL_LW:
    kinship_rel_lw = strtod(optarg, &endptr);
    if (errno != 0 || *endptr != '\0') {
      fprintf(stderr, "ERROR: unable to parse --kinship_rel_lw argument as floating point value\n");
      if (errno != 0)
        perror("strtod");
      exit(2);
    }

    if (kinship_rel_lw < 0 || kinship_rel_lw > 1) {
      fprintf(stderr, "ERROR: --kinship_rel_up value must be between 0 and 1\n");
      exit(5);
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
  //fprintf(out, "  --ibd_thres <#>\tthreshold in cM to cluster relatives \n");
  fprintf(out, "  --pc \tallows to analyze parent-child pairs \n");
  fprintf(out, "  --ibd0 <#>\tibd0 proportion to classify parent-child relationship, only useful when --pc exists (default 0.01)\n");
  fprintf(out, "  --kinship_lw <#>\tupper bound of kinship to choose second degree relatives \n");
  fprintf(out, "  --kinship_up <#>\tlower bound of kinship to choose second degree relatives \n");
  fprintf(out, "  --kinship_rel_lw <#>\tlower bound of kinship to choose mutual relatives \n");
  fprintf(out, "  --kinship_rel_up <#>\tupper bound of kinship to choose mutual relatives \n");
  fprintf(out, "\n");
}


////////////////////////////////////////////////////////////////////////////

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

// structure to store the segments of one pair
struct dataBlock {
    public:
    double total=0;
    double max=0;
    vector<segment> segmentVec;
};

//used to store mutual relatives for clustering
struct MutRelative {
    public:
    int id;
    double coverage;
    dataBlock intersect;
    dataBlock ibd_x1;
    dataBlock ibd_x2;
};

bool sortByCoverage(const MutRelative &lhs, const MutRelative &rhs) { return lhs.coverage > rhs.coverage; };

struct hap {
    public:
    vector<int> rel_list;
    dataBlock ibd_union1;
    dataBlock ibd_union2;
    dataBlock intersect;
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
    int is_pc;
    double coverage;
    double ratio1;
    double ratio2;

    bool operator==(const secondPair& other) { return id1 == other.id1 && id2 == other.id2; }
    
};

// read in relatives
istream& operator>>(istream& is, pairRelated& pair)
{

    is >> pair.id1 >> pair.id2 >> pair.kinship >> pair.ibd2 >> pair.segCount >> pair.degree;
    return is;
};

//writing the output to a file
ofstream & operator << (ofstream& ofs, secondPair &entry)
{
    ofs << entry.id1 << "," << entry.id2 <<  "," << entry.is_pc << ","
    << entry.coverage << "," << entry.ratio1 << "," << entry.ratio2 <<"\n" ;
    return ofs;
};

//*************************************************************
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
            
            int m = 0;
            int n = 0;
            int num;
            vector<segment> newresult;
            segIndex1 = segByChrom1[chr];
            segIndex2 = segByChrom2[chr];
            while(m < segByChrom1[chr].size() && n < segByChrom2[chr].size()){
        
                seg1 = xy.segmentVec[segIndex1[m]];
                seg2 = tmp.segmentVec[segIndex2[n]];
                num = newresult.size();

                if(num != 0){
                    segment& seg = newresult[num-1];
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
            for (int i = 0; i < newresult.size();i++){
                result.segmentVec.push_back(newresult[i]);
                result.total += abs(newresult[i].value2 - newresult[i].value1);
                result.max = max(result.max,abs(newresult[i].value2 - newresult[i].value1));
            }
        }
        return result;
    }
};// end combine


dataBlock find_ibd011(dataBlock ibd11, dataBlock ibd_all){
    dataBlock result;
    // no IBD between ys
    if(ibd_all.total == 0){
        result=ibd11;
        return result;
    }else if(ibd11.total == 0){
        result=ibd11;
        return result;
    }else{
        array<vector<int>,22> segByChrom1;
        array<vector<int>,22> segByChrom2;
        
        segment s1,s2;
        
        // store segments index for different chromosomes
        for(int i=0; i < ibd11.segmentVec.size(); i++){
            s1=ibd11.segmentVec[i];
            segByChrom1[s1.chrom-1].push_back(i);
        }
        
        for(int j=0; j < ibd_all.segmentVec.size(); j++){
            s2=ibd_all.segmentVec[j];
            segByChrom2[s2.chrom-1].push_back(j);
        }
        
        segment seg1;
        segment seg2;
        
        vector<int> segIndex1, segIndex2;
        dataBlock complement;
        // find remaining intervals of ibd_all
        for(int chr=0; chr<22;chr++){
            int m = segByChrom1[chr].size();
            if (m!=0){
                int n = 0;
                vector<segment> newresult;
                segIndex1 = segByChrom1[chr];
                segIndex2 = segByChrom2[chr];
                
                segment currentSeg;
                currentSeg.chrom = chr+1;
                if (segByChrom2[chr].size() == 0){
                    currentSeg.value1 =ibd11.segmentVec[segIndex1[0]].value1;
                    currentSeg.value2 = ibd11.segmentVec[segIndex1[m-1]].value2;
                    newresult.push_back(currentSeg);
                }else{
                    if(ibd11.segmentVec[segIndex1[0]].value1 <ibd_all.segmentVec[segIndex2[0]].value1){
                        currentSeg.value1 =ibd11.segmentVec[segIndex1[0]].value1;
                        currentSeg.value2 = ibd_all.segmentVec[segIndex2[0]].value1;
                        newresult.push_back(currentSeg);
                    }
                    while(n+1 < segByChrom2[chr].size()){
                        currentSeg.value1 = ibd_all.segmentVec[segIndex2[n]].value2;
                        currentSeg.value2 = ibd_all.segmentVec[segIndex2[n+1]].value1;
                        newresult.push_back(currentSeg);
                        n++;
                    }
                    if(ibd_all.segmentVec[segIndex2[n]].value2 <  ibd11.segmentVec[segIndex1[m-1]].value2){
                        currentSeg.value1 =ibd_all.segmentVec[segIndex2[n]].value2;
                        currentSeg.value2 = ibd11.segmentVec[segIndex1[m-1]].value2;
                        newresult.push_back(currentSeg);
                    }
                }// end else
                    
                for (int i = 0; i < newresult.size();i++){
                    complement.segmentVec.push_back(newresult[i]);
                    complement.total += abs(newresult[i].value2 - newresult[i].value1);
                    complement.max = max(complement.max,abs(newresult[i].value2 - newresult[i].value1));
                }
            }// end if
        }//end chr
        result = intersection(ibd11, complement);
        return result;
    }
};// end finding ibd011

//*************************************************************

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
  char *outRelFile = new char[ outPrefixLen + 4 + 1 ];
    
    double ibd2_thres = CmdLineOpts::ibd2;
    int pc = CmdLineOpts::pcflag;
    double ibd0_thres = CmdLineOpts::ibd0;
    int deg_cutoff = CmdLineOpts::maxDegree;
    double kin_up = CmdLineOpts::kinship_up;
    double kin_rel_up = CmdLineOpts::kinship_rel_up;
    double kin_lw = CmdLineOpts::kinship_lw;
    double kin_rel_lw = CmdLineOpts::kinship_rel_lw;
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
    vector<secondPair> setOfPair;
    vector<pairRelated> setOfRelatives;
    //vector<secondPair> setOfPC;
    secondPair pair;
    // total number of samples in second degree pairs
    int NUM = 0;
    
    vector<int> dc_list;
    while(getline(idfile,id)){
        stringstream istring(id);
        
        istring >> pair_kinship;
        // second degree into ID
        if(ID.find(pair_kinship.id1) == ID.end()){
              ID[pair_kinship.id1] = ID.size();
              NUM++;
        }// end if id1
        if(ID.find(pair_kinship.id2) == ID.end()){
            ID[pair_kinship.id2] = ID.size();
            NUM++;
        }// end if id2
        // get the list of possible dc pairs
        if (pair_kinship.degree == 2 && pair_kinship.ibd2>ibd2_thres){
            dc_list.push_back(ID[pair_kinship.id1]);
            dc_list.push_back(ID[pair_kinship.id2]);
        }
        
        int pairFlag = (pair_kinship.kinship >= kin_lw && pair_kinship.kinship <= kin_up && pair_kinship.ibd2<=ibd2_thres );
        
        if (pairFlag){
            if(find(dc_list.begin(), dc_list.end(), ID[pair_kinship.id1]) == dc_list.end() && find(dc_list.begin(), dc_list.end(), ID[pair_kinship.id2]) == dc_list.end()){
                    pair.id1 = pair_kinship.id1;
                    pair.id2 = pair_kinship.id2;
                    pair.kinship = pair_kinship.kinship;
                    pair.is_pc = 0;
                setOfPair.push_back(pair);
            }
        }
        
        if (pc && pair_kinship.degree == 1 ){
            double ibd0 = 1 - pair_kinship.kinship*4 + pair_kinship.ibd2;
        
            if(ibd0 < ibd0_thres){
                pair.id1 = pair_kinship.id1;
                    pair.id2 = pair_kinship.id2;
                    pair.kinship = pair_kinship.kinship;
                    pair.is_pc = 1;
                setOfPair.push_back(pair);
            }
        }
        int relativeFlag;
        if (deg_cutoff!=6){
            kin_rel_lw = pow(2.0, -(deg_cutoff+1.5));
        }
        relativeFlag = (pair_kinship.kinship>=kin_rel_lw && pair_kinship.kinship<=kin_rel_up);

        if(relativeFlag){
            setOfRelatives.push_back(pair_kinship);
        }
    }// end while read in second pairs list

   //************************************************************************

   // const int NUM1=ID.size();
    unordered_map<int,set<int>> relatives;
   
    for(vector<pairRelated>::iterator it3 = setOfRelatives.begin(); it3 != setOfRelatives.end(); ++it3){
        // id1 in the ID and is second degree samples
        pair.id1 = (*it3).id1;
        pair.id2 = (*it3).id2;
        pair.kinship = (*it3).kinship;

        //if pair in second degree set
        if(find(setOfPair.begin(), setOfPair.end(), pair) != setOfPair.end()){
            continue;
        }
        
        relatives[ID[(*it3).id1]].insert(ID[(*it3).id2]);
        relatives[ID[(*it3).id2]].insert(ID[(*it3).id1]);
        
    }// end for loop

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
        assert(sample.value1 <= sample.value2  && "The IBD segment must have the start position in front of the end position");
        // only read in segments for samples in ID
        if (ID.find(sample.id1) == ID.end() || ID.find(sample.id2) == ID.end()){
            continue;
        }
        int row = ID[sample.id1]-1;
        int col = ID[sample.id2]-1;

        double length;
        segment oneSegment;
        
        oneSegment.chrom = sample.chrom;
        oneSegment.value1 = sample.value1;
        oneSegment.value2 = sample.value2;

        
        IBDsegment[row][col].segmentVec.push_back(oneSegment);
        length = abs(oneSegment.value2 - oneSegment.value1);
  
        IBDsegment[row][col].total += length;
        
    }//end read in lines

    dataBlock x1y;
    dataBlock x2y;
    
    
    vector<secondPair>::iterator it;
    for(it = setOfPair.begin(); it != setOfPair.end(); ++it){
        int i = ID[(*it).id1]-1;
        int j = ID[(*it).id2]-1;

        vector<int> commonRelatives;
        
        set_intersection(relatives[ID[(*it).id1]].begin(), relatives[ID[(*it).id1]].end(), relatives[ID[(*it).id2]].begin(), relatives[ID[(*it).id2]].end(), back_inserter(commonRelatives));
        
        vector<hap> haplotypes;
        dataBlock ibd_rel;
        vector<MutRelative> relativeVec;
        
        if (commonRelatives.size() !=0){
            // sort relatives
            
            for (vector<int>::iterator it1 = commonRelatives.begin(); it1!=commonRelatives.end(); ++it1){
                // the index of relative
                int k = (*it1)-1;
                //get coverage of each relative
                MutRelative rel;
                dataBlock x1y;
                dataBlock x2y;
            
                rel.id=k+1;
    
                if (IBDsegment[k][i].total == 0){
                    x1y = IBDsegment[i][k];
                }else{
                    x1y = IBDsegment[k][i];
                }
                if (IBDsegment[k][j].total == 0){
                    x2y = IBDsegment[j][k];
                }else {
                    x2y = IBDsegment[k][j];
                }
                rel.intersect = intersection(x1y,x2y);
                rel.coverage = max(x1y.total, x2y.total);
                rel.ibd_x1 = x1y;
                rel.ibd_x2 = x2y;
               
              //  if (max(x1y.total, x2y.total)/ min(x1y.total, x2y.total) <= 8){
                    relativeVec.push_back(rel);
             //   }
            }
            sort (relativeVec.begin(), relativeVec.end(), sortByCoverage);
            //take the union of relatives
            
            dataBlock tmp_intersect;
            for (vector<MutRelative>::iterator it2 = relativeVec.begin(); it2!=relativeVec.end(); ++it2){
                
                int h = (*it2).id-1;
                tmp_intersect = combine((*it2).intersect, tmp_intersect);
                // initialize a temporary haplotype
                if (haplotypes.empty()){
                    hap temp_hap1;
                    temp_hap1.rel_list.push_back(h+1);
                   // temp_hap1.intersect=(*it2).intersect;
                    temp_hap1.ibd_union1=(*it2).ibd_x1;
                    temp_hap1.ibd_union2=(*it2).ibd_x2;
                    haplotypes.push_back(temp_hap1);
    //                (*it2).indicator = 1;
                }else{
                
                int indicator = 0;
               
                dataBlock tmp1;
                dataBlock tmp2;
                dataBlock ibd_check1;
                dataBlock ibd_check2;
                vector<int> list;
                dataBlock ibd_x1_update = (*it2).ibd_x1;
                dataBlock ibd_x2_update = (*it2).ibd_x2;
                for (vector<hap>::iterator it1 = haplotypes.begin(); it1 !=haplotypes.end(); ++it1){
                    tmp1 = intersection(ibd_x1_update, (*it1).ibd_union1);
                    tmp2 = intersection(ibd_x2_update, (*it1).ibd_union2);
        
                        list = (*it1).rel_list;
                        dataBlock tmp;
                        for (vector<int>::iterator r = list.begin(); r != list.end(); ++r){
                            int l = (*r) -1;
             
                            dataBlock sharing;
                            if (IBDsegment[l][h].total == 0){
                                 sharing = IBDsegment[h][l];
                            }else{
                                sharing = IBDsegment[l][h];
                            }
                            tmp = combine(tmp, sharing);
                        } // end for loop of relatives in this haplotype
                    // check whether ibd 011
                    ibd_check1 = intersection(tmp1, tmp);
                    ibd_check2 = intersection(tmp2, tmp);
                    
                    ibd_x1_update=find_ibd011(ibd_x1_update,ibd_check1 );
                    ibd_x2_update=find_ibd011(ibd_x2_update, ibd_check2);
                    
                        if(max(ibd_x1_update.total, ibd_x2_update.total) <=10){
                        (*it1).rel_list.push_back(h+1);
                        (*it1).ibd_union1=combine( (*it1).ibd_union1, ibd_x1_update);
                        (*it1).ibd_union2=combine( (*it1).ibd_union2, ibd_x2_update);
                        indicator = 1;
                        break;
                    } // end if check the intersection
                    
                }// end for loop of haplotypes
                // ibd 011 case
                if (indicator==0){
                        hap temp_hap1;
                        temp_hap1.rel_list.push_back(h+1);
                        temp_hap1.ibd_union1=ibd_x1_update;
                        temp_hap1.ibd_union2=ibd_x2_update;
                        haplotypes.push_back(temp_hap1);
                        
                    }
                }// end else
            }// end for loop of common relatives
        double n = tmp_intersect.total;
        double d1 = 0;
        double d2 = 0;
        
    for (vector<hap>::iterator it3 = haplotypes.begin(); it3 !=haplotypes.end(); ++it3){
          
            d1 += (*it3).ibd_union1.total;
            d2 += (*it3).ibd_union2.total;
               
       } //end for loop of haplotypes
        (*it).coverage = max(d1,d2);
        (*it).ratio1 = n / d1;
        (*it).ratio2 = n / d2;
        }else{
            (*it).coverage = 0;
            (*it).ratio1 = 0 ;
            (*it).ratio2 = 0;
        }// end if relatives not empty
//}
        out << (*it);
    } // end for loop over 2nd degree pairs
   
    IBDdata.close();
    out.close();
    
    return 0;
}
