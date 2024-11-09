#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

using std::FILE;
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ifstream;
using std::stringstream;


static const unsigned int UNINITIALIZED = -1;



class Cache_line{
public:
    unsigned last_used = UNINITIALIZED;
    bool dirty_line = 0;
    unsigned tag = UNINITIALIZED;
    unsigned address = UNINITIALIZED;


};

class Cache_Level{
public:
    unsigned int size;
    unsigned int waysNum;
    unsigned int Cyc;
    unsigned int linesNum;

    std::vector<std::vector<Cache_line>> cache_lines;

    //Cache_Level() ;
    Cache_Level(unsigned int Size, unsigned int Cyc, unsigned int waysNum)  {
        this->size = Size;
        this->Cyc = Cyc;
        this->waysNum = waysNum;
        this->linesNum = 0;

    }
};


class Cache{
public:
    Cache_Level* l1;
    Cache_Level* l2;
    unsigned int bSize;
    unsigned int memCyc;
    bool writeAlloc;
    unsigned int l1_waysNum;
    unsigned int l2_waysNum;
    double accTime;
    int l1_hit;
    unsigned int l2_hit;
    unsigned int l1_misses;
    unsigned int l2_misses;
    unsigned int count;

    Cache(unsigned int memCyc, unsigned int bSize, bool writeAlloc, unsigned int l1_size, unsigned int l1_assoc, unsigned int l1_cyc,
          unsigned int l2_size, unsigned int l2_assoc, unsigned int l2_cyc) {
        this->memCyc = memCyc;
        this->count = 0;
        // Convert block size from log2 value to actual size
        int actualBlockSize = pow(2, bSize);
        this->bSize = actualBlockSize;
        this->writeAlloc = writeAlloc;

        // Convert sizes and associativities from log2 values to actual values
        int l1_actualSize = pow(2, l1_size);
        int l1_actualAssoc = pow(2, l1_assoc);
        int l2_actualSize = pow(2, l2_size);
        int l2_actualAssoc = pow(2, l2_assoc);

        // Calculate the number of sets for each cache level
        int l1_numSets = l1_actualSize / (actualBlockSize * l1_actualAssoc);
        int l2_numSets = l2_actualSize / (actualBlockSize * l2_actualAssoc);

        // Initialize cache levels
        l1 = new Cache_Level(l1_actualSize, l1_cyc, l1_actualAssoc);
        l2 = new Cache_Level(l2_actualSize, l2_cyc, l2_actualAssoc);

        // Set ways and other stats
        l1_waysNum = l1_actualAssoc;
        l2_waysNum = l2_actualAssoc;
        accTime = 0;
        l1_misses = 0;
        l2_misses = 0;
        l1_hit = 0;
        l2_hit = 0;

        // Initialize cache lines for L1 and L2 caches
        l1->linesNum = l1_numSets;
        l2->linesNum = l2_numSets;
        for (int i = 0; i < l1_actualAssoc; ++i) {
            l1->cache_lines.push_back(std::vector<Cache_line>(l1->linesNum));
        }

        for (int i = 0; i < l2_actualAssoc; ++i) {
            l2->cache_lines.push_back(std::vector<Cache_line>(l2->linesNum));
        }
    }


    ~Cache() {
        delete l1;
        delete l2;
    }
    unsigned int Log2n(unsigned int n)
    {
        return (n > 1) ? 1 + Log2n(n / 2) : 0;
    }
    unsigned int extractSet(Cache_Level* l, unsigned int input){

        unsigned int offsetBits = Log2n(this->bSize);
        unsigned int numberOfSets = l->linesNum;
        unsigned int setBits = Log2n(numberOfSets);

        unsigned int mask = (1 << setBits) - 1;
        unsigned int set = (input >> offsetBits) & mask;

        return set;
    }
    unsigned int extractTag(Cache_Level* l,unsigned int input)
    {
        unsigned int tag;
        tag = input>>(Log2n(this->bSize) + Log2n(l->linesNum));
        return tag;
    }
    void cacheReadUpdate(unsigned int op)
    {
        count++;
        unsigned int tag1 = extractTag(l1,op);
        unsigned int tag2 = extractTag(l2,op);
        unsigned int set1 = extractSet(l1,op);
        unsigned int set2 = extractSet(l2,op);


        unsigned int l1_lru = 0; //we only need the index of the lru way
        unsigned int l2_lru = 0;

        accTime += l1->Cyc;
        for (int i=0; i< this->l1_waysNum; i++)
        {

            if(this->l1->cache_lines[i][set1].tag == tag1)
            {
                //hit in l1
                l1_hit++;
                l1->cache_lines[i][set1].last_used = accTime;
                return;
            }

            if( (l1->cache_lines[i][set1].address != UNINITIALIZED) && (l1->cache_lines[i][set1].last_used < l1->cache_lines[l1_lru][set1].last_used))
            {
                l1_lru = i;
            }

        }
        //miss in l1
        l1_misses++;

        accTime +=l2->Cyc;
        unsigned int count_l2_updates = 0;
        for (int i=0; i< this->l2_waysNum; i++)
        {
            if(this->l2->cache_lines[i][set2].tag == tag2)
            {
                //miss in l1 hit in l2
                l2_hit++;
                l2->cache_lines[i][set2].last_used = accTime;
                count_l2_updates++;

                //bring over to l1
                for (int j = 0; j < l1_waysNum; j++){
                    if (l1->cache_lines[j][set1].address == UNINITIALIZED){

                        l1->cache_lines[j][set1].address = op;
                        l1->cache_lines[j][set1].tag = extractTag(l1, op);
                        l1->cache_lines[j][set1].dirty_line = false;
                        l1->cache_lines[j][set1].last_used = accTime;
                        return;
                    }
                }


                // if we reached this point that means we must evict and in case of dirty we need to write into higher level cache before we evict
                if (l1->cache_lines[l1_lru][set1].dirty_line){
                    unsigned int l2_dirty_set = extractSet(l2, l1->cache_lines[l1_lru][set1].address);

                    for (unsigned int k = 0; k < l2_waysNum; k++){
                        if ((l2->cache_lines[i][set2].address != UNINITIALIZED) && (l2->cache_lines[k][l2_dirty_set].address == l1->cache_lines[l1_lru][set1].address) ){
                            //found the line that we need to update in l2 due to dirty evicted line in l1
                            l2->cache_lines[k][l2_dirty_set].dirty_line = true;
                            l2->cache_lines[k][l2_dirty_set].last_used = accTime+count_l2_updates;
                            count_l2_updates++;

                        }

                    }

                }
                l1->cache_lines[l1_lru][set1].address = op;
                l1->cache_lines[l1_lru][set1].tag = extractTag(l1, op);
                l1->cache_lines[l1_lru][set1].dirty_line = false;
                l1->cache_lines[l1_lru][set1].last_used = accTime;


                return;

            }


            if( (l2->cache_lines[i][set2].address != UNINITIALIZED) && (l2->cache_lines[i][set2].last_used < l2->cache_lines[l2_lru][set2].last_used))
            {
                l2_lru = i;
            }

        }

        l2_misses++;
        accTime += memCyc;

        // First, attempt to bring the line into L2
        bool l2_line_allocated = false;
        for (int i = 0; i < this->l2_waysNum; i++)
        {
            if (l2->cache_lines[i][set2].address == UNINITIALIZED)
            {
                // Allocate the new line in L2
                l2->cache_lines[i][set2].address = op;
                l2->cache_lines[i][set2].tag = tag2;
                l2->cache_lines[i][set2].dirty_line = false;
                l2->cache_lines[i][set2].last_used = accTime;
                count_l2_updates++;
                l2_line_allocated = true;
                break;
            }
        }

        // If L2 is fully occupied, use the LRU line
        if (!l2_line_allocated)
        {

            //snoop L1 to check if the evacuated address is empty.
            for (unsigned int k = 0; k < l1_waysNum; k++)
            {
                unsigned int l1_invalid_set = extractSet(l1, l2->cache_lines[l2_lru][set2].address);

                if ((l1->cache_lines[k][l1_invalid_set].address != UNINITIALIZED) && (l1->cache_lines[k][l1_invalid_set].address == l2->cache_lines[l2_lru][set2].address))
                {
                    //found the relevant line in L1. We need to invalidate it.
                    l1->cache_lines[k][l1_invalid_set].dirty_line = false;
                    l1->cache_lines[k][l1_invalid_set].last_used = UNINITIALIZED;
                    l1->cache_lines[k][l1_invalid_set].address = UNINITIALIZED;
                    l1->cache_lines[k][l1_invalid_set].tag = UNINITIALIZED;
                    break;
                }
            }

            l2->cache_lines[l2_lru][set2].address = op;
            l2->cache_lines[l2_lru][set2].tag = tag2;
            l2->cache_lines[l2_lru][set2].dirty_line = false;
            l2->cache_lines[l2_lru][set2].last_used = accTime;
            count_l2_updates++;
        }

        // Then, try to bring the line into L1
        bool l1_line_allocated = false;
        for (int j = 0; j < l1_waysNum; j++)
        {
            if (l1->cache_lines[j][set1].address == UNINITIALIZED)
            {
                l1->cache_lines[j][set1].address = op;
                l1->cache_lines[j][set1].tag = tag1;
                l1->cache_lines[j][set1].dirty_line = false;
                l1->cache_lines[j][set1].last_used = accTime;
                l1_line_allocated = true;
                break;
            }
        }

        // If L1 is fully occupied, use the LRU line
        if (!l1_line_allocated)
        {
            // Check and handle dirty line before eviction if necessary
            if (l1->cache_lines[l1_lru][set1].dirty_line)
            {
                unsigned int l2_dirty_set = extractSet(l2, l1->cache_lines[l1_lru][set1].address);
                for (unsigned int k = 0; k < l2_waysNum; k++)
                {
                    if ((l2->cache_lines[k][l2_dirty_set].address != UNINITIALIZED) && (l2->cache_lines[k][l2_dirty_set].address == l1->cache_lines[l1_lru][set1].address))
                    {
                        l2->cache_lines[k][l2_dirty_set].dirty_line = true;
                        l2->cache_lines[k][l2_dirty_set].last_used = accTime + count_l2_updates;
                        count_l2_updates++;
                        break;
                    }
                }
            }

            l1->cache_lines[l1_lru][set1].address = op;
            l1->cache_lines[l1_lru][set1].tag = tag1;
            l1->cache_lines[l1_lru][set1].dirty_line = false;
            l1->cache_lines[l1_lru][set1].last_used = accTime;
        }

    }
    void cacheWriteUpdate(unsigned int op)
    {

        count++;
        unsigned int tag1 = extractTag(l1,op);
        unsigned int tag2 = extractTag(l2,op);
        unsigned int set1 = extractSet(l1,op);
        unsigned int set2 = extractSet(l2,op);


        unsigned int l1_lru = 0; //we only need the index of the lru way
        unsigned int l2_lru = 0;

        accTime += l1->Cyc;
        for (int i=0; i< this->l1_waysNum; i++)
        {

            if(this->l1->cache_lines[i][set1].tag == tag1)
            {
                //hit in l1
                l1_hit++;
                l1->cache_lines[i][set1].dirty_line = true;
                l1->cache_lines[i][set1].last_used = accTime;
                return;
            }

            if( (l1->cache_lines[i][set1].address != UNINITIALIZED) && (l1->cache_lines[i][set1].last_used < l1->cache_lines[l1_lru][set1].last_used))
            {
                l1_lru = i;
            }

        }
        //miss in l1
        l1_misses++;

        unsigned int count_l2_updates = 0;
        accTime +=l2->Cyc;
        for (int i=0; i< this->l2_waysNum; i++)
        {
            if(this->l2->cache_lines[i][set2].tag == tag2)
            {
                //miss in l1 hit in l2
                l2_hit++;
                l2->cache_lines[i][set2].last_used = accTime;
                if (!writeAlloc) l2->cache_lines[i][set2].dirty_line = true;
                count_l2_updates++;


                //in case we're in writeAlloc we need to bring the line into l1
                if (writeAlloc)
                {
                    for (int j = 0; j < l1_waysNum; j++){
                        if (l1->cache_lines[j][set1].address == UNINITIALIZED){

                            l1->cache_lines[j][set1].address = l2->cache_lines[i][set2].address;
                            l1->cache_lines[j][set1].tag = extractTag(l1, l2->cache_lines[i][set2].address);
                            l1->cache_lines[j][set1].dirty_line = true;
                            l1->cache_lines[j][set1].last_used = accTime;
                            return;
                        }
                    }


                    // if we reached this point that means we must evict and in case of dirty we need to write into higher level cache before we evict
                    if (l1->cache_lines[l1_lru][set1].dirty_line){
                        unsigned int l2_dirty_set = extractSet(l2, l1->cache_lines[l1_lru][set1].address);

                        for (unsigned int k = 0; k < l2_waysNum; k++){
                            if ((l2->cache_lines[i][set2].address != UNINITIALIZED) && (l2->cache_lines[k][l2_dirty_set].address == l1->cache_lines[l1_lru][set1].address) ){
                                //found the line that we need to update in l2 due to dirty evicted line in l1
                                l2->cache_lines[k][l2_dirty_set].dirty_line = true;
                                l2->cache_lines[k][l2_dirty_set].last_used = accTime+count_l2_updates;
                                count_l2_updates++;

                            }

                        }

                    }
                    l1->cache_lines[l1_lru][set1].address = l2->cache_lines[i][set2].address;
                    l1->cache_lines[l1_lru][set1].tag = extractTag(l1, l2->cache_lines[i][set2].address);
                    l1->cache_lines[l1_lru][set1].dirty_line = true;
                    l1->cache_lines[l1_lru][set1].last_used = accTime;

                }
                return;

            }


            if( (l2->cache_lines[i][set2].address != UNINITIALIZED) && (l2->cache_lines[i][set2].last_used < l2->cache_lines[l2_lru][set2].last_used))
            {
                l2_lru = i;
            }

        }

        l2_misses++;
        accTime += memCyc;

        // if we had 2 misses and it's a writeAlloc we need to bring the line into L2 and L1
        if(writeAlloc)
        {
            // First, attempt to bring the line into L2
            bool l2_line_allocated = false;
            for (int i = 0; i < this->l2_waysNum; i++)
            {
                if (l2->cache_lines[i][set2].address == UNINITIALIZED)
                {
                    // Allocate the new line in L2
                    l2->cache_lines[i][set2].address = op;
                    l2->cache_lines[i][set2].tag = tag2;
                    l2->cache_lines[i][set2].dirty_line = false;
                    l2->cache_lines[i][set2].last_used = accTime;
                    count_l2_updates++;
                    l2_line_allocated = true;
                    break;
                }
            }

            // If L2 is fully occupied, use the LRU line
            if (!l2_line_allocated)
            {

                //snoop L1 to check if the evacuated address is empty.
                for (unsigned int k = 0; k < l1_waysNum; k++)
                {
                    unsigned int l1_invalid_set = extractSet(l1, l2->cache_lines[l2_lru][set2].address);

                    if ((l1->cache_lines[k][l1_invalid_set].address != UNINITIALIZED) && (l1->cache_lines[k][l1_invalid_set].address == l2->cache_lines[l2_lru][set2].address))
                    {
                        //found the relevant line in L1. We need to invalidate it.
                        l1->cache_lines[k][l1_invalid_set].dirty_line = false;
                        l1->cache_lines[k][l1_invalid_set].last_used = UNINITIALIZED;
                        l1->cache_lines[k][l1_invalid_set].address = UNINITIALIZED;
                        l1->cache_lines[k][l1_invalid_set].tag = UNINITIALIZED;
                        break;
                    }
                }

                l2->cache_lines[l2_lru][set2].address = op;
                l2->cache_lines[l2_lru][set2].tag = tag2;
                l2->cache_lines[l2_lru][set2].dirty_line = false;
                l2->cache_lines[l2_lru][set2].last_used = accTime+count_l2_updates;
                count_l2_updates++;
            }

            // Then, try to bring the line into L1
            bool l1_line_allocated = false;
            for (int j = 0; j < l1_waysNum; j++)
            {
                if (l1->cache_lines[j][set1].address == UNINITIALIZED)
                {
                    l1->cache_lines[j][set1].address = op;
                    l1->cache_lines[j][set1].tag = tag1;
                    l1->cache_lines[j][set1].dirty_line = true;
                    l1->cache_lines[j][set1].last_used = accTime;
                    l1_line_allocated = true;
                    break;
                }
            }

            // If L1 is fully occupied, use the LRU line
            if (!l1_line_allocated)
            {
                // Check and handle dirty line before eviction if necessary
                if (l1->cache_lines[l1_lru][set1].dirty_line)
                {
                    unsigned int l2_dirty_set = extractSet(l2, l1->cache_lines[l1_lru][set1].address);
                    for (unsigned int k = 0; k < l2_waysNum; k++)
                    {
                        if ((l2->cache_lines[k][l2_dirty_set].address != UNINITIALIZED) && (l2->cache_lines[k][l2_dirty_set].address == l1->cache_lines[l1_lru][set1].address))
                        {
                            l2->cache_lines[k][l2_dirty_set].dirty_line = true;
                            l2->cache_lines[k][l2_dirty_set].last_used = accTime+count_l2_updates;
                            count_l2_updates++;
                            break;
                        }
                    }
                }

                l1->cache_lines[l1_lru][set1].address = op;
                l1->cache_lines[l1_lru][set1].tag = tag1;
                l1->cache_lines[l1_lru][set1].dirty_line = true;
                l1->cache_lines[l1_lru][set1].last_used = accTime;
            }
        }

    }





    void cache_func(unsigned int op, char operation){
        double total_time1;
        double total_time2;
        if(operation == 'w')
        {
            cacheWriteUpdate(op);
        }
        else
        {
            cacheReadUpdate(op);
        }
        //printf("\n misses = ******************8 %d , %d \n", l1_misses , l2_misses);
        total_time1 = (l1_hit+l1_misses);
        total_time2 = (l2_hit+l2_misses);

        if(total_time1 > total_time2)
        {
            return;
        }

    }
};
int main(int argc, char **argv) {

    if (argc < 19) {
        cerr << "Not enough arguments" << endl;
        return 0;
    }

    // Get input arguments

    // File
    // Assuming it is the first argument
    char* fileString = argv[1];
    ifstream file(fileString); //input file stream
    string line;
    if (!file || !file.good()) {
        // File doesn't exist or some other error
        cerr << "File not found" << endl;
        return 0;
    }

    unsigned MemCyc = 0, BSize = 0, L1Size = 0, L2Size = 0, L1Assoc = 0,
            L2Assoc = 0, L1Cyc = 0, L2Cyc = 0, WrAlloc = 0;

    for (int i = 2; i < 19; i += 2) {
        string s(argv[i]);
        if (s == "--mem-cyc") {
            MemCyc = atoi(argv[i + 1]);
        } else if (s == "--bsize") {
            BSize = atoi(argv[i + 1]);
        } else if (s == "--l1-size") {
            L1Size = atoi(argv[i + 1]);
        } else if (s == "--l2-size") {
            L2Size = atoi(argv[i + 1]);
        } else if (s == "--l1-cyc") {
            L1Cyc = atoi(argv[i + 1]);
        } else if (s == "--l2-cyc") {
            L2Cyc = atoi(argv[i + 1]);
        } else if (s == "--l1-assoc") {
            L1Assoc = atoi(argv[i + 1]);
        } else if (s == "--l2-assoc") {
            L2Assoc = atoi(argv[i + 1]);
        } else if (s == "--wr-alloc") {
            WrAlloc = atoi(argv[i + 1]);
        } else {
            cerr << "Error in arguments" << endl;
            return 0;
        }
    }
    Cache* cache = new Cache(MemCyc, BSize, WrAlloc, L1Size, L1Assoc, L1Cyc, L2Size, L2Assoc, L2Cyc);
    double total_time_l1 = 0;
    double total_time_l2 = 0;
    double total_time_mem = 0;
    double total_time = 0;
    while (getline(file, line)) {

        stringstream ss(line);
        string address;
        char operation = 0; // read (R) or write (W)
        if (!(ss >> operation >> address)) {
            // Operation appears in an Invalid format
            cout << "Command Format error" << endl;
            return 0;
        }

        // DEBUG - remove this line
        //cout << "operation: " << operation;

        string cutAddress = address.substr(2); // Removing the "0x" part of the address

        // DEBUG - remove this line
        //cout << ", address (hex)" << cutAddress;

        unsigned long int num;
        num = strtoul(cutAddress.c_str(), NULL, 16);
        cache->cache_func(num, operation);
        // DEBUG - remove this line
        //cout << " (dec) " << num << endl;

    }

    total_time_l1 = (cache->l1_hit + cache->l1_misses);
    total_time_l2 = (cache->l2_hit + cache->l2_misses);
    total_time_mem = cache->l2_misses;

    double L1MissRate = cache->l1_misses/total_time_l1;
    double L2MissRate = cache->l2_misses/total_time_l2;
    double avgAccTime = cache->accTime/cache->count;

    printf("L1miss=%.03f ", L1MissRate);
    printf("L2miss=%.03f ", L2MissRate);
    printf("AccTimeAvg=%.03f\n", avgAccTime);

    return 0;
}


