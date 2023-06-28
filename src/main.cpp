//=================================================================================================
// sensor_frame_gen
//
// Author: Doug Wolf
//
// Command line options:
//
//   -config <filename>      : specifies the name of a configuration file
//
//   -trace <cell_number>    : instead of creating an output file, traces a cell in an existing 
//                             file.
//
//   -dict                   : instead of creating an output file, display data dictionary
//
//   -load <filename> <addr> <size_limit>
//                           : instead of creating output file, loads a file into the specified
//                             RAM physical address
//
//=================================================================================================

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <exception>
#include <fcntl.h>
#include <cstdarg>
#include <cstring>
#include <iostream>
#include <memory>
#include <map>
#include <vector>
#include <fstream>
#include "config_file.h"
#include "PhysMem.h"
#include "changelog.h"

using namespace std;
 
void     execute(const char** argv);
void     loadNucleotides();
void     loadFragments();
void     loadDistribution();
uint32_t findLongestSequence();
uint32_t verifyDistributionIsValid();
void     writeOutputFile(uint32_t frameGroupCount);
void     parseCommandLine(const char** argv);
void     trace(uint32_t cellNumber);
void     readConfigurationFile(string filename);
void     loadFile(string filename, string address);
void     printDictionary();
uint64_t stringTo64(const string& str);

// Define a convenient type to encapsulate a vector of strings
typedef vector<string> strvec_t;

// Contains nucleic acid fragement definitions
map<string, vector<string>> fragment;

// Contains nucleotide definitions
map<string, vector<int>> nucleotide;

// This list defines each fragment distribution in the distribution definitions file
struct distribution_t
{
    int      first, last, step;
    strvec_t cellValue;
};
vector<distribution_t> distributionList;

// This object maps physical RAM address into userspace.
PhysMem RAM;

// This is the number of cells in a single data row on the chip
const int ROW_SIZE = 2048;


//=================================================================================================
// Command line options
//=================================================================================================
struct cmdline_t
{
    bool     load;
    string   filename;
    string   address;
    string   sizeLimit;
    
    bool     trace;
    uint32_t cellNumber;
    
    bool     dict;
    
    string   config;
} cmdLine;
//=================================================================================================


//=================================================================================================
// Variable names in this structure should exactly match the configuration file
//=================================================================================================
struct config_t
{
    uint32_t         adc_per_nucleotide;
    uint64_t         random_seed;
    uint32_t         cells_per_frame;
    uint64_t         ring_buffer_size;
    uint32_t         data_frames;
    uint8_t          filler_value;
    string           nucleotide_file;
    string           fragment_file;
    string           distribution_file;
    string           output_file;

} config;
//=================================================================================================



//=================================================================================================
// main() - Execution starts here.
//=================================================================================================
int main(int argc, const char** argv)
{
    // Tell the world what software version we are
    printf("Version %s\n", VERSION_REV);

    try
    {
        execute(argv);
    }
    catch(const std::exception& e)
    {
        cerr << e.what() << '\n';
    }
}
//=================================================================================================


//=================================================================================================
// concatVec() - Helper functions that concatenate vectors together
//=================================================================================================
void concatVec(vector<int>&    v1, vector<int>&    v2) {v1.insert(v1.end(), v2.begin(), v2.end());}
void concatVec(vector<string>& v1, vector<string>& v2) {v1.insert(v1.end(), v2.begin(), v2.end());}
//=================================================================================================

//=================================================================================================
// throwRuntime() - Throws a runtime exception
//=================================================================================================
static void throwRuntime(const char* fmt, ...)
{
    char buffer[1024];
    va_list ap;
    va_start(ap, fmt);
    vsprintf(buffer, fmt, ap);
    va_end(ap);

    throw runtime_error(buffer);
}
//=================================================================================================



//=================================================================================================
// radix() - Examine an input string and return a 16 if the first two characters are "0x" or "0X",
//           otherwise it returns 10.
//=================================================================================================
int radix(const string& str)
{
    // Get a pointer to the characters
    const char* p = str.c_str(); 

    // Skip past any spaces or tabs
    while (*p == 32 || *p == 9) ++p;

    // Tell the caller whether this string is radix 10 or radix 16
    return (p[0] == '0' && (p[1] == 'x' || p[1] == 'X')) ? 16 : 10;  
}
//=================================================================================================


//=================================================================================================
// to_int() - Converts an ASCII string to an integer
//=================================================================================================
int to_int(const char* str)
{
    while (*str == 32 || *str == 9) ++str;
    if (*str == 0) return 0;
    return stoi(str, nullptr, radix(str));
}
//=================================================================================================


//=================================================================================================
// showHelp() - Display help-text and quit
//=================================================================================================
void showHelp()
{
    printf
    (
        "Usage:\n"
        "  sfg [-config <filename>]\n"
        "  sfg -trace <cell_number>\n"
        "  sfg -dict\n"
        "  sfg -load <filename> <address> <size_limit>\n"
        "\n"
        "  <address> and <size_limit> may be expressed in either decimal or hex, and may\n"
        "  include optional K, M, or G suffixes.   Verilog-style underscores are allowed\n"
        "  in hex values.\n"
    );

    // Terminate the program
    exit(0);
}
//=================================================================================================


//=================================================================================================
// parseCommandLine() - Parse the command line parameters into the -cmdLine structure
//=================================================================================================
void parseCommandLine(const char** argv)
{
    int i=0;

    // Loop through each command line parameter
    while (argv[++i])
    {
        // Fetch this parameter
        string token = argv[i];

        // Handle the "-help" command line switch
        if (token == "-help" || token == "-h" || token == "?")
            showHelp();

        // Handle "-load" command line switch
        if (token == "-load")
        {
            cmdLine.load = true;

            if (argv[i+1])
                cmdLine.filename = argv[++i];                
            else
                throwRuntime("Missing filename on -load");                

            if (argv[i+1])
                cmdLine.address = argv[++i];                
            else
                throwRuntime("Missing address on -load"); 

            if (argv[i+1])
                cmdLine.sizeLimit = argv[++i];                
            else
                throwRuntime("Missing size limit on -load"); 

            continue;
        }

        // Handle the "-trace" command line switch
        if (token == "-trace")
        {
            cmdLine.trace = true;
            if (argv[i+1])
                cmdLine.cellNumber = atoi(argv[++i]);                
            else
                throwRuntime("Missing parameter on -trace");                
            continue;
        }

        // Handle the "-dict" command line switch
        if (token == "-dict")
        {
            cmdLine.dict = true;
            continue;            
        }

        // Handle the "-config" command line switch
        if (token == "-config")
        {
            if (argv[i+1])
                cmdLine.config = argv[++i];
            else
                throwRuntime("Missing parameter on -config");
            continue;
        }

        printf("Illegal command line parameter '%s'\n", token.c_str());
        exit(1);
    }
}
//=================================================================================================




//=================================================================================================
// execute() - Top level code for program logic
//=================================================================================================
void execute(const char** argv)
{
    // Ensure that comma-separators get printed for numbers
    setlocale(LC_ALL, "");

    // Parse the command line
    parseCommandLine(argv);

    // If we're just loading a data-file into a contig-buffer, make it so
    if (cmdLine.load)
    {
        loadFile(cmdLine.filename, cmdLine.address);
        exit(0);
    }

    // Fetch the configuration values from the file and populate the global "config" structure
    readConfigurationFile(cmdLine.config);

    // Seed the random number generator
    srand(config.random_seed);

    // If we're supposed to trace a single cell, make it so
    if (cmdLine.trace)
    {
        trace(cmdLine.cellNumber);
        exit(0);
    }

    // Load the nucleotide definitions
    loadNucleotides();

    // Load the fragment definitions
    loadFragments();

    // Load the fragment sequence distribution definitions
    loadDistribution();

    // Find out how many frame groups we need to write to the output file
    uint32_t frameGroupCount = verifyDistributionIsValid();

    // Either print out the data dictionary, or create the output file
    if (cmdLine.dict)
        printDictionary();
    else
        writeOutputFile(frameGroupCount);
}
//=================================================================================================



//=================================================================================================
// getNextCommaSeparatedToken() - Fetches the next comma separated token from a line of text
//
// Passed:  p (by REFERENCE!) = A pointer to the text being parsed
//          token             = A pointer to where the extracted token should be stored
//
// Returns: true if a token was extracted, false if no more tokens available on the line
//
// Note: This routine will accomodate lines containing optional carriage-returns
//=================================================================================================
bool getNextCommaSeparatedToken(const char*& p, char* token)
{
    // Clear the caller's 'token' field in case we can't find a token
    *token = 0;

    // Skip over white-space
    while (*p == 32 || *p == 9) ++p;

    // If we've hit the end of the input line, tell the caller
    if (*p == 0 || *p == 10 || *p == 13) return false;

    // Extract the token into the buffer
    while (!(*p == 32 || *p == 9 || *p == 10 || *p == 13 || *p == 0 || *p == ',' || *p == '='))
    {
        *token++ = *p++;
    }

    // Nul-terminate the extracted token
    *token = 0;

    // Skip over any trailing whitespace
    while (*p == 32 || *p == 9) ++p;

    // If there is a comma or equal-sign, skip over it
    if (*p == ',' || *p == '=') ++p;

    // Tell the caller that they have extracted a token
    return true;
}
//=================================================================================================


//=================================================================================================
// getNextCommaSeparatedInt() - Similar to "getNextCommaSeparatedToken()", but converts the token
//                              to an integer
//=================================================================================================
bool getNextCommaSeparatedInt(const char*& p, int* pValue)
{
    char token[1000];

    // Fetch the next token
    bool status = getNextCommaSeparatedToken(p, token);    

    // Convert that token to an integer
    *pValue = stringTo64(token);

    // And tell the caller whether or not there was a token available
    return status;
}
//=================================================================================================



//=================================================================================================
// loadNucleotides() - Load nucleotide definitions into RAM
//
// On Exit: the global "nucleotide" object contains nucleotide definitions
//=================================================================================================
void loadNucleotides()
{
    char name[1000], buffer[1000];
    vector<int> v;
    string line;

    // Fetch the filename of the nucleotide definiton file
    const char* filename = config.nucleotide_file.c_str();

    // Open the input file
    ifstream file(filename);   

    // If we can't open the input file, complain
    if (!file.is_open()) throwRuntime("%s not found", filename);

    // Loop through each line of the input file
    while (getline(file, line))
    {
        // Get a pointer to the line of text we just read
        const char* p = line.c_str();

        // Skip over whitespace
        while (*p == 32 || *p == 9) ++p;

        // If the line is blank, skip it
        if (*p == 0 || *p == 10 || *p == 13) continue;

        // Any line starting with '#' is a comment
        if (*p == '#') continue;

        // Any line starting with '//' is a comment
        if (p[0] == '/' && p[1] == '/') continue;

        // Clear the fragment value vector
        v.clear();

        // Fetch the nucleotide name
        getNextCommaSeparatedToken(p, name);

        // If the name is blank, skip this line
        if (name[0] == 0) continue;

        // If the name is more than a single character, complain
        if (strlen(name) != 1) throwRuntime("Illegal nucleotide: %s", name);

        // Fetch every integer value after the name and stuff them into a vector
        while (getNextCommaSeparatedToken(p, buffer))
        {
            v.push_back(to_int(buffer));
        }

        // Save this nucleotide data into our global variable
        nucleotide[name] = v;
    }
}
//=================================================================================================




//=================================================================================================
// displayFragment() - Just prints out the specified fragement definition.   This routine is 
//                     useful when debugging this code.
//=================================================================================================
void displayFragment(const char* name)
{
    auto it = fragment.find(name); 
    if (it == fragment.end()) throwRuntime("Unknown fragment '%s'", name);
    auto& v= it->second;
    printf("%s: ", name);
    for (auto s : v) printf(" %s", s.c_str());
    printf("\n");
}
//=================================================================================================



//=================================================================================================
// readFragmentFromFile() - Reads a file, converts every byte of that file into an ASCII Hex value
//                          appends those ASCII hex values into a string vector, and returns the
//                          string vector
//=================================================================================================
strvec_t readFragmentFromFile(const char* filename)
{
    strvec_t retval;
    char buffer[20];      

    // Open the input file
    int fd = open(filename, O_RDONLY);

    // If we can't open the file, complain about it
    if (fd < 0) throwRuntime("Can't open fragment file '%s'", filename);

    // Find out how large the file is in bytes
    auto fileSize = lseek64(fd, 0, SEEK_END);

    // Allocate sufficient RAM to contain the entire file
    unique_ptr<uint8_t[]> fragData(new uint8_t[fileSize]);

    // Rewind the file, read the entire file into RAM, and close it
    lseek64(fd, 0, SEEK_SET);
    read(fd, fragData.get(), fileSize);
    close(fd);

    // Loop through every byte of data in the file...
    for (uint32_t i = 0; i<fileSize; ++i)
    {
        // Convert this byte to ASCII hex
        sprintf(buffer, "0x%02X", fragData[i]);

        // And add the ASCII hex version of this byte to our return vector
        retval.push_back(buffer);
    }

    // Hand the vector of ASCII data to the caller
    return retval;
}
//=================================================================================================


//=================================================================================================
// tokenToStringVec() - Breaks a token into a vector of strings
//
// If the token starts with a digit, the token is returned as is.
// Any nucleotide name in the token is its own element in the output vector
// Any fragment name in the token is turned into nucleotides in the output vector
//=================================================================================================
strvec_t tokenToStringVec(const char* token)
{
    char name[1000];
    vector<string> retval;

    // If this string is an integer, just append it to retval
    if (token[0] >= '0' && token[0] <= '9')
    {
        retval.push_back(token);
        return retval;
    }

    // If this string begins with a '@', remainder of the string is the
    // the filename of a binary file that contains the actual ADC data
    if (token[0] == '@')
    {
        return readFragmentFromFile(token+1);
    }

    // Point to the start of the token
    const char *in = token;

    // We're going to loop through each character of the input token
    while (*in)
    {
        // We'll be copying characters into array 'name'
        char* out = name;

        // Fetch the next input character
        int c = *in++;

        // If that input character was an open-paren, keep extracting characters
        // until we encounter the matching close-paren
        if (c == '(') while (true)
        {
            c = *in++;
            if (c == ')') break;
            if (c == 0) throwRuntime("Unbalanced parenthesis in fragment file");
            *out++ = c;
        } 

        // Otherwise, simply place that input character into 'name'
        else *out++ = c;

        // nul-terminate 'name'
        *out++ = 0;

        // Was the name we just extracted a nucleotide name?
        auto it1 = nucleotide.find(name);

        // If it was, then append the name of the nucleotide to 'retval'
        if (it1 != nucleotide.end())
        {
            for (int ii=0; ii<config.adc_per_nucleotide; ++ii)
                retval.push_back(name);
            continue;
        }

        // Was the name we just extracted a fragment name?
        auto it2 = fragment.find(name);

        // If it was, then append the fragment definition to 'retval'
        if (it2 != fragment.end())
        {
            concatVec(retval, it2->second);
            continue;                
        }

        // If we get here, we've encountered a fragment definition that
        // contains a symbol that is neither a nucleotide name, nor a fragment
        // name, nor an integer literal
        throwRuntime("Unknown fragment/nucleotide %s", name);
    }    

    // Hand the resulting vector of strings to the caller
    return retval;
}
//=================================================================================================


//=================================================================================================
// loadFragments() - Load fragment definitions into RAM
//
// On Exit: the global "fragment" object contains fragment definitions
//=================================================================================================
void loadFragments()
{
    char fragmentName[1000], buffer[1000];
    vector<string> v;
    string line;

    // Fetch the filename of the fragment definiton file
    const char* filename = config.fragment_file.c_str();

    // Open the input file
    ifstream file(filename);   

    // If we can't open the input file, complain
    if (!file.is_open()) throwRuntime("%s not found", filename);

    // Loop through each line of the input file
    while (getline(file, line))
    {
        // Get a pointer to the line of text we just read
        const char* p = line.c_str();

        // Skip over whitespace
        while (*p == 32 || *p == 9) ++p;

        // If the line is blank, skip it
        if (*p == 0 || *p == 10 || *p == 13) continue;

        // Any line starting with '#' is a comment
        if (*p == '#') continue;

        // Any line starting with '//' is a comment
        if (p[0] == '/' && p[1] == '/') continue;

        // Clear the fragment value vector
        v.clear();

        // Fetch the fragment name
        getNextCommaSeparatedToken(p, fragmentName);

        // If the fragment name is blank, skip this line
        if (fragmentName[0] == 0) continue;

        // Fragments are not allowed to share a name with a nucleotide
        if (nucleotide.find(fragmentName) != nucleotide.end())
        {
            throwRuntime("Fragment '%s' shares name with nucleotide", fragmentName);            
        }
    
        // Fetch every token on the line, and turn it into a list of strings
        // representing either nucleotides or integer literals
        while (getNextCommaSeparatedToken(p, buffer))
        {
            vector<string> ntides = tokenToStringVec(buffer);
            concatVec(v, ntides);
        }

        // Save this fragment data into our global variable
        fragment[fragmentName] = v;
    }
}
//=================================================================================================





//=================================================================================================
// dumpDistributionList() - Displays the distribution list for debugging purposes
//=================================================================================================
void dumpDistributionList()
{
    for (auto& r : distributionList)
    {
        printf("%i : %i : %i  *** ", r.first, r.last, r.step);
        for (auto s : r.cellValue) printf("%s  ", s.c_str());
        printf("\n");

    }
}
//=================================================================================================



//=================================================================================================
// loadDistribution() - Loads the fragment distribution definitions into RAM
//
// On Exit: the global "distribuitionList" object contains distribution definitions
//=================================================================================================
void loadDistribution()
{
    char fragmentName[1000];
    distribution_t distRecord;
    string line;

    // Get a handy reference to the vector of cell values in a distribution record
    auto& drcv = distRecord.cellValue;

    // Fetch the filename of the fragment distribiution definiton file
    const char* filename = config.distribution_file.c_str();

    // Open the input file
    ifstream file(filename);   

    // If we can't open the input file, complain
    if (!file.is_open()) throwRuntime("%s not found", filename);

     // Loop through each line of the input file
    while (getline(file, line))
    {
        // Get a pointer to the line of text we just read
        const char* p = line.c_str();

        // Skip over whitespace
        while (*p == 32 || *p == 9) ++p;

        // If the line is blank, skip it
        if (*p == 0 || *p == 10 || *p == 13) continue;

        // Any line starting with '#' is a comment
        if (*p == '#') continue;

        // Any line starting with '//' is a comment
        if (p[0] == '/' && p[1] == '/') continue;

        // Look for the '$' delimeter that begins a list of fragment IDs
        const char* delimeter = strchr((char*)p, '$');

        // If that delimeter doesn't exist, this isn't a valid distribution definition
        if (delimeter == nullptr) continue;

        // Replace the '$' with a carriage return to fake an "end of line"
        *(char*)delimeter++ = 13;

        // Just in case the user added a comma after the '$', consume it 
        while (*delimeter == 32 || *delimeter == 9) ++delimeter;
        if (*delimeter == ',') ++delimeter;

        // Get the first cell number, the last cell number, and the step-size
        getNextCommaSeparatedInt(p, &distRecord.first);
        getNextCommaSeparatedInt(p, &distRecord.last );
        getNextCommaSeparatedInt(p, &distRecord.step );

        // Ensure that the first cell number in the distribution is valid
        if (distRecord.first < 1 || distRecord.first > config.cells_per_frame)
        {
            throwRuntime("Invalid cell number %i", distRecord.first);
        }

        // If no "last cell" was specified, this distribution is just for the first cell
        if (distRecord.last == 0) distRecord.last = distRecord.first;

        // If no 'step' is specified, we're defining every cell from 'first' to 'last'
        if (distRecord.step == 0) distRecord.step = 1;

        // Clear the vector that will hold fragment data values
        drcv.clear();

        // Point to the comma separated fragement ids that come after the '$' delimeter
        p = delimeter;

        // Loop through every fragment name in the comma separated list...
        while (getNextCommaSeparatedToken(p, fragmentName))
        {
            // If we don't recognize this fragment name, complain
            if (fragment.find(fragmentName) == fragment.end())
            {
                throwRuntime("Undefined fragment name '%s'", fragmentName);
            }

            // Get a reference to the cell values for this fragment
            auto& fragcv = fragment[fragmentName];

            // Append the cell values for this fragment to the distribution record
            concatVec(drcv, fragcv);
        }

        // And add this distribution record to the distribution list
        distributionList.push_back(distRecord);
    }
}
//=================================================================================================



//=================================================================================================
// nucleotideToADC() - Return an ADC value that is valid for the specified nucleotide
//=================================================================================================
int nucleotideToADC(string name)
{
    // Get the name of the nucleotide as a c-string
    const char* cname = name.c_str();

    // If the name is an integer literal, return its value
    if (name[0] >= '0' && name[0] <= '9') return to_int(cname);
    
    // Does this nucleotide name exist?
    auto it = nucleotide.find(name);

    // If it doesn't return an empty result
    if (it == nucleotide.end()) throwRuntime("Unknown nucleotide '%s'", cname);

    // Get a handy reference to the integers that define this nucleotide
    auto& adc_values = it->second;

    // Select a random index into the vector 'adc_values'
    int idx = rand() % adc_values.size();

    // And return the ADC value at that random index
    return adc_values[idx];
}
//=================================================================================================



//=================================================================================================
// findLongestSequence() - Finds and returns the number of frames requires by the longest sequence
//                         in the distributionList
//=================================================================================================
uint32_t findLongestSequence()
{
    uint32_t longestLength = 0;

    // Loop through every record in the distribution list and keep
    // track of the length of the longest sequence of fragments we find    
    for (auto& distRec : distributionList)
    {
        // Keep track of the length of the longest sequence of fragments we find
        if (distRec.cellValue.size() > longestLength) longestLength = distRec.cellValue.size();
    };

    // Hand the caller the length of the longest sequence of fragments
    return longestLength;
}
//=================================================================================================


//=================================================================================================
// verifyDistributionIsValid() - Checks to make sure that number of frame groups implied by the
//                               longest fragement sequence will fit into the contiguous buffer.
//
// Returns: The number of frame groups that will be written to the output file
//=================================================================================================
uint32_t verifyDistributionIsValid()
{

    // Ensure that the number of cells in a single frame is a multiple of the row size
    if (config.cells_per_frame % ROW_SIZE != 0)
    {
        printf("\nConfig value 'cells_per_frame' must a multiple of %i\n", ROW_SIZE);
        exit(1);        
    }

    // What's the maximum number of frames that will fit into the contig buffer?
    uint32_t maxFrames = config.ring_buffer_size / config.cells_per_frame;

    // What is the maximum number of frames required by any fragment sequence?
    uint32_t longestSequence = findLongestSequence();

    // A "frame group" is a set of data frames.
    uint32_t frameGroupLength = config.data_frames;

    // How many frames groups are required to express our longest sequence?
    uint32_t frameGroupCount = longestSequence / config.data_frames + 1;

    // How many data frames are in 'frameGroupCount' frameg groups?
    uint32_t totalReqdFrames = frameGroupCount * frameGroupLength;

    // How many bytes will that number of frames occupy in the contiguous buffer?
    uint64_t totalContigReqd = (uint64_t)totalReqdFrames * (uint64_t)config.cells_per_frame;

    // Tell the user basic statistics about this run
    printf("%'16u Frames in the longest fragment sequence\n", longestSequence);
    printf("%'16u Frames in a frame group\n", frameGroupLength);
    printf("%'16u Frame group(s) required\n", frameGroupCount);
    printf("%'16u Frames will fit into the contiguous buffer\n", maxFrames);
    printf("%'16u Frames required in total\n", totalReqdFrames);
    printf("%'16lu Bytes required in total\n", totalContigReqd);

    // If the longest fragment sequence is too long to fit into the contiguous buffer,
    // complain and drop dead
    if (totalReqdFrames > maxFrames)
    {
        printf("\nThe specified fragment distribution won't fit into the contiguous buffer!\n");
        exit(1);
    }

    // Tell the caller how many frame groups we're going to output
    return frameGroupCount;
}
//=================================================================================================


//=================================================================================================
// buildDataFrame() - Uses the fragment-sequence distribution list to create a data frame
//=================================================================================================
void buildDataFrame(uint8_t* frame, uint32_t frameNumber)
{
    // Every cell in the frame starts out quiescient
    memset(frame, config.filler_value, config.cells_per_frame);

    // Loop through every distribution record in the distribution list
    for (auto& dr : distributionList)
    {
        // If this fragment sequence contains a value for this frame number...
        if (frameNumber < dr.cellValue.size())
        {
            // Populate the appropriate cells with the data value for this frame
            for (uint32_t cellNumber = dr.first-1; cellNumber < dr.last; cellNumber += dr.step)
            {
                string nucleotideName = dr.cellValue[frameNumber];
                frame[cellNumber] = nucleotideToADC(nucleotideName);
            }
        }
    }
}
//=================================================================================================


//=================================================================================================
// writeOutputFile() - Creates the output file
//=================================================================================================
void writeOutputFile(uint32_t frameGroupCount)
{
    uint32_t i, frameNumber = 0;

    // Fetch the name of the file we're going to create
    const char* filename = config.output_file.c_str();
   
    // Open the file we're going to write, and complain if we can't
    FILE* ofile = fopen(filename, "w");
    if (ofile == nullptr) throwRuntime("Can't create %s", filename);

    // Allocate sufficient RAM to contain an entire raw data frame
    unique_ptr<uint8_t> framePtr(new uint8_t[config.cells_per_frame]);

    // Get a pointer to the frame data
    uint8_t* frame  = framePtr.get();

    // Loop through each frame group
    for (int32_t frameGroup = 0; frameGroup < frameGroupCount; ++frameGroup)
    {
        // For each data frame in this frame group...
        for (i=0; i<config.data_frames; ++i)
        {
            // Build the raw data frame for this frame number
            buildDataFrame(frame, frameNumber++);
            
            // And write the resulting frame to the output file
            fwrite(frame, 1, config.cells_per_frame, ofile);
        }
    }

    // We're done with the output file
    fclose(ofile);
}
//=================================================================================================




//=================================================================================================
// trace() - Displays the value of a single cell for every frame in the output file
//=================================================================================================
void trace(uint32_t cellNumber)
{
    bool first = true;

    // Fetch the name of the file we're going to open
    const char* filename = config.output_file.c_str();

    // Open the file we're going to read, and complain if we can't
    FILE* ifile = fopen(filename, "r");
    if (ifile == nullptr) throwRuntime("Can't create %s", filename);

    // Allocate sufficient RAM to contain an entire data frame
    unique_ptr<uint8_t> framePtr(new uint8_t[config.cells_per_frame]);

    // Get a pointer to the frame data
    uint8_t* frame = framePtr.get();

    // Loop through each frame of the file...
    while (fread(frame, 1, config.cells_per_frame, ifile) == config.cells_per_frame)
    {
        
        // If this isn't the first value we've output, print a comma separator
        #if 0
        if (!first) printf(", ");
        first = false;
        #endif
        
        // And display the value of the cell number that was specified on the command line
        printf("%d\n", frame[cellNumber]);
    }
    
    // Terminate the line of text in the output
    printf("\n");
}
//=================================================================================================



//=================================================================================================
// readConfigurationFile() - Reads in the configuration file and populates the global "config"
//                           structure.
//=================================================================================================
void readConfigurationFile(string filename)
{
    CConfigFile cf;
    string cells_per_frame, ring_buffer_size;

    // Declare a default filename
    const char* cfilename = "sensor_frame_gen.conf";

    // If the filename passed by the caller isn't blank, that's our filename
    if (!filename.empty()) cfilename = filename.c_str();

    // Read and parse the configuration file and complain if we can't
    if (!cf.read(cfilename, false)) throwRuntime("Can't read %s", cfilename);

    // Fetch each configuration
    cf.get("cells_per_frame",     &cells_per_frame          );
    cf.get("ring_buffer_size",    &ring_buffer_size         );
    cf.get("adc_per_nucleotide",  &config.adc_per_nucleotide);
    cf.get("random_seed",         &config.random_seed       );
    cf.get("data_frames",         &config.data_frames       );
    cf.get("filler_value",        &config.filler_value      );
    cf.get("nucleotide_file",     &config.nucleotide_file   );
    cf.get("fragment_file",       &config.fragment_file     );
    cf.get("distribution_file",   &config.distribution_file );
    cf.get("output_file",         &config.output_file       );

    // Convert the scaled integer strings into binary values
    config.cells_per_frame = stringTo64(cells_per_frame);
    config.ring_buffer_size     = stringTo64(ring_buffer_size);
}
//=================================================================================================


//=================================================================================================
// getFileSize() - Returns the size (in bytes) of the input file
//
// Passed:  descriptor = The file descriptor of the file we want the size of
//
// Returns: The size of the file in bytes
//=================================================================================================
size_t getFileSize(int descriptor)
{
    // Find out how big the file is
    off64_t result = lseek64(descriptor, 0, SEEK_END);

    // Rewind back to the start of the file
    lseek(descriptor, 0, SEEK_SET);

    // And hand the size of the input-file to the caller
    return result;
}
//=================================================================================================


//=================================================================================================
// fillBuffer() - This stuffs some data into a DMA buffer
//
//                          <<< THIS ROUTINE IS A KLUDGE >>>
//
// Because of yet unresolved issues with very slow-writes to the DMA buffer, we are reading the
// file into a local user-space buffer then copying it into the DMA buffer.    For reasons we don't
// yet understand, the MMU allows us to copy a user-space buffer into the DMA space buffer faster
// than it allows us to write to it directly.
//
//                               <<< THIS IS A HACK >>>
//
// The hack will be fixed when we figure out how to write a device driver that can allocate
// very large contiguous blocks.
//  
//=================================================================================================
void fillBuffer(int fd, size_t fileSize)
{
    // We will load the file into the buffer in blocks of data this size (i.e, 1 GB)
    const uint32_t FRAME_SIZE = 0x40000000;

    // This is the number of bytes that have been loaded into the buffer
    uint64_t bytesLoaded = 0;

    // Get a pointer to the start of the contiguous buffer
    uint8_t* ptr = RAM.bptr();

    // Allocate a RAM buffer in userspace
    uint8_t* localBuffer = new uint8_t[FRAME_SIZE];

    // Compute how many bytes of data to load...
    uint64_t bytesRemaining = fileSize;

    // Display the completion percentage
    printf("Percent loaded =   0");
    fflush(stdout);

    // While there is still data to load from the file...
    while (bytesRemaining)
    {
        // We'd like to load the entire remainder of the file
        size_t blockSize = bytesRemaining;

        // We're going to load this file in chunks of no more than 1 GB
        if (blockSize > FRAME_SIZE) blockSize = FRAME_SIZE;

        // Load this chunk of the file into our local user-space buffer
        size_t rc = read(fd, localBuffer, blockSize);
        if (rc != blockSize)
        {
            perror("\nread");
            exit(1);
        }

        // Copy the userspace buffer into the contiguous block of physical RAM
        memcpy(ptr, localBuffer, blockSize);

        // Bump the pointer to where the next chunk will be stored
        ptr += blockSize;

        // And keep track of how many bytes are left to load
        bytesRemaining -= blockSize;

        // Compute and display the completion percentage
        bytesLoaded += blockSize;
        int pct = 100 * bytesLoaded / fileSize;
        printf("\b\b\b%3i", pct);
        fflush(stdout);
    }

    // Finish the "percent complete" display
    printf("\b\b\b100\n");

    // Free up the localBuffer so we don't leak memory
    delete[] localBuffer;
}
//=================================================================================================


//=================================================================================================
// stringTo64() - Converts a character string to a 64-bit integer after stripping out any
//                underscore characters from the input string.  Also scales the return value
//                according to any K, M, or G suffix on the string
//=================================================================================================
uint64_t stringTo64(const string& str)
{
    char buffer[100], *out = buffer;
    uint64_t multiplier = 0;

    // Point to the beginning of the input string
    const char* in = str.c_str();

    // This loop is going to strip underscores from the input string
    while (*in)
    {
        // If we've reached the end of the input string, we're done
        if (*in == 0) break;

        // If we're at the last character of the output buffer, we're done
        if ((out - buffer) == sizeof(buffer)-1) break;

        // If the input character is an underscore, skip it
        if (*in == '_')
        {
            ++in;
            continue;
        }

        // Move the input character to the output buffer
        *out++ = *in++;
    }

    // Nul-terminate the output string
    *out = 0;

    // If the buffer is empty, the result is 0
    if (buffer[0] == 0) return 0;

    // Fetch the final character
    char letter = *(out-1);

    // Use the suffix (if any) to determine the multiplier
         if (letter >= '0' && letter <= '9') multiplier = 1;
    else if (letter >= 'a' && letter <= 'f') multiplier = 1;
    else if (letter >= 'A' && letter <= 'F') multiplier = 1;
    else if (letter == 'K') multiplier = 1024;
    else if (letter == 'M') multiplier = 1024 * 1024;
    else if (letter == 'G') multiplier = 1024 * 1024 * 1024;

    // If an invalid suffix was specified, whine about it
    if (multiplier == 0) throwRuntime("Invalid suffix on %s", str.c_str());

    // If there was a suffix on the string, remove it
    if (multiplier > 1) *(out-1) = 0;

    // Convert the ASCII string to a numeric value
    int64_t value = stoull(buffer, 0, radix(buffer));

    // Hand the resulting value to the caller
    return value * multiplier;
}
//=================================================================================================



//=================================================================================================
// loadFile() - Loads a data-file into RAM at a defined physical address
//=================================================================================================
void loadFile(string filename, string address)
{
    // Ensure that we're running as the root user
    if (geteuid() != 0) throw runtime_error("Must be root to run.  Use sudo.");

    // Open the data file
    int fd = open(filename.c_str(), O_RDONLY);
    if (fd < 0) throwRuntime("Can't open '%s'", filename.c_str());

    // Find out how big the input file is
    size_t fileSize = getFileSize(fd);

    // Find out how large our contiguous buffer is
    size_t sizeLimit = stringTo64(cmdLine.sizeLimit);

    // Ensure that the file doesn't exceed the size of the buffer it's being loaded in to
    if (fileSize > sizeLimit) throwRuntime("%s is too big to fit into buffer", filename.c_str());

    // Convert the string-form of the address to binary
    uint64_t physAddr = stringTo64(address);

    // Make at least minimal effort to ensure the user doesn't blow away their system
    if (physAddr == 0) throwRuntime("Loading to RAM address 0 not permitted");

    // Tell the user what we're doing...
    printf("Mapping RAM...\n");

    // Map the physical RAM space
    RAM.map(physAddr, fileSize);

    // Tell the user what's taking so long...
    printf("Loading %s into RAM at address %s\n", filename.c_str(), address.c_str());

    // Load the data file into the RAM buffer
    fillBuffer(fd, fileSize);

    // Close the input file, we're done
    close(fd);
}
//=================================================================================================



//=================================================================================================
// printDictionary() - Display the name of each fragment along with its length (in frames), then
//                     display every fragment sequence along with its length (in frames)
//=================================================================================================
void printDictionary()
{
    const char* name;
    int length;
    char buffer[100];

    // Display a legend for the fragment dictionary
    printf("\n");
    printf("%30s    Size\n", "Fragment Name");
    printf("------------------------------------------\n");

    // Display the name and length of each fragment
    for (auto it=fragment.begin(); it != fragment.end(); ++it)
    {
        name = it->first.c_str();
        length = it->second.size();
        printf("%30s %7i\n", name, length);
    }

    // Leave a couple of blank lines between the fragements and the sequences
    printf("\n\n");
    printf("%30s    Size\n", "Distribution Name");
    printf("------------------------------------------\n");


    // The name of the distribution will be stuffed into "buffer"
    name = buffer;
    
    // Display the name and length of each distribution sequence
    for (auto d : distributionList)
    {
        sprintf(buffer, "%i,%i,%i", d.first, d.last, d.step);
        length = d.cellValue.size();
        printf("%30s %7i\n", name, length);        
    }
}
//=================================================================================================




