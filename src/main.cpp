
#include <iostream>

#ifdef  _WIN32
#include "win\dirent.h"
#else
#include "dirent.h"
#endif

#include "WfipsGrid.h"
#include "ogrsf_frmts.h"
#include "gdal_alg.h"

#include <algorithm>
#include <omp.h>
#include <set>
#include <vector>
#include <unordered_map>
#include <cstdio>
#include <chrono>
#include <gdal_utils.h>

#include "sqlite3.h"

#include "WfipsGridData.h"

#include "MyPolygon.h"

#include "FileGDBAPI.h"

#ifndef EQUAL
#define EQUAL(a,b) (strcmp(a,b)==0)
#endif

using namespace FileGDBAPI;
using std::multimap;
using std::string;
using std::vector;
using std::chrono::steady_clock;

/***************************************************************************
// Global Constants
***************************************************************************/
static const double epsilon = 0.000000000000000001;
static const double cellWidthInMeters = 2000.0; // 2 km
static const double numMicrosecondsPerSecond = 1000000.0;

/***************************************************************************
// Data Structures
***************************************************************************/

struct WfipsData
{
    OGRSpatialReference spatialReference;
    vector<SBoundingBox> cellBoundingBoxes;
    CWfipsGridData gridData;

    double GeoTransform[6];

    int numRows = -1;
    int numCols = -1;
};

struct FireshedData
{
    vector<unordered_multimap<int, int>> wfipscellsToFireOriginsForSingleGDB;
    vector<vector<int>> originCellsForWfipscell;
    unordered_map<int, int> finalIndexToWfipsCellMap;
    vector<vector<int>> numWfipscellOriginPairs;
    map<int, int> finalOriginCellToTotalPairCountMap;
};

SBoundingBox GetBoundingBox(vector<MyPoint> theRingString);

int FindCellOriginPairsInRAM(const int numThreads, const bool verbose, const steady_clock::time_point startClock, const int numEdgeCellDivisions, string& gdbPath, vector<string>& gdbNameList, FireshedData& fireshedData, WfipsData& wfipsData);
void ConsolidateFinalData(const int num_shape_files, FireshedData& fireshedData);
int FillWfipsData(WfipsData& wfipsData, std::string dataPath);
int CreateFireShedDB(const bool verbose, sqlite3* db, FireshedData& fireshedData, WfipsData& wfipsData);
void Usage();

fgdbError ReadGeometryFromGDBTable(Table& pyromeTable, WfipsData& wfipsData, vector<MyMultiPolygon>& fireMultipolygons, vector<SBoundingBox>& fireBoundingBoxes, vector<int>& fireNumbers, vector<int>& originCells);

bool AreClose(double a, double b);
bool IsDirectory(const string& path);

enum ESRIExtent
{
    minX = 0,
    minY = 1,
    maxX = 2,
    maxY = 3
};

int main(int argc, char *argv[])
{
    DIR *dir;
    struct dirent *ent;

    string fileOrDirectoryName = "";
    string extension = "";
    size_t pos = -1;
    int type = 0;
    string dataPath = "";
    string outPath = "";
    string gridPath = "";
    string numThreadsString = "";
    int rc = 0;
    bool verboseParameter = false;
    const int SUCCESS = 0;

    const int max_argument_index = argc - 1;
    int argIndex = 1;
    bool isOutPathSpecified = false;
    bool isGridPathSpecified = false;
    bool isNumThreadsSpecified = false;

    // Parse commandline arguments
    if (argc > 1)
    {
        while (argIndex < argc)
        {
            if (EQUAL(argv[argIndex], "--gdb"))
            {
                if ((argIndex + 1) > max_argument_index) // An error has occurred
                {
                    // Report error
                    printf("ERROR: No path to geodatabase entered\n");
                    Usage(); // Exits program
                }
                dataPath = argv[++argIndex];
            }
            else if (EQUAL(argv[argIndex], "--out"))
            {
                if ((argIndex + 1) > max_argument_index) // An error has occurred
                {
                    // Report error
                    printf("ERROR: No output path entered\n");
                    Usage(); // Exits program
                }
                isOutPathSpecified = true;
                outPath = argv[++argIndex];
            }
            else if (EQUAL(argv[argIndex], "--grid"))
            {
                if ((argIndex + 1) > max_argument_index) // An error has occurred
                {
                    // Report error
                    printf("ERROR: No grid path entered\n");
                    Usage(); // Exits program
                }
                isGridPathSpecified = true;
                gridPath = argv[++argIndex];
            }
            else if (EQUAL(argv[argIndex], "--t"))
            {
                if ((argIndex + 1) > max_argument_index) // An error has occurred
                {
                    // Report error
                    printf("ERROR: No thread number entered\n");
                    Usage(); // Exits program
                }
                isNumThreadsSpecified = true;
                numThreadsString = argv[++argIndex];
            }
            else if (EQUAL(argv[argIndex], "--verbose"))
            {
                verboseParameter = true;
            }
            else
            {
                printf("ERROR: %s is an invalid argument\n", argv[argIndex]);
                Usage(); // Exits program
            }
            argIndex++;
        }
    }
    else
    {
        printf("ERROR: no shapefile path given\n", argv[argIndex]);
        Usage(); // Exits program
    }

    const bool verbose = verboseParameter;

    std::replace(dataPath.begin(), dataPath.end(), '\\', '/'); // replace all '\\' with '/'
    if ((dataPath.back() != '/') && (dataPath.back() != '\\'))
    {
        dataPath.push_back('/');
    }

    if (isGridPathSpecified)
    {
        std::replace(gridPath.begin(), gridPath.end(), '\\', '/'); // replace all '\\' with '/'
    }
    if (!isGridPathSpecified)
    {
        gridPath = dataPath;
    }
    else if ((gridPath.back() != '/') && (gridPath.back() != '\\'))
    {
        gridPath.push_back('/');
    }

    if (isOutPathSpecified)
    {
        std::replace(outPath.begin(), outPath.end(), '\\', '/'); // replace all '\\' with '/'
    }
    if (!isOutPathSpecified)
    {
        outPath = dataPath;
    }
    else if ((outPath.back() != '/') && (outPath.back() != '\\'))
    {
        outPath.push_back('/');
    }

    int numThreadsArg = -1;
    if (isNumThreadsSpecified)
    {
        numThreadsArg = atoi(numThreadsString.c_str());
    }

    if (isNumThreadsSpecified && (numThreadsArg < 1))
    {
        printf("ERROR: Number of threads must be greater than zero\n", argv[argIndex]);
        Usage(); // Exits program
    }

    vector<string> gdbNameList;
    vector<string> gdbFullPathList;

    if ((dir = opendir(dataPath.c_str())) != nullptr)
    {
        /* read all the files and directories within directory */
        while ((ent = readdir(dir)) != nullptr)
        {
            fileOrDirectoryName = ent->d_name;

            string directoryName = fileOrDirectoryName;
            pos = directoryName.find_last_of('.');
            if (pos != string::npos)
            {
                extension = directoryName.substr(pos, directoryName.size());
                std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
                if (extension == ".gdb")
                {
                    //printf("%s\n", fileName.c_str());
                    gdbNameList.push_back(directoryName);
                    gdbFullPathList.push_back(dataPath + directoryName);
                }
            }            
        }
        closedir(dir);
    }
    else
    {
        /* could not open directory */
        perror("");
        return EXIT_FAILURE;
    }

    if (gdbFullPathList.size() < 1)
    {
        printf("Error: No geodatabase directories found, exiting program\n");
        return EXIT_FAILURE;
    }

    const int gdbListSize = gdbFullPathList.size();

    sqlite3 *db = nullptr;

    string outPutDbFullPath = outPath + "firesheds.db";
    rc = sqlite3_open_v2(outPutDbFullPath.c_str(), &db,
        SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE,
        NULL);

    if (rc != SQLITE_OK || db == nullptr)
    {
        printf("Error: Could not create firesheds.db (does the output directory exist?)\nexiting program\n");
        return EXIT_FAILURE;
    }

    GDALAllRegister();

    char *MapESRIProjStrings[] =
    {
        "",
        "PROJCS[\"Albers\",GEOGCS[\"GCS_North_American_1983\",DATUM[\"D_North_American_1983\",SPHEROID[\"GRS_1980\",6378137,298.257222101]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]],PROJECTION[\"Albers\"],PARAMETER[\"standard_parallel_1\",29.5],PARAMETER[\"standard_parallel_2\",45.5],PARAMETER[\"latitude_of_origin\",23],PARAMETER[\"central_meridian\",-96],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"METERS\",1]]"
    };

    CPLStringList papszPrj5070;
    papszPrj5070.AddString(MapESRIProjStrings[1]);

    WfipsData wfipsData;
    rc = wfipsData.spatialReference.importFromESRI(papszPrj5070);

    if (verbose)
    {
        printf("Loading WFIPS data\n");
    }

    rc = FillWfipsData(wfipsData, dataPath);
    if (rc != SUCCESS)
    {
        printf("Error: WFIPS data loading failed\n");
        printf("    Make sure the file \"WFIPSGrid.tif\" exists in\n    %s\n\n", dataPath.c_str());
        return EXIT_FAILURE;
    }
    else if (verbose)
    {
        printf("WFIPS data populated\n\n");
    }

    printf("Processing all geodatabases in\n    %s\n\nPlease wait...\n\n", dataPath.c_str());

    FireshedData fireshedData;
    fireshedData.wfipscellsToFireOriginsForSingleGDB.resize(gdbListSize);

    steady_clock::time_point startClock = std::chrono::steady_clock::now();
    const int numCellEdgeDivisions = 16; // Number of times to subdivide WFIPS cells (currently 2000m, makes 125m subcells) 
                                         // to get closer to the resolution used in FSim (apparently ~135m)

    FindCellOriginPairsInRAM(numThreadsArg, verbose, startClock, numCellEdgeDivisions, dataPath, gdbNameList, fireshedData, wfipsData);

    ConsolidateFinalData(gdbListSize, fireshedData);
    CreateFireShedDB(verbose, db, fireshedData, wfipsData);

    long long int elapsedTimeInMicroseconds = std::chrono::microseconds(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startClock)).count();
    printf("Successfully processed all geodatabases in\n    %s\n", dataPath.c_str());
    printf("    and created \"firesheds.db\" in\n    %s\n\n", outPath.c_str());
    printf("Total time elapsed is %4.2f seconds\n\n", elapsedTimeInMicroseconds / numMicrosecondsPerSecond);

    return SUCCESS;
}

/***************************************************************************
// Shapefile Reading Function
***************************************************************************/
int FindCellOriginPairsInRAM(const int numThreadsArg, const bool verbose, const steady_clock::time_point startClock, const int numEdgeCellDivisions, string& gdbPath, vector<string>& gdbNameList, FireshedData& fireshedData, WfipsData& wfipsData)
{
    fgdbError hr;

    int numThreads = 1;
    const int availableThreads = omp_get_num_procs();
    if (numThreadsArg == -1)
    {
        numThreads = availableThreads;
    }
    else
    {
        if (numThreadsArg > availableThreads)
        {
            numThreads = availableThreads;
        }
        else
        {
            numThreads = numThreadsArg;
        }
    }

    omp_set_num_threads(numThreads);

    volatile int numPyromesProcessed = 0;

    const int gdbListSize = gdbNameList.size();
    int numPyromes = 0;
    // Get total count of pyromes
    for (int gdbIndex = 0; gdbIndex < gdbListSize; gdbIndex++)
    {
        // Open the geodatabase.
        Geodatabase geodatabase;
        std::replace(gdbPath.begin(), gdbPath.end(), '\\', '/'); // replace all '\\' with '/'
        wstring wstrDataPath(gdbPath.begin(), gdbPath.end());
        wstring wstrGdbName(gdbNameList[gdbIndex].begin(), gdbNameList[gdbIndex].end());
        wstring wstrFullGdbPath = wstrDataPath + wstrGdbName;

        if ((hr = OpenGeodatabase(wstrFullGdbPath, geodatabase)) != S_OK)
        {
            printf("An error occurred while opening the geodatabase\n");
            printf("ESRI Error Code ( %d )\n", hr);
        }
        else
        {
            // Get all layer names in current gdb
            vector<wstring> pyromeNames;
            hr = geodatabase.GetChildDatasets(L"\\", L"", pyromeNames); // root path is "\\", empty second parameter returns all
            numPyromes += pyromeNames.size();
        }

        // Close the geodatabase
        if ((hr = CloseGeodatabase(geodatabase)) != S_OK)
        {
            printf("An error occurred while closing the geodatabase\n");
        }
    }

    for (int gdbIndex = 0; gdbIndex < gdbListSize; gdbIndex++)
    {
        // Open the geodatabase.
        Geodatabase geodatabase;
        std::replace(gdbPath.begin(), gdbPath.end(), '\\', '/'); // replace all '\\' with '/'
        wstring wstrDataPath(gdbPath.begin(), gdbPath.end());
        wstring wstrGdbName(gdbNameList[gdbIndex].begin(), gdbNameList[gdbIndex].end());
        wstring wstrFullGdbPath = wstrDataPath + wstrGdbName;

        if ((hr = OpenGeodatabase(wstrFullGdbPath, geodatabase)) != S_OK)
        {
            printf("An error occurred while opening the geodatabase\n");
            printf("ESRI Error Code ( %d )\n", hr);
        }
        else
        {
            // Get all layer names in current gdb
            vector<wstring> pyromeNames;
            hr = geodatabase.GetChildDatasets(L"\\", L"", pyromeNames); // root path is "\\", empty second parameter returns all

            for (int pyromeIndex = 0; pyromeIndex < pyromeNames.size(); pyromeIndex++)
            {
                vector<MyMultiPolygon> fireMultiPolygons;
                vector<SBoundingBox> fireBoundingBoxes;
                vector<int> fireOriginCells;
                vector<int> fireNumbers;
                Table pyromeTable;
                // Begin load only mode. This shuts off the update of all indexes. Set a write lock.
                 // Setting the write lock will improve performance when bulk loading.
                pyromeTable.LoadOnlyMode(true);
                pyromeTable.SetWriteLock();

                string currentPyromeName(pyromeNames[pyromeIndex].begin(), pyromeNames[pyromeIndex].end());
                std::remove(currentPyromeName.begin(), currentPyromeName.end(), '\\');
                currentPyromeName.pop_back(); // weird issue with repeated last character from wstring

                // Open the Pyrome table.
                if ((hr = geodatabase.OpenTable(pyromeNames[pyromeIndex], pyromeTable)) != S_OK)
                {
                    printf("An error occurred while opening geodatabase table %s\n", currentPyromeName.c_str());
                }
                else
                {
                    if(verbose)
                    {
                        printf("Reading geometry data from table %s\n", currentPyromeName.c_str());
                    }
                    ReadGeometryFromGDBTable(pyromeTable, wfipsData, fireMultiPolygons, fireBoundingBoxes, fireNumbers, fireOriginCells);
                    // Close the table
                    if ((hr = geodatabase.CloseTable(pyromeTable)) != S_OK)
                    {
                        printf("An error occurred while closing the table\n");
                    }

                    unordered_multimap<int, int> tempTotalWfipscellsToFireOrigins;

                    #pragma omp parallel for schedule(dynamic, 1) shared(numPyromesProcessed)
                    for (int fireIndex = 0; fireIndex < fireMultiPolygons.size(); fireIndex++)
                    {
                        int fireNumber = fireNumbers[fireIndex];
                        int origin = fireOriginCells[fireIndex];

                        unordered_map<int, int> tempSingleFireWfipscellsToFireOrigins; // used to check if cell is already burned by fire

                            // Get the bounding box for the current fire
                        SBoundingBox fireBoundingBox = fireBoundingBoxes[fireIndex];

                        // Always assume cell containing the origin is burned by the fire
                        tempSingleFireWfipscellsToFireOrigins.insert(std::make_pair(origin, origin));
                        for (int subPolygonIndex = 0; subPolygonIndex < fireMultiPolygons[fireIndex].polygonPart.size(); subPolygonIndex++)
                        {
                            int upperLeftCell = wfipsData.gridData.WG_GetCellIndex(fireBoundingBox.minX, fireBoundingBox.maxY);
                            int lowerRightCell = wfipsData.gridData.WG_GetCellIndex(fireBoundingBox.maxX, fireBoundingBox.minY);

                            int minRow = 0;
                            int maxRow = 0;
                            int minCol = 0;
                            int maxCol = 0;

                            // Find rows and columns that correspond to the upper-left and lower-right corners of fire bounding box
                            wfipsData.gridData.GetRowColFromWfipsCell(upperLeftCell, &minRow, &minCol);
                            wfipsData.gridData.GetRowColFromWfipsCell(lowerRightCell, &maxRow, &maxCol);

                            int cellIndex = 0;
                            for (int row = minRow; row <= maxRow; row++)
                            {
                                for (int col = minCol; col <= maxCol; col++)
                                {
                                    cellIndex = wfipsData.gridData.GetWfipsCellIndex(row, col);

                                    // Below each cell within the current fire's bounding box not already marked as burned by that fire
                                    // is subdivided into numEdgeCellDivisions subcells. The nodes along the outer edges of the subcells
                                    // are then used as test points in Mark's point in poly function IsOverlapping(). We test only points
                                    // along the edge as we only need to check where the fire crosses into a WFIPS cell from another cell
                                    // as origin cells are always assumed to be burned by a fire. If any node on the edge of the subcells
                                    // for the current WFIPS cell is found to be within the fire, we mark that cell as burned and move on
                                    // to the next cell within the fire's bounding box.

                                    // Need to check if cell has already been burned by one of the current fire's subpolygons so as not to
                                    // count cells as being burned by the same fire more than once  
                                    bool isCellAlreadyBurned = !(tempSingleFireWfipscellsToFireOrigins.find(cellIndex) == tempSingleFireWfipscellsToFireOrigins.end());
                                    if (!isCellAlreadyBurned)
                                    {
                                        SBoundingBox cellBoundingBox = wfipsData.cellBoundingBoxes[cellIndex];

                                        const double xMin = cellBoundingBox.minX;
                                        const double yMin = cellBoundingBox.minY;
                                        const double xMax = cellBoundingBox.maxX;
                                        const double yMax = cellBoundingBox.maxY;

                                        MyPoint testPoint;
                                        testPoint.x = xMin;
                                        testPoint.y = yMin;
                                        const int numEdgeNodes = numEdgeCellDivisions + 1;
                                        const double nodeStepSizeInMeters = cellWidthInMeters * (1.0 / numEdgeCellDivisions);

                                        bool isCellInFirePolygon = false;
                                        double xEdgeNodeCoord = 0;
                                        double yEdgeNodeCoord = 0;

                                        const int bottomEdge = 0; // yNodeIndex value for top edge of cell
                                        const int topEdge = numEdgeCellDivisions; // yNodeIndex value for bottom edge of WFIPS cell
                                        const int leftEdge = 0; // xNodeIndex value for left edge of cell
                                        const int rightEdge = numEdgeCellDivisions; // xNodeIndex value for right edge of WFIPS cell

                                        for (int xNodeIndex = 0; xNodeIndex < numEdgeNodes; xNodeIndex++)
                                        {
                                            if (!isCellInFirePolygon)
                                            {
                                                xEdgeNodeCoord = xMin + (nodeStepSizeInMeters * xNodeIndex);
                                                testPoint.x = xEdgeNodeCoord;
                                                for (int yNodeIndex = 0; yNodeIndex < numEdgeNodes; yNodeIndex++)
                                                {
                                                    yEdgeNodeCoord = yMin + (nodeStepSizeInMeters * yNodeIndex);
                                                    testPoint.y = yEdgeNodeCoord;

                                                    if ((yNodeIndex == topEdge) || (yNodeIndex == bottomEdge))
                                                    {
                                                        // Check all subcell nodes along the cell top or bottom egdes of WFIPS cell  
                                                        isCellInFirePolygon = MyPolygonUtility::IsOverlapping(testPoint,
                                                            fireMultiPolygons[fireIndex].polygonPart[subPolygonIndex].pointString);
                                                        if (isCellInFirePolygon)
                                                        {
                                                            break; // Leave yNode loop
                                                        }
                                                    }
                                                    else
                                                    {
                                                        // If not along top or bottom edge of WFIPS cell, check only subcell nodes on left and right edges
                                                        if ((xNodeIndex == leftEdge) || (xNodeIndex == rightEdge))
                                                        {
                                                            isCellInFirePolygon = MyPolygonUtility::IsOverlapping(testPoint,
                                                                fireMultiPolygons[fireIndex].polygonPart[subPolygonIndex].pointString);
                                                            if (isCellInFirePolygon)
                                                            {
                                                                break; // Leave yNode loop
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            else
                                            {
                                                break;  // Leave xNode loop
                                            }
                                        }
                                        if (isCellInFirePolygon)
                                        {
                                            // Mark current cell as being burned by the fire
                                            tempSingleFireWfipscellsToFireOrigins.insert(std::make_pair(cellIndex, origin));
                                        }
                                    }
                                }
                            }
                        }
                  
                        // Update shared wfipscells to origin map for current geodatabase
                        #pragma omp critical
                        {
                            if (verbose & (fireIndex % 1000 == 0))
                            {
                                double currentPyromeProgress = (fireIndex / (fireMultiPolygons.size() * 1.0)) * 100.00;
                                double totalProgress = (numPyromesProcessed / (numPyromes * 1.0)) * 100.00;
                                long long int elapsedTimeInMicroseconds = std::chrono::microseconds(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startClock)).count();

                                printf("%4.2f percent of pyrome %s\n    is complete\n", currentPyromeProgress, currentPyromeName.c_str());
                                printf("Processed %d pyromes out of %d\n", numPyromesProcessed, numPyromes);
                                printf("%4.2f percent of all pyromes is complete\n", totalProgress);
                                printf("Total time elapsed is %4.2f seconds\n\n", elapsedTimeInMicroseconds / numMicrosecondsPerSecond);
                            }

                            int wfipscell;
                            int origin;
                            for (auto iterator = tempSingleFireWfipscellsToFireOrigins.begin(); iterator != tempSingleFireWfipscellsToFireOrigins.end(); iterator++)
                            {
                                wfipscell = iterator->first;
                                origin = iterator->second;
                                fireshedData.wfipscellsToFireOriginsForSingleGDB[gdbIndex].insert(std::make_pair(wfipscell, origin));
                            }
                        }
                        // Clear list of cells burned by current fire for next loop iteration
                        tempSingleFireWfipscellsToFireOrigins.clear();
                    }

                    #pragma omp atomic
                    numPyromesProcessed++;

                    if (verbose)
                    {
                        double totalProgress = (numPyromesProcessed / (numPyromes * 1.0)) * 100.00;
                        long long int elapsedTimeInMicroseconds = std::chrono::microseconds(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startClock)).count();
                        printf("Processed %d pyromes out of %d in\n", numPyromesProcessed, numPyromes);
                        printf("    %s\n    %4.2f percent of all pyromes to be processed are complete\n", gdbPath.c_str(), totalProgress);
                        printf(" Total time elapsed is %4.2f seconds\n\n", elapsedTimeInMicroseconds / numMicrosecondsPerSecond);
                    }


                    fireNumbers.clear();
                    fireOriginCells.clear();
                    tempTotalWfipscellsToFireOrigins.clear();
                }
            }

            // Close the geodatabase
            if ((hr = CloseGeodatabase(geodatabase)) != S_OK)
            {
                printf("An error occurred while closing the geodatabase\n");
            }
        }
    }
 
    return 0;
}

int FillWfipsData(WfipsData& wfipsData, std::string dataPath)
{
    int rc = wfipsData.gridData.LoadData((const char*)dataPath.c_str());

    wfipsData.numRows = wfipsData.gridData.GetNumY();
    wfipsData.numCols = wfipsData.gridData.GetNumX();

    wfipsData.gridData.GetWfipsGridGeotransform(wfipsData.GeoTransform);

    SBoundingBox currentCellBoundingBox;

    int cellIndex = 0;
    for (int row = 0; row < wfipsData.numRows; row++)
    {
        for (int col = 0; col < wfipsData.numCols; col++)
        {
            cellIndex = wfipsData.gridData.GetWfipsCellIndex(row, col);
            wfipsData.gridData.WG_GetCellCoords(cellIndex,
                &currentCellBoundingBox.minX,
                &currentCellBoundingBox.minY,
                &currentCellBoundingBox.maxX,
                &currentCellBoundingBox.maxY);
            wfipsData.cellBoundingBoxes.push_back(currentCellBoundingBox);
        }
    }

    return rc;
}

void ConsolidateFinalData(const int num_gdbs, FireshedData& fireshedData)
{
    multiset<pair<int, int>> totalWfipscellsToFireOriginPairs;
    vector<pair<int, int>> totalPairsForSingleOrigin;
    int wfipscell = -1;
    int origin = -1;
    int wfipscellPrevious = -1;
    int originPrevious = -1;
    int wfipsCellIndex = -1;
    int numPairs = -1;
    int numTotalPairs = -1;

    for (int gdbIndex = 0; gdbIndex < num_gdbs; gdbIndex++)
    {
        for (auto iterator = fireshedData.wfipscellsToFireOriginsForSingleGDB[gdbIndex].begin(); iterator != fireshedData.wfipscellsToFireOriginsForSingleGDB[gdbIndex].end(); iterator++)
        {
            wfipscell = iterator->first;
            origin = iterator->second;
            totalWfipscellsToFireOriginPairs.insert(std::make_pair(wfipscell, origin));
        }

        fireshedData.wfipscellsToFireOriginsForSingleGDB[gdbIndex].clear();
    }

    // Get the data for num_pairs field
    for (auto iterator = totalWfipscellsToFireOriginPairs.begin(); iterator != totalWfipscellsToFireOriginPairs.end(); iterator++)
    {
        wfipscell = iterator->first;
        origin = iterator->second;
        bool wfipsCellChanged = false;
        bool originChanged = false;

        if (wfipscell != wfipscellPrevious) // Wfispcell changed, add new vector rows
        {
            wfipsCellChanged = true;
            wfipscellPrevious = wfipscell;
            wfipsCellIndex++;
            // Add current wfipscell index to map so it can be retrieved by row index later
            fireshedData.finalIndexToWfipsCellMap.insert(std::make_pair(wfipsCellIndex, wfipscell));
            vector<int> originVectorRow;
            fireshedData.originCellsForWfipscell.push_back(originVectorRow);
            vector<int> numPairsVectorRow;
            fireshedData.numWfipscellOriginPairs.push_back(numPairsVectorRow);
        }
        if (origin != originPrevious) // Origin changed
        {
            originChanged = true;
            originPrevious = origin;
        }

        if (originChanged || wfipsCellChanged) // Add element to current vector row
        {
            numPairs = 1;
            fireshedData.originCellsForWfipscell[wfipsCellIndex].push_back(origin);
            fireshedData.numWfipscellOriginPairs[wfipsCellIndex].push_back(numPairs);
        }
        else // Increment count of pairs and overwrite element containing count for current wfipscell
        {
            numPairs++;
            fireshedData.numWfipscellOriginPairs[wfipsCellIndex].back() = numPairs;
        }
    }

    // Get the data for total_for_origin field
    for (int wfipsCellIndex = 0; wfipsCellIndex < fireshedData.finalIndexToWfipsCellMap.size(); wfipsCellIndex++)
    {
        int wfipscell = fireshedData.finalIndexToWfipsCellMap.at(wfipsCellIndex);
        for (int originCellIndex = 0; originCellIndex < fireshedData.originCellsForWfipscell[wfipsCellIndex].size(); originCellIndex++)
        {
            origin = fireshedData.originCellsForWfipscell[wfipsCellIndex][originCellIndex];
            if (origin == wfipscell) // Fire started in that cell
            {
                int numTotalPairs = fireshedData.numWfipscellOriginPairs[wfipsCellIndex][originCellIndex];
                fireshedData.finalOriginCellToTotalPairCountMap.insert(std::make_pair(origin, numTotalPairs));
            }
        }
    }

    totalWfipscellsToFireOriginPairs.clear();
}

int CreateFireShedDB(const bool verbose, sqlite3* db, FireshedData& fireshedData, WfipsData& wfipsData)
{
    int rc = 0;
    sqlite3_stmt *fireshedsStmt;
    sqlite3_stmt *totalsStmt;
    char *sqlErrMsg = 0;

    string createFireshedDBSQLString = "CREATE TABLE IF NOT EXISTS firesheds " \
        "(wfipscell INTEGER, " \
        "origin INTEGER, " \
        "num_pairs INTEGER)";
    
    string insertFireshedsSQLString = "INSERT INTO firesheds(" \
        "wfipscell, " \
        "origin, " \
        "num_pairs) " \
        "VALUES " \
        "(:wfipscell, " \
        ":origin, " \
        ":num_pairs)";

    string createTotalFiresDBSQLString = "CREATE TABLE IF NOT EXISTS totals_for_origins " \
        "(origin INTEGER, " \
        "total_fires INTEGER)";

    string insertTotalsFiresSQLString = "INSERT INTO totals_for_origins " \
        "(origin, " \
        "total_fires) " \
        "VALUES " \
        "(:origin, " \
        ":total_fires)";

    int wfipscell = -1;
    int bindColumnIndex = -1;
    int numPairs = -1;
    int origin = -1;


    rc = sqlite3_exec(db, "DROP TABLE IF EXISTS 'firesheds'", NULL, NULL, &sqlErrMsg);
    rc = sqlite3_exec(db, "DROP TABLE IF EXISTS 'totals_for_origins'", NULL, NULL, &sqlErrMsg);

    rc = sqlite3_exec(db, createFireshedDBSQLString.c_str(), NULL, NULL, &sqlErrMsg);
    rc = sqlite3_exec(db, createTotalFiresDBSQLString.c_str(), NULL, NULL, &sqlErrMsg);

    if (rc != SQLITE_OK)
    {
        printf("Error: Could not create schema on firesheds.db\n\n");
        return rc;
    }

    rc = sqlite3_prepare_v2(db, insertFireshedsSQLString.c_str(), -1, &fireshedsStmt, NULL); // Prepare SQL statement
    if (rc != SQLITE_OK)
    {
        printf("Error: Could not create prepare SQL statement\n\n");
        return rc;
    }

    rc = sqlite3_prepare_v2(db, insertTotalsFiresSQLString.c_str(), -1, &totalsStmt, NULL); // Prepare SQL statement
    if (rc != SQLITE_OK)
    {
        printf("Error: Could not create prepare total fires SQL statement\n\n");
        return rc;
    }

    rc = sqlite3_exec(db, "PRAGMA SYNCHRONOUS=OFF", NULL, NULL, NULL);
    rc = sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &sqlErrMsg);

    if (verbose)
    {
        printf("Creating \"firesheds.db\", please wait...\n\n");
    }

    for (int wfipsCellIndex = 0; wfipsCellIndex < fireshedData.originCellsForWfipscell.size(); wfipsCellIndex++)
    {
        for (int originCellIndex = 0; originCellIndex < fireshedData.originCellsForWfipscell[wfipsCellIndex].size(); originCellIndex++)
        {
            wfipscell = fireshedData.finalIndexToWfipsCellMap.at(wfipsCellIndex);
            bindColumnIndex = sqlite3_bind_parameter_index(fireshedsStmt, ":wfipscell");
            rc = sqlite3_bind_int(fireshedsStmt, bindColumnIndex, wfipscell);

            origin = fireshedData.originCellsForWfipscell[wfipsCellIndex][originCellIndex];
            bindColumnIndex = sqlite3_bind_parameter_index(fireshedsStmt, ":origin");
            rc = sqlite3_bind_int(fireshedsStmt, bindColumnIndex, origin);

            numPairs = fireshedData.numWfipscellOriginPairs[wfipsCellIndex][originCellIndex];
            bindColumnIndex = sqlite3_bind_parameter_index(fireshedsStmt, ":num_pairs");
            rc = sqlite3_bind_int(fireshedsStmt, bindColumnIndex, numPairs);

            rc = sqlite3_step(fireshedsStmt);
            rc = sqlite3_reset(fireshedsStmt);
        }
    }

    fireshedData.finalIndexToWfipsCellMap.clear();
    fireshedData.originCellsForWfipscell.clear();
    fireshedData.numWfipscellOriginPairs.clear();

    rc = sqlite3_reset(fireshedsStmt);
    rc = sqlite3_finalize(fireshedsStmt);

    int totalForOrigin = -1;
    for (auto iterator = fireshedData.finalOriginCellToTotalPairCountMap.begin(); iterator != fireshedData.finalOriginCellToTotalPairCountMap.end(); iterator++)
    {
        origin = iterator->first;
        bindColumnIndex = sqlite3_bind_parameter_index(totalsStmt, ":origin");
        rc = sqlite3_bind_int(totalsStmt, bindColumnIndex, origin);

        totalForOrigin = iterator->second;
        bindColumnIndex = sqlite3_bind_parameter_index(totalsStmt, ":total_fires");
        rc = sqlite3_bind_int(totalsStmt, bindColumnIndex, totalForOrigin);

        rc = sqlite3_step(totalsStmt);
        rc = sqlite3_reset(totalsStmt);
    }

    rc = sqlite3_reset(totalsStmt);
    rc = sqlite3_finalize(totalsStmt);

    rc = sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &sqlErrMsg);
    rc = sqlite3_exec(db, "PRAGMA SYNCHRONOUS=ON", NULL, NULL, NULL);

    //create indices, index naming convention: idx_table_name_field_name
    printf("Creating indices on tables...\n");
    rc = sqlite3_exec(db, "CREATE INDEX idx_firesheds_wfipscell ON firesheds(wfipscell ASC)",
        NULL, NULL, &sqlErrMsg);

    rc = sqlite3_exec(db, "CREATE INDEX idx_firesheds_origin ON firesheds(origin ASC)",
        NULL, NULL, &sqlErrMsg);

    rc = sqlite3_exec(db, "CREATE UNIQUE INDEX idx_totals_for_origins_origin ON totals_for_origins(origin ASC)",
        NULL, NULL, &sqlErrMsg);

    // "Vacuum" the database to free unused memory
    if (verbose)
    {
        printf("Optimizing database file size, please wait...\n\n");
    }
    rc = sqlite3_exec(db, "VACUUM", NULL, NULL, NULL);

    rc = sqlite3_close(db);

    sqlite3_free(sqlErrMsg);

    return rc;
}

bool AreClose(double a, double b)
{
    return fabs(a - b) < epsilon;
}

SBoundingBox GetBoundingBox(vector<MyPoint> theRingString)
{
    SBoundingBox theBoundingBox;
    theBoundingBox.minX = theRingString[0].x;
    theBoundingBox.maxX = theRingString[0].x;
    theBoundingBox.minY = theRingString[0].y;
    theBoundingBox.maxY = theRingString[0].y;

    for (int i = 0; i < theRingString.size(); i++)
    {
        if (theRingString[i].x < theBoundingBox.minX)
        {
            theBoundingBox.minX = theRingString[i].x;
        }
        if (theRingString[i].x > theBoundingBox.maxX)
        {
            theBoundingBox.maxX = theRingString[i].x;
        }
        if (theRingString[i].y < theBoundingBox.minY)
        {
            theBoundingBox.minY = theRingString[i].y;
        }
        if (theRingString[i].y > theBoundingBox.maxY)
        {
            theBoundingBox.maxY = theRingString[i].y;
        }
    }

    return theBoundingBox;
}

void Usage()
{
    printf("\nUsage:\n");
    printf("firesheds <--gdb path>  Required\n");
    printf("          [--grid path] Optional\n");
    printf("          [--out path]  Optional\n");
    printf("          [--t number]  Optional\n");
    printf("          [--verbose]   Optional\n");
    printf("--gdb <path>          Required: Specifies path to input geodatabase\n");
    printf("--grid <path>           Optional: Specifies path to WFIPSgrid.tif\n");
    printf("                        default: grid path is same as geodatabase path\n");
    printf("--out <path>            Optional: Specifies output path\n");
    printf("                        default: output path is same as geodatabase path\n");
    printf("--t <number>            Optional: Specifies number of threads\n");
    printf("                        default: number of available processors\n");
    printf("--verbose               Optional: Provides detailed progress information\n");
    printf("\n\nA valid path to FSIM shapefile data must exist\n");
    printf("\nA valid path to WFIPSgrid.tif must exist (same path as geodatabase by default)\n");
    printf("\nIf out path is specified, it must exist (same path as geodatabase by default)\n");
    printf("\nIf thread number is specified, it must be greater than zero\n");
    printf("Example usage:\n");
    printf("firesheds --gdb C:/my_geodatabases --grid C:/grid --out C:/output --t 4 --verbose\n\n");

    exit(1); // Exit with error code 1
}

fgdbError ReadGeometryFromGDBTable(Table& pyromeTable, WfipsData& wfipsData, vector<MyMultiPolygon>& fireMultipolygons, vector<SBoundingBox>& fireBoundingBoxes, vector<int>& fireNumbers, vector<int>& originCells)
{
    fgdbError hr;

    EnumRows enumRows;
    int rowCount;
    pyromeTable.GetRowCount(rowCount);

    enum field_number // retrieved from reading field information using GetFields()
    {
        fire_number = 2,
        origin_x = 9,
        origin_y = 10
    };

    if ((hr = pyromeTable.Search(L"*", L"GIS_SizeAc >= 150", true, enumRows)) != S_OK)
    {
        wcout << "An error occurred while performing the attribute query." << endl;
        return E_FAIL;
    }

    Row row;
    
    // Read in all results from query
    while ((hr = enumRows.Next(row)) == S_OK)
    {
        double fireNumber; // fire number is stored as double in the gdb
        double x;
        double y;

        std::vector<FieldDef> fieldDefs;

        row.GetFields(fieldDefs);

        row.GetDouble(field_number::fire_number, fireNumber);
        row.GetDouble(field_number::origin_x, x);
        row.GetDouble(field_number::origin_y, y);

        fireNumbers.push_back(static_cast<int>(fireNumber));
        originCells.push_back(wfipsData.gridData.WG_GetCellIndex(x, y));

        MultiPartShapeBuffer fireGeometry;
        row.GetGeometry(fireGeometry);

        int numPoints;
        fireGeometry.GetNumPoints(numPoints);

        int numParts;
        fireGeometry.GetNumParts(numParts);

        double xMin, xMax, yMin, yMax;

        double* xyExtent;
        fireGeometry.GetExtent(xyExtent);
        
        SBoundingBox boundingBox;
        boundingBox.minX = xyExtent[ESRIExtent::minX];
        boundingBox.minY = xyExtent[ESRIExtent::minY];
        boundingBox.maxX = xyExtent[ESRIExtent::maxX];
        boundingBox.maxY = xyExtent[ESRIExtent::maxY];

        fireBoundingBoxes.push_back(boundingBox);

        int* partsArray;
        fireGeometry.GetParts(partsArray);

        Point* points;
        fireGeometry.GetPoints(points);

        MyMultiPolygon multiPolygon;

        int polygonPartIndex = -1;
        int partCount = 0;
        for (int i = 0; i < numPoints; i++)
        {
            if ((partCount < numParts) && (i == partsArray[partCount]))
            {
                MyPolygonPart polygonPart;
                multiPolygon.polygonPart.push_back(polygonPart);
                partCount++;
                polygonPartIndex++;
            }
            MyPoint myPoint;
            myPoint.x = points[i].x;
            myPoint.y = points[i].y;
            multiPolygon.polygonPart[polygonPartIndex].pointString.push_back(myPoint);
        }
        fireMultipolygons.push_back(multiPolygon);
    }
    enumRows.Close();

    return S_OK;
}

bool IsDirectory(const string& path)
{
    struct stat statbuf;
    if (stat(path.c_str(), &statbuf) != 0)
    {
        return 0;
    }
    return S_ISDIR(statbuf.st_mode);
}
